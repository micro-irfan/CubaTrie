#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <zlib.h>

#include "kseq.h"         
KSEQ_INIT(gzFile, gzread)  

#include "trie.h" 
#include "utils.h"
#include "kmer.h"

static char *revcomp_new_n(const char *s, size_t n) {
    char *rc = (char*)malloc(n + 1);
    if (!rc) return NULL;
    for (size_t i = 0; i < n; ++i) {
        char c = s[n - 1 - i], r = 'N';
        if (c=='A') r='T'; else if (c=='C') r='G';
        else if (c=='G') r='C'; else if (c=='T') r='A';
        rc[i] = r;
    }
    rc[n] = '\0';
    return rc;
}

const size_t STEP = 1000000;
static const size_t MT_TASK_QUEUE_CAPACITY = 1024;
static const size_t MT_SAM_QUEUE_CAPACITY = 16384;
static const size_t MT_READ_BATCH_TARGET_BASES = 4u * 1024u * 1024u;
static const size_t MT_READ_BATCH_MAX_READS = 4096;

static void sam_write_unmapped(FILE *sam_fp,
                               const char *read_name,
                               const char *read_seq,
                               const char *read_qual,
                               size_t read_len) {
    if (!sam_fp || !read_name || !read_seq) return;
    const char *qual = (read_qual && strlen(read_qual) == read_len) ? read_qual : "*";
    size_t qname_len = name_len_no_rc_suffix(read_name);
    fprintf(sam_fp, "%.*s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
            (int)qname_len, read_name, read_seq, qual);
}

void total_hits (const kh_counter_t *m, const size_t nreads) {
    size_t total = 0;
    for (khint_t i = kh_begin(m); i != kh_end(m); ++i) 
        if (kh_exist(m, i)) total += kh_val(m, i);

    fprintf(stderr, "Finished Processing %zu reads...\n", nreads);

    double pct = nreads ? 100.0 * (double)total / (double)nreads : 0.0;
    fprintf(stderr, "Reads with matches: %zu (%.2f%%)\n", total, pct);
}

typedef struct {
    char *name;
    char *seq;
    char *qual;
    size_t len;
} FastqTask;

typedef struct {
    FastqTask *tasks;
    size_t n;
    size_t bases;
} FastqBatch;

typedef struct {
    FastqBatch *buf;
    size_t cap;
    size_t head;
    size_t tail;
    size_t size;
    int done;
    int aborted;
    pthread_mutex_t mtx;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} FastqTaskQueue;

typedef struct {
    char *data;
    size_t len;
} SamChunk;

typedef struct {
    SamChunk *buf;
    size_t cap;
    size_t head;
    size_t tail;
    size_t size;
    int done;
    int aborted;
    pthread_mutex_t mtx;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} SamChunkQueue;

typedef struct {
    SamChunkQueue *queue;
    FILE *sam_fp;
    int status;
} SamWriterCtx;

typedef struct {
    FastqTaskQueue *queue;
    SamChunkQueue *sam_queue;
    const KmerBitset *seed_index;
    uint32_t min_len;
    uint32_t max_len;
    int seed_mm;
    int k_mm;
    int exclude_multihit;
    int sam_soft_clip;
    int sam_emit_unmapped;
    FILE *sam_fp;
    kh_counter_t *local_counts;
    int status;
} FastqWorkerCtx;

static char *dup_cstr_n(const char *s, size_t n) {
    char *p = (char*)malloc(n + 1);
    if (!p) return NULL;
    memcpy(p, s, n);
    p[n] = '\0';
    return p;
}

static void free_fastq_task(FastqTask *task) {
    if (!task) return;
    free(task->name);
    free(task->seq);
    free(task->qual);
    task->name = NULL;
    task->seq = NULL;
    task->qual = NULL;
    task->len = 0;
}

static void free_fastq_batch(FastqBatch *batch) {
    if (!batch || !batch->tasks) return;
    for (size_t i = 0; i < batch->n; ++i) {
        free_fastq_task(&batch->tasks[i]);
    }
    free(batch->tasks);
    batch->tasks = NULL;
    batch->n = 0;
    batch->bases = 0;
}

static void free_sam_chunk(SamChunk *chunk) {
    if (!chunk) return;
    free(chunk->data);
    chunk->data = NULL;
    chunk->len = 0;
}

static void strset_free_with_keys(khash_t(strset) *set) {
    if (!set) return;
    for (khint_t i = kh_begin(set); i != kh_end(set); ++i) {
        if (kh_exist(set, i)) free((char*)kh_key(set, i));
    }
    kh_destroy(strset, set);
}

static int merge_counter_into(kh_counter_t *dst, const kh_counter_t *src) {
    for (khint_t i = kh_begin(src); i != kh_end(src); ++i) {
        if (!kh_exist(src, i)) continue;
        size_t val = kh_val(src, i);
        if (counter_add_with_init(dst, kh_key(src, i), val, val) != 0) return -1;
    }
    return 0;
}

static int add_matches_to_counter_hits(khash_t(strset) *found_sequences,
                                       kh_counter_t *map) {
    for (khint_t i = kh_begin(found_sequences); i != kh_end(found_sequences); ++i) {
        if (!kh_exist(found_sequences, i)) continue;
        if (counter_add_with_init(map, kh_key(found_sequences, i), 1, 1) != 0) return -1;
    }
    return 0;
}

static int sam_chunk_queue_push(SamChunkQueue *q, SamChunk *chunk);

static int process_one_read(const char *read_name,
                            const char *read_seq,
                            const char *read_qual,
                            size_t read_len,
                            kh_counter_t *counts,
                            const KmerBitset *seed_index,
                            int seed_mm,
                            size_t min_len,
                            size_t max_len,
                            int k_mm,
                            int exclude_multihit,
                            int sam_soft_clip,
                            int sam_emit_unmapped,
                            FILE *sam_fp,
                            SamChunkQueue *sam_queue,
                            int counts_preseeded) {
    int status = 0;
    FILE *sam_out = sam_fp;
    FILE *sam_mem = NULL;
    char *sam_buf = NULL;
    size_t sam_len = 0;
    kvec_t(uint32_t) hits;
    kv_init(hits);
    khash_t(strset) *matches = NULL;

    if (sam_fp && sam_queue) {
        sam_mem = open_memstream(&sam_buf, &sam_len);
        if (!sam_mem) {
            kv_destroy(hits);
            return 1;
        }
        sam_out = sam_mem;
    }

    find_kmer_bitset(read_seq, read_len, seed_index, seed_mm, &hits);

    if (hits.n == 0) {
        if (sam_out) {
            if (sam_emit_unmapped) {
                if (sam_queue) {
                    sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len);
                } else {
                    flockfile(sam_out);
                    sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len);
                    funlockfile(sam_out);
                }
            }
        }
        goto cleanup;
    }

    matches = kh_init(strset);
    if (!matches) {
        status = 1;
        goto cleanup;
    }

    find_matches_seeded(read_seq,
                        read_len,
                        &hits,
                        (uint32_t)min_len,
                        (uint32_t)max_len,
                        seed_index,
                        seed_mm,
                        matches,
                        k_mm,
                        read_name,
                        read_qual,
                        sam_out,
                        sam_soft_clip);

    if (sam_out && kh_size(matches) == 0) {
        if (sam_emit_unmapped) {
            if (sam_queue) {
                sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len);
            } else {
                flockfile(sam_out);
                sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len);
                funlockfile(sam_out);
            }
        }
    }

    if (!(exclude_multihit && kh_size(matches) > 1)) {
        if (counts_preseeded) {
            add_to_counter(matches, counts);
        } else if (add_matches_to_counter_hits(matches, counts) != 0) {
            status = 1;
            goto cleanup;
        }
    }

cleanup:
    if (matches) strset_free_with_keys(matches);
    kv_destroy(hits);

    if (sam_mem) {
        if (fflush(sam_mem) != 0) status = 1;
        fclose(sam_mem);
        sam_mem = NULL;

        if (status != 0) {
            free(sam_buf);
            return 1;
        }

        if (sam_len > 0) {
            SamChunk chunk = { sam_buf, sam_len };
            if (sam_chunk_queue_push(sam_queue, &chunk) != 0) {
                free(sam_buf);
                return 1;
            }
        } else {
            free(sam_buf);
        }
    }
    return status;
}

// Shared ring-queue implementation used by both SAM chunks and FASTQ batches.
#define DEFINE_RING_QUEUE_FUNCS(PREFIX, QueueType, ItemType, FreeItemFn)                  \
static int PREFIX##_init(QueueType *q, size_t cap) {                                       \
    if (!q || cap == 0) return -1;                                                         \
    memset(q, 0, sizeof(*q));                                                               \
    q->buf = (ItemType*)calloc(cap, sizeof(*q->buf));                                      \
    if (!q->buf) return -1;                                                                 \
    q->cap = cap;                                                                           \
    if (pthread_mutex_init(&q->mtx, NULL) != 0) {                                           \
        free(q->buf);                                                                       \
        q->buf = NULL;                                                                      \
        return -1;                                                                          \
    }                                                                                       \
    if (pthread_cond_init(&q->not_empty, NULL) != 0) {                                     \
        pthread_mutex_destroy(&q->mtx);                                                     \
        free(q->buf);                                                                       \
        q->buf = NULL;                                                                      \
        return -1;                                                                          \
    }                                                                                       \
    if (pthread_cond_init(&q->not_full, NULL) != 0) {                                      \
        pthread_cond_destroy(&q->not_empty);                                                \
        pthread_mutex_destroy(&q->mtx);                                                     \
        free(q->buf);                                                                       \
        q->buf = NULL;                                                                      \
        return -1;                                                                          \
    }                                                                                       \
    return 0;                                                                               \
}                                                                                           \
                                                                                            \
static void PREFIX##_request_abort(QueueType *q) {                                          \
    pthread_mutex_lock(&q->mtx);                                                            \
    q->aborted = 1;                                                                         \
    pthread_cond_broadcast(&q->not_empty);                                                  \
    pthread_cond_broadcast(&q->not_full);                                                   \
    pthread_mutex_unlock(&q->mtx);                                                          \
}                                                                                           \
                                                                                            \
static void PREFIX##_finish(QueueType *q) {                                                 \
    pthread_mutex_lock(&q->mtx);                                                            \
    q->done = 1;                                                                            \
    pthread_cond_broadcast(&q->not_empty);                                                  \
    pthread_mutex_unlock(&q->mtx);                                                          \
}                                                                                           \
                                                                                            \
/* Returns 0 on success, -1 on abort. */                                                    \
static int PREFIX##_push(QueueType *q, ItemType *item) {                                   \
    pthread_mutex_lock(&q->mtx);                                                            \
    while (q->size == q->cap && !q->aborted) {                                              \
        pthread_cond_wait(&q->not_full, &q->mtx);                                           \
    }                                                                                       \
    if (q->aborted) {                                                                       \
        pthread_mutex_unlock(&q->mtx);                                                      \
        return -1;                                                                          \
    }                                                                                       \
                                                                                            \
    q->buf[q->tail] = *item;                                                                \
    q->tail = (q->tail + 1) % q->cap;                                                       \
    q->size++;                                                                              \
    pthread_cond_signal(&q->not_empty);                                                     \
    pthread_mutex_unlock(&q->mtx);                                                          \
    return 0;                                                                               \
}                                                                                           \
                                                                                            \
/* Returns 0 on success, 1 on completion, -1 on abort. */                                   \
static int PREFIX##_pop(QueueType *q, ItemType *item_out) {                                \
    pthread_mutex_lock(&q->mtx);                                                            \
    while (q->size == 0 && !q->done && !q->aborted) {                                       \
        pthread_cond_wait(&q->not_empty, &q->mtx);                                          \
    }                                                                                       \
                                                                                            \
    if (q->aborted) {                                                                       \
        pthread_mutex_unlock(&q->mtx);                                                      \
        return -1;                                                                          \
    }                                                                                       \
                                                                                            \
    if (q->size == 0 && q->done) {                                                          \
        pthread_mutex_unlock(&q->mtx);                                                      \
        return 1;                                                                           \
    }                                                                                       \
                                                                                            \
    *item_out = q->buf[q->head];                                                            \
    memset(&q->buf[q->head], 0, sizeof(q->buf[q->head]));                                  \
    q->head = (q->head + 1) % q->cap;                                                       \
    q->size--;                                                                              \
    pthread_cond_signal(&q->not_full);                                                      \
    pthread_mutex_unlock(&q->mtx);                                                          \
    return 0;                                                                               \
}                                                                                           \
                                                                                            \
static void PREFIX##_destroy(QueueType *q) {                                                \
    if (!q) return;                                                                         \
    for (size_t i = 0; i < q->size; ++i) {                                                  \
        size_t idx = (q->head + i) % q->cap;                                                \
        FreeItemFn(&q->buf[idx]);                                                           \
    }                                                                                       \
    free(q->buf);                                                                           \
    pthread_cond_destroy(&q->not_empty);                                                    \
    pthread_cond_destroy(&q->not_full);                                                     \
    pthread_mutex_destroy(&q->mtx);                                                         \
}

DEFINE_RING_QUEUE_FUNCS(sam_chunk_queue, SamChunkQueue, SamChunk, free_sam_chunk)
DEFINE_RING_QUEUE_FUNCS(fastq_queue, FastqTaskQueue, FastqBatch, free_fastq_batch)

static void *sam_writer_main(void *arg) {
    SamWriterCtx *ctx = (SamWriterCtx*)arg;
    for (;;) {
        SamChunk chunk = {0};
        int pop_status = sam_chunk_queue_pop(ctx->queue, &chunk);
        if (pop_status == 1) break;
        if (pop_status < 0) break;

        if (chunk.data && chunk.len > 0) {
            size_t wrote = fwrite(chunk.data, 1, chunk.len, ctx->sam_fp);
            if (wrote != chunk.len) {
                ctx->status = 1;
                free_sam_chunk(&chunk);
                sam_chunk_queue_request_abort(ctx->queue);
                break;
            }
        }
        free_sam_chunk(&chunk);
    }
    return NULL;
}


static void *fastq_worker_main(void *arg) {
    FastqWorkerCtx *ctx = (FastqWorkerCtx*)arg;
    ctx->local_counts = kh_init(counter);
    if (!ctx->local_counts) {
        ctx->status = 1;
        fastq_queue_request_abort(ctx->queue);
        if (ctx->sam_queue) sam_chunk_queue_request_abort(ctx->sam_queue);
        return NULL;
    }

    for (;;) {
        FastqBatch batch = {0};
        int pop_status = fastq_queue_pop(ctx->queue, &batch);
        if (pop_status == 1) break;
        if (pop_status < 0) break;

        /*
         * In MT SAM mode, aggregate a whole FASTQ batch into one memory stream
         * before enqueueing to reduce per-read stream/queue overhead.
         */
        FILE *sam_batch_fp = NULL;
        char *sam_batch_buf = NULL;
        size_t sam_batch_len = 0;
        FILE *sam_target_fp = ctx->sam_fp;
        SamChunkQueue *sam_target_queue = ctx->sam_queue;
        if (ctx->sam_queue && ctx->sam_fp) {
            sam_batch_fp = open_memstream(&sam_batch_buf, &sam_batch_len);
            if (!sam_batch_fp) {
                ctx->status = 1;
                fastq_queue_request_abort(ctx->queue);
                sam_chunk_queue_request_abort(ctx->sam_queue);
                free_fastq_batch(&batch);
                break;
            }
            sam_target_fp = sam_batch_fp;
            sam_target_queue = NULL;
        }

        for (size_t i = 0; i < batch.n; ++i) {
            FastqTask *task = &batch.tasks[i];
            if (process_one_read(task->name,
                                 task->seq,
                                 task->qual,
                                 task->len,
                                 ctx->local_counts,
                                 ctx->seed_index,
                                 ctx->seed_mm,
                                 ctx->min_len,
                                 ctx->max_len,
                                 ctx->k_mm,
                                 ctx->exclude_multihit,
                                 ctx->sam_soft_clip,
                                 ctx->sam_emit_unmapped,
                                 sam_target_fp,
                                 sam_target_queue,
                                 0) != 0) {
                ctx->status = 1;
                fastq_queue_request_abort(ctx->queue);
                if (ctx->sam_queue) sam_chunk_queue_request_abort(ctx->sam_queue);
                break;
            }
        }

        if (sam_batch_fp) {
            if (ctx->status == 0 && fflush(sam_batch_fp) != 0) {
                ctx->status = 1;
                fastq_queue_request_abort(ctx->queue);
                sam_chunk_queue_request_abort(ctx->sam_queue);
            }
            if (fclose(sam_batch_fp) != 0 && ctx->status == 0) {
                ctx->status = 1;
                fastq_queue_request_abort(ctx->queue);
                sam_chunk_queue_request_abort(ctx->sam_queue);
            }
            sam_batch_fp = NULL;

            if (sam_batch_buf) {
                if (ctx->status == 0 && sam_batch_len > 0) {
                    SamChunk chunk = { sam_batch_buf, sam_batch_len };
                    if (sam_chunk_queue_push(ctx->sam_queue, &chunk) != 0) {
                        free(sam_batch_buf);
                        ctx->status = 1;
                        fastq_queue_request_abort(ctx->queue);
                        sam_chunk_queue_request_abort(ctx->sam_queue);
                    }
                } else {
                    free(sam_batch_buf);
                }
                sam_batch_buf = NULL;
                sam_batch_len = 0;
            }
        }

        free_fastq_batch(&batch);
        if (ctx->status != 0) {
            break;
        }
    }
    return NULL;
}

static int load_fastq_single(const char *path,
                             TrieNode *root,
                             kh_counter_t *counts,
                             int kmerlen,
                             int seed_mm,
                             size_t *min_len,
                             size_t *max_len,
                             int k_mm,
                             int exclude_multihit,
                             int sam_soft_clip,
                             int sam_emit_unmapped,
                             FILE *sam_fp) {
    gzFile fp = gzopen(path, "rb");
    if (fp == 0) { perror("gzopen"); return 1; }

    kseq_t *ks = kseq_init(fp);
    if (!ks) { gzclose(fp); return 1; }

    size_t nreads = 0;
    int l = 0;
    int status = 0;
    KmerBitset *seed_index = kmer_bitset_from_trie(root, (size_t)kmerlen);
    if (!seed_index) {
        fprintf(stderr, "Failed to build seed bitset index.\n");
        kseq_destroy(ks);
        gzclose(fp);
        return 1;
    }

    while ((l = kseq_read(ks)) >= 0) {
        if (++nreads % STEP == 0) {
            fprintf(stderr, "processing %zu reads...\r", nreads);
            fflush(stderr);
        }

        if (process_one_read(ks->name.s,
                             ks->seq.s,
                             ks->qual.l ? ks->qual.s : NULL,
                             ks->seq.l,
                             counts,
                             seed_index,
                             seed_mm,
                             *min_len,
                             *max_len,
                             k_mm,
                             exclude_multihit,
                             sam_soft_clip,
                             sam_emit_unmapped,
                             sam_fp,
                             NULL,
                             1) != 0) {
            status = 1;
            break;
        }
    }

    if (l < -1) status = 1;

    if (status == 0) total_hits(counts, nreads);

    kseq_destroy(ks);
    gzclose(fp);
    kmer_bitset_destroy(seed_index);
    return status;
}

static int load_fastq_mt(const char *path,
                         TrieNode *root,
                         kh_counter_t *counts,
                         int kmerlen,
                         int seed_mm,
                         size_t *min_len,
                         size_t *max_len,
                         int k_mm,
                         int exclude_multihit,
                         int sam_soft_clip,
                         int sam_emit_unmapped,
                         FILE *sam_fp,
                         unsigned threads) {
    SamChunkQueue sam_queue;
    SamWriterCtx writer_ctx;
    pthread_t writer_tid;
    int writer_started = 0;
    int use_sam_writer = (sam_fp != NULL);

    gzFile fp = gzopen(path, "rb");
    if (fp == 0) { perror("gzopen"); return 1; }

    kseq_t *ks = kseq_init(fp);
    if (!ks) { gzclose(fp); return 1; }

    FastqTaskQueue queue;
    if (fastq_queue_init(&queue, MT_TASK_QUEUE_CAPACITY) != 0) {
        kseq_destroy(ks);
        gzclose(fp);
        return 1;
    }
    KmerBitset *seed_index = kmer_bitset_from_trie(root, (size_t)kmerlen);
    if (!seed_index) {
        fprintf(stderr, "Failed to build seed bitset index.\n");
        fastq_queue_destroy(&queue);
        kseq_destroy(ks);
        gzclose(fp);
        return 1;
    }

    if (use_sam_writer) {
        if (sam_chunk_queue_init(&sam_queue, MT_SAM_QUEUE_CAPACITY) != 0) {
            kmer_bitset_destroy(seed_index);
            fastq_queue_destroy(&queue);
            kseq_destroy(ks);
            gzclose(fp);
            return 1;
        }
        memset(&writer_ctx, 0, sizeof(writer_ctx));
        writer_ctx.queue = &sam_queue;
        writer_ctx.sam_fp = sam_fp;
        if (pthread_create(&writer_tid, NULL, sam_writer_main, &writer_ctx) != 0) {
            sam_chunk_queue_destroy(&sam_queue);
            kmer_bitset_destroy(seed_index);
            fastq_queue_destroy(&queue);
            kseq_destroy(ks);
            gzclose(fp);
            return 1;
        }
        writer_started = 1;
    } else {
        memset(&writer_ctx, 0, sizeof(writer_ctx));
    }

    pthread_t *tids = (pthread_t*)calloc(threads, sizeof(*tids));
    FastqWorkerCtx *ctxs = (FastqWorkerCtx*)calloc(threads, sizeof(*ctxs));
    if (!tids || !ctxs) {
        free(tids);
        free(ctxs);
        if (writer_started) {
            sam_chunk_queue_request_abort(&sam_queue);
            pthread_join(writer_tid, NULL);
            sam_chunk_queue_destroy(&sam_queue);
        }
        kmer_bitset_destroy(seed_index);
        fastq_queue_destroy(&queue);
        kseq_destroy(ks);
        gzclose(fp);
        return 1;
    }

    size_t started = 0;
    int status = 0;
    for (size_t i = 0; i < threads; ++i) {
        ctxs[i].queue = &queue;
        ctxs[i].sam_queue = use_sam_writer ? &sam_queue : NULL;
        ctxs[i].seed_index = seed_index;
        ctxs[i].min_len = (uint32_t)(*min_len);
        ctxs[i].max_len = (uint32_t)(*max_len);
        ctxs[i].seed_mm = seed_mm;
        ctxs[i].k_mm = k_mm;
        ctxs[i].exclude_multihit = exclude_multihit;
        ctxs[i].sam_soft_clip = sam_soft_clip;
        ctxs[i].sam_emit_unmapped = sam_emit_unmapped;
        ctxs[i].sam_fp = sam_fp;
        ctxs[i].local_counts = NULL;
        ctxs[i].status = 0;
        if (pthread_create(&tids[i], NULL, fastq_worker_main, &ctxs[i]) != 0) {
            status = 1;
            fastq_queue_request_abort(&queue);
            if (use_sam_writer) sam_chunk_queue_request_abort(&sam_queue);
            break;
        }
        started++;
    }

    size_t nreads = 0;
    int l = 0;
    FastqBatch pending = {0};
    pending.tasks = (FastqTask*)calloc(MT_READ_BATCH_MAX_READS, sizeof(*pending.tasks));
    if (!pending.tasks) {
        status = 1;
        fastq_queue_request_abort(&queue);
        if (use_sam_writer) sam_chunk_queue_request_abort(&sam_queue);
    }
    if (status == 0) {
        while ((l = kseq_read(ks)) >= 0) {
            FastqTask task = {0};
            task.seq = dup_cstr_n(ks->seq.s, ks->seq.l);
            task.len = ks->seq.l;
            if (use_sam_writer) {
                task.name = dup_cstr_n(ks->name.s, ks->name.l);
                task.qual = ks->qual.l ? dup_cstr_n(ks->qual.s, ks->qual.l) : NULL;
            }

            if (!task.seq || (use_sam_writer && (!task.name || (ks->qual.l && !task.qual)))) {
                free_fastq_task(&task);
                status = 1;
                fastq_queue_request_abort(&queue);
                if (use_sam_writer) sam_chunk_queue_request_abort(&sam_queue);
                break;
            }

            pending.tasks[pending.n++] = task;
            pending.bases += task.len;
            if (pending.bases >= MT_READ_BATCH_TARGET_BASES || pending.n == MT_READ_BATCH_MAX_READS) {
                FastqBatch out = pending;
                if (fastq_queue_push(&queue, &out) != 0) {
                    free_fastq_batch(&out);
                    status = 1;
                    fastq_queue_request_abort(&queue);
                    if (use_sam_writer) sam_chunk_queue_request_abort(&sam_queue);
                    pending.tasks = NULL;
                    pending.n = 0;
                    pending.bases = 0;
                    break;
                }
                pending.tasks = (FastqTask*)calloc(MT_READ_BATCH_MAX_READS, sizeof(*pending.tasks));
                pending.n = 0;
                pending.bases = 0;
                if (!pending.tasks) {
                    status = 1;
                    fastq_queue_request_abort(&queue);
                    if (use_sam_writer) sam_chunk_queue_request_abort(&sam_queue);
                    break;
                }
            }

            if (++nreads % STEP == 0) {
                fprintf(stderr, "processing %zu reads...\r", nreads);
                fflush(stderr);
            }
        }
        if (l < -1) status = 1;
    }

    if (pending.tasks) {
        if (status == 0 && pending.n > 0) {
            FastqBatch out = pending;
            if (fastq_queue_push(&queue, &out) != 0) {
                free_fastq_batch(&out);
                status = 1;
                fastq_queue_request_abort(&queue);
                if (use_sam_writer) sam_chunk_queue_request_abort(&sam_queue);
            }
            pending.tasks = NULL;
            pending.n = 0;
            pending.bases = 0;
        }
        free_fastq_batch(&pending);
    }

    fastq_queue_finish(&queue);
    for (size_t i = 0; i < started; ++i) {
        pthread_join(tids[i], NULL);
    }

    if (writer_started) {
        if (status == 0) {
            sam_chunk_queue_finish(&sam_queue);
        } else {
            sam_chunk_queue_request_abort(&sam_queue);
        }
        pthread_join(writer_tid, NULL);
        if (writer_ctx.status != 0) status = 1;
    }

    if (status == 0) {
        for (size_t i = 0; i < started; ++i) {
            if (ctxs[i].status != 0) {
                status = 1;
                break;
            }
        }
    }

    if (status == 0) {
        for (size_t i = 0; i < started; ++i) {
            if (ctxs[i].local_counts && merge_counter_into(counts, ctxs[i].local_counts) != 0) {
                status = 1;
                break;
            }
        }
    }

    if (status == 0) total_hits(counts, nreads);

    for (size_t i = 0; i < started; ++i) {
        if (ctxs[i].local_counts) counter_free(ctxs[i].local_counts);
    }

    free(tids);
    free(ctxs);
    if (writer_started) sam_chunk_queue_destroy(&sam_queue);
    kmer_bitset_destroy(seed_index);
    fastq_queue_destroy(&queue);
    kseq_destroy(ks);
    gzclose(fp);
    return status;
}

int load_fastq(const char *path, TrieNode *root, kh_counter_t *counts,
               int kmerlen, int seed_mm, size_t *min_len, size_t *max_len, int k_mm,
               int exclude_multihit, FILE *sam_fp, int sam_soft_clip,
               int sam_emit_unmapped,
               unsigned threads) {
    if (threads <= 1) {
        return load_fastq_single(path, root, counts, kmerlen, seed_mm, min_len, max_len, k_mm,
                                 exclude_multihit, sam_soft_clip, sam_emit_unmapped, sam_fp);
    }
    return load_fastq_mt(path, root, counts, kmerlen, seed_mm, min_len, max_len, k_mm,
                         exclude_multihit, sam_soft_clip, sam_emit_unmapped, sam_fp, threads);
}



/* Loads reference sequences into trie.
   Returns 0 on success, 1 on I/O/memory error, 2 on duplicate when dup-policy=error.
   Works for FASTA or FASTQ, plain or .gz. */
int load_reference(const char *path, TrieNode *root, kh_counter_t *map, 
                   size_t *min_out, size_t *max_out, size_t *n_out, 
                   int add_revcomp, size_t *kmer_len, TrieDupPolicy dup_policy) 
{
    gzFile fp = gzopen(path, "rb");
    if (fp == 0) { perror("gzopen"); return 1; }

    kseq_t *ks = kseq_init(fp);
    if (!ks) { gzclose(fp); return 1; }

    size_t inserted = 0, min_len = (size_t)-1, max_len = 0;
    size_t k = kmer_len ? *kmer_len : 0;
    int l;
    while ((l = kseq_read(ks)) >= 0) {
        /* ks->name.s    : sequence name (kstring)
           ks->comment.s : optional comment (may be empty)
           ks->seq.s     : sequence string
           ks->qual.s    : FASTQ qualities (empty for FASTA)
           ks->seq.l     : sequence length
        */

        if (k > 0 && ks->seq.l < k) continue;

        char *seq = malloc(ks->seq.l + 1);
        if (!seq) {
            perror("malloc");
            kseq_destroy(ks);
            gzclose(fp);
            return 1;
        }
        memcpy(seq, ks->seq.s, ks->seq.l);
        seq[ks->seq.l] = '\0';

        normalize_acgt(seq);

        size_t L = ks->seq.l; // length of this record's sequence
        if (L < min_len) min_len = L;
        if (L > max_len) max_len = L;

        /* Use name up to first whitespace (common FASTA behavior) */
        char name_buf[256];
        size_t i = 0;
        while (ks->name.s[i] && !isspace((unsigned char)ks->name.s[i]) && i < sizeof(name_buf)-1) {
            name_buf[i] = ks->name.s[i];
            ++i;
        }
        name_buf[i] = '\0';

        const TrieNode *existing = trie_find_exact(root, seq);
        if (existing) {
            if (dup_policy != TRIE_DUP_IGNORE) {
                fprintf(stderr,
                        "Duplicate sequence detected: new seqID '%s' already exists as '%s'%s. "
                        "Possible reverse-complement or accidental duplicate.\n",
                        name_buf,
                        existing->name ? existing->name : "(unknown)",
                        existing->rev ? " [rev]" : "");
            }

            free(seq);
            if (dup_policy == TRIE_DUP_ERROR) {
                kseq_destroy(ks);
                gzclose(fp);
                return 2;
            }
            continue;
        }

        TrieInsertStatus st = trie_insert(root, seq, name_buf, 0);
        if (st == TRIE_INSERT_OOM) {
            free(seq);
            kseq_destroy(ks);
            gzclose(fp);
            return 1;
        }
        if (st == TRIE_INSERT_INVALID_BASE) {
            fprintf(stderr, "Skipping seqID '%s': invalid base found (expected A/C/G/T only).\n", name_buf);
            free(seq);
            continue;
        }
        if (st == TRIE_INSERT_DUP) {
            free(seq);
            if (dup_policy == TRIE_DUP_ERROR) {
                kseq_destroy(ks);
                gzclose(fp);
                return 2;
            }
            continue;
        }

        counter_inc(map, name_buf);
        inserted++;

        if (add_revcomp) {
            char *rc = revcomp_new_n(seq, ks->seq.l);
            if (rc) {
                char rcname[256];
                snprintf(rcname, sizeof(rcname), "%s/rc", name_buf);
                const TrieNode *existing_rc = trie_find_exact(root, rc);
                if (existing_rc) {
                    if (dup_policy != TRIE_DUP_IGNORE) {
                        fprintf(stderr,
                                "Duplicate sequence detected: new seqID '%s' already exists as '%s'%s. "
                                "Possible reverse-complement or accidental duplicate.\n",
                                rcname,
                                existing_rc->name ? existing_rc->name : "(unknown)",
                                existing_rc->rev ? " [rev]" : "");
                    }
                    if (dup_policy == TRIE_DUP_ERROR) {
                        free(rc);
                        free(seq);
                        kseq_destroy(ks);
                        gzclose(fp);
                        return 2;
                    }
                } else {
                    TrieInsertStatus rc_st = trie_insert(root, rc, rcname, 1);
                    if (rc_st == TRIE_INSERT_OOM) {
                        free(rc);
                        free(seq);
                        kseq_destroy(ks);
                        gzclose(fp);
                        return 1;
                    }
                    if (rc_st == TRIE_INSERT_INVALID_BASE) {
                        fprintf(stderr, "Skipping seqID '%s': invalid base found (expected A/C/G/T only).\n", rcname);
                    } else if (rc_st == TRIE_INSERT_DUP) {
                        if (dup_policy == TRIE_DUP_ERROR) {
                            free(rc);
                            free(seq);
                            kseq_destroy(ks);
                            gzclose(fp);
                            return 2;
                        }
                    } else {
                        counter_inc(map, rcname);
                        inserted++;
                    }
                }
                free(rc);
            }
        }
        free(seq);
    }

    if (min_out) *min_out = min_len;
    if (max_out) *max_out = max_len;
    if (n_out)   *n_out  = inserted;

    kseq_destroy(ks);
    gzclose(fp);
    return 0;
}
