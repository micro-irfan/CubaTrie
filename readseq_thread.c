#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "readseq_internal.h"

static const size_t STEP = 1000000;
static const size_t MT_TASK_QUEUE_CAPACITY = 1024;
static const size_t MT_READ_BATCH_TARGET_BASES = 16u * 1024u * 1024u;
static const size_t MT_READ_BATCH_MAX_READS = 4096;

static int anchor_runtime_init_for_count_mode(AnchorRuntime *out,
                                              const AnchorConfig *cfg,
                                              size_t ref_len) {
    if (!out) return -1;
    if (!cfg || !cfg->enabled) return anchor_runtime_init(out, cfg, ref_len);
    // In count mode, enforce insert length to match reference length so anchor
    // distance stays consistent with the mapped short-reference length.
    return anchor_runtime_init(out, cfg, ref_len);
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
    FastqTaskQueue *queue;
    const KmerBitset *seed_index;
    uint32_t min_len;
    uint32_t max_len;
    const AnchorRuntime *anchor_runtime;
    int seed_mm;
    int k_mm;
    int exclude_multihit;
    int sam_soft_clip;
    int sam_emit_unmapped;
    FILE *sam_fp;
    pthread_mutex_t *sam_write_mtx;
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

static int merge_counter_into(kh_counter_t *dst, const kh_counter_t *src) {
    for (khint_t i = kh_begin(src); i != kh_end(src); ++i) {
        if (!kh_exist(src, i)) continue;
        size_t val = kh_val(src, i);
        if (counter_add_with_init(dst, kh_key(src, i), val, val) != 0) return -1;
    }
    return 0;
}

static void total_hits(const kh_counter_t *m, const size_t nreads) {
    size_t total = 0;
    for (khint_t i = kh_begin(m); i != kh_end(m); ++i)
        if (kh_exist(m, i)) total += kh_val(m, i);

    fprintf(stderr, "Finished Processing %zu reads...\n", nreads);

    double pct = nreads ? 100.0 * (double)total / (double)nreads : 0.0;
    fprintf(stderr, "Total counted assignments: %zu (%.2f%% of reads)\n", total, pct);
}

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
static int PREFIX##_push(QueueType *q, ItemType *item) {                                   \
    pthread_mutex_lock(&q->mtx);                                                            \
    while (q->size == q->cap && !q->aborted) {                                              \
        pthread_cond_wait(&q->not_full, &q->mtx);                                           \
    }                                                                                       \
    if (q->aborted) {                                                                       \
        pthread_mutex_unlock(&q->mtx);                                                      \
        return -1;                                                                          \
    }                                                                                       \
    q->buf[q->tail] = *item;                                                                \
    q->tail = (q->tail + 1) % q->cap;                                                       \
    q->size++;                                                                              \
    pthread_cond_signal(&q->not_empty);                                                     \
    pthread_mutex_unlock(&q->mtx);                                                          \
    return 0;                                                                               \
}                                                                                           \
                                                                                            \
static int PREFIX##_pop(QueueType *q, ItemType *item_out) {                                \
    pthread_mutex_lock(&q->mtx);                                                            \
    while (q->size == 0 && !q->done && !q->aborted) {                                       \
        pthread_cond_wait(&q->not_empty, &q->mtx);                                          \
    }                                                                                       \
    if (q->aborted) {                                                                       \
        pthread_mutex_unlock(&q->mtx);                                                      \
        return -1;                                                                          \
    }                                                                                       \
    if (q->size == 0 && q->done) {                                                          \
        pthread_mutex_unlock(&q->mtx);                                                      \
        return 1;                                                                           \
    }                                                                                       \
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

DEFINE_RING_QUEUE_FUNCS(fastq_queue, FastqTaskQueue, FastqBatch, free_fastq_batch)

static void *fastq_worker_main(void *arg) {
    FastqWorkerCtx *ctx = (FastqWorkerCtx*)arg;
    ctx->local_counts = kh_init(counter);
    if (!ctx->local_counts) {
        ctx->status = 1;
        fastq_queue_request_abort(ctx->queue);
        return NULL;
    }

    for (;;) {
        FastqBatch batch = {0};
        int pop_status = fastq_queue_pop(ctx->queue, &batch);
        if (pop_status == 1) break;
        if (pop_status < 0) break;

        FILE *sam_batch_fp = NULL;
        char *sam_batch_buf = NULL;
        size_t sam_batch_len = 0;
        FILE *sam_target_fp = ctx->sam_fp;
        if (ctx->sam_fp && ctx->sam_write_mtx) {
            sam_batch_fp = open_memstream(&sam_batch_buf, &sam_batch_len);
            if (!sam_batch_fp) {
                ctx->status = 1;
                fastq_queue_request_abort(ctx->queue);
                free_fastq_batch(&batch);
                break;
            }
            sam_target_fp = sam_batch_fp;
        }

        for (size_t i = 0; i < batch.n; ++i) {
            FastqTask *task = &batch.tasks[i];
            if (readseq_process_one_read(task->name,
                                         task->seq,
                                         task->qual,
                                         task->len,
                                         ctx->local_counts,
                                         ctx->seed_index,
                                         ctx->anchor_runtime,
                                         ctx->seed_mm,
                                         ctx->min_len,
                                         ctx->max_len,
                                         ctx->k_mm,
                                         ctx->exclude_multihit,
                                         ctx->sam_soft_clip,
                                         ctx->sam_emit_unmapped,
                                         sam_target_fp,
                                         0) != 0) {
                ctx->status = 1;
                fastq_queue_request_abort(ctx->queue);
                break;
            }
        }

        if (sam_batch_fp) {
            if (ctx->status == 0 && fflush(sam_batch_fp) != 0) {
                ctx->status = 1;
                fastq_queue_request_abort(ctx->queue);
            }
            if (fclose(sam_batch_fp) != 0 && ctx->status == 0) {
                ctx->status = 1;
                fastq_queue_request_abort(ctx->queue);
            }
            sam_batch_fp = NULL;
            if (sam_batch_buf) {
                if (ctx->status == 0 && sam_batch_len > 0) {
                    pthread_mutex_lock(ctx->sam_write_mtx);
                    size_t wrote = fwrite(sam_batch_buf, 1, sam_batch_len, ctx->sam_fp);
                    pthread_mutex_unlock(ctx->sam_write_mtx);
                    if (wrote != sam_batch_len) {
                        ctx->status = 1;
                        fastq_queue_request_abort(ctx->queue);
                    }
                }
                free(sam_batch_buf);
                sam_batch_buf = NULL;
                sam_batch_len = 0;
            }
        }

        free_fastq_batch(&batch);
        if (ctx->status != 0) break;
    }
    return NULL;
}

int readseq_load_fastq_mt(const char *path,
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
                          unsigned threads,
                          const AnchorConfig *anchor_cfg) {
    int use_sam_output = (sam_fp != NULL);
    pthread_mutex_t sam_write_mtx;
    int sam_write_mtx_init = 0;
    AnchorRuntime anchor_runtime;

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
    if (anchor_runtime_init_for_count_mode(&anchor_runtime, anchor_cfg, *min_len) != 0) {
        fprintf(stderr, "Failed to initialize anchor matcher. Check adapter sequences.\n");
        kmer_bitset_destroy(seed_index);
        fastq_queue_destroy(&queue);
        kseq_destroy(ks);
        gzclose(fp);
        return 1;
    }

    if (use_sam_output && pthread_mutex_init(&sam_write_mtx, NULL) != 0) {
        kmer_bitset_destroy(seed_index);
        fastq_queue_destroy(&queue);
        kseq_destroy(ks);
        gzclose(fp);
        return 1;
    }
    sam_write_mtx_init = use_sam_output ? 1 : 0;

    pthread_t *tids = (pthread_t*)calloc(threads, sizeof(*tids));
    FastqWorkerCtx *ctxs = (FastqWorkerCtx*)calloc(threads, sizeof(*ctxs));
    if (!tids || !ctxs) {
        free(tids);
        free(ctxs);
        if (sam_write_mtx_init) pthread_mutex_destroy(&sam_write_mtx);
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
        ctxs[i].seed_index = seed_index;
        ctxs[i].min_len = (uint32_t)(*min_len);
        ctxs[i].max_len = (uint32_t)(*max_len);
        ctxs[i].anchor_runtime = &anchor_runtime;
        ctxs[i].seed_mm = seed_mm;
        ctxs[i].k_mm = k_mm;
        ctxs[i].exclude_multihit = exclude_multihit;
        ctxs[i].sam_soft_clip = sam_soft_clip;
        ctxs[i].sam_emit_unmapped = sam_emit_unmapped;
        ctxs[i].sam_fp = use_sam_output ? sam_fp : NULL;
        ctxs[i].sam_write_mtx = use_sam_output ? &sam_write_mtx : NULL;
        ctxs[i].local_counts = NULL;
        ctxs[i].status = 0;
        if (pthread_create(&tids[i], NULL, fastq_worker_main, &ctxs[i]) != 0) {
            status = 1;
            fastq_queue_request_abort(&queue);
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
    }
    if (status == 0) {
        while ((l = kseq_read(ks)) >= 0) {
            FastqTask task = {0};
            task.seq = dup_cstr_n(ks->seq.s, ks->seq.l);
            task.len = ks->seq.l;
            if (use_sam_output) {
                task.name = dup_cstr_n(ks->name.s, ks->name.l);
                task.qual = ks->qual.l ? dup_cstr_n(ks->qual.s, ks->qual.l) : NULL;
            }

            if (!task.seq || (use_sam_output && (!task.name || (ks->qual.l && !task.qual)))) {
                free_fastq_task(&task);
                status = 1;
                fastq_queue_request_abort(&queue);
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

    if (status == 0 && use_sam_output && fflush(sam_fp) != 0) {
        status = 1;
    }

    if (status == 0) total_hits(counts, nreads);

    for (size_t i = 0; i < started; ++i) {
        if (ctxs[i].local_counts) counter_free(ctxs[i].local_counts);
    }

    free(tids);
    free(ctxs);
    if (sam_write_mtx_init) pthread_mutex_destroy(&sam_write_mtx);
    kmer_bitset_destroy(seed_index);
    fastq_queue_destroy(&queue);
    kseq_destroy(ks);
    gzclose(fp);
    return status;
}
