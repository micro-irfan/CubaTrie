#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "readseq.h"
#include "readseq_anchor.h"

typedef struct {
    int use_gzip;
    FILE *fp;
    gzFile gzf;
} FastqOut;

typedef struct {
    size_t total;
    size_t kept;
    size_t discarded;
    size_t kept_rc;
} CutStats;

typedef struct {
    char *name;
    char *seq;
    char *qual;
    size_t len;
    size_t qual_len;
} CutTask;

typedef struct {
    CutTask *buf;
    size_t cap;
    size_t head;
    size_t tail;
    size_t size;
    int done;
    int aborted;
    pthread_mutex_t mtx;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} CutTaskQueue;

typedef struct {
    CutTaskQueue *queue;
    const AnchorRuntime *anchor_runtime;
    int check_revcomp;
    FastqOut *out;
    pthread_mutex_t *out_mtx;
    CutStats stats;
    int status;
    kh_counter_t *counts;
} CutWorkerCtx;

typedef struct {
    const char *sequence;
    size_t count;
} CutCountEntry;

static const size_t CUT_TASK_QUEUE_CAPACITY = 4096;

static int has_gz_suffix(const char *path) {
    if (!path) return 0;
    size_t n = strlen(path);
    return n >= 3 && path[n - 3] == '.' && path[n - 2] == 'g' && path[n - 1] == 'z';
}

static int fastq_out_open(FastqOut *out, const char *path) {
    if (!out || !path) return -1;
    memset(out, 0, sizeof(*out));
    if (strcmp(path, "-") == 0) {
        out->fp = stdout;
        out->use_gzip = 0;
        return 0;
    }
    if (has_gz_suffix(path)) {
        out->gzf = gzopen(path, "wb");
        if (!out->gzf) return -1;
        out->use_gzip = 1;
        return 0;
    }
    out->fp = fopen(path, "w");
    if (!out->fp) return -1;
    out->use_gzip = 0;
    return 0;
}

static int fastq_out_write(FastqOut *out,
                           const char *name,
                           const char *seq,
                           const char *qual,
                           size_t len) {
    if (!out || !name || !seq || !qual) return -1;
    if (out->use_gzip) {
        if (gzprintf(out->gzf, "@%s\n", name) <= 0) return -1;
        if ((size_t)gzwrite(out->gzf, seq, (unsigned)len) != len) return -1;
        if (gzprintf(out->gzf, "\n+\n") <= 0) return -1;
        if ((size_t)gzwrite(out->gzf, qual, (unsigned)len) != len) return -1;
        if (gzprintf(out->gzf, "\n") <= 0) return -1;
        return 0;
    }
    if (fprintf(out->fp, "@%s\n", name) < 0) return -1;
    if (fwrite(seq, 1, len, out->fp) != len) return -1;
    if (fputs("\n+\n", out->fp) == EOF) return -1;
    if (fwrite(qual, 1, len, out->fp) != len) return -1;
    if (fputc('\n', out->fp) == EOF) return -1;
    return 0;
}

static void fastq_out_close(FastqOut *out) {
    if (!out) return;
    if (out->use_gzip) {
        if (out->gzf) gzclose(out->gzf);
    } else if (out->fp && out->fp != stdout) {
        fclose(out->fp);
    }
}

static char rc_base(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 'T';
        case 'c': return 'G';
        case 'g': return 'C';
        case 't': return 'A';
        case 'u': return 'A';
        case 'U': return 'A';
        default:  return 'N';
    }
}

static char *revcomp_new_n(const char *s, size_t n) {
    char *rc = (char*)malloc(n + 1);
    if (!rc) return NULL;
    for (size_t i = 0; i < n; ++i) rc[i] = rc_base(s[n - 1 - i]);
    rc[n] = '\0';
    return rc;
}

static char *reverse_new_n(const char *s, size_t n) {
    char *r = (char*)malloc(n + 1);
    if (!r) return NULL;
    for (size_t i = 0; i < n; ++i) r[i] = s[n - 1 - i];
    r[n] = '\0';
    return r;
}

static char *dup_cstr_n(const char *s, size_t n) {
    if (!s) return NULL;
    char *p = (char*)malloc(n + 1);
    if (!p) return NULL;
    memcpy(p, s, n);
    p[n] = '\0';
    return p;
}

static int cut_counter_inc_slice(kh_counter_t *m, const char *seq, size_t len) {
    if (!m || !seq) return -1;
    char *owned = dup_cstr_n(seq, len);
    if (!owned) return -1;

    int ret = 0;
    khiter_t it = kh_put(counter, m, owned, &ret);
    if (ret < 0) {
        free(owned);
        return -1;
    }
    if (ret == 0) {
        kh_val(m, it) += 1;
        free(owned);
    } else {
        kh_key(m, it) = owned;
        kh_val(m, it) = 1;
    }
    return 0;
}

static int cut_counter_merge_into(kh_counter_t *dst, const kh_counter_t *src) {
    if (!dst || !src) return -1;
    for (khint_t i = kh_begin(src); i != kh_end(src); ++i) {
        if (!kh_exist(src, i)) continue;
        if (counter_add_with_init(dst, kh_key(src, i), kh_val(src, i), kh_val(src, i)) != 0) {
            return -1;
        }
    }
    return 0;
}

static int cmp_cut_count_desc_then_seq_asc(const void *pa, const void *pb) {
    const CutCountEntry *a = (const CutCountEntry*)pa;
    const CutCountEntry *b = (const CutCountEntry*)pb;
    if (a->count != b->count) return (a->count > b->count) ? -1 : 1;
    return strcmp(a->sequence, b->sequence);
}

static int cut_write_count_csv_sorted(const kh_counter_t *counts, const char *path) {
    if (!counts || !path) return -1;
    FILE *fp = (strcmp(path, "-") == 0) ? stdout : fopen(path, "w");
    if (!fp) {
        perror("fopen");
        return -1;
    }

    size_t n = 0;
    for (khint_t i = kh_begin(counts); i != kh_end(counts); ++i) {
        if (kh_exist(counts, i)) ++n;
    }

    if (fprintf(fp, "sequence,count\n") < 0) {
        if (fp != stdout) fclose(fp);
        return -1;
    }
    if (n == 0) {
        if (fp != stdout) fclose(fp);
        return 0;
    }

    CutCountEntry *arr = (CutCountEntry*)malloc(n * sizeof(*arr));
    if (!arr) {
        if (fp != stdout) fclose(fp);
        return -1;
    }

    size_t j = 0;
    for (khint_t i = kh_begin(counts); i != kh_end(counts); ++i) {
        if (!kh_exist(counts, i)) continue;
        arr[j].sequence = kh_key(counts, i);
        arr[j].count = kh_val(counts, i);
        ++j;
    }

    qsort(arr, n, sizeof(*arr), cmp_cut_count_desc_then_seq_asc);

    for (size_t i = 0; i < n; ++i) {
        if (fprintf(fp, "%s,%zu\n", arr[i].sequence, arr[i].count) < 0) {
            free(arr);
            if (fp != stdout) fclose(fp);
            return -1;
        }
    }

    free(arr);
    if (fp != stdout) fclose(fp);
    return 0;
}

static void cut_task_free(CutTask *task) {
    if (!task) return;
    free(task->name);
    free(task->seq);
    free(task->qual);
    task->name = NULL;
    task->seq = NULL;
    task->qual = NULL;
    task->len = 0;
    task->qual_len = 0;
}

static int cut_queue_init(CutTaskQueue *q, size_t cap) {
    if (!q || cap == 0) return -1;
    memset(q, 0, sizeof(*q));
    q->buf = (CutTask*)calloc(cap, sizeof(*q->buf));
    if (!q->buf) return -1;
    q->cap = cap;
    if (pthread_mutex_init(&q->mtx, NULL) != 0) {
        free(q->buf);
        q->buf = NULL;
        return -1;
    }
    if (pthread_cond_init(&q->not_empty, NULL) != 0) {
        pthread_mutex_destroy(&q->mtx);
        free(q->buf);
        q->buf = NULL;
        return -1;
    }
    if (pthread_cond_init(&q->not_full, NULL) != 0) {
        pthread_cond_destroy(&q->not_empty);
        pthread_mutex_destroy(&q->mtx);
        free(q->buf);
        q->buf = NULL;
        return -1;
    }
    return 0;
}

static void cut_queue_request_abort(CutTaskQueue *q) {
    if (!q) return;
    pthread_mutex_lock(&q->mtx);
    q->aborted = 1;
    pthread_cond_broadcast(&q->not_empty);
    pthread_cond_broadcast(&q->not_full);
    pthread_mutex_unlock(&q->mtx);
}

static void cut_queue_finish(CutTaskQueue *q) {
    if (!q) return;
    pthread_mutex_lock(&q->mtx);
    q->done = 1;
    pthread_cond_broadcast(&q->not_empty);
    pthread_mutex_unlock(&q->mtx);
}

static int cut_queue_push(CutTaskQueue *q, CutTask *task) {
    if (!q || !task) return -1;
    pthread_mutex_lock(&q->mtx);
    while (q->size == q->cap && !q->aborted) {
        pthread_cond_wait(&q->not_full, &q->mtx);
    }
    if (q->aborted) {
        pthread_mutex_unlock(&q->mtx);
        return -1;
    }

    q->buf[q->tail] = *task;
    memset(task, 0, sizeof(*task));
    q->tail = (q->tail + 1) % q->cap;
    q->size++;
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->mtx);
    return 0;
}

static int cut_queue_pop(CutTaskQueue *q, CutTask *task_out) {
    if (!q || !task_out) return -1;
    pthread_mutex_lock(&q->mtx);
    while (q->size == 0 && !q->done && !q->aborted) {
        pthread_cond_wait(&q->not_empty, &q->mtx);
    }

    if (q->aborted) {
        pthread_mutex_unlock(&q->mtx);
        return -1;
    }
    if (q->size == 0 && q->done) {
        pthread_mutex_unlock(&q->mtx);
        return 1;
    }

    *task_out = q->buf[q->head];
    memset(&q->buf[q->head], 0, sizeof(q->buf[q->head]));
    q->head = (q->head + 1) % q->cap;
    q->size--;
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->mtx);
    return 0;
}

static void cut_queue_destroy(CutTaskQueue *q) {
    if (!q) return;
    if (q->buf) {
        for (size_t i = 0; i < q->size; ++i) {
            size_t idx = (q->head + i) % q->cap;
            cut_task_free(&q->buf[idx]);
        }
    }
    free(q->buf);
    q->buf = NULL;
    pthread_cond_destroy(&q->not_empty);
    pthread_cond_destroy(&q->not_full);
    pthread_mutex_destroy(&q->mtx);
}

static int cut_process_record(const AnchorRuntime *ar,
                              int check_revcomp,
                              const char *name,
                              const char *seq,
                              const char *qual,
                              size_t len,
                              size_t qual_len,
                              FastqOut *out,
                              pthread_mutex_t *out_mtx,
                              CutStats *stats,
                              kh_counter_t *counts) {
    if (!ar || !name || !seq || !stats) return -1;

    stats->total++;
    if (!qual || qual_len != len) {
        stats->discarded++;
        return 0;
    }

    size_t start = 0, trim_len = 0;
    int ok = anchor_extract_window_range(ar, seq, len, &start, &trim_len);
    const char *src_seq = seq;
    const char *src_qual = qual;
    char *rc_seq = NULL;
    char *rc_qual = NULL;

    if (ok != 0 && check_revcomp) {
        rc_seq = revcomp_new_n(seq, len);
        rc_qual = reverse_new_n(qual, qual_len);
        if (!rc_seq || !rc_qual) {
            free(rc_seq);
            free(rc_qual);
            return -1;
        }
        ok = anchor_extract_window_range(ar, rc_seq, len, &start, &trim_len);
        if (ok == 0) {
            src_seq = rc_seq;
            src_qual = rc_qual;
            stats->kept_rc++;
        }
    }

    if (ok != 0 || start + trim_len > len) {
        stats->discarded++;
        free(rc_seq);
        free(rc_qual);
        return 0;
    }

    if (out) {
        int wr = 0;
        if (out_mtx) pthread_mutex_lock(out_mtx);
        wr = fastq_out_write(out, name, src_seq + start, src_qual + start, trim_len);
        if (out_mtx) pthread_mutex_unlock(out_mtx);
        if (wr != 0) {
            free(rc_seq);
            free(rc_qual);
            return -1;
        }
    }

    if (counts && cut_counter_inc_slice(counts, src_seq + start, trim_len) != 0) {
        free(rc_seq);
        free(rc_qual);
        return -1;
    }

    free(rc_seq);
    free(rc_qual);
    stats->kept++;
    return 0;
}

static void *cut_worker_main(void *arg) {
    CutWorkerCtx *ctx = (CutWorkerCtx*)arg;

    for (;;) {
        CutTask task = {0};
        int pop_rc = cut_queue_pop(ctx->queue, &task);
        if (pop_rc == 1) break;
        if (pop_rc < 0) break;

        if (cut_process_record(ctx->anchor_runtime,
                               ctx->check_revcomp,
                               task.name,
                               task.seq,
                               task.qual,
                               task.len,
                               task.qual_len,
                               ctx->out,
                               ctx->out_mtx,
                               &ctx->stats,
                               ctx->counts) != 0) {
            ctx->status = 1;
            cut_queue_request_abort(ctx->queue);
            cut_task_free(&task);
            break;
        }

        cut_task_free(&task);
    }

    return NULL;
}

static int cut_fastq_by_anchors_single(const char *in_path,
                                        const char *out_path,
                                        const AnchorRuntime *ar,
                                        int check_revcomp,
                                        const char *count_out_path) {
    gzFile in_fp = gzopen(in_path, "rb");
    if (!in_fp) {
        perror("gzopen");
        return 1;
    }
    kseq_t *ks = kseq_init(in_fp);
    if (!ks) {
        gzclose(in_fp);
        return 1;
    }

    int out_opened = 0;
    FastqOut out = {0};
    if (out_path) {
        if (fastq_out_open(&out, out_path) != 0) {
            perror("open output");
            kseq_destroy(ks);
            gzclose(in_fp);
            return 1;
        }
        out_opened = 1;
    }

    kh_counter_t *counts = NULL;
    if (count_out_path) {
        counts = kh_init(counter);
        if (!counts) {
            if (out_opened) fastq_out_close(&out);
            kseq_destroy(ks);
            gzclose(in_fp);
            return 1;
        }
    }

    CutStats stats = {0};
    int l = 0;
    int status = 0;
    while ((l = kseq_read(ks)) >= 0) {
        if (cut_process_record(ar,
                               check_revcomp,
                               ks->name.s,
                               ks->seq.s,
                               ks->qual.l ? ks->qual.s : NULL,
                               ks->seq.l,
                               ks->qual.l,
                               out_opened ? &out : NULL,
                               NULL,
                               &stats,
                               counts) != 0) {
            status = 1;
            break;
        }
    }

    if (l < -1) status = 1;

    if (status == 0 && counts && cut_write_count_csv_sorted(counts, count_out_path) != 0) {
        status = 1;
    }

    if (counts) counter_free(counts);
    if (out_opened) fastq_out_close(&out);
    kseq_destroy(ks);
    gzclose(in_fp);

    fprintf(stderr,
            "cut summary: total=%zu kept=%zu discarded=%zu kept_revcomp=%zu\n",
            stats.total, stats.kept, stats.discarded, stats.kept_rc);
    return status;
}

static int cut_fastq_by_anchors_mt(const char *in_path,
                                    const char *out_path,
                                    const AnchorRuntime *ar,
                                    int check_revcomp,
                                    unsigned threads,
                                    const char *count_out_path) {
    int status = 0;
    int queue_inited = 0;
    int out_mtx_inited = 0;
    pthread_t *tids = NULL;
    CutWorkerCtx *ctxs = NULL;
    kh_counter_t *all_counts = NULL;

    gzFile in_fp = gzopen(in_path, "rb");
    if (!in_fp) {
        perror("gzopen");
        return 1;
    }
    kseq_t *ks = kseq_init(in_fp);
    if (!ks) {
        gzclose(in_fp);
        return 1;
    }

    int out_opened = 0;
    FastqOut out = {0};
    if (out_path) {
        if (fastq_out_open(&out, out_path) != 0) {
            perror("open output");
            kseq_destroy(ks);
            gzclose(in_fp);
            return 1;
        }
        out_opened = 1;
    }

    if (count_out_path) {
        all_counts = kh_init(counter);
        if (!all_counts) {
            if (out_opened) fastq_out_close(&out);
            kseq_destroy(ks);
            gzclose(in_fp);
            return 1;
        }
    }

    CutTaskQueue queue;
    if (cut_queue_init(&queue, CUT_TASK_QUEUE_CAPACITY) != 0) {
        if (all_counts) counter_free(all_counts);
        if (out_opened) fastq_out_close(&out);
        kseq_destroy(ks);
        gzclose(in_fp);
        return 1;
    }
    queue_inited = 1;

    pthread_mutex_t out_mtx;
    if (out_opened) {
        if (pthread_mutex_init(&out_mtx, NULL) != 0) {
            status = 1;
            goto cleanup;
        }
        out_mtx_inited = 1;
    }

    tids = (pthread_t*)calloc(threads, sizeof(*tids));
    ctxs = (CutWorkerCtx*)calloc(threads, sizeof(*ctxs));
    if (!tids || !ctxs) {
        status = 1;
        goto cleanup;
    }

    size_t started = 0;
    for (size_t i = 0; i < threads; ++i) {
        ctxs[i].queue = &queue;
        ctxs[i].anchor_runtime = ar;
        ctxs[i].check_revcomp = check_revcomp;
        ctxs[i].out = out_opened ? &out : NULL;
        ctxs[i].out_mtx = out_opened ? &out_mtx : NULL;
        memset(&ctxs[i].stats, 0, sizeof(ctxs[i].stats));
        ctxs[i].status = 0;
        ctxs[i].counts = NULL;
        if (count_out_path) {
            ctxs[i].counts = kh_init(counter);
            if (!ctxs[i].counts) {
                status = 1;
                cut_queue_request_abort(&queue);
                break;
            }
        }
        if (pthread_create(&tids[i], NULL, cut_worker_main, &ctxs[i]) != 0) {
            status = 1;
            cut_queue_request_abort(&queue);
            break;
        }
        started++;
    }

    int l = 0;
    if (status == 0) {
        while ((l = kseq_read(ks)) >= 0) {
            CutTask task = {0};
            task.name = dup_cstr_n(ks->name.s, ks->name.l);
            task.seq = dup_cstr_n(ks->seq.s, ks->seq.l);
            task.len = ks->seq.l;
            task.qual_len = ks->qual.l;
            task.qual = ks->qual.l ? dup_cstr_n(ks->qual.s, ks->qual.l) : NULL;

            if (!task.name || !task.seq || (ks->qual.l && !task.qual)) {
                cut_task_free(&task);
                status = 1;
                cut_queue_request_abort(&queue);
                break;
            }

            if (cut_queue_push(&queue, &task) != 0) {
                cut_task_free(&task);
                status = 1;
                cut_queue_request_abort(&queue);
                break;
            }
        }
        if (l < -1) status = 1;
    }

    cut_queue_finish(&queue);
    for (size_t i = 0; i < started; ++i) {
        pthread_join(tids[i], NULL);
    }

    CutStats stats = {0};
    for (size_t i = 0; i < started; ++i) {
        stats.total += ctxs[i].stats.total;
        stats.kept += ctxs[i].stats.kept;
        stats.discarded += ctxs[i].stats.discarded;
        stats.kept_rc += ctxs[i].stats.kept_rc;
        if (ctxs[i].status != 0) status = 1;
        if (status == 0 && all_counts && ctxs[i].counts &&
            cut_counter_merge_into(all_counts, ctxs[i].counts) != 0) {
            status = 1;
        }
    }

    if (status == 0 && all_counts && cut_write_count_csv_sorted(all_counts, count_out_path) != 0) {
        status = 1;
    }

    fprintf(stderr,
            "cut summary: total=%zu kept=%zu discarded=%zu kept_revcomp=%zu\n",
            stats.total, stats.kept, stats.discarded, stats.kept_rc);

cleanup:
    if (ctxs) {
        for (size_t i = 0; i < threads; ++i) {
            if (ctxs[i].counts) counter_free(ctxs[i].counts);
        }
    }
    if (all_counts) counter_free(all_counts);
    free(tids);
    free(ctxs);
    if (out_mtx_inited) pthread_mutex_destroy(&out_mtx);
    if (queue_inited) cut_queue_destroy(&queue);
    if (out_opened) fastq_out_close(&out);
    kseq_destroy(ks);
    gzclose(in_fp);
    return status;
}

int cut_fastq_by_anchors(const char *in_path,
                         const char *out_path,
                         const AnchorConfig *anchor_cfg,
                         size_t min_len,
                         size_t max_len,
                         int check_revcomp,
                         unsigned threads,
                         const char *count_out_path) {
    if (!in_path || !anchor_cfg || !anchor_cfg->enabled) return 1;
    if (!out_path && !count_out_path) {
        fprintf(stderr, "cut mode requires at least one output: FASTQ (-o/--output) or counts CSV (-c/--count).\n");
        return 1;
    }
    if (min_len == 0 || max_len == 0 || min_len > max_len) return 1;

    AnchorRuntime ar;
    if (anchor_runtime_init_range(&ar, anchor_cfg, min_len, max_len) != 0) {
        fprintf(stderr, "Failed to initialize anchor matcher for cut mode.\n");
        return 1;
    }

    if (threads <= 1) {
        return cut_fastq_by_anchors_single(in_path, out_path, &ar, check_revcomp, count_out_path);
    }
    return cut_fastq_by_anchors_mt(in_path, out_path, &ar, check_revcomp, threads, count_out_path);
}
