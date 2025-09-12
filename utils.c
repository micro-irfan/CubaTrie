#include "utils.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

void mv_init(kFoundVec *mv) { mv->a=NULL; mv->n=mv->m=0; }
void mv_push(kFoundVec *mv, size_t pos, const char *name, const char *seq) {
    if (mv->n == mv->m) {
        size_t nm = mv->m ? mv->m<<1 : 8;
        Match *p = (Match*)realloc(mv->a, nm*sizeof *p);
        if (!p) { fprintf(stderr,"OOM\n"); exit(1); }
        mv->a = p; mv->m = nm;
    }
    mv->a[mv->n].pos = pos;
    mv->a[mv->n].name = strndup(name, strlen(name)); 
    mv->a[mv->n].seq = strndup(seq, strlen(seq));
    mv->n++;
}
void mv_free(kFoundVec *mv) {
    for (size_t i=0;i<mv->n;i++) {
        free(mv->a[i].name);
        free(mv->a[i].seq);
    }
    free(mv->a);
}

static char *dup_cstr_n(const char *s, size_t n) {
    char *p = (char*)malloc(n + 1);
    if (!p) return NULL;
    memcpy(p, s, n);
    p[n] = '\0';
    return p;
}


/* -------- COUNTER CONTAINER -------- */
static char *dup_cstr(const char *s) {
    size_t n = strlen(s);
    char *p = (char*)malloc(n + 1);
    if (!p) return NULL;
    memcpy(p, s, n + 1);
    return p;
}

void counter_inc(kh_counter_t *m, const char *key) {
    int ret;
    khiter_t k = kh_put(counter, m, key, &ret);
    if (ret < 0) { perror("kh_put"); exit(1); }
    if (ret == 0) {
        // existed
        kh_val(m, k)++;
    } else {
        // new slot — own a copy of the key
        char *owned = dup_cstr(key);
        if (!owned) { perror("malloc"); exit(1); }
        kh_key(m, k) = owned;
        kh_val(m, k) = 0;
    }
}



void counter_free(kh_counter_t *m) {
    if (!m) return;
    for (khint_t k = kh_begin(m); k != kh_end(m); ++k)
        if (kh_exist(m, k)) free((void*)kh_key(m, k));
    kh_destroy(counter, m);
}


/* helper: strdup(s without trailing "/rc") */
static char *dup_base_no_rc(const char *s) {
    size_t len = strlen(s);
    if (len >= 3 && s[len-3] == '/' && s[len-2] == 'r' && s[len-1] == 'c') {
        len -= 3;
    }
    char *out = (char*)malloc(len + 1);
    if (!out) return NULL;
    memcpy(out, s, len);
    out[len] = '\0';
    return out;
}

/* your comparator (if you don’t already have it) */
static int cmp_count_desc_then_key_asc(const void *pa, const void *pb) {
    const counter_entry_t *a = (const counter_entry_t*)pa;
    const counter_entry_t *b = (const counter_entry_t*)pb;
    if (a->count != b->count) return (a->count > b->count) ? -1 : 1;
    return strcmp(a->key, b->key);
}

/* write a kh_counter_t to CSV, sorted; keys are NOT freed here */
static int write_map_sorted_csv_(const kh_counter_t *m, const char *path, size_t topk) {
    /* count */
    size_t n = 0;
    for (khint_t i = kh_begin(m); i != kh_end(m); ++i)
        if (kh_exist(m, i)) ++n;

    FILE *fp = (path && strcmp(path, "-")==0) ? stdout : fopen(path ? path : "counts.csv", "w");
    if (!fp) { perror("fopen"); return 1; }

    fprintf(fp, "key,count\n");
    if (n == 0) {
        if (fp != stdout) fclose(fp);
        return 0;
    }

    /* collect */
    counter_entry_t *arr = (counter_entry_t*)malloc(n * sizeof(*arr));
    if (!arr) { perror("malloc"); if (fp != stdout) fclose(fp); return 1; }
    size_t j = 0;
    for (khint_t i = kh_begin(m); i != kh_end(m); ++i)
        if (kh_exist(m, i)) {
            arr[j].key   = kh_key(m, i);
            arr[j].count = kh_val(m, i);
            ++j;
        }

    /* sort & write */
    qsort(arr, n, sizeof(*arr), cmp_count_desc_then_key_asc);
    size_t limit = (topk > 0 && topk < n) ? topk : n;
    size_t hit = 0;
    for (size_t i = 0; i < limit; ++i) {
        if (arr[i].count > 0) ++hit;
        fprintf(fp, "%s,%zu\n", arr[i].key, arr[i].count);
    }
    
    fprintf(stderr, "\nSummary:\nTotal sequences searched: %zu \n", n);
    fprintf(stderr, "Sequences with hits: %zu \n", hit);
    fprintf(stderr, "Results Written to %s \n", path);
    
    

    free(arr);
    if (fp != stdout) fclose(fp);
    return 0;
}

/* PUBLIC: collapse "/rc" suffix into its base key before writing */
int counter_write_csv_sorted_collapse_rc(const kh_counter_t *m, const char *path, size_t topk)
{
    /* 1) aggregate into a new map keyed by the base (suffix removed) */
    kh_counter_t *agg = kh_init(counter);
    if (!agg) { fprintf(stderr, "OOM\n"); return 1; }

    for (khint_t i = kh_begin(m); i != kh_end(m); ++i) {
        if (!kh_exist(m, i)) continue;

        const char *orig = kh_key(m, i);
        size_t val = kh_val(m, i);

        char *base = dup_base_no_rc(orig);       /* malloc copy of canonical key */
        if (!base) { fprintf(stderr, "OOM\n"); kh_destroy(counter, agg); return 1; }

        int ret;
        khiter_t it = kh_put(counter, agg, base, &ret);
        if (ret < 0) {                           /* OOM inserting */
            free(base);
            kh_destroy(counter, agg);
            fprintf(stderr, "OOM\n");
            return 1;
        } else if (ret == 0) {                   /* existed: add and free our temp key */
            kh_val(agg, it) += val;
            free(base);
        } else {                                 /* new: table takes ownership of base */
            kh_key(agg, it) = base;              /* important: store our malloc'd key */
            kh_val(agg, it) = val;
        }
    }

    /* 2) write sorted CSV from the aggregated map */
    int rc = write_map_sorted_csv_(agg, path, topk);

    /* 3) free keys in agg and destroy */
    for (khint_t i = kh_begin(agg); i != kh_end(agg); ++i)
        if (kh_exist(agg, i)) free((char*)kh_key(agg, i));
    kh_destroy(counter, agg);

    return rc;
}





// Extract 2-bit base at position i (0..3). Assumes same packing as above.
static inline uint32_t get2(const uint64_t *text, size_t i) {
    return (text[i >> 5] >> ((i & 31) * 2)) & 3U;
}

// Check if base i is ambiguous (N). Only if you built a mask.
static inline int isN(const uint64_t *mask, size_t i) {
    return mask ? (int)((mask[i >> 6] >> (i & 63)) & 1U) : 0;
}