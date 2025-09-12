#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "kmer.h"
#include "utils.h"

// Trie API
typedef struct TrieNode TrieNode;
bool trie_prefix_search(const TrieNode *root,
                        const uint64_t *b2, 
                        size_t n);

// comparator for qsort (ascending)
static int cmp_u32_asc(const void *a, const void *b) {
    uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

// small strdup (portable)
static char *dup_cstr_n(const char *s, size_t n) {
    char *p = (char*)malloc(n + 1);
    if (!p) return NULL;
    memcpy(p, s, n);
    p[n] = '\0';
    return p;
}



// reverse-complement a 2-bit base
static inline uint64_t rc2(uint64_t b){ return (~b) & 3; }

static inline uint64_t hash64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    return x ^ (x >> 33);
}



// Pack text[0..n) into 2-bit words (LSB-first). Also fill ambig-mask (1 bit per base) if provided.
static inline int encode_kmer(const char *s, size_t n, uint64_t *out) {
    uint64_t code = 0;
    for (size_t i = 0; i < n; ++i) {
        int b = nt2bits(s[i]);
        if (b < 0) return -1;
        code = (code << 2) | (uint64_t)b;

    }

    *out = code;
    return 0;
}

void find_kmer(const char *s, const size_t s_len, 
               const TrieNode *trie,
               khash_t(kset64) *kmer_hit,
               khash_t(kset64) *kmer_search,
               size_t k,
               u32vec_t *hits) {

    if (!s || k <= 0 || (size_t)k > s_len) return;                    

    for (size_t i = 0; i + (size_t)k <= s_len; ++i) {
        // Build a temp key buffer for lookup (NUL-terminated)
        char *kmer = dup_cstr_n(s + i, (size_t)k);
        if (!kmer) { /* OOM */ break; }

        printf("%s %d\n", kmer, i);

        uint64_t code;
        int n_found = encode_kmer(kmer, k, &code);
        free(kmer);

        if (n_found < 0) continue;

        if (kh_get(kset64, kmer_hit, (khint64_t)code) != kh_end(kmer_hit)) {
            kv_push(uint32_t, 0, *hits, i);
            continue;
        }

        // if kmer in searched cache; continue
        if (kh_get(kset64, kmer_search, (khint64_t)code) != kh_end(kmer_search)) {
            continue;
        }

        if (trie_prefix_search(trie, &code, k)) {
            kv_push(uint32_t, 0, *hits, i);
            int ret;
            khiter_t sk = kh_put(kset64, kmer_hit, (khint64_t)code, &ret);
            if (ret > 0) kh_key(kmer_hit, sk) = (khint64_t)code;
        }

        int ret;
        khiter_t sk = kh_put(kset64, kmer_search, (khint64_t)code, &ret);
        if (ret > 0) kh_key(kmer_search, sk) = (khint64_t)code;
    }

    // sort hits (ascending)
    if (hits->n > 1) qsort(hits->a, hits->n, sizeof(uint32_t), cmp_u32_asc);

    return;

}

// ---------- core routine ----------
/* sequence/seq_len : the read bases
   hit              : vector of candidate start offsets
   min_len, max_len : same semantics as in Python
   trie             : trie root
   found_sequences  : output set of names (non-owning; don't free keys here)
*/

void add_to_counter(khash_t(strset) *found_sequences,
                    kh_counter_t *map)
{
    for (khint_t i = kh_begin(found_sequences);
         i != kh_end(found_sequences); ++i)
    {
        if (!kh_exist(found_sequences, i)) continue;
        const char *name = kh_key(found_sequences, i);  // non-owning pointer
        counter_inc(map, name);
    }
}

void find_matches(const char *sequence, size_t seq_len,
                  const u32vec_t *hit,
                  uint32_t min_len, uint32_t max_len,
                  const TrieNode *trie,
                  khash_t(strset) *found_sequences,
                  int k_mm)
{
    khash_t(posset) *start_pos_cache = kh_init(posset);   // dedupe absolute start positions
    kFoundVec matches; kv_init(matches);

    for (size_t idx = 0; idx < hit->n; ++idx) {
        uint32_t h = hit->a[idx];

        if ((size_t)h + (size_t)min_len > seq_len) continue;

        size_t qlen = (size_t)max_len;
        if ((size_t)h + qlen > seq_len) 
            continue;

        // call trie, gather matches for this slice
        matches.n = 0;  // reuse buffer
        
        char *slice = malloc(max_len + 1);
        if (!slice) /* handle OOM */;
        memcpy(slice, sequence + h, max_len);
        slice[max_len] = '\0';

        normalize_acgt(slice);

        // printf("%s at position %d \n", slice, h);

        trie_search_exact(trie, slice, min_len, k_mm, &matches);

        for (size_t m = 0; m < matches.n; ++m) {
            uint32_t abs_start = matches.a[m].pos + h;

            // if abs_start in start_pos_cache: continue
            khiter_t it = kh_get(posset, start_pos_cache, (khint_t)abs_start);
            if (it != kh_end(start_pos_cache)) continue;

            // mark visited start position
            int ret;
            it = kh_put(posset, start_pos_cache, (khint_t)abs_start, &ret);
            (void)ret;

            // found_sequences.add(m.sequence_name) (non-owning)
            khiter_t jt = kh_put(strset, found_sequences,
                                 (char*)matches.a[m].name, &ret);

            if (ret > 0) {
                // insert; reuse pointer owned by trie (don’t free in this set)
                kh_key(found_sequences, jt) = strdup(matches.a[m].name);
            }
        }

        free(slice);
    }

    kv_destroy(matches);
    kh_destroy(posset, start_pos_cache);
}