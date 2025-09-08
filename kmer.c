#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "kmer.h"

// Trie API
typedef struct TrieNode TrieNode;
bool trie_prefix_search(const TrieNode *root, const char *kmer);

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


void find_kmer(const char *s, const size_t s_len, 
               const TrieNode *trie,
               khash_t(strset) *kmer_hit,
               khash_t(strset) *kmer_search,
               size_t k,
               u32vec_t *hits) {

    if (!s || k <= 0 || (size_t)k > s_len) return;                    

    for (size_t i = 0; i + (size_t)k <= s_len; ++i) {
        // Build a temp key buffer for lookup (NUL-terminated)
        char *kmer = dup_cstr_n(s + i, (size_t)k);
        if (!kmer) { /* OOM */ break; }
        
        // if Kmer in kmer_hit; add pos to check
        if (kh_get(strset, kmer_hit, kmer) != kh_end(kmer_hit)) {
            kv_push(uint32_t, 0, *hits, i);
            continue;
        }

        // if kmer in searched cache; continue
        if (kh_get(strset, kmer_search, kmer) != kh_end(kmer_search)) {
            continue;
        }

        if (trie_prefix_search(trie, kmer)) {
            kv_push(uint32_t, 0, *hits, i);
            int ret;
            khiter_t hk = kh_put(strset, kmer_hit, kmer, &ret);
            if (ret > 0) kh_key(kmer_hit, hk) = strdup(kmer); // non-owning: reuse pointer
        }

        int ret;
        khiter_t sk = kh_put(strset, kmer_search, kmer, &ret);
        if (ret > 0) kh_key(kmer_search, sk) = strdup(kmer);
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
                  khash_t(strset) *found_sequences)
{
    khash_t(posset) *start_pos_cache = kh_init(posset);   // dedupe absolute start positions
    kFoundVec matches; kv_init(matches);

    for (size_t idx = 0; idx < hit->n; ++idx) {
        uint32_t h = hit->a[idx];

        //fprintf(stderr, "1 Looking For At Pos %u,%u,%u,%u\n", h, h + max_len, h + min_len, seq_len);

        if ((size_t)h + (size_t)min_len > seq_len) continue;

        // fprintf(stderr, "2 Looking For At Pos %u,%u,%u,%u\n", h, h + max_len, h + min_len, seq_len);

        size_t qlen = (size_t)max_len;
        if ((size_t)h + qlen > seq_len) 
            // qlen = seq_len - (size_t)h;
            continue;

        // call trie, gather matches for this slice
        matches.n = 0;  // reuse buffer
        
        char *slice = malloc(max_len + 1);
        if (!slice) /* handle OOM */;
        memcpy(slice, sequence + h, max_len);
        slice[max_len] = '\0';

        trie_search_exact(trie, slice, min_len, &matches);

        //fprintf(stderr, "Number of Matches: %u\n", matches.n);
        
        for (size_t m = 0; m < matches.n; ++m) {
            uint32_t abs_start = matches.a[m].pos + h;
            //fprintf(stderr, "Matches Start Pos: %u\n", abs_start);

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