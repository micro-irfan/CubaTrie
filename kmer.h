#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "kvec.h"

typedef kvec_t(uint32_t) u32vec_t;

#include "khash.h"

// Map: char* (k-mer) -> vector of uint32_t positions
KHASH_MAP_INIT_STR(k2pos, u32vec_t)
KHASH_MAP_INIT_INT64(k2pos64, u32vec_t)   // map: uint64_t -> positions

// Sets of strings
KHASH_SET_INIT_INT(posset)    // set<uint32_t> 
KHASH_SET_INIT_STR(strset)    // set<char*>    

#include "trie.h"

typedef struct TrieCursorState TrieCursorState;

typedef struct {
    size_t k;
    size_t nbits;
    size_t nwords;
    uint64_t *words;
    int32_t *state_idx;
    TrieCursorState *states;
    size_t nstates;
    size_t mstates;
} KmerBitset;

KmerBitset *kmer_bitset_from_trie(const TrieNode *root, size_t k);
void kmer_bitset_destroy(KmerBitset *index);
int kmer_bitset_test(const KmerBitset *index, uint64_t code);

void find_kmer_bitset(const char *s,
                      size_t s_len,
                      const KmerBitset *index,
                      int seed_mm,
                      u32vec_t *hits);


void add_to_counter(khash_t(strset) *found_sequences,
                    kh_counter_t *map);

void find_matches(const char *sequence, size_t seq_len,
                  const u32vec_t *hit,
                  uint32_t min_len, uint32_t max_len,
                  const TrieNode *trie,
                  khash_t(strset) *found_sequences, 
                  int k_mm,
                  const char *read_name,
                  const char *read_qual,
                  FILE *sam_fp);

void find_matches_seeded(const char *sequence, size_t seq_len,
                         const u32vec_t *hit,
                         uint32_t min_len, uint32_t max_len,
                         const TrieNode *trie,
                         const KmerBitset *seed_index,
                         int seed_mm,
                         khash_t(strset) *found_sequences,
                         int k_mm,
                         const char *read_name,
                         const char *read_qual,
                         FILE *sam_fp);
