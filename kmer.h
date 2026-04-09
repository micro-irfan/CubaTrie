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
KHASH_SET_INIT_INT64(kset64)         // set<uint64_t>

// Sets of strings
KHASH_SET_INIT_INT(posset)    // set<uint32_t> 
KHASH_SET_INIT_STR(strset)    // set<char*>    

#include "trie.h"

void find_kmer(const char *s, const size_t s_len, 
               const TrieNode *trie,
               khash_t(kset64) *kmer_hit,
               khash_t(kset64) *kmer_search,
               size_t k,
               u32vec_t *hits);                    

khash_t(kset64) *kmer_set_from_trie(const TrieNode *root, size_t k);


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
