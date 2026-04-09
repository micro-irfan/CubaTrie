#ifndef TRIE_H
#define TRIE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "utils.h"

/* Alphabet size for DNA (A,C,G,T) */
#define ALPHABET_SIZE 4

typedef struct TrieNode {
    // For each nucleotide, we store an edge label (suffix string) and a child node
    char *edge_label[ALPHABET_SIZE];   // malloc'd strings (may be NULL)
    struct TrieNode *child[ALPHABET_SIZE];   // child pointers (may be NULL)

    // Terminal info (sequence ends at this node)
    bool end, rev;
    char *name, *seq, *name_alt;
    size_t seq_len;
} TrieNode;

typedef enum TrieDupPolicy {
    TRIE_DUP_ERROR = 0,  // stop loading on duplicate sequence
    TRIE_DUP_WARN  = 1,  // warn and skip duplicate sequence
    TRIE_DUP_IGNORE = 2  // silently skip duplicate sequence
} TrieDupPolicy;

typedef enum TrieInsertStatus {
    TRIE_INSERT_OK = 0,
    TRIE_INSERT_DUP = 1,
    TRIE_INSERT_INVALID_BASE = 2,
    TRIE_INSERT_OOM = 3
} TrieInsertStatus;


/* --------- Trie API --------- */
/* Create an empty trie node */
TrieNode *trie_create_node(void);

/* Insert a DNA sequence with a name into the trie */
TrieInsertStatus trie_insert(TrieNode *root, const char *seq, const char *sequence_name, int rev);

/* Find terminal node for exact sequence; NULL if not present */
const TrieNode *trie_find_exact(const TrieNode *root, const char *seq);

/* Write minimal SAM header (@HD + @SQ) from terminal trie entries */
void trie_write_sam_header(FILE *fp, const TrieNode *root);

/* Free the entire trie (all nodes + strings) */
void trie_free_node(TrieNode *root);

/* Check if `text` is a prefix of some stored sequence in the trie */
bool trie_prefix_search(const TrieNode *root,
                        const uint64_t *b2, 
                        size_t n);

void trie_search_exact(const TrieNode *root, 
                       const char *text, 
                       size_t min_len, int k_mm,  
                       kFoundVec *out);


#endif /* TRIE_H */
