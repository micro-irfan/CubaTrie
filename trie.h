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


/* --------- Trie API --------- */
/* Create an empty trie node */
TrieNode *trie_create_node(void);

/* Insert a DNA sequence with a name into the trie */
void trie_insert(TrieNode *root, const char *seq, const char *sequence_name);

/* Free the entire trie (all nodes + strings) */
void trie_free_node(TrieNode *root);

/* Check if `text` is a prefix of some stored sequence in the trie */
bool trie_prefix_search(const TrieNode *root, const char *text);

void trie_search_exact(const TrieNode *root, const char *text, size_t min_len, kFoundVec *out);

#endif /* TRIE_H */