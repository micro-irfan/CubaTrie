#include <stdio.h>
#include <stdlib.h>

#define _POSIX_C_SOURCE 200809L
#include <string.h>
#include <stdbool.h>
#include "trie.h"
#include "utils.h"

#define ALPHABET_SIZE 4  

TrieNode *trie_create_node(void) {
    TrieNode *root = (TrieNode*)calloc(1, sizeof *root);
    return root;
}
    
// shallow recursive free
void trie_free_node(TrieNode *root) { 
    if (!root) return;
    for (int i = 0; i < ALPHABET_SIZE; ++i) {
        free(root->edge_label[i]);
        trie_free_node(root->child[i]);
    }
    free(root->name);
    free(root->seq);
    free(root);
}


// longest common prefix
static size_t lcp(const char *a, const char *b) { 
    size_t i = 0;
    while (a[i] && b[i] && a[i] == b[i]) ++i;
    return i;
}

void trie_insert(TrieNode *root, const char *seq, const char *sequence_name, int rev) {
    TrieNode *node = root;
    size_t pos = 0, len = strlen(seq);

    while (pos < len) {
        int idx = nt2bits(seq[pos]);
        if (idx < 0) { ++pos; continue; } // skip invalids

        if (node->edge_label[idx] == NULL) {
            node->edge_label[idx] = strndup(seq + pos, len - pos);
            node->child[idx] = trie_create_node();
            node->child[idx]->end = true;
            node->child[idx]->name = strndup(sequence_name, strlen(sequence_name));
            node->child[idx]->seq = strndup(seq, len);
            node->child[idx]->seq_len = len;
            node->child[idx]->rev = rev;
            return;
        }

        char *label = node->edge_label[idx];
        size_t common = lcp(seq + pos, label);

        if (common == strlen(label)) {          // full label match → descend
            pos += common;
            node = node->child[idx];
            continue;
        }

        // split edge
        TrieNode *mid = trie_create_node();
        TrieNode *old_child = node->child[idx];

        size_t len_label = strlen(label);

        if (common > len_label) {
            // lcp must never exceed the shorter string; treat as logic error.
            // handle as you prefer: abort, clamp, or return error.
            common = len_label;
        }

        const char *suffix_ptr = label + common;
        size_t suffix_len = len_label - common;

        char *old_suffix = strndup(suffix_ptr, suffix_len);

        char *shared = (char*)malloc(common + 1);
        memcpy(shared, label, common);
        shared[common] = '\0';

        free(node->edge_label[idx]);
        node->edge_label[idx] = shared;
        node->child[idx] = mid;

        int old_idx = nt2bits(old_suffix[0]);
        mid->edge_label[old_idx] = old_suffix;
        mid->child[old_idx] = old_child;

        const char *new_suffix = seq + pos + common;
        if (*new_suffix == '\0') {
            mid->end = true;
            mid->name = strndup(sequence_name, strlen(sequence_name));
            mid->seq = strndup(seq, len);
            mid->seq_len = len;
            mid->rev = rev;
        } else {
            int new_idx = nt2bits(*new_suffix);
            if (new_idx >= 0) {
                mid->edge_label[new_idx] = strndup(new_suffix, strlen(new_suffix));
                mid->child[new_idx] = trie_create_node();
                mid->child[new_idx]->end = true;
                mid->child[new_idx]->name = strndup(sequence_name, strlen(sequence_name));
                mid->child[new_idx]->seq = strndup(seq, len);
                mid->child[new_idx]->seq_len = len;
                mid->child[new_idx]->rev = rev;
            }
        }
        return;
    }

    // sequence ended exactly at current node
    node->end = true;
    free(node->name);
    node->name = strndup(sequence_name, strlen(sequence_name));
    free(node->seq);
    node->seq = strndup(seq, len);
    node->seq_len = len;
}


/* -------- EXACT SEARCH --------
   Scans `text` and, at every start, tries to follow trie edges
   by matching each edge label exactly against the text.
   When a terminal node is reached, we report a match.
*/

void trie_search_exact(const TrieNode *root, 
                       const char *text, size_t min_len, int k_mm,
                       kFoundVec *out) 
{
    size_t n = strlen(text);

    for (size_t start = 0; start + min_len <= n; ++start) {
        const TrieNode *node = root;
        size_t cur = start;
        int used = 0;

        while (cur < n) {
            int idx = nt2bits(text[cur]);
            if (idx < 0) break;                            

            const char *label = node->edge_label[idx];
            if (!label) break;                            
            
            size_t lablen = strlen(label);
            // text too short for label
            if (cur + lablen > n) break;
            
            // compare label with mismatches allowed
            int mm = 0;
            for (size_t j = 0; j < lablen; ++j) {
                if (text[cur + j] != label[j]) {
                    if (++mm > k_mm - used) { mm = k_mm + 1; break; }
                }
            }
            if (mm > k_mm - used) break;

            used += mm;
            cur  += lablen;
            node  = node->child[idx];

            // full-sequence check not needed anymore (we already matched along the path),
            // but keep it if your trie can have shared labels/aliases.
            if (node->end && node->seq_len >= min_len) 
                mv_push(out, (uint32_t)start, node->name, node->seq); 
                
            
            // edge label must match exactly
            // if (memcmp(text + cur, label, lablen) != 0)   
            //     break;
            
            // consume edge and descend
            // cur += lablen;
            // node = node->child[idx];

            // if this node is terminal, we have an exact match of its full sequence
            // if (node->end && node->seq_len >= min_len) {
            //     if (start + node->seq_len <= n &&
            //         memcmp(text + start, node->seq, node->seq_len) == 0) {
            //         mv_push(out, start, node->name, node->seq);
            //     }
            //     // continue; there might be a longer sequence sharing this path
            // }
        }
    }
}


bool trie_prefix_search(const TrieNode *root,
                        const uint64_t *b2, 
                        size_t n) {

    const TrieNode *node = root;
    size_t p = 0;

    while (p < n) {
        uint32_t idx = (uint32_t)((*b2 >> (2 * (n - 1 - (int)p))) & 3u);
        
        if (idx < 0 || idx > 3) return false;                 

        const char *label = node->edge_label[idx];   // invalid char in query
     
        // no edge starting with this nt
        if (!label) return false;                  

        size_t lablen = strlen(label);
        size_t i = 0;

        // Compare query with the edge label 
        while (i < lablen && p + i < n && kmer2bit(*b2, n, p + i) == label[i]) {
            ++i;
        }

        // whole query is consumed (possibly mid-edge) → query is a prefix
        if (p + i == n) return true;

        // Mismatch before finishing the edge and query still has chars → not a prefix
        if (i < lablen) return false;

        // Full edge matched; descend and continue with remaining query
        p += lablen;
        node = node->child[idx];
    }

    return true; // empty string is trivially a prefix
}




