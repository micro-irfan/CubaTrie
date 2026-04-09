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

static bool seq_is_acgt_only(const char *seq) {
    if (!seq) return false;
    for (size_t i = 0; seq[i]; ++i) {
        if (nt2bits(seq[i]) < 0) return false;
    }
    return true;
}

static size_t name_len_no_rc(const char *name) {
    if (!name) return 0;
    size_t n = strlen(name);
    if (n >= 3 && name[n-3] == '/' && name[n-2] == 'r' && name[n-1] == 'c') return n - 3;
    return n;
}

const TrieNode *trie_find_exact(const TrieNode *root, const char *seq) {
    if (!root || !seq) return NULL;

    const TrieNode *node = root;
    size_t pos = 0, len = strlen(seq);

    while (pos < len) {
        int idx = nt2bits(seq[pos]);
        if (idx < 0) return NULL;

        const char *label = node->edge_label[idx];
        if (!label) return NULL;

        size_t lablen = strlen(label);
        if (pos + lablen > len) return NULL;
        if (strncmp(seq + pos, label, lablen) != 0) return NULL;

        node = node->child[idx];
        if (!node) return NULL;
        pos += lablen;
    }

    return node->end ? node : NULL;
}

TrieInsertStatus trie_insert(TrieNode *root, const char *seq, const char *sequence_name, int rev) {
    if (!root || !seq || !sequence_name) return TRIE_INSERT_INVALID_BASE;
    if (!seq_is_acgt_only(seq)) return TRIE_INSERT_INVALID_BASE;
    if (trie_find_exact(root, seq)) return TRIE_INSERT_DUP;

    TrieNode *node = root;
    size_t pos = 0, len = strlen(seq), name_len = strlen(sequence_name);

    while (pos < len) {
        int idx = nt2bits(seq[pos]);
        if (idx < 0) return TRIE_INSERT_INVALID_BASE;

        if (node->edge_label[idx] == NULL) {
            char *edge = strndup(seq + pos, len - pos);
            TrieNode *child = trie_create_node();
            char *name_copy = strndup(sequence_name, name_len);
            char *seq_copy = strndup(seq, len);

            if (!edge || !child || !name_copy || !seq_copy) {
                free(edge);
                free(name_copy);
                free(seq_copy);
                if (child) trie_free_node(child);
                return TRIE_INSERT_OOM;
            }

            node->edge_label[idx] = edge;
            node->child[idx] = child;
            child->end = true;
            child->name = name_copy;
            child->seq = seq_copy;
            child->seq_len = len;
            child->rev = rev;
            return TRIE_INSERT_OK;
        }

        char *label = node->edge_label[idx];
        size_t common = lcp(seq + pos, label);

        if (common == strlen(label)) {
            pos += common;
            node = node->child[idx];
            continue;
        }

        TrieNode *mid = NULL;
        TrieNode *new_child = NULL;
        TrieNode *old_child = node->child[idx];
        char *shared = NULL;
        char *old_suffix = NULL;
        char *new_edge = NULL;
        char *new_name = NULL;
        char *new_seq = NULL;

        size_t len_label = strlen(label);
        if (common > len_label) common = len_label;

        const char *suffix_ptr = label + common;
        size_t suffix_len = len_label - common;
        const char *new_suffix = seq + pos + common;

        mid = trie_create_node();
        old_suffix = strndup(suffix_ptr, suffix_len);
        shared = (char*)malloc(common + 1);
        if (!mid || !old_suffix || !shared) {
            free(old_suffix);
            free(shared);
            if (mid) trie_free_node(mid);
            return TRIE_INSERT_OOM;
        }
        memcpy(shared, label, common);
        shared[common] = '\0';

        int old_idx = nt2bits(old_suffix[0]);
        if (old_idx < 0) {
            free(old_suffix);
            free(shared);
            trie_free_node(mid);
            return TRIE_INSERT_INVALID_BASE;
        }

        int new_idx = -1;
        if (*new_suffix == '\0') {
            new_name = strndup(sequence_name, name_len);
            new_seq = strndup(seq, len);
            if (!new_name || !new_seq) {
                free(old_suffix);
                free(shared);
                free(new_name);
                free(new_seq);
                trie_free_node(mid);
                return TRIE_INSERT_OOM;
            }
        } else {
            new_idx = nt2bits(*new_suffix);
            if (new_idx < 0) {
                free(old_suffix);
                free(shared);
                trie_free_node(mid);
                return TRIE_INSERT_INVALID_BASE;
            }
            new_edge = strndup(new_suffix, strlen(new_suffix));
            new_child = trie_create_node();
            new_name = strndup(sequence_name, name_len);
            new_seq = strndup(seq, len);
            if (!new_edge || !new_child || !new_name || !new_seq) {
                free(old_suffix);
                free(shared);
                free(new_edge);
                free(new_name);
                free(new_seq);
                if (new_child) trie_free_node(new_child);
                trie_free_node(mid);
                return TRIE_INSERT_OOM;
            }
        }

        // commit split only after all allocations succeed
        char *old_label = node->edge_label[idx];
        node->edge_label[idx] = shared;
        node->child[idx] = mid;

        mid->edge_label[old_idx] = old_suffix;
        mid->child[old_idx] = old_child;

        if (*new_suffix == '\0') {
            mid->end = true;
            mid->name = new_name;
            mid->seq = new_seq;
            mid->seq_len = len;
            mid->rev = rev;
        } else {
            mid->edge_label[new_idx] = new_edge;
            mid->child[new_idx] = new_child;
            mid->child[new_idx]->end = true;
            mid->child[new_idx]->name = new_name;
            mid->child[new_idx]->seq = new_seq;
            mid->child[new_idx]->seq_len = len;
            mid->child[new_idx]->rev = rev;
        }

        free(old_label);
        return TRIE_INSERT_OK;
    }

    // sequence ended exactly at current node
    char *new_name = strndup(sequence_name, name_len);
    char *new_seq = strndup(seq, len);
    if (!new_name || !new_seq) {
        free(new_name);
        free(new_seq);
        return TRIE_INSERT_OOM;
    }

    free(node->name);
    free(node->seq);
    node->name = new_name;
    node->seq = new_seq;
    node->end = true;
    node->seq_len = len;
    node->rev = rev;
    return TRIE_INSERT_OK;
}

typedef struct SamNameSeen {
    char *name;
    struct SamNameSeen *next;
} SamNameSeen;

static int sam_seen_contains(SamNameSeen *head, const char *name, size_t len) {
    for (SamNameSeen *p = head; p; p = p->next) {
        if (strncmp(p->name, name, len) == 0 && p->name[len] == '\0') return 1;
    }
    return 0;
}

static int sam_seen_add(SamNameSeen **head, const char *name, size_t len) {
    SamNameSeen *node = (SamNameSeen*)calloc(1, sizeof(*node));
    if (!node) return -1;
    node->name = (char*)malloc(len + 1);
    if (!node->name) {
        free(node);
        return -1;
    }
    memcpy(node->name, name, len);
    node->name[len] = '\0';
    node->next = *head;
    *head = node;
    return 0;
}

static void sam_seen_free(SamNameSeen *head) {
    while (head) {
        SamNameSeen *next = head->next;
        free(head->name);
        free(head);
        head = next;
    }
}

static void trie_write_sam_sq_rec(FILE *fp, const TrieNode *node, SamNameSeen **seen) {
    if (!fp || !node) return;
    if (node->end && node->name && node->seq_len > 0) {
        size_t name_len = name_len_no_rc(node->name);
        if (name_len > 0 && !sam_seen_contains(*seen, node->name, name_len)) {
            if (sam_seen_add(seen, node->name, name_len) == 0) {
                fprintf(fp, "@SQ\tSN:%.*s\tLN:%zu\n", (int)name_len, node->name, node->seq_len);
            } else {
                fprintf(fp, "@SQ\tSN:%.*s\tLN:%zu\n", (int)name_len, node->name, node->seq_len);
            }
        }
    }
    for (int i = 0; i < ALPHABET_SIZE; ++i) {
        trie_write_sam_sq_rec(fp, node->child[i], seen);
    }
}

void trie_write_sam_header(FILE *fp, const TrieNode *root) {
    if (!fp || !root) return;
    fprintf(fp, "@HD\tVN:1.6\tSO:unknown\n");
    SamNameSeen *seen = NULL;
    trie_write_sam_sq_rec(fp, root, &seen);
    sam_seen_free(seen);
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

            if (node->end && node->seq_len >= min_len)
                mv_push(out, (uint32_t)start, node->name, node->seq, used);
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
        if (idx > 3) return false;

        const char *label = node->edge_label[idx];
        if (!label) return false;

        size_t lablen = strlen(label);
        size_t i = 0;

        while (i < lablen && p + i < n && kmer2bit(*b2, n, p + i) == label[i]) {
            ++i;
        }

        if (p + i == n) return true;
        if (i < lablen) return false;

        p += lablen;
        node = node->child[idx];
    }

    return true;
}
