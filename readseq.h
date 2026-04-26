#ifndef READSEQ_H
#define READSEQ_H

#include <stddef.h>
#include <stdio.h>

#include "trie.h"
#include "utils.h"

typedef struct {
    int enabled;
    const char *anchor5;
    const char *anchor3;
    int max_error;
} AnchorConfig;

int load_reference(const char *path,
                   TrieNode *root,
                   kh_counter_t *map,
                   size_t *min_out,
                   size_t *max_out,
                   size_t *n_out,
                   int add_revcomp,
                   size_t *kmer_len,
                   TrieDupPolicy dup_policy);

int load_fastq(const char *path,
               TrieNode *root,
               kh_counter_t *counts,
               int kmerlen,
               int seed_mm,
               size_t *min_len,
               size_t *max_len,
               int k_mm,
               int exclude_multihit,
               FILE *sam_fp,
               int sam_soft_clip,
               int sam_emit_unmapped,
               unsigned threads,
               const AnchorConfig *anchor_cfg);

int cut_fastq_by_anchors(const char *in_path,
                         const char *out_path,
                         const AnchorConfig *anchor_cfg,
                         size_t min_len,
                         size_t max_len,
                         int check_revcomp,
                         unsigned threads);

#endif
