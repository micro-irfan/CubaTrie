#ifndef READSEQ_INTERNAL_H
#define READSEQ_INTERNAL_H

#include <stddef.h>
#include <stdio.h>

#include "readseq.h"
#include "readseq_anchor.h"
#include "kmer.h"

int readseq_process_one_read(const char *read_name,
                             const char *read_seq,
                             const char *read_qual,
                             size_t read_len,
                             kh_counter_t *counts,
                             const KmerBitset *seed_index,
                             const AnchorRuntime *anchor_runtime,
                             int seed_mm,
                             size_t min_len,
                             size_t max_len,
                             int k_mm,
                             int exclude_multihit,
                             int sam_soft_clip,
                             int sam_emit_unmapped,
                             FILE *sam_fp,
                             int counts_preseeded);

int readseq_load_fastq_mt(const char *path,
                          TrieNode *root,
                          kh_counter_t *counts,
                          int kmerlen,
                          int seed_mm,
                          size_t *min_len,
                          size_t *max_len,
                          int k_mm,
                          int exclude_multihit,
                          int sam_soft_clip,
                          int sam_emit_unmapped,
                          FILE *sam_fp,
                          unsigned threads,
                          const AnchorConfig *anchor_cfg);

#endif
