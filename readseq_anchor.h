#ifndef READSEQ_ANCHOR_H
#define READSEQ_ANCHOR_H

#include <stddef.h>
#include <stdint.h>

#include "readseq.h"

typedef struct {
    size_t len;
    char seq[64]; // uppercase A/C/G/T, null-terminated (max len 63)
    uint64_t mask[256];
} BitapPattern;

typedef struct {
    size_t start;
    size_t end; // exclusive
    int errors;
} AnchorMatch;

typedef struct {
    int enabled;
    int max_error;
    int has_anchor5;
    int has_anchor3;
    size_t expected_insert_len;
    size_t min_insert_len;
    size_t max_insert_len;
    BitapPattern a5;
    BitapPattern a3;
    BitapPattern a5_rc;
    BitapPattern a3_rc;
} AnchorRuntime;

int anchor_runtime_init(AnchorRuntime *out,
                        const AnchorConfig *cfg,
                        size_t expected_insert_len);

int anchor_runtime_init_range(AnchorRuntime *out,
                              const AnchorConfig *cfg,
                              size_t min_insert_len,
                              size_t max_insert_len);

int anchor_extract_window(const AnchorRuntime *ar,
                          const char *read_seq,
                          size_t read_len,
                          size_t *insert_start_out,
                          size_t *insert_len_out);

int anchor_extract_window_range(const AnchorRuntime *ar,
                                const char *read_seq,
                                size_t read_len,
                                size_t *insert_start_out,
                                size_t *insert_len_out);

#endif
