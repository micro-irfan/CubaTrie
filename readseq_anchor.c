#include "readseq_anchor.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

static char *anchor_revcomp_new_n(const char *s, size_t n) {
    char *rc = (char*)malloc(n + 1);
    if (!rc) return NULL;
    for (size_t i = 0; i < n; ++i) {
        char c = s[n - 1 - i];
        char r = 'N';
        if (c == 'A') r = 'T';
        else if (c == 'C') r = 'G';
        else if (c == 'G') r = 'C';
        else if (c == 'T') r = 'A';
        rc[i] = r;
    }
    rc[n] = '\0';
    return rc;
}

static int bitap_pattern_build(BitapPattern *out, const char *pattern) {
    if (!out || !pattern) return -1;
    size_t m = strlen(pattern);
    if (m == 0 || m > 63) return -1;

    uint64_t empty = ((uint64_t)1 << m) - 1ULL;
    out->len = m;
    for (size_t i = 0; i < 256; ++i) out->mask[i] = empty;

    for (size_t i = 0; i < m; ++i) {
        unsigned char c = (unsigned char)toupper((unsigned char)pattern[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return -1;
        out->seq[i] = (char)c;
        out->mask[c] &= ~(1ULL << i);
    }
    out->seq[m] = '\0';
    return 0;
}

// Semi-global edit-distance scan: find earliest-ending substring match under one error budget.
static int bitap_find_first(const BitapPattern *pat,
                            const char *text,
                            size_t text_len,
                            size_t from_pos,
                            int max_error,
                            AnchorMatch *out_hit) {
    if (!pat || !text || !out_hit) return -1;
    if (from_pos >= text_len || pat->len == 0) return 1;

    size_t m = pat->len;
    size_t n = text_len - from_pos;
    if (m > 63) return -1; // guarded by pattern build, keep here for safety

    int prev_buf[64];
    int curr_buf[64];
    size_t prev_start_buf[64];
    size_t curr_start_buf[64];
    int *prev = prev_buf;
    int *curr = curr_buf;
    size_t *prev_start = prev_start_buf;
    size_t *curr_start = curr_start_buf;

    int mm = max_error;
    if (mm < 0) mm = 0;
    if ((size_t)mm > m) mm = (int)m;

    for (size_t i = 0; i <= m; ++i) {
        prev[i] = (int)i;
        prev_start[i] = 0;
    }

    for (size_t col = 1; col <= n; ++col) {
        curr[0] = 0;
        curr_start[0] = col; // empty pattern can start at current boundary

        unsigned char tc = (unsigned char)toupper((unsigned char)text[from_pos + col - 1]);

        for (size_t row = 1; row <= m; ++row) {
            int cost = (tc == (unsigned char)pat->seq[row - 1]) ? 0 : 1;

            int sub = prev[row - 1] + cost; // match/substitute
            int ins = prev[row] + 1;        // insertion in text
            int del = curr[row - 1] + 1;    // deletion in text

            int best = sub;
            size_t best_start = prev_start[row - 1];

            if (ins < best || (ins == best && prev_start[row] < best_start)) {
                best = ins;
                best_start = prev_start[row];
            }
            if (del < best || (del == best && curr_start[row - 1] < best_start)) {
                best = del;
                best_start = curr_start[row - 1];
            }

            curr[row] = best;
            curr_start[row] = best_start;
        }

        if (curr[m] <= mm) {
            out_hit->start = from_pos + curr_start[m];
            out_hit->end = from_pos + col;
            out_hit->errors = curr[m];
            return 0;
        }

        int *tmp_i = prev;
        prev = curr;
        curr = tmp_i;
        size_t *tmp_s = prev_start;
        prev_start = curr_start;
        curr_start = tmp_s;
    }

    return 1;
}

static size_t anchor_collect_matches(const BitapPattern *pat,
                                     const char *read_seq,
                                     size_t read_len,
                                     int max_error,
                                     size_t from_pos,
                                     AnchorMatch *hits,
                                     size_t hits_cap) {
    if (!pat || !read_seq || !hits || hits_cap == 0) return 0;
    size_t n = 0;
    size_t cursor = from_pos;
    while (n < hits_cap && cursor < read_len) {
        AnchorMatch h = {0};
        int rc = bitap_find_first(pat, read_seq, read_len, cursor, max_error, &h);
        if (rc != 0) break;
        hits[n++] = h;
        cursor = (h.start + 1 > cursor) ? (h.start + 1) : (cursor + 1);
    }
    return n;
}

static size_t anchor_collect_matches_exact(const BitapPattern *pat,
                                           const char *read_seq,
                                           size_t read_len,
                                           size_t from_pos,
                                           AnchorMatch *hits,
                                           size_t hits_cap) {
    if (!pat || !read_seq || !hits || hits_cap == 0) return 0;
    size_t m = pat->len;
    if (m == 0 || from_pos >= read_len || read_len - from_pos < m) return 0;

    size_t n = 0;
    for (size_t s = from_pos; s + m <= read_len && n < hits_cap; ++s) {
        size_t j = 0;
        while (j < m) {
            unsigned char tc = (unsigned char)toupper((unsigned char)read_seq[s + j]);
            if (tc != (unsigned char)pat->seq[j]) break;
            ++j;
        }
        if (j == m) {
            hits[n].start = s;
            hits[n].end = s + m;
            hits[n].errors = 0;
            ++n;
        }
    }
    return n;
}

typedef struct {
    int valid;
    int ambiguous;
    int errors;
    size_t first_start;
    size_t insert_start;
    size_t insert_len;
} AnchorCandidate;

static int candidate_score_cmp(const AnchorCandidate *a, const AnchorCandidate *b) {
    if (a->errors != b->errors) return (a->errors < b->errors) ? -1 : 1;
    if (a->first_start != b->first_start) return (a->first_start < b->first_start) ? -1 : 1;
    return 0;
}

static int candidate_same_window(const AnchorCandidate *a, const AnchorCandidate *b) {
    return a->insert_start == b->insert_start && a->insert_len == b->insert_len;
}

static void candidate_consider(AnchorCandidate *best, const AnchorCandidate *cand) {
    if (!best || !cand || !cand->valid) return;
    if (!best->valid) {
        *best = *cand;
        best->ambiguous = 0;
        return;
    }

    int cmp = candidate_score_cmp(cand, best);
    if (cmp < 0) {
        *best = *cand;
        best->ambiguous = 0;
        return;
    }
    if (cmp == 0 && !candidate_same_window(cand, best)) {
        best->ambiguous = 1;
    }
}

static void candidate_from_pair(const AnchorMatch *start_hit,
                                const AnchorMatch *end_hit,
                                size_t min_len,
                                size_t max_len,
                                AnchorCandidate *out) {
    if (!start_hit || !end_hit || !out) return;
    if (end_hit->start <= start_hit->end) return;
    size_t payload_len = end_hit->start - start_hit->end;
    if (payload_len < min_len || payload_len > max_len) return;

    out->valid = 1;
    out->errors = start_hit->errors + end_hit->errors;
    out->first_start = start_hit->start;
    out->insert_start = start_hit->end;
    out->insert_len = payload_len;
}

static AnchorCandidate anchor_find_best_pair_exact(const BitapPattern *start_pat,
                                                   const BitapPattern *end_pat,
                                                   const char *read_seq,
                                                   size_t read_len,
                                                   size_t min_len,
                                                   size_t max_len) {
    enum { START_CAP = 8, END_CAP = 8 };
    AnchorCandidate best = {0};
    AnchorMatch start_hits[START_CAP];
    size_t nstart = anchor_collect_matches_exact(start_pat, read_seq, read_len, 0, start_hits, START_CAP);

    for (size_t i = 0; i < nstart; ++i) {
        AnchorMatch s = start_hits[i];
        if (s.end >= read_len) continue;

        AnchorMatch end_hits[END_CAP];
        size_t nend = anchor_collect_matches_exact(end_pat, read_seq, read_len, s.end, end_hits, END_CAP);
        for (size_t j = 0; j < nend; ++j) {
            AnchorCandidate cand = {0};
            candidate_from_pair(&s, &end_hits[j], min_len, max_len, &cand);
            candidate_consider(&best, &cand);
        }

        if (best.valid && !best.ambiguous && best.errors == 0 && best.first_start == s.start) {
            return best;
        }
    }
    return best;
}

static AnchorCandidate anchor_find_best_5only_exact(const BitapPattern *start_pat,
                                                    const char *read_seq,
                                                    size_t read_len,
                                                    size_t min_len,
                                                    size_t max_len) {
    enum { START_CAP = 8 };
    AnchorCandidate best = {0};
    AnchorMatch hits[START_CAP];
    size_t n = anchor_collect_matches_exact(start_pat, read_seq, read_len, 0, hits, START_CAP);
    int fixed_len_mode = (min_len == max_len);
    size_t fixed_len = min_len;

    for (size_t i = 0; i < n; ++i) {
        if (hits[i].end > read_len) continue;
        size_t payload_len = 0;
        size_t insert_start = hits[i].end;
        if (fixed_len_mode) {
            if (insert_start + fixed_len > read_len) continue;
            payload_len = fixed_len;
        } else {
            payload_len = read_len - insert_start;
            if (payload_len < min_len || payload_len > max_len) continue;
        }

        AnchorCandidate cand = {0};
        cand.valid = 1;
        cand.errors = 0;
        cand.first_start = hits[i].start;
        cand.insert_start = insert_start;
        cand.insert_len = payload_len;
        candidate_consider(&best, &cand);
        if (best.valid && !best.ambiguous) return best;
    }
    return best;
}

static AnchorCandidate anchor_find_best_3only_exact(const BitapPattern *end_pat,
                                                    const char *read_seq,
                                                    size_t read_len,
                                                    size_t min_len,
                                                    size_t max_len) {
    enum { END_CAP = 8 };
    AnchorCandidate best = {0};
    AnchorMatch hits[END_CAP];
    size_t n = anchor_collect_matches_exact(end_pat, read_seq, read_len, 0, hits, END_CAP);
    int fixed_len_mode = (min_len == max_len);
    size_t fixed_len = min_len;

    for (size_t i = 0; i < n; ++i) {
        size_t payload_len = 0;
        size_t insert_start = 0;
        if (fixed_len_mode) {
            if (hits[i].start < fixed_len) continue;
            payload_len = fixed_len;
            insert_start = hits[i].start - fixed_len;
        } else {
            payload_len = hits[i].start;
            if (payload_len < min_len || payload_len > max_len) continue;
            insert_start = 0;
        }

        AnchorCandidate cand = {0};
        cand.valid = 1;
        cand.errors = 0;
        cand.first_start = hits[i].start;
        cand.insert_start = insert_start;
        cand.insert_len = payload_len;
        candidate_consider(&best, &cand);
        if (best.valid && !best.ambiguous) return best;
    }
    return best;
}

static AnchorCandidate anchor_find_best_pair(const BitapPattern *start_pat,
                                             const BitapPattern *end_pat,
                                             const char *read_seq,
                                             size_t read_len,
                                             int max_error,
                                             size_t min_len,
                                             size_t max_len) {
    enum { START_CAP = 8, END_CAP = 8 };
    AnchorCandidate best = {0};
    AnchorMatch start_hits[START_CAP];
    size_t nstart = anchor_collect_matches(start_pat, read_seq, read_len, max_error, 0, start_hits, START_CAP);

    for (size_t i = 0; i < nstart; ++i) {
        AnchorMatch s = start_hits[i];
        if (s.end >= read_len) continue;

        AnchorMatch end_hits[END_CAP];
        size_t nend = anchor_collect_matches(end_pat, read_seq, read_len, max_error, s.end, end_hits, END_CAP);
        for (size_t j = 0; j < nend; ++j) {
            AnchorCandidate cand = {0};
            candidate_from_pair(&s, &end_hits[j], min_len, max_len, &cand);
            candidate_consider(&best, &cand);
        }

        // For this orientation, after finishing one start position, a 0-error best is unbeatable.
        if (best.valid && !best.ambiguous && best.errors == 0 && best.first_start == s.start) {
            return best;
        }
    }
    return best;
}

static AnchorCandidate anchor_find_best_5only(const BitapPattern *start_pat,
                                              const char *read_seq,
                                              size_t read_len,
                                              int max_error,
                                              size_t min_len,
                                              size_t max_len) {
    enum { START_CAP = 8 };
    AnchorCandidate best = {0};
    AnchorMatch hits[START_CAP];
    size_t n = anchor_collect_matches(start_pat, read_seq, read_len, max_error, 0, hits, START_CAP);
    int fixed_len_mode = (min_len == max_len);
    size_t fixed_len = min_len;

    for (size_t i = 0; i < n; ++i) {
        if (hits[i].end > read_len) continue;
        size_t payload_len = 0;
        size_t insert_start = hits[i].end;
        if (fixed_len_mode) {
            if (insert_start + fixed_len > read_len) continue;
            payload_len = fixed_len;
        } else {
            payload_len = read_len - insert_start;
            if (payload_len < min_len || payload_len > max_len) continue;
        }

        AnchorCandidate cand = {0};
        cand.valid = 1;
        cand.errors = hits[i].errors;
        cand.first_start = hits[i].start;
        cand.insert_start = insert_start;
        cand.insert_len = payload_len;
        candidate_consider(&best, &cand);

        // For one-sided 5' mode, the first valid 0-error hit is optimal.
        if (best.valid && !best.ambiguous && best.errors == 0) {
            return best;
        }
    }
    return best;
}

static AnchorCandidate anchor_find_best_3only(const BitapPattern *end_pat,
                                              const char *read_seq,
                                              size_t read_len,
                                              int max_error,
                                              size_t min_len,
                                              size_t max_len) {
    enum { END_CAP = 8 };
    AnchorCandidate best = {0};
    AnchorMatch hits[END_CAP];
    size_t n = anchor_collect_matches(end_pat, read_seq, read_len, max_error, 0, hits, END_CAP);
    int fixed_len_mode = (min_len == max_len);
    size_t fixed_len = min_len;

    for (size_t i = 0; i < n; ++i) {
        size_t payload_len = 0;
        size_t insert_start = 0;
        if (fixed_len_mode) {
            if (hits[i].start < fixed_len) continue;
            payload_len = fixed_len;
            insert_start = hits[i].start - fixed_len;
        } else {
            payload_len = hits[i].start;
            if (payload_len < min_len || payload_len > max_len) continue;
            insert_start = 0;
        }

        AnchorCandidate cand = {0};
        cand.valid = 1;
        cand.errors = hits[i].errors;
        cand.first_start = hits[i].start;
        cand.insert_start = insert_start;
        cand.insert_len = payload_len;
        candidate_consider(&best, &cand);

        // For one-sided 3' mode, the first valid 0-error hit is optimal.
        if (best.valid && !best.ambiguous && best.errors == 0) {
            return best;
        }
    }
    return best;
}

int anchor_runtime_init(AnchorRuntime *out,
                        const AnchorConfig *cfg,
                        size_t expected_insert_len) {
    return anchor_runtime_init_range(out, cfg, expected_insert_len, expected_insert_len);
}

int anchor_runtime_init_range(AnchorRuntime *out,
                              const AnchorConfig *cfg,
                              size_t min_insert_len,
                              size_t max_insert_len) {
    if (!out) return -1;
    memset(out, 0, sizeof(*out));
    if (!cfg || !cfg->enabled) return 0;
    int has5 = (cfg->anchor5 && cfg->anchor5[0] != '\0');
    int has3 = (cfg->anchor3 && cfg->anchor3[0] != '\0');
    if (!has5 && !has3) return -1;
    if (min_insert_len > max_insert_len) return -1;

    out->enabled = 1;
    out->max_error = (cfg->max_error < 0) ? 0 : cfg->max_error;
    out->has_anchor5 = has5;
    out->has_anchor3 = has3;
    out->expected_insert_len = min_insert_len;
    out->min_insert_len = min_insert_len;
    out->max_insert_len = max_insert_len;

    char *a5_rc = NULL;
    char *a3_rc = NULL;
    if (has5) {
        a5_rc = anchor_revcomp_new_n(cfg->anchor5, strlen(cfg->anchor5));
        if (!a5_rc) return -1;
    }
    if (has3) {
        a3_rc = anchor_revcomp_new_n(cfg->anchor3, strlen(cfg->anchor3));
        if (!a3_rc) {
            free(a5_rc);
            return -1;
        }
    }

    int rc = 0;
    if (has5) {
        rc |= bitap_pattern_build(&out->a5, cfg->anchor5);
        rc |= bitap_pattern_build(&out->a5_rc, a5_rc);
    }
    if (has3) {
        rc |= bitap_pattern_build(&out->a3, cfg->anchor3);
        rc |= bitap_pattern_build(&out->a3_rc, a3_rc);
    }
    free(a5_rc);
    free(a3_rc);
    if (rc != 0) return -1;

    return 0;
}

int anchor_extract_window(const AnchorRuntime *ar,
                          const char *read_seq,
                          size_t read_len,
                          size_t *insert_start_out,
                          size_t *insert_len_out) {
    return anchor_extract_window_range(ar, read_seq, read_len, insert_start_out, insert_len_out);
}

int anchor_extract_window_range(const AnchorRuntime *ar,
                                const char *read_seq,
                                size_t read_len,
                                size_t *insert_start_out,
                                size_t *insert_len_out) {
    if (!ar || !ar->enabled || !read_seq || !insert_start_out || !insert_len_out) return -1;
    size_t min_len = ar->min_insert_len;
    size_t max_len = ar->max_insert_len;
    if (min_len > max_len) return -1;

    // Fast path: try exact matching first for all error budgets (including --anchor-error 0).
    if (ar->has_anchor5 && ar->has_anchor3) {
        AnchorCandidate best = {0};
        AnchorCandidate fw = anchor_find_best_pair_exact(&ar->a5, &ar->a3, read_seq, read_len, min_len, max_len);
        AnchorCandidate rv = anchor_find_best_pair_exact(&ar->a3_rc, &ar->a5_rc, read_seq, read_len, min_len, max_len);
        candidate_consider(&best, &fw);
        candidate_consider(&best, &rv);
        if (best.valid) {
            if (best.ambiguous) return 1;
            *insert_start_out = best.insert_start;
            *insert_len_out = best.insert_len;
            return 0;
        }
    } else if (ar->has_anchor5) {
        AnchorCandidate best = anchor_find_best_5only_exact(&ar->a5, read_seq, read_len, min_len, max_len);
        if (best.valid) {
            if (best.ambiguous) return 1;
            *insert_start_out = best.insert_start;
            *insert_len_out = best.insert_len;
            return 0;
        }
    } else if (ar->has_anchor3) {
        AnchorCandidate best = anchor_find_best_3only_exact(&ar->a3, read_seq, read_len, min_len, max_len);
        if (best.valid) {
            if (best.ambiguous) return 1;
            *insert_start_out = best.insert_start;
            *insert_len_out = best.insert_len;
            return 0;
        }
    }

    if (ar->has_anchor5 && ar->has_anchor3) {
        AnchorCandidate best = {0};
        AnchorCandidate fw = anchor_find_best_pair(&ar->a5, &ar->a3,
                                                   read_seq, read_len,
                                                   ar->max_error, min_len, max_len);
        AnchorCandidate rv = anchor_find_best_pair(&ar->a3_rc, &ar->a5_rc,
                                                   read_seq, read_len,
                                                   ar->max_error, min_len, max_len);

        candidate_consider(&best, &fw);
        candidate_consider(&best, &rv);

        if (!best.valid || best.ambiguous) return 1;
        *insert_start_out = best.insert_start;
        *insert_len_out = best.insert_len;
        return 0;
    }

    if (ar->has_anchor5) {
        AnchorCandidate best = anchor_find_best_5only(&ar->a5, read_seq, read_len,
                                                      ar->max_error, min_len, max_len);
        if (!best.valid || best.ambiguous) return 1;
        *insert_start_out = best.insert_start;
        *insert_len_out = best.insert_len;
        return 0;
    }

    if (ar->has_anchor3) {
        AnchorCandidate best = anchor_find_best_3only(&ar->a3, read_seq, read_len,
                                                      ar->max_error, min_len, max_len);
        if (!best.valid || best.ambiguous) return 1;
        *insert_start_out = best.insert_start;
        *insert_len_out = best.insert_len;
        return 0;
    }

    return 1;
}
