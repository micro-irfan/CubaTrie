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

static size_t usize_abs_diff(size_t a, size_t b) {
    return (a >= b) ? (a - b) : (b - a);
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

    int found = 0;
    int best_err = 0;
    size_t best_start = 0;
    size_t best_end = 0;
    size_t best_span_diff = 0;

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
            size_t cand_start = from_pos + curr_start[m];
            size_t cand_end = from_pos + col;
            size_t cand_len = (cand_end > cand_start) ? (cand_end - cand_start) : 0;
            size_t cand_span_diff = usize_abs_diff(cand_len, m);
            int cand_err = curr[m];

            if (!found ||
                cand_start < best_start ||
                (cand_start == best_start && cand_err < best_err) ||
                (cand_start == best_start && cand_err == best_err && cand_span_diff < best_span_diff) ||
                (cand_start == best_start && cand_err == best_err && cand_span_diff == best_span_diff && cand_end < best_end)) {
                found = 1;
                best_err = cand_err;
                best_start = cand_start;
                best_end = cand_end;
                best_span_diff = cand_span_diff;
            }
        }

        int *tmp_i = prev;
        prev = curr;
        curr = tmp_i;
        size_t *tmp_s = prev_start;
        prev_start = curr_start;
        curr_start = tmp_s;
    }

    if (!found) return 1;
    out_hit->start = best_start;
    out_hit->end = best_end;
    out_hit->errors = best_err;
    return 0;
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

static int bitap_find_first_exact(const BitapPattern *pat,
                                  const char *text,
                                  size_t text_len,
                                  size_t from_pos,
                                  AnchorMatch *out_hit) {
    if (!pat || !text || !out_hit) return -1;
    if (from_pos >= text_len || pat->len == 0) return 1;
    if (pat->len > 63) return -1;

    size_t m = pat->len;
    uint64_t full_mask = ((uint64_t)1 << m) - 1ULL;
    uint64_t accept_bit = (uint64_t)1 << (m - 1);
    uint64_t state = 0ULL;

    for (size_t i = from_pos; i < text_len; ++i) {
        unsigned char c = (unsigned char)toupper((unsigned char)text[i]);
        uint64_t eq = (~pat->mask[c]) & full_mask;
        state = ((state << 1) | 1ULL) & eq;
        if (state & accept_bit) {
            size_t end = i + 1;
            out_hit->start = end - m;
            out_hit->end = end;
            out_hit->errors = 0;
            return 0;
        }
    }
    return 1;
}

static size_t anchor_collect_matches_exact(const BitapPattern *pat,
                                           const char *read_seq,
                                           size_t read_len,
                                           size_t from_pos,
                                           AnchorMatch *hits,
                                           size_t hits_cap) {
    if (!pat || !read_seq || !hits || hits_cap == 0 || from_pos >= read_len) return 0;
    size_t n = 0;
    size_t cursor = from_pos;
    while (n < hits_cap && cursor < read_len) {
        AnchorMatch h = {0};
        int rc = bitap_find_first_exact(pat, read_seq, read_len, cursor, &h);
        if (rc != 0) break;
        hits[n++] = h;
        cursor = (h.start + 1 > cursor) ? (h.start + 1) : (cursor + 1);
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
    AnchorOrientation orientation;
    int has_anchor5;
    size_t anchor5_start;
    size_t anchor5_end;
    int anchor5_errors;
    int has_anchor3;
    size_t anchor3_start;
    size_t anchor3_end;
    int anchor3_errors;
} AnchorCandidate;

static int candidate_score_cmp(const AnchorCandidate *a, const AnchorCandidate *b) {
    if (a->errors != b->errors) return (a->errors < b->errors) ? -1 : 1;
    if (a->first_start != b->first_start) return (a->first_start < b->first_start) ? -1 : 1;
    if (a->insert_len != b->insert_len) return (a->insert_len < b->insert_len) ? -1 : 1;
    if (a->insert_start != b->insert_start) return (a->insert_start < b->insert_start) ? -1 : 1;
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
                                AnchorOrientation orientation,
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
    out->orientation = orientation;
    out->has_anchor5 = 1;
    out->has_anchor3 = 1;
    if (orientation == ANCHOR_ORIENT_RC) {
        out->anchor3_start = start_hit->start;
        out->anchor3_end = start_hit->end;
        out->anchor3_errors = start_hit->errors;
        out->anchor5_start = end_hit->start;
        out->anchor5_end = end_hit->end;
        out->anchor5_errors = end_hit->errors;
    } else {
        out->anchor5_start = start_hit->start;
        out->anchor5_end = start_hit->end;
        out->anchor5_errors = start_hit->errors;
        out->anchor3_start = end_hit->start;
        out->anchor3_end = end_hit->end;
        out->anchor3_errors = end_hit->errors;
    }
}

static AnchorCandidate anchor_find_best_pair_exact(const BitapPattern *start_pat,
                                                   const BitapPattern *end_pat,
                                                   AnchorOrientation orientation,
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
            candidate_from_pair(&s, &end_hits[j], orientation, min_len, max_len, &cand);
            candidate_consider(&best, &cand);
        }

        if (best.valid && !best.ambiguous && best.errors == 0 && best.first_start == s.start) {
            return best;
        }
    }
    return best;
}

static AnchorCandidate anchor_find_best_5only_exact(const BitapPattern *start_pat,
                                                    AnchorOrientation orientation,
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
        cand.orientation = orientation;
        cand.has_anchor5 = 1;
        cand.anchor5_start = hits[i].start;
        cand.anchor5_end = hits[i].end;
        cand.anchor5_errors = hits[i].errors;
        candidate_consider(&best, &cand);
        if (best.valid && !best.ambiguous) return best;
    }
    return best;
}

static AnchorCandidate anchor_find_best_3only_exact(const BitapPattern *end_pat,
                                                    AnchorOrientation orientation,
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
        cand.orientation = orientation;
        cand.has_anchor3 = 1;
        cand.anchor3_start = hits[i].start;
        cand.anchor3_end = hits[i].end;
        cand.anchor3_errors = hits[i].errors;
        candidate_consider(&best, &cand);
        if (best.valid && !best.ambiguous) return best;
    }
    return best;
}

static AnchorCandidate anchor_find_best_pair(const BitapPattern *start_pat,
                                             const BitapPattern *end_pat,
                                             AnchorOrientation orientation,
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
            candidate_from_pair(&s, &end_hits[j], orientation, min_len, max_len, &cand);
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
                                              AnchorOrientation orientation,
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
        cand.orientation = orientation;
        cand.has_anchor5 = 1;
        cand.anchor5_start = hits[i].start;
        cand.anchor5_end = hits[i].end;
        cand.anchor5_errors = hits[i].errors;
        candidate_consider(&best, &cand);

        // For one-sided 5' mode, the first valid 0-error hit is optimal.
        if (best.valid && !best.ambiguous && best.errors == 0) {
            return best;
        }
    }
    return best;
}

static AnchorCandidate anchor_find_best_3only(const BitapPattern *end_pat,
                                              AnchorOrientation orientation,
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
        cand.orientation = orientation;
        cand.has_anchor3 = 1;
        cand.anchor3_start = hits[i].start;
        cand.anchor3_end = hits[i].end;
        cand.anchor3_errors = hits[i].errors;
        candidate_consider(&best, &cand);

        // For one-sided 3' mode, the first valid 0-error hit is optimal.
        if (best.valid && !best.ambiguous && best.errors == 0) {
            return best;
        }
    }
    return best;
}

static AnchorCandidate anchor_find_best_start_only_exact(const BitapPattern *start_pat,
                                                         AnchorOrientation orientation,
                                                         int is_a5_anchor,
                                                         const char *read_seq,
                                                         size_t read_len) {
    enum { START_CAP = 8 };
    AnchorCandidate best = {0};
    AnchorMatch hits[START_CAP];
    size_t n = anchor_collect_matches_exact(start_pat, read_seq, read_len, 0, hits, START_CAP);
    for (size_t i = 0; i < n; ++i) {
        AnchorCandidate cand = {0};
        cand.valid = 1;
        cand.errors = hits[i].errors;
        cand.first_start = hits[i].start;
        cand.insert_start = hits[i].start;
        cand.insert_len = start_pat->len;
        cand.orientation = orientation;
        if (is_a5_anchor) {
            cand.has_anchor5 = 1;
            cand.anchor5_start = hits[i].start;
            cand.anchor5_end = hits[i].end;
            cand.anchor5_errors = hits[i].errors;
        } else {
            cand.has_anchor3 = 1;
            cand.anchor3_start = hits[i].start;
            cand.anchor3_end = hits[i].end;
            cand.anchor3_errors = hits[i].errors;
        }
        candidate_consider(&best, &cand);
    }
    return best;
}

static AnchorCandidate anchor_find_best_start_only(const BitapPattern *start_pat,
                                                   AnchorOrientation orientation,
                                                   int is_a5_anchor,
                                                   const char *read_seq,
                                                   size_t read_len,
                                                   int max_error) {
    enum { START_CAP = 8 };
    AnchorCandidate best = {0};
    AnchorMatch hits[START_CAP];
    size_t n = anchor_collect_matches(start_pat, read_seq, read_len, max_error, 0, hits, START_CAP);
    for (size_t i = 0; i < n; ++i) {
        AnchorCandidate cand = {0};
        cand.valid = 1;
        cand.errors = hits[i].errors;
        cand.first_start = hits[i].start;
        cand.insert_start = hits[i].start;
        cand.insert_len = start_pat->len;
        cand.orientation = orientation;
        if (is_a5_anchor) {
            cand.has_anchor5 = 1;
            cand.anchor5_start = hits[i].start;
            cand.anchor5_end = hits[i].end;
            cand.anchor5_errors = hits[i].errors;
        } else {
            cand.has_anchor3 = 1;
            cand.anchor3_start = hits[i].start;
            cand.anchor3_end = hits[i].end;
            cand.anchor3_errors = hits[i].errors;
        }
        candidate_consider(&best, &cand);
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
    return anchor_extract_window_range_info(ar, read_seq, read_len, insert_start_out, insert_len_out, NULL);
}

int anchor_extract_window_range(const AnchorRuntime *ar,
                                const char *read_seq,
                                size_t read_len,
                                size_t *insert_start_out,
                                size_t *insert_len_out) {
    return anchor_extract_window_range_info(ar, read_seq, read_len, insert_start_out, insert_len_out, NULL);
}

static void anchor_window_info_fill(const AnchorCandidate *best,
                                    AnchorWindowInfo *info_out) {
    if (!best || !info_out) return;
    memset(info_out, 0, sizeof(*info_out));
    info_out->orientation = best->orientation;
    info_out->insert_start = best->insert_start;
    info_out->insert_len = best->insert_len;
    info_out->has_anchor5 = best->has_anchor5;
    info_out->anchor5_start = best->anchor5_start;
    info_out->anchor5_end = best->anchor5_end;
    info_out->anchor5_errors = best->anchor5_errors;
    info_out->has_anchor3 = best->has_anchor3;
    info_out->anchor3_start = best->anchor3_start;
    info_out->anchor3_end = best->anchor3_end;
    info_out->anchor3_errors = best->anchor3_errors;
}

int anchor_extract_window_range_info(const AnchorRuntime *ar,
                                     const char *read_seq,
                                     size_t read_len,
                                     size_t *insert_start_out,
                                     size_t *insert_len_out,
                                     AnchorWindowInfo *info_out) {
    if (!ar || !ar->enabled || !read_seq || !insert_start_out || !insert_len_out) return -1;
    size_t min_len = ar->min_insert_len;
    size_t max_len = ar->max_insert_len;
    if (min_len > max_len) return -1;

    // Fast path: try exact matching first for all error budgets (including --anchor-error 0).
    if (ar->has_anchor5 && ar->has_anchor3) {
        AnchorCandidate best = {0};
        AnchorCandidate fw = anchor_find_best_pair_exact(&ar->a5, &ar->a3, ANCHOR_ORIENT_FWD,
                                                         read_seq, read_len, min_len, max_len);
        AnchorCandidate rv = anchor_find_best_pair_exact(&ar->a3_rc, &ar->a5_rc, ANCHOR_ORIENT_RC,
                                                         read_seq, read_len, min_len, max_len);
        candidate_consider(&best, &fw);
        candidate_consider(&best, &rv);
        if (best.valid) {
            if (best.ambiguous) return 1;
            *insert_start_out = best.insert_start;
            *insert_len_out = best.insert_len;
            anchor_window_info_fill(&best, info_out);
            return 0;
        }
    } else if (ar->has_anchor5) {
        AnchorCandidate best = anchor_find_best_5only_exact(&ar->a5, ANCHOR_ORIENT_FWD,
                                                             read_seq, read_len, min_len, max_len);
        if (best.valid) {
            if (best.ambiguous) return 1;
            *insert_start_out = best.insert_start;
            *insert_len_out = best.insert_len;
            anchor_window_info_fill(&best, info_out);
            return 0;
        }
    } else if (ar->has_anchor3) {
        AnchorCandidate best = anchor_find_best_3only_exact(&ar->a3, ANCHOR_ORIENT_FWD,
                                                             read_seq, read_len, min_len, max_len);
        if (best.valid) {
            if (best.ambiguous) return 1;
            *insert_start_out = best.insert_start;
            *insert_len_out = best.insert_len;
            anchor_window_info_fill(&best, info_out);
            return 0;
        }
    }

    if (ar->has_anchor5 && ar->has_anchor3) {
        AnchorCandidate best = {0};
        AnchorCandidate fw = anchor_find_best_pair(&ar->a5, &ar->a3,
                                                   ANCHOR_ORIENT_FWD,
                                                   read_seq, read_len,
                                                   ar->max_error, min_len, max_len);
        AnchorCandidate rv = anchor_find_best_pair(&ar->a3_rc, &ar->a5_rc,
                                                   ANCHOR_ORIENT_RC,
                                                   read_seq, read_len,
                                                   ar->max_error, min_len, max_len);

        candidate_consider(&best, &fw);
        candidate_consider(&best, &rv);

        if (!best.valid || best.ambiguous) return 1;
        *insert_start_out = best.insert_start;
        *insert_len_out = best.insert_len;
        anchor_window_info_fill(&best, info_out);
        return 0;
    }

    if (ar->has_anchor5) {
        AnchorCandidate best = anchor_find_best_5only(&ar->a5, ANCHOR_ORIENT_FWD,
                                                      read_seq, read_len,
                                                      ar->max_error, min_len, max_len);
        if (!best.valid || best.ambiguous) return 1;
        *insert_start_out = best.insert_start;
        *insert_len_out = best.insert_len;
        anchor_window_info_fill(&best, info_out);
        return 0;
    }

    if (ar->has_anchor3) {
        AnchorCandidate best = anchor_find_best_3only(&ar->a3, ANCHOR_ORIENT_FWD,
                                                      read_seq, read_len,
                                                      ar->max_error, min_len, max_len);
        if (!best.valid || best.ambiguous) return 1;
        *insert_start_out = best.insert_start;
        *insert_len_out = best.insert_len;
        anchor_window_info_fill(&best, info_out);
        return 0;
    }

    return 1;
}

int anchor_extract_two_sided_partial_start_info(const AnchorRuntime *ar,
                                                const char *read_seq,
                                                size_t read_len,
                                                AnchorWindowInfo *info_out) {
    if (!ar || !ar->enabled || !read_seq || !info_out) return -1;
    if (!ar->has_anchor5 || !ar->has_anchor3) return -1;

    AnchorCandidate best = {0};
    AnchorCandidate fw = anchor_find_best_start_only_exact(&ar->a5, ANCHOR_ORIENT_FWD, 1,
                                                           read_seq, read_len);
    AnchorCandidate rv = anchor_find_best_start_only_exact(&ar->a3_rc, ANCHOR_ORIENT_RC, 0,
                                                           read_seq, read_len);
    candidate_consider(&best, &fw);
    candidate_consider(&best, &rv);
    if (best.valid) {
        if (best.ambiguous) return 1;
        anchor_window_info_fill(&best, info_out);
        return 0;
    }

    best.valid = 0;
    best.ambiguous = 0;
    fw = anchor_find_best_start_only(&ar->a5, ANCHOR_ORIENT_FWD, 1, read_seq, read_len, ar->max_error);
    rv = anchor_find_best_start_only(&ar->a3_rc, ANCHOR_ORIENT_RC, 0, read_seq, read_len, ar->max_error);
    candidate_consider(&best, &fw);
    candidate_consider(&best, &rv);
    if (!best.valid || best.ambiguous) return 1;
    anchor_window_info_fill(&best, info_out);
    return 0;
}

int anchor_extract_two_sided_best_pair_info(const AnchorRuntime *ar,
                                            const char *read_seq,
                                            size_t read_len,
                                            AnchorWindowInfo *info_out,
                                            int *ambiguous_out) {
    if (!ar || !ar->enabled || !read_seq) return -1;
    if (!ar->has_anchor5 || !ar->has_anchor3) return -1;

    if (ambiguous_out) *ambiguous_out = 0;
    if (info_out) memset(info_out, 0, sizeof(*info_out));

    // Diagnostic search intentionally ignores insert-length constraints.
    // This lets callers explain why constrained extraction failed.
    AnchorCandidate best = {0};
    AnchorCandidate fw = anchor_find_best_pair_exact(&ar->a5, &ar->a3,
                                                     ANCHOR_ORIENT_FWD,
                                                     read_seq, read_len,
                                                     1, (size_t)-1);
    AnchorCandidate rv = anchor_find_best_pair_exact(&ar->a3_rc, &ar->a5_rc,
                                                     ANCHOR_ORIENT_RC,
                                                     read_seq, read_len,
                                                     1, (size_t)-1);
    candidate_consider(&best, &fw);
    candidate_consider(&best, &rv);

    if (!best.valid) {
        best.valid = 0;
        best.ambiguous = 0;
        fw = anchor_find_best_pair(&ar->a5, &ar->a3,
                                   ANCHOR_ORIENT_FWD,
                                   read_seq, read_len,
                                   ar->max_error, 1, (size_t)-1);
        rv = anchor_find_best_pair(&ar->a3_rc, &ar->a5_rc,
                                   ANCHOR_ORIENT_RC,
                                   read_seq, read_len,
                                   ar->max_error, 1, (size_t)-1);
        candidate_consider(&best, &fw);
        candidate_consider(&best, &rv);
    }

    if (!best.valid) return 1;
    if (ambiguous_out) *ambiguous_out = best.ambiguous ? 1 : 0;
    if (info_out) anchor_window_info_fill(&best, info_out);
    return 0;
}
