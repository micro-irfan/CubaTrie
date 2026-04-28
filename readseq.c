#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "kseq.h"         
KSEQ_INIT(gzFile, gzread)  

#include "readseq_internal.h"

const size_t STEP = 1000000;

static void sam_write_unmapped(FILE *sam_fp,
                               const char *read_name,
                               const char *read_seq,
                               const char *read_qual,
                               size_t read_len,
                               const char *extra_tag) {
    if (!sam_fp || !read_name || !read_seq) return;
    const char *qual = (read_qual && strlen(read_qual) == read_len) ? read_qual : "*";
    size_t qname_len = name_len_no_rc_suffix(read_name);
    if (extra_tag && extra_tag[0] != '\0') {
        fprintf(sam_fp, "%.*s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s%s\n",
                (int)qname_len, read_name, read_seq, qual, extra_tag);
    } else {
        fprintf(sam_fp, "%.*s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
                (int)qname_len, read_name, read_seq, qual);
    }
}

void total_hits (const kh_counter_t *m, const size_t nreads) {
    size_t total = 0;
    for (khint_t i = kh_begin(m); i != kh_end(m); ++i) 
        if (kh_exist(m, i)) total += kh_val(m, i);

    fprintf(stderr, "Finished Processing %zu reads...\n", nreads);

    double pct = nreads ? 100.0 * (double)total / (double)nreads : 0.0;
    fprintf(stderr, "Reads with matches: %zu (%.2f%%)\n", total, pct);
}

static _Thread_local khash_t(strset) *tls_matches = NULL;

typedef struct {
    char op;
    char ref_base;
    char read_base;
} AnchorEditOp;

static int anchor_runtime_init_for_count_mode(AnchorRuntime *out,
                                              const AnchorConfig *cfg,
                                              size_t ref_len) {
    if (!out) return -1;
    if (!cfg || !cfg->enabled) return anchor_runtime_init(out, cfg, ref_len);
    int has5 = (cfg->anchor5 && cfg->anchor5[0] != '\0');
    int has3 = (cfg->anchor3 && cfg->anchor3[0] != '\0');
    if (has5 && has3) {
        // In count mode, paired anchors define the insert boundaries directly.
        // Keep one-sided anchor behavior unchanged (fixed reference length).
        return anchor_runtime_init_range(out, cfg, 1, (size_t)-1);
    }
    return anchor_runtime_init(out, cfg, ref_len);
}

static void strset_clear_with_keys(khash_t(strset) *set) {
    if (!set) return;
    for (khint_t i = kh_begin(set); i != kh_end(set); ++i) {
        if (kh_exist(set, i)) free((char*)kh_key(set, i));
    }
    kh_clear(strset, set);
}

static khash_t(strset) *tls_get_matches(void) {
    if (!tls_matches) tls_matches = kh_init(strset);
    return tls_matches;
}

static int anchor_tag_buf_append_char(char *buf, size_t cap, size_t *len, char c) {
    if (!buf || !len || *len + 1 >= cap) return -1;
    buf[*len] = c;
    *len += 1;
    buf[*len] = '\0';
    return 0;
}

static int anchor_tag_buf_append_cstr(char *buf, size_t cap, size_t *len, const char *s) {
    if (!buf || !len || !s) return -1;
    while (*s) {
        if (anchor_tag_buf_append_char(buf, cap, len, *s++) != 0) return -1;
    }
    return 0;
}

static int anchor_tag_buf_append_size(char *buf, size_t cap, size_t *len, size_t v) {
    char tmp[32];
    size_t n = 0;
    do {
        tmp[n++] = (char)('0' + (v % 10));
        v /= 10;
    } while (v > 0);
    if (!buf || !len || *len + n >= cap) return -1;
    for (size_t i = 0; i < n; ++i) {
        buf[*len + i] = tmp[n - 1 - i];
    }
    *len += n;
    buf[*len] = '\0';
    return 0;
}

static int anchor_compute_md(const char *ref_seq,
                             size_t ref_len,
                             const char *read_seq,
                             size_t read_len,
                             char *md_out,
                             size_t md_cap) {
    if (!ref_seq || !read_seq || !md_out || md_cap == 0) return -1;
    if (ref_len > 63 || read_len > 63) return -1;

    int dp[64][64];
    char bt[64][64];
    for (size_t i = 0; i <= ref_len; ++i) {
        dp[i][0] = (int)i;
        bt[i][0] = 'D';
    }
    for (size_t j = 0; j <= read_len; ++j) {
        dp[0][j] = (int)j;
        bt[0][j] = 'I';
    }
    bt[0][0] = 'M';

    for (size_t i = 1; i <= ref_len; ++i) {
        for (size_t j = 1; j <= read_len; ++j) {
            int sub_cost = (ref_seq[i - 1] == read_seq[j - 1]) ? 0 : 1;
            int best = dp[i - 1][j - 1] + sub_cost;
            char op = 'M';

            int ins = dp[i][j - 1] + 1;
            if (ins < best) {
                best = ins;
                op = 'I';
            }
            int del = dp[i - 1][j] + 1;
            if (del < best) {
                best = del;
                op = 'D';
            }
            dp[i][j] = best;
            bt[i][j] = op;
        }
    }

    AnchorEditOp rev_ops[128];
    size_t nops = 0;
    size_t i = ref_len;
    size_t j = read_len;
    while (i > 0 || j > 0) {
        char op = bt[i][j];
        AnchorEditOp e = {0};
        if (i > 0 && j > 0 && op == 'M') {
            e.op = 'M';
            e.ref_base = ref_seq[i - 1];
            e.read_base = read_seq[j - 1];
            --i;
            --j;
        } else if (j > 0 && (i == 0 || op == 'I')) {
            e.op = 'I';
            e.read_base = read_seq[j - 1];
            --j;
        } else if (i > 0) {
            e.op = 'D';
            e.ref_base = ref_seq[i - 1];
            --i;
        } else {
            return -1;
        }
        if (nops >= sizeof(rev_ops) / sizeof(rev_ops[0])) return -1;
        rev_ops[nops++] = e;
    }

    AnchorEditOp ops[128];
    for (size_t k = 0; k < nops; ++k) ops[k] = rev_ops[nops - 1 - k];

    size_t md_len = 0;
    size_t run = 0;
    int in_del = 0;
    md_out[0] = '\0';
    for (size_t k = 0; k < nops; ++k) {
        AnchorEditOp e = ops[k];
        if (e.op == 'M') {
            if (e.ref_base == e.read_base) {
                run++;
            } else {
                if (anchor_tag_buf_append_size(md_out, md_cap, &md_len, run) != 0) return -1;
                if (anchor_tag_buf_append_char(md_out, md_cap, &md_len, e.ref_base) != 0) return -1;
                run = 0;
            }
            in_del = 0;
            continue;
        }
        if (e.op == 'D') {
            if (!in_del) {
                if (anchor_tag_buf_append_size(md_out, md_cap, &md_len, run) != 0) return -1;
                run = 0;
                if (anchor_tag_buf_append_char(md_out, md_cap, &md_len, '^') != 0) return -1;
                in_del = 1;
            }
            if (anchor_tag_buf_append_char(md_out, md_cap, &md_len, e.ref_base) != 0) return -1;
            continue;
        }
        in_del = 0;
    }
    if (anchor_tag_buf_append_size(md_out, md_cap, &md_len, run) != 0) return -1;
    return 0;
}

static int anchor_build_sam_tag(const AnchorRuntime *ar,
                                const char *read_seq,
                                size_t read_len,
                                const AnchorWindowInfo *info,
                                char *tag_out,
                                size_t tag_cap) {
    if (!ar || !read_seq || !info || !tag_out || tag_cap == 0) return -1;
    if (info->insert_start + info->insert_len > read_len) return -1;

    size_t len = 0;
    tag_out[0] = '\0';
    if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, "\tZA:Z:ori=") != 0) return -1;
    if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len,
                                   info->orientation == ANCHOR_ORIENT_RC ? "RC" : "FWD") != 0) return -1;
    if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ";ins=") != 0) return -1;
    if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->insert_start + 1) != 0) return -1;
    if (anchor_tag_buf_append_char(tag_out, tag_cap, &len, ',') != 0) return -1;
    if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->insert_len) != 0) return -1;

    if (info->has_anchor5) {
        if (info->anchor5_end < info->anchor5_start) return -1;
        size_t seg_len = info->anchor5_end - info->anchor5_start;
        char md5[256];
        const char *a5_ref = (info->orientation == ANCHOR_ORIENT_RC) ? ar->a5_rc.seq : ar->a5.seq;
        if (info->anchor5_start >= read_len || info->anchor5_end > read_len) return -1;
        if (anchor_compute_md(a5_ref, strlen(a5_ref), read_seq + info->anchor5_start, seg_len, md5, sizeof(md5)) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ";a5=") != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->anchor5_start + 1) != 0) return -1;
        if (anchor_tag_buf_append_char(tag_out, tag_cap, &len, '-') != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->anchor5_end) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ",ed=") != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, (size_t)((info->anchor5_errors < 0) ? 0 : info->anchor5_errors)) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ",md=") != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, md5) != 0) return -1;
    }

    if (info->has_anchor3) {
        if (info->anchor3_end < info->anchor3_start) return -1;
        size_t seg_len = info->anchor3_end - info->anchor3_start;
        char md3[256];
        const char *a3_ref = (info->orientation == ANCHOR_ORIENT_RC) ? ar->a3_rc.seq : ar->a3.seq;
        if (info->anchor3_start >= read_len || info->anchor3_end > read_len) return -1;
        if (anchor_compute_md(a3_ref, strlen(a3_ref), read_seq + info->anchor3_start, seg_len, md3, sizeof(md3)) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ";a3=") != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->anchor3_start + 1) != 0) return -1;
        if (anchor_tag_buf_append_char(tag_out, tag_cap, &len, '-') != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->anchor3_end) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ",ed=") != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, (size_t)((info->anchor3_errors < 0) ? 0 : info->anchor3_errors)) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ",md=") != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, md3) != 0) return -1;
    }
    return 0;
}

static int anchor_build_partial_two_sided_start_sam_tag(const AnchorRuntime *ar,
                                                        const char *read_seq,
                                                        size_t read_len,
                                                        const AnchorWindowInfo *info,
                                                        char *tag_out,
                                                        size_t tag_cap) {
    if (!ar || !read_seq || !info || !tag_out || tag_cap == 0) return -1;
    size_t len = 0;
    tag_out[0] = '\0';
    if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, "\tZA:Z:ori=") != 0) return -1;
    if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len,
                                   info->orientation == ANCHOR_ORIENT_RC ? "RC" : "FWD") != 0) return -1;
    if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ";partial=1;") != 0) return -1;

    if (info->orientation == ANCHOR_ORIENT_FWD && info->has_anchor5) {
        if (info->anchor5_end < info->anchor5_start || info->anchor5_end > read_len) return -1;
        size_t seg_len = info->anchor5_end - info->anchor5_start;
        char md[256];
        if (anchor_compute_md(ar->a5.seq, strlen(ar->a5.seq),
                              read_seq + info->anchor5_start, seg_len,
                              md, sizeof(md)) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, "a5=") != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->anchor5_start + 1) != 0) return -1;
        if (anchor_tag_buf_append_char(tag_out, tag_cap, &len, '-') != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->anchor5_end) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ",ed=") != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, (size_t)((info->anchor5_errors < 0) ? 0 : info->anchor5_errors)) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ",md=") != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, md) != 0) return -1;
        return 0;
    }

    if (info->orientation == ANCHOR_ORIENT_RC && info->has_anchor3) {
        if (info->anchor3_end < info->anchor3_start || info->anchor3_end > read_len) return -1;
        size_t seg_len = info->anchor3_end - info->anchor3_start;
        char md[256];
        if (anchor_compute_md(ar->a3_rc.seq, strlen(ar->a3_rc.seq),
                              read_seq + info->anchor3_start, seg_len,
                              md, sizeof(md)) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, "a3rc=") != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->anchor3_start + 1) != 0) return -1;
        if (anchor_tag_buf_append_char(tag_out, tag_cap, &len, '-') != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, info->anchor3_end) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ",ed=") != 0) return -1;
        if (anchor_tag_buf_append_size(tag_out, tag_cap, &len, (size_t)((info->anchor3_errors < 0) ? 0 : info->anchor3_errors)) != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, ",md=") != 0) return -1;
        if (anchor_tag_buf_append_cstr(tag_out, tag_cap, &len, md) != 0) return -1;
        return 0;
    }

    return -1;
}

static int add_matches_to_counter_hits(khash_t(strset) *found_sequences,
                                       kh_counter_t *map) {
    for (khint_t i = kh_begin(found_sequences); i != kh_end(found_sequences); ++i) {
        if (!kh_exist(found_sequences, i)) continue;
        if (counter_add_with_init(map, kh_key(found_sequences, i), 1, 1) != 0) return -1;
    }
    return 0;
}

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
                             int counts_preseeded) {
    int status = 0;
    FILE *sam_out = sam_fp;
    kvec_t(uint32_t) hits;
    kv_init(hits);
    khash_t(strset) *matches = NULL;
    khash_t(strset) *sam_matches = NULL;
    const char *search_seq = read_seq;
    const char *search_qual = read_qual;
    size_t search_len = read_len;
    size_t anchor_insert_start = 0;
    size_t anchor_insert_len = 0;
    AnchorWindowInfo anchor_info = {0};
    int have_anchor_info = 0;
    char anchor_sam_tag[1024];
    int have_anchor_sam_tag = 0;

    kmer_clear_sam_read_override();
    kmer_clear_sam_optional_tag_override();

    if (anchor_runtime && anchor_runtime->enabled) {
        if (anchor_extract_window_range_info(anchor_runtime,
                                             read_seq,
                                             read_len,
                                             &anchor_insert_start,
                                             &anchor_insert_len,
                                             &anchor_info) != 0) {
            if (sam_out && sam_emit_unmapped) {
                const char *partial_tag = NULL;
                char partial_anchor_tag[1024];
                AnchorWindowInfo partial_info = {0};
                if (anchor_runtime->has_anchor5 &&
                    anchor_runtime->has_anchor3 &&
                    anchor_extract_two_sided_partial_start_info(anchor_runtime,
                                                                read_seq,
                                                                read_len,
                                                                &partial_info) == 0 &&
                    anchor_build_partial_two_sided_start_sam_tag(anchor_runtime,
                                                                 read_seq,
                                                                 read_len,
                                                                 &partial_info,
                                                                 partial_anchor_tag,
                                                                 sizeof(partial_anchor_tag)) == 0) {
                    partial_tag = partial_anchor_tag;
                }
                sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len, partial_tag);
            }
            goto cleanup;
        }
        have_anchor_info = 1;
        if (anchor_build_sam_tag(anchor_runtime,
                                 read_seq,
                                 read_len,
                                 &anchor_info,
                                 anchor_sam_tag,
                                 sizeof(anchor_sam_tag)) == 0) {
            have_anchor_sam_tag = 1;
        }
        search_seq = read_seq + anchor_insert_start;
        search_len = anchor_insert_len;
        if (read_qual && strlen(read_qual) == read_len &&
            anchor_insert_start + anchor_insert_len <= read_len) {
            search_qual = read_qual + anchor_insert_start;
        } else {
            search_qual = NULL;
        }
    }

    matches = tls_get_matches();
    if (!matches) {
        status = 1;
        goto cleanup;
    }
    strset_clear_with_keys(matches);

    if (anchor_runtime && anchor_runtime->enabled) {
        // 1) Find hits in the anchor-extracted segment (required).
        find_kmer_bitset(search_seq, search_len, seed_index, seed_mm, &hits);
        if (hits.n == 0) {
            if (sam_out && sam_emit_unmapped) {
                sam_write_unmapped(sam_out,
                                   read_name,
                                   read_seq,
                                   read_qual,
                                   read_len,
                                   have_anchor_sam_tag ? anchor_sam_tag : NULL);
            }
            goto cleanup;
        }

        find_matches_seeded(search_seq,
                            search_len,
                            &hits,
                            (uint32_t)min_len,
                            (uint32_t)max_len,
                            seed_index,
                            seed_mm,
                            matches,
                            k_mm,
                            read_name,
                            NULL,
                            NULL,
                            sam_soft_clip,
                            0);

        size_t anchor_hits = kh_size(matches);
        if (anchor_hits == 0) {
            if (sam_out && sam_emit_unmapped) {
                sam_write_unmapped(sam_out,
                                   read_name,
                                   read_seq,
                                   read_qual,
                                   read_len,
                                   have_anchor_sam_tag ? anchor_sam_tag : NULL);
            }
            goto cleanup;
        }

        // 2) Scan the read remainder for additional mappings (multi-hit classification).
        if (anchor_insert_start >= min_len) {
            hits.n = 0;
            find_kmer_bitset(read_seq, anchor_insert_start, seed_index, seed_mm, &hits);
            if (hits.n > 0) {
                find_matches_seeded(read_seq,
                                    anchor_insert_start,
                                    &hits,
                                    (uint32_t)min_len,
                                    (uint32_t)max_len,
                                    seed_index,
                                    seed_mm,
                                    matches,
                                    k_mm,
                                    read_name,
                                    NULL,
                                    NULL,
                                    sam_soft_clip,
                                    0);
            }
        }

        size_t tail_start = anchor_insert_start + anchor_insert_len;
        if (tail_start < read_len && (read_len - tail_start) >= min_len) {
            const char *tail_seq = read_seq + tail_start;
            hits.n = 0;
            find_kmer_bitset(tail_seq, read_len - tail_start, seed_index, seed_mm, &hits);
            if (hits.n > 0) {
                find_matches_seeded(tail_seq,
                                    read_len - tail_start,
                                    &hits,
                                    (uint32_t)min_len,
                                    (uint32_t)max_len,
                                    seed_index,
                                    seed_mm,
                                    matches,
                                    k_mm,
                                    read_name,
                                    NULL,
                                    NULL,
                                    sam_soft_clip,
                                    0);
            }
        }

        // 3) Emit SAM for anchor-extracted hits only.
        // NH should reflect records emitted to SAM (not suppressed remainder hits).
        if (sam_out && anchor_hits > 0) {
            sam_matches = kh_init(strset);
            if (!sam_matches) {
                status = 1;
                goto cleanup;
            }
            hits.n = 0;
            find_kmer_bitset(search_seq, search_len, seed_index, seed_mm, &hits);
            if (hits.n > 0) {
                const char *full_read_qual =
                    (read_qual && strlen(read_qual) == read_len) ? read_qual : NULL;
                if (have_anchor_info && have_anchor_sam_tag) {
                    kmer_set_sam_optional_tag_override(anchor_sam_tag);
                }
                kmer_set_sam_read_override(read_seq,
                                           full_read_qual,
                                           (int)read_len,
                                           (int)anchor_insert_start);
                find_matches_seeded(search_seq,
                                    search_len,
                                    &hits,
                                    (uint32_t)min_len,
                                    (uint32_t)max_len,
                                    seed_index,
                                    seed_mm,
                                    sam_matches,
                                    k_mm,
                                    read_name,
                                    search_qual,
                                    sam_out,
                                    sam_soft_clip,
                                    0);
                if (have_anchor_info && have_anchor_sam_tag) kmer_clear_sam_optional_tag_override();
                kmer_clear_sam_read_override();
            }
        }
    } else {
        // Standard path (no anchor gating).
        find_kmer_bitset(search_seq, search_len, seed_index, seed_mm, &hits);

        if (hits.n == 0) {
            if (sam_out && sam_emit_unmapped) {
                sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len, NULL);
            }
            goto cleanup;
        }

        find_matches_seeded(search_seq,
                            search_len,
                            &hits,
                            (uint32_t)min_len,
                            (uint32_t)max_len,
                            seed_index,
                            seed_mm,
                            matches,
                            k_mm,
                            read_name,
                            search_qual,
                            sam_out,
                            sam_soft_clip,
                            0);

        if (sam_out && kh_size(matches) == 0 && sam_emit_unmapped) {
            sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len, NULL);
        }
    }

    if (!(exclude_multihit && kh_size(matches) > 1)) {
        if (counts_preseeded) {
            add_to_counter(matches, counts);
        } else if (add_matches_to_counter_hits(matches, counts) != 0) {
            status = 1;
            goto cleanup;
        }
    }

cleanup:
    kmer_clear_sam_read_override();
    kmer_clear_sam_optional_tag_override();
    if (sam_matches) {
        strset_clear_with_keys(sam_matches);
        kh_destroy(strset, sam_matches);
    }
    if (matches) strset_clear_with_keys(matches);
    kv_destroy(hits);
    return status;
}

static int load_fastq_single(const char *path,
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
                             const AnchorConfig *anchor_cfg) {
    gzFile fp = gzopen(path, "rb");
    if (fp == 0) { perror("gzopen"); return 1; }

    kseq_t *ks = kseq_init(fp);
    if (!ks) { gzclose(fp); return 1; }

    size_t nreads = 0;
    int l = 0;
    int status = 0;
    AnchorRuntime anchor_runtime;
    KmerBitset *seed_index = kmer_bitset_from_trie(root, (size_t)kmerlen);
    if (!seed_index) {
        fprintf(stderr, "Failed to build seed bitset index.\n");
        kseq_destroy(ks);
        gzclose(fp);
        return 1;
    }
    if (anchor_runtime_init_for_count_mode(&anchor_runtime, anchor_cfg, *min_len) != 0) {
        fprintf(stderr, "Failed to initialize anchor matcher. Check adapter sequences.\n");
        kseq_destroy(ks);
        gzclose(fp);
        kmer_bitset_destroy(seed_index);
        return 1;
    }

    while ((l = kseq_read(ks)) >= 0) {
        if (++nreads % STEP == 0) {
            fprintf(stderr, "processing %zu reads...\r", nreads);
            fflush(stderr);
        }

        if (readseq_process_one_read(ks->name.s,
                                     ks->seq.s,
                                     ks->qual.l ? ks->qual.s : NULL,
                                     ks->seq.l,
                                     counts,
                                     seed_index,
                                     &anchor_runtime,
                                     seed_mm,
                                     *min_len,
                                     *max_len,
                                     k_mm,
                                     exclude_multihit,
                                     sam_soft_clip,
                                     sam_emit_unmapped,
                                     sam_fp,
                                     1) != 0) {
            status = 1;
            break;
        }
    }

    if (l < -1) status = 1;

    if (status == 0) total_hits(counts, nreads);

    kseq_destroy(ks);
    gzclose(fp);
    kmer_bitset_destroy(seed_index);
    return status;
}

int load_fastq(const char *path, TrieNode *root, kh_counter_t *counts,
               int kmerlen, int seed_mm, size_t *min_len, size_t *max_len, int k_mm,
               int exclude_multihit, FILE *sam_fp, int sam_soft_clip,
               int sam_emit_unmapped,
               unsigned threads,
               const AnchorConfig *anchor_cfg) {
    if (threads <= 1) {
        return load_fastq_single(path, root, counts, kmerlen, seed_mm, min_len, max_len, k_mm,
                                 exclude_multihit, sam_soft_clip, sam_emit_unmapped, sam_fp, anchor_cfg);
    }
    return readseq_load_fastq_mt(path, root, counts, kmerlen, seed_mm, min_len, max_len, k_mm,
                                 exclude_multihit, sam_soft_clip, sam_emit_unmapped, sam_fp, threads, anchor_cfg);
}
