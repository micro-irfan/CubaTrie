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
                               size_t read_len) {
    if (!sam_fp || !read_name || !read_seq) return;
    const char *qual = (read_qual && strlen(read_qual) == read_len) ? read_qual : "*";
    size_t qname_len = name_len_no_rc_suffix(read_name);
    fprintf(sam_fp, "%.*s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
            (int)qname_len, read_name, read_seq, qual);
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

    if (anchor_runtime && anchor_runtime->enabled) {
        if (anchor_extract_window(anchor_runtime,
                                  read_seq,
                                  read_len,
                                  &anchor_insert_start,
                                  &anchor_insert_len) != 0) {
            if (sam_out && sam_emit_unmapped) {
                sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len);
            }
            goto cleanup;
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
                sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len);
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
                sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len);
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

        // 3) Emit SAM for anchor-extracted hits only, with NH based on total combined hits.
        if (sam_out && anchor_hits > 0) {
            sam_matches = kh_init(strset);
            if (!sam_matches) {
                status = 1;
                goto cleanup;
            }
            hits.n = 0;
            find_kmer_bitset(search_seq, search_len, seed_index, seed_mm, &hits);
            if (hits.n > 0) {
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
                                    (int)kh_size(matches));
            }
        }
    } else {
        // Standard path (no anchor gating).
        find_kmer_bitset(search_seq, search_len, seed_index, seed_mm, &hits);

        if (hits.n == 0) {
            if (sam_out && sam_emit_unmapped) {
                sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len);
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
            sam_write_unmapped(sam_out, read_name, read_seq, read_qual, read_len);
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
    if (anchor_runtime_init(&anchor_runtime, anchor_cfg, *min_len) != 0) {
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
