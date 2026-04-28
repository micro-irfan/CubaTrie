#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "../trie.h"
#include "../utils.h"
#include "../readseq.h"

static int file_exists(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return 0;
    fclose(fp);
    return 1;
}

static char *read_text_normalized(const char *path, size_t *out_len) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return NULL;
    if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); return NULL; }
    long sz = ftell(fp);
    if (sz < 0) { fclose(fp); return NULL; }
    rewind(fp);

    char *raw = (char*)malloc((size_t)sz + 1);
    if (!raw) { fclose(fp); return NULL; }
    size_t got = fread(raw, 1, (size_t)sz, fp);
    fclose(fp);
    raw[got] = '\0';

    char *norm = (char*)malloc(got + 1);
    if (!norm) { free(raw); return NULL; }
    size_t j = 0;
    for (size_t i = 0; i < got; ++i) {
        if (raw[i] == '\r') {
            if (i + 1 < got && raw[i + 1] == '\n') continue;
            norm[j++] = '\n';
            continue;
        }
        norm[j++] = raw[i];
    }
    norm[j] = '\0';
    free(raw);
    if (out_len) *out_len = j;
    return norm;
}

static int compare_text_files(const char *expected, const char *actual, const char *label) {
    size_t e_len = 0, a_len = 0;
    char *e = read_text_normalized(expected, &e_len);
    char *a = read_text_normalized(actual, &a_len);
    if (!e || !a) {
        fprintf(stderr, "ERROR: failed to read %s files (%s, %s)\n", label, expected, actual);
        free(e); free(a);
        return 1;
    }
    int diff = (e_len != a_len) || (memcmp(e, a, e_len) != 0);
    if (diff) {
        fprintf(stderr, "MISMATCH: %s differs.\nExpected: %s\nActual:   %s\n", label, expected, actual);
    }
    free(e);
    free(a);
    return diff;
}

static size_t counter_get_or_zero(const kh_counter_t *m, const char *key) {
    if (!m || !key) return 0;
    khiter_t it = kh_get(counter, m, key);
    if (it == kh_end(m)) return 0;
    return kh_val(m, it);
}

static int write_test_fastq_gz(const char *path) {
    gzFile fp = gzopen(path, "wb");
    if (!fp) return 1;
    const char *fastq =
        "@r_multi\n"
        "ACGTA\n"
        "+\n"
        "IIIII\n"
        "@r_single\n"
        "AACGT\n"
        "+\n"
        "IIIII\n";
    int ok = (gzputs(fp, fastq) >= 0) ? 0 : 1;
    gzclose(fp);
    return ok;
}

static int write_test_fastq_with_unmapped_gz(const char *path) {
    gzFile fp = gzopen(path, "wb");
    if (!fp) return 1;
    const char *fastq =
        "@r_map\n"
        "ACGTA\n"
        "+\n"
        "IIIII\n"
        "@r_unmapped\n"
        "TTTTT\n"
        "+\n"
        "IIIII\n";
    int ok = (gzputs(fp, fastq) >= 0) ? 0 : 1;
    gzclose(fp);
    return ok;
}

static int write_test_fastq_anchor_nh_case_gz(const char *path) {
    gzFile fp = gzopen(path, "wb");
    if (!fp) return 1;
    const char *fastq =
        "@r_anchor_nh\n"
        "TGCAGGGGTTACGTAAAA\n"
        "+\n"
        "IIIIIIIIIIIIIIIIII\n";
    int ok = (gzputs(fp, fastq) >= 0) ? 0 : 1;
    gzclose(fp);
    return ok;
}

static int write_test_fastq_anchor_unmapped_case_gz(const char *path) {
    gzFile fp = gzopen(path, "wb");
    if (!fp) return 1;
    const char *fastq =
        "@r_anchor_unmapped\n"
        "GGGGTTTT\n"
        "+\n"
        "IIIIIIII\n";
    int ok = (gzputs(fp, fastq) >= 0) ? 0 : 1;
    gzclose(fp);
    return ok;
}

static int write_test_fastq_anchor_partial_a5_case_gz(const char *path) {
    gzFile fp = gzopen(path, "wb");
    if (!fp) return 1;
    const char *fastq =
        "@r_anchor_partial_a5\n"
        "GGGGTTTTAAAA\n"
        "+\n"
        "IIIIIIIIIIII\n";
    int ok = (gzputs(fp, fastq) >= 0) ? 0 : 1;
    gzclose(fp);
    return ok;
}

static int run_sam_unmapped_case(const char *reads_path,
                                 const char *sam_path,
                                 int sam_emit_unmapped,
                                 int *has_mapped_out,
                                 int *has_unmapped_out) {
    TrieNode *root = trie_create_node();
    kh_counter_t *map = kh_init(counter);
    if (!root || !map) {
        if (root) trie_free_node(root);
        if (map) kh_destroy(counter, map);
        return 1;
    }
    if (trie_insert(root, "ACGT", "refA", 0) != TRIE_INSERT_OK) {
        counter_free(map);
        trie_free_node(root);
        return 1;
    }
    if (counter_add_with_init(map, "refA", 0, 0) != 0) {
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    FILE *sam_fp = fopen(sam_path, "w");
    if (!sam_fp) {
        counter_free(map);
        trie_free_node(root);
        return 1;
    }
    trie_write_sam_header(sam_fp, root);

    size_t min_len = 4, max_len = 4;
    int rc = load_fastq(reads_path, root, map, 4, 0, &min_len, &max_len, 0, 0,
                        sam_fp, 0, sam_emit_unmapped, 1, NULL);
    fclose(sam_fp);

    if (rc != 0) {
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    size_t sam_len = 0;
    char *sam_txt = read_text_normalized(sam_path, &sam_len);
    (void)sam_len;
    if (!sam_txt) {
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    *has_mapped_out = (strstr(sam_txt, "\t0\trefA\t1\t255\t") != NULL);
    *has_unmapped_out = (strstr(sam_txt, "\t4\t*\t0\t0\t*\t*\t0\t0\t") != NULL);

    free(sam_txt);
    counter_free(map);
    trie_free_node(root);
    return 0;
}

static int run_multihit_counting_case(const char *reads_path, int exclude_multihit,
                                      size_t *refA_out, size_t *refB_out) {
    TrieNode *root = trie_create_node();
    kh_counter_t *map = kh_init(counter);
    if (!root || !map) {
        if (root) trie_free_node(root);
        if (map) kh_destroy(counter, map);
        return 1;
    }
    if (trie_insert(root, "ACGT", "refA", 0) != TRIE_INSERT_OK ||
        trie_insert(root, "CGTA", "refB", 0) != TRIE_INSERT_OK) {
        counter_free(map);
        trie_free_node(root);
        return 1;
    }
    if (counter_add_with_init(map, "refA", 0, 0) != 0 ||
        counter_add_with_init(map, "refB", 0, 0) != 0) {
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    size_t min_len = 4, max_len = 4;
    int rc = load_fastq(reads_path, root, map, 4, 0, &min_len, &max_len, 0, exclude_multihit, NULL, 0, 1, 1, NULL);
    if (rc != 0) {
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    *refA_out = counter_get_or_zero(map, "refA");
    *refB_out = counter_get_or_zero(map, "refB");
    counter_free(map);
    trie_free_node(root);
    return 0;
}

static int test_exclude_multihit_counting(void) {
    const char *reads_path = "tests/golden/tmp.multihit.fastq.gz";
    if (write_test_fastq_gz(reads_path) != 0) {
        fprintf(stderr, "ERROR: failed to write test FASTQ fixture: %s\n", reads_path);
        return 1;
    }

    size_t refA_keep = 0, refB_keep = 0;
    size_t refA_excl = 0, refB_excl = 0;
    int failed = 0;

    if (run_multihit_counting_case(reads_path, 0, &refA_keep, &refB_keep) != 0) failed = 1;
    if (run_multihit_counting_case(reads_path, 1, &refA_excl, &refB_excl) != 0) failed = 1;

    remove(reads_path);

    if (failed) {
        fprintf(stderr, "ERROR: multihit counting test setup failed.\n");
        return 1;
    }

    // Without exclusion: r_multi contributes to refA and refB; r_single contributes to refA.
    if (refA_keep != 2 || refB_keep != 1) {
        fprintf(stderr, "MISMATCH: expected keep counts refA=2 refB=1, got refA=%zu refB=%zu\n",
                refA_keep, refB_keep);
        return 1;
    }
    // With exclusion: r_multi is ignored for counting; only r_single counts.
    if (refA_excl != 1 || refB_excl != 0) {
        fprintf(stderr, "MISMATCH: expected exclude counts refA=1 refB=0, got refA=%zu refB=%zu\n",
                refA_excl, refB_excl);
        return 1;
    }
    return 0;
}

static int test_no_sam_unmapped_output(void) {
    const char *reads_path = "tests/golden/tmp.unmapped.fastq.gz";
    const char *sam_keep = "tests/golden/tmp.unmapped.keep.sam";
    const char *sam_drop = "tests/golden/tmp.unmapped.drop.sam";

    if (write_test_fastq_with_unmapped_gz(reads_path) != 0) {
        fprintf(stderr, "ERROR: failed to write unmapped FASTQ fixture: %s\n", reads_path);
        return 1;
    }

    int keep_has_mapped = 0, keep_has_unmapped = 0;
    int drop_has_mapped = 0, drop_has_unmapped = 0;
    int failed = 0;

    if (run_sam_unmapped_case(reads_path, sam_keep, 1, &keep_has_mapped, &keep_has_unmapped) != 0) failed = 1;
    if (run_sam_unmapped_case(reads_path, sam_drop, 0, &drop_has_mapped, &drop_has_unmapped) != 0) failed = 1;

    remove(reads_path);
    remove(sam_keep);
    remove(sam_drop);

    if (failed) {
        fprintf(stderr, "ERROR: unmapped SAM test setup failed.\n");
        return 1;
    }
    if (!keep_has_mapped || !drop_has_mapped) {
        fprintf(stderr, "MISMATCH: expected mapped SAM records in both modes.\n");
        return 1;
    }
    if (!keep_has_unmapped) {
        fprintf(stderr, "MISMATCH: expected unmapped SAM record when sam_emit_unmapped=1.\n");
        return 1;
    }
    if (drop_has_unmapped) {
        fprintf(stderr, "MISMATCH: unexpected unmapped SAM record when sam_emit_unmapped=0.\n");
        return 1;
    }
    return 0;
}

static int test_anchor_sam_nh_consistency(void) {
    const char *reads_path = "tests/golden/tmp.anchor_nh.fastq.gz";
    const char *sam_path = "tests/golden/tmp.anchor_nh.sam";
    int failed = 0;

    if (write_test_fastq_anchor_nh_case_gz(reads_path) != 0) {
        fprintf(stderr, "ERROR: failed to write anchor NH FASTQ fixture: %s\n", reads_path);
        return 1;
    }

    TrieNode *root = trie_create_node();
    kh_counter_t *map = kh_init(counter);
    if (!root || !map) {
        if (root) trie_free_node(root);
        if (map) kh_destroy(counter, map);
        remove(reads_path);
        return 1;
    }
    if (trie_insert(root, "ACGT", "refA", 0) != TRIE_INSERT_OK ||
        trie_insert(root, "TGCA", "refB", 0) != TRIE_INSERT_OK) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        return 1;
    }
    if (counter_add_with_init(map, "refA", 0, 0) != 0 ||
        counter_add_with_init(map, "refB", 0, 0) != 0) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        return 1;
    }

    FILE *sam_fp = fopen(sam_path, "w");
    if (!sam_fp) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        return 1;
    }
    trie_write_sam_header(sam_fp, root);

    AnchorConfig anchor_cfg = {0};
    anchor_cfg.enabled = 1;
    anchor_cfg.anchor5 = "GGGGTT";
    anchor_cfg.anchor3 = NULL;
    anchor_cfg.max_error = 0;

    size_t min_len = 4, max_len = 4;
    int rc = load_fastq(reads_path, root, map, 4, 0, &min_len, &max_len, 0, 0,
                        sam_fp, 0, 0, 1, &anchor_cfg);
    fclose(sam_fp);
    if (rc != 0) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        remove(sam_path);
        return 1;
    }

    size_t sam_len = 0;
    char *sam_txt = read_text_normalized(sam_path, &sam_len);
    (void)sam_len;
    if (!sam_txt) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        remove(sam_path);
        return 1;
    }

    if (strstr(sam_txt, "r_anchor_nh\t0\trefA\t1\t255\t") == NULL) {
        fprintf(stderr, "MISMATCH: expected one mapped anchor SAM record for refA.\n");
        failed = 1;
    }
    if (strstr(sam_txt, "r_anchor_nh\t0\trefB\t1\t255\t") != NULL) {
        fprintf(stderr, "MISMATCH: unexpected non-anchor SAM record for refB.\n");
        failed = 1;
    }
    if (strstr(sam_txt, "\tNH:i:2") != NULL) {
        fprintf(stderr, "MISMATCH: NH:i:2 present despite only one emitted alignment.\n");
        failed = 1;
    }

    free(sam_txt);
    counter_free(map);
    trie_free_node(root);
    remove(reads_path);
    remove(sam_path);
    return failed;
}

static int test_anchor_unmapped_includes_za_tag(void) {
    const char *reads_path = "tests/golden/tmp.anchor_unmapped.fastq.gz";
    const char *sam_path = "tests/golden/tmp.anchor_unmapped.sam";
    int failed = 0;

    if (write_test_fastq_anchor_unmapped_case_gz(reads_path) != 0) {
        fprintf(stderr, "ERROR: failed to write anchor-unmapped FASTQ fixture: %s\n", reads_path);
        return 1;
    }

    TrieNode *root = trie_create_node();
    kh_counter_t *map = kh_init(counter);
    if (!root || !map) {
        if (root) trie_free_node(root);
        if (map) kh_destroy(counter, map);
        remove(reads_path);
        return 1;
    }
    if (trie_insert(root, "ACGT", "refA", 0) != TRIE_INSERT_OK) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        return 1;
    }
    if (counter_add_with_init(map, "refA", 0, 0) != 0) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        return 1;
    }

    FILE *sam_fp = fopen(sam_path, "w");
    if (!sam_fp) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        return 1;
    }
    trie_write_sam_header(sam_fp, root);

    AnchorConfig anchor_cfg = {0};
    anchor_cfg.enabled = 1;
    anchor_cfg.anchor5 = "GGGG";
    anchor_cfg.anchor3 = NULL;
    anchor_cfg.max_error = 0;

    size_t min_len = 4, max_len = 4;
    int rc = load_fastq(reads_path, root, map, 4, 0, &min_len, &max_len, 0, 0,
                        sam_fp, 0, 1, 1, &anchor_cfg);
    fclose(sam_fp);
    if (rc != 0) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        remove(sam_path);
        return 1;
    }

    size_t sam_len = 0;
    char *sam_txt = read_text_normalized(sam_path, &sam_len);
    (void)sam_len;
    if (!sam_txt) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        remove(sam_path);
        return 1;
    }

    if (strstr(sam_txt, "r_anchor_unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGTTTT\tIIIIIIII\tZA:Z:ori=FWD;ins=5,4;a5=1-4,ed=0,md=4\n") == NULL) {
        fprintf(stderr, "MISMATCH: expected unmapped anchor SAM record with ZA tag.\n");
        failed = 1;
    }

    free(sam_txt);
    counter_free(map);
    trie_free_node(root);
    remove(reads_path);
    remove(sam_path);
    return failed;
}

static int test_two_sided_partial_start_anchor_tag(void) {
    const char *reads_path = "tests/golden/tmp.anchor_partial_a5.fastq.gz";
    const char *sam_path = "tests/golden/tmp.anchor_partial_a5.sam";
    int failed = 0;

    if (write_test_fastq_anchor_partial_a5_case_gz(reads_path) != 0) {
        fprintf(stderr, "ERROR: failed to write anchor-partial FASTQ fixture: %s\n", reads_path);
        return 1;
    }

    TrieNode *root = trie_create_node();
    kh_counter_t *map = kh_init(counter);
    if (!root || !map) {
        if (root) trie_free_node(root);
        if (map) kh_destroy(counter, map);
        remove(reads_path);
        return 1;
    }
    if (trie_insert(root, "ACGT", "refA", 0) != TRIE_INSERT_OK) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        return 1;
    }
    if (counter_add_with_init(map, "refA", 0, 0) != 0) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        return 1;
    }

    FILE *sam_fp = fopen(sam_path, "w");
    if (!sam_fp) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        return 1;
    }
    trie_write_sam_header(sam_fp, root);

    AnchorConfig anchor_cfg = {0};
    anchor_cfg.enabled = 1;
    anchor_cfg.anchor5 = "GGGG";
    anchor_cfg.anchor3 = "CCCC";
    anchor_cfg.max_error = 0;

    size_t min_len = 4, max_len = 4;
    int rc = load_fastq(reads_path, root, map, 4, 0, &min_len, &max_len, 0, 0,
                        sam_fp, 0, 1, 1, &anchor_cfg);
    fclose(sam_fp);
    if (rc != 0) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        remove(sam_path);
        return 1;
    }

    size_t sam_len = 0;
    char *sam_txt = read_text_normalized(sam_path, &sam_len);
    (void)sam_len;
    if (!sam_txt) {
        counter_free(map);
        trie_free_node(root);
        remove(reads_path);
        remove(sam_path);
        return 1;
    }

    if (strstr(sam_txt, "r_anchor_partial_a5\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGTTTTAAAA\tIIIIIIIIIIII\tZA:Z:ori=FWD;partial=1;a5=1-4,ed=0,md=4\n") == NULL) {
        fprintf(stderr, "MISMATCH: expected two-sided partial start anchor ZA tag on unmapped SAM.\n");
        failed = 1;
    }

    free(sam_txt);
    counter_free(map);
    trie_free_node(root);
    remove(reads_path);
    remove(sam_path);
    return failed;
}

int main(void) {
    if (test_exclude_multihit_counting() != 0) return 1;
    if (test_no_sam_unmapped_output() != 0) return 1;
    if (test_anchor_sam_nh_consistency() != 0) return 1;
    if (test_anchor_unmapped_includes_za_tag() != 0) return 1;
    if (test_two_sided_partial_start_anchor_tag() != 0) return 1;

    const char *ref_path = "tests/golden/ref.fa";
    const char *reads_path = "tests/golden/reads.fastq.gz";
    const char *expected_csv = "tests/golden/expected.csv";
    const char *expected_sam = "tests/golden/expected.sam";
    const char *actual_csv = "tests/golden/actual.csv";
    const char *actual_sam = "tests/golden/actual.sam";

    if (!file_exists(ref_path) || !file_exists(reads_path) ||
        !file_exists(expected_csv) || !file_exists(expected_sam)) {
        fprintf(stderr,
                "SKIP golden test: add fixtures at\n"
                "  %s\n  %s\n  %s\n  %s\n",
                ref_path, reads_path, expected_csv, expected_sam);
        return 0;
    }

    TrieNode *root = trie_create_node();
    kh_counter_t *map = kh_init(counter);
    if (!root || !map) {
        fprintf(stderr, "ERROR: allocation failed\n");
        if (root) trie_free_node(root);
        if (map) kh_destroy(counter, map);
        return 1;
    }

    size_t min_len = 0, max_len = 0, n_loaded = 0;
    size_t kmer_len = 4; // tiny fixture default
    int rc = load_reference(ref_path, root, map, &min_len, &max_len, &n_loaded,
                            1, &kmer_len, TRIE_DUP_ERROR);
    if (rc != 0) {
        fprintf(stderr, "ERROR: load_reference failed (rc=%d)\n", rc);
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    FILE *sam_fp = fopen(actual_sam, "w");
    if (!sam_fp) {
        perror("fopen");
        counter_free(map);
        trie_free_node(root);
        return 1;
    }
    trie_write_sam_header(sam_fp, root);

    rc = load_fastq(reads_path, root, map, (int)kmer_len, 0, &min_len, &max_len, 0, 0, sam_fp, 0, 1, 1, NULL);
    fclose(sam_fp);
    if (rc != 0) {
        fprintf(stderr, "ERROR: load_fastq failed (rc=%d)\n", rc);
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    if (counter_write_csv_sorted_collapse_rc(map, actual_csv, 0) != 0) {
        fprintf(stderr, "ERROR: failed to write actual csv\n");
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    int failed = 0;
    failed |= compare_text_files(expected_csv, actual_csv, "CSV");
    failed |= compare_text_files(expected_sam, actual_sam, "SAM");

    counter_free(map);
    trie_free_node(root);

    if (failed) return 1;
    fprintf(stderr, "Golden test passed.\n");
    return 0;
}
