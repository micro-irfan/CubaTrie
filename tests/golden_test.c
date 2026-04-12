#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "../trie.h"
#include "../utils.h"

// readseq.c entry points
int load_reference(const char *path, TrieNode *root, kh_counter_t *map,
                   size_t *min_out, size_t *max_out, size_t *n_out,
                   int add_revcomp, size_t *kmer_len, TrieDupPolicy dup_policy);
int load_fastq(const char *path, TrieNode *root, kh_counter_t *counts,
               int kmerlen, int seed_mm, size_t *min_len, size_t *max_len, int k_mm,
               int exclude_multihit, FILE *sam_fp, int sam_soft_clip, unsigned threads);

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
    int rc = load_fastq(reads_path, root, map, 4, 0, &min_len, &max_len, 0, exclude_multihit, NULL, 0, 1);
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

int main(void) {
    if (test_exclude_multihit_counting() != 0) return 1;

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

    rc = load_fastq(reads_path, root, map, (int)kmer_len, 0, &min_len, &max_len, 0, 0, sam_fp, 0, 1);
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
