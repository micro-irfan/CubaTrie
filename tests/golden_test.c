#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../trie.h"
#include "../utils.h"

// readseq.c entry points
int load_reference(const char *path, TrieNode *root, kh_counter_t *map,
                   size_t *min_out, size_t *max_out, size_t *n_out,
                   int add_revcomp, size_t *kmer_len, TrieDupPolicy dup_policy);
int load_fastq(const char *path, TrieNode *root, kh_counter_t *counts,
               int kmerlen, int seed_mm, size_t *min_len, size_t *max_len, int k_mm, FILE *sam_fp,
               unsigned threads);

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

int main(void) {
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

    rc = load_fastq(reads_path, root, map, (int)kmer_len, 0, &min_len, &max_len, 0, sam_fp, 1);
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
