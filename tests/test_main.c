#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../trie.h"
#include "../kmer.h"
#include "../readseq_anchor.h"
#include "../utils.h"

static char *read_tmpfile_all(FILE *fp) {
    if (!fp) return NULL;
    fflush(fp);
    if (fseek(fp, 0, SEEK_END) != 0) return NULL;
    long n = ftell(fp);
    if (n < 0) return NULL;
    rewind(fp);

    char *buf = (char*)malloc((size_t)n + 1);
    if (!buf) return NULL;
    size_t got = fread(buf, 1, (size_t)n, fp);
    buf[got] = '\0';
    return buf;
}

static int count_substr(const char *s, const char *needle) {
    int c = 0;
    size_t nlen = strlen(needle);
    if (!nlen) return 0;
    const char *p = s;
    while ((p = strstr(p, needle)) != NULL) {
        ++c;
        p += nlen;
    }
    return c;
}

static void free_strset_keys_and_destroy(khash_t(strset) *set) {
    if (!set) return;
    for (khint_t i = kh_begin(set); i != kh_end(set); ++i) {
        if (kh_exist(set, i)) free((char*)kh_key(set, i));
    }
    kh_destroy(strset, set);
}

static void test_trie_insert_statuses(void) {
    TrieNode *root = trie_create_node();
    assert(root != NULL);

    assert(trie_insert(root, "ACGT", "id1", 0) == TRIE_INSERT_OK);
    assert(trie_insert(root, "ACGT", "id2", 0) == TRIE_INSERT_DUP);
    assert(trie_insert(root, "ACNT", "bad", 0) == TRIE_INSERT_INVALID_BASE);

    const TrieNode *n = trie_find_exact(root, "ACGT");
    assert(n != NULL);
    assert(n->name != NULL);
    assert(strcmp(n->name, "id1") == 0); // duplicate insert must not overwrite existing payload

    assert(trie_find_exact(root, "ACG") == NULL);
    trie_free_node(root);
}

static void test_trie_search_exact_mm(void) {
    TrieNode *root = trie_create_node();
    assert(root != NULL);
    assert(trie_insert(root, "ACGT", "refA", 0) == TRIE_INSERT_OK);

    kFoundVec out;
    kv_init(out);

    trie_search_exact(root, "AACCTAA", 4, 1, &out);
    assert(out.n >= 1);
    assert(strcmp(out.a[0].name, "refA") == 0);
    assert(out.a[0].mm == 1);
    mv_free(&out);

    kv_init(out);
    trie_search_exact(root, "AACCTAA", 4, 0, &out);
    assert(out.n == 0);
    mv_free(&out);

    trie_free_node(root);
}

static void test_find_kmer_bitset_seed_mm(void) {
    TrieNode *root = trie_create_node();
    assert(root != NULL);
    assert(trie_insert(root, "ACGT", "seedRef", 0) == TRIE_INSERT_OK);

    KmerBitset *index = kmer_bitset_from_trie(root, 4);
    assert(index != NULL);

    u32vec_t hits;
    kv_init(hits);

    find_kmer_bitset("AACGTAA", 7, index, 0, &hits);
    assert(hits.n == 1);
    assert(hits.a[0] == 1);
    kv_destroy(hits);

    kv_init(hits);
    find_kmer_bitset("AACCTAA", 7, index, 0, &hits);
    assert(hits.n == 0);
    kv_destroy(hits);

    kv_init(hits);
    find_kmer_bitset("AACCTAA", 7, index, 1, &hits);
    assert(hits.n == 1);
    assert(hits.a[0] == 1);
    kv_destroy(hits);

    kmer_bitset_destroy(index);
    trie_free_node(root);
}

static void test_find_matches_sam_line(void) {
    TrieNode *root = trie_create_node();
    assert(root != NULL);
    assert(trie_insert(root, "ACGT", "ref1/rc", 1) == TRIE_INSERT_OK);

    u32vec_t hit;
    kv_init(hit);
    kv_push(uint32_t, 0, hit, 2);

    khash_t(strset) *found = kh_init(strset);
    assert(found != NULL);

    FILE *fp = tmpfile();
    assert(fp != NULL);

    find_matches("TTACCTAA", 8, &hit, 4, 4, root, found, 1, "readX/rc", "HHHXHHHH", fp, 0);

    char *sam = read_tmpfile_all(fp);
    assert(sam != NULL);

    assert(strstr(sam, "readX\t16\tref1\t1\t255\t2H4M2H\t*\t0\t0\tACCT\tHHXH\tNM:i:1") != NULL);
    assert(strstr(sam, "/rc\t") == NULL);
    assert(strstr(sam, "\tNH:i:") == NULL);

    free(sam);
    fclose(fp);
    kv_destroy(hit);
    free_strset_keys_and_destroy(found);
    trie_free_node(root);
}

static void test_find_matches_sam_line_soft_clip(void) {
    TrieNode *root = trie_create_node();
    assert(root != NULL);
    assert(trie_insert(root, "ACGT", "ref1/rc", 1) == TRIE_INSERT_OK);

    u32vec_t hit;
    kv_init(hit);
    kv_push(uint32_t, 0, hit, 2);

    khash_t(strset) *found = kh_init(strset);
    assert(found != NULL);

    FILE *fp = tmpfile();
    assert(fp != NULL);

    find_matches("TTACCTAA", 8, &hit, 4, 4, root, found, 1, "readX/rc", "HHHXHHHH", fp, 1);

    char *sam = read_tmpfile_all(fp);
    assert(sam != NULL);

    assert(strstr(sam, "readX\t16\tref1\t1\t255\t2S4M2S\t*\t0\t0\tTTACCTAA\tHHHXHHHH\tNM:i:1") != NULL);
    assert(strstr(sam, "/rc\t") == NULL);
    assert(strstr(sam, "\tNH:i:") == NULL);

    free(sam);
    fclose(fp);
    kv_destroy(hit);
    free_strset_keys_and_destroy(found);
    trie_free_node(root);
}

static void test_find_matches_multihit_qname_tag(void) {
    TrieNode *root = trie_create_node();
    assert(root != NULL);
    assert(trie_insert(root, "ACGT", "refA", 0) == TRIE_INSERT_OK);
    assert(trie_insert(root, "CGTA", "refB", 0) == TRIE_INSERT_OK);

    u32vec_t hit;
    kv_init(hit);
    kv_push(uint32_t, 0, hit, 0);
    kv_push(uint32_t, 0, hit, 1);

    khash_t(strset) *found = kh_init(strset);
    assert(found != NULL);

    FILE *fp = tmpfile();
    assert(fp != NULL);

    find_matches("ACGTA", 5, &hit, 4, 4, root, found, 0, "readMulti", NULL, fp, 0);

    char *sam = read_tmpfile_all(fp);
    assert(sam != NULL);

    assert(strstr(sam, "readMulti\t0\trefA\t1\t255\t4M1H\t*\t0\t0\tACGT\t*\tNM:i:0\tNH:i:2") != NULL);
    assert(strstr(sam, "readMulti\t0\trefB\t1\t255\t1H4M\t*\t0\t0\tCGTA\t*\tNM:i:0\tNH:i:2") != NULL);
    assert(count_substr(sam, "\tNH:i:2") == 2);

    free(sam);
    fclose(fp);
    kv_destroy(hit);
    free_strset_keys_and_destroy(found);
    trie_free_node(root);
}

static void test_find_matches_md_tag(void) {
    TrieNode *root = trie_create_node();
    assert(root != NULL);
    assert(trie_insert(root, "ACGT", "ref1/rc", 1) == TRIE_INSERT_OK);

    u32vec_t hit;
    kv_init(hit);
    kv_push(uint32_t, 0, hit, 2);

    khash_t(strset) *found = kh_init(strset);
    assert(found != NULL);

    FILE *fp = tmpfile();
    assert(fp != NULL);

    find_matches("TTACCTAA", 8, &hit, 4, 4, root, found, 1, "readMD/rc", "HHHXHHHH", fp, 0);

    char *sam = read_tmpfile_all(fp);
    assert(sam != NULL);

    assert(strstr(sam, "readMD\t16\tref1\t1\t255\t2H4M2H\t*\t0\t0\tACCT\tHHXH\tNM:i:1\tMD:Z:2G1\n") != NULL);
    assert(count_substr(sam, "\tMD:Z:") == 1);

    free(sam);
    fclose(fp);
    kv_destroy(hit);
    free_strset_keys_and_destroy(found);
    trie_free_node(root);
}

static void test_sam_header_strip_and_dedupe(void) {
    TrieNode *root = trie_create_node();
    assert(root != NULL);

    // Same base name after stripping /rc; should appear once in @SQ.
    assert(trie_insert(root, "ACGT", "geneA", 0) == TRIE_INSERT_OK);
    assert(trie_insert(root, "TGCA", "geneA/rc", 1) == TRIE_INSERT_OK);
    assert(trie_insert(root, "CCCC", "geneB", 0) == TRIE_INSERT_OK);

    FILE *fp = tmpfile();
    assert(fp != NULL);
    trie_write_sam_header(fp, root);

    char *hdr = read_tmpfile_all(fp);
    assert(hdr != NULL);

    assert(strstr(hdr, "@HD\tVN:1.6\tSO:unknown\n") != NULL);
    assert(count_substr(hdr, "@SQ\tSN:geneA\tLN:") == 1);
    assert(count_substr(hdr, "@SQ\tSN:geneB\tLN:") == 1);
    assert(strstr(hdr, "/rc") == NULL);

    free(hdr);
    fclose(fp);
    trie_free_node(root);
}

static void test_anchor_extract_window_variants(void) {
    AnchorConfig cfg = {0};
    AnchorRuntime ar = {0};
    size_t start = 0, len = 0;

    assert(anchor_runtime_init_range(&ar, &cfg, 3, 3) == 0); // disabled config is a no-op
    assert(anchor_extract_window_range(&ar, "ACGT", 4, &start, &len) == -1);

    cfg.enabled = 1;
    cfg.anchor5 = "AGTC";
    cfg.anchor3 = "TCCA";
    cfg.max_error = 0;

    assert(anchor_runtime_init_range(&ar, &cfg, 3, 3) == 0);

    // Forward orientation: anchor5 ... payload ... anchor3.
    assert(anchor_extract_window_range(&ar, "CCAGTCGATTCCAGG", 15, &start, &len) == 0);
    assert(start == 6);
    assert(len == 3);
    assert(strncmp("CCAGTCGATTCCAGG" + start, "GAT", len) == 0);
    assert(anchor_extract_window(&ar, "CCAGTCGATTCCAGG", 15, &start, &len) == 0);
    assert(start == 6);
    assert(len == 3);

    // Reverse-complement orientation: anchor3_rc ... payload ... anchor5_rc.
    assert(anchor_extract_window_range(&ar, "AATGGACATGACTCC", 15, &start, &len) == 0);
    assert(start == 6);
    assert(len == 3);
    assert(strncmp("AATGGACATGACTCC" + start, "CAT", len) == 0);

    // Ambiguous: repeated 5' anchor should be rejected.
    assert(anchor_extract_window_range(&ar, "AGTCAAAGTCTCCA", 14, &start, &len) == 1);

    cfg.anchor3 = NULL;
    assert(anchor_runtime_init_range(&ar, &cfg, 4, 4) == 0);
    assert(anchor_extract_window_range(&ar, "TTAGTCGGGG", 10, &start, &len) == 0);
    assert(start == 6);
    assert(len == 4);
    assert(strncmp("TTAGTCGGGG" + start, "GGGG", len) == 0);
    // With fixed-length mode (min==max), one-sided 5' takes exactly that many nt after anchor.
    assert(anchor_extract_window_range(&ar, "TTAGTCGGGGAAAA", 14, &start, &len) == 0);
    assert(start == 6);
    assert(len == 4);
    assert(strncmp("TTAGTCGGGGAAAA" + start, "GGGG", len) == 0);

    cfg.anchor5 = NULL;
    cfg.anchor3 = "TCCA";
    assert(anchor_runtime_init_range(&ar, &cfg, 4, 4) == 0);
    assert(anchor_extract_window_range(&ar, "GGGGTCCAAA", 10, &start, &len) == 0);
    assert(start == 0);
    assert(len == 4);
    assert(strncmp("GGGGTCCAAA" + start, "GGGG", len) == 0);
    // With fixed-length mode (min==max), one-sided 3' takes exactly that many nt before anchor.
    assert(anchor_extract_window_range(&ar, "AAAAGGGGTCCAAA", 14, &start, &len) == 0);
    assert(start == 4);
    assert(len == 4);
    assert(strncmp("AAAAGGGGTCCAAA" + start, "GGGG", len) == 0);

    // Real-world regression: exact anchors should recover this 20nt guide.
    cfg.anchor5 = "ATTTTCAATTTAACGTCG";
    cfg.anchor3 = "GTTTTAGAGCTAGAAATA";
    cfg.max_error = 0;
    const char *real_read =
        "CTCCGTGACCTATTTTCAATTTAACGTCGTCCGATGCAGGCTCCAGAGGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGCCTACCTGGAGCCTGAGA";
    assert(anchor_runtime_init_range(&ar, &cfg, 20, 20) == 0);
    assert(anchor_extract_window_range(&ar, real_read, strlen(real_read), &start, &len) == 0);
    assert(len == 20);
    assert(strncmp(real_read + start, "TCCGATGCAGGCTCCAGAGG", len) == 0);
}

static void test_anchor_extract_window_realworld_mm_indel(void) {
    AnchorConfig cfg = {0};
    AnchorRuntime ar = {0};
    size_t start = 0, len = 0;

    cfg.enabled = 1;
    cfg.anchor5 = "ATTTTCAATTTAACGTCG";
    cfg.anchor3 = "GTTTTAGAGCTAGAAATA";

    const char *guide = "TCCGATGCAGGCTCCAGAGG";
    const char *real_read =
        "CTCCGTGACCTATTTTCAATTTAACGTCGTCCGATGCAGGCTCCAGAGGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGCCTACCTGGAGCCTGAGA";
    size_t real_len = strlen(real_read);

    const char *a5_pos = strstr(real_read, cfg.anchor5);
    const char *a3_pos = strstr(real_read, cfg.anchor3);
    assert(a5_pos != NULL);
    assert(a3_pos != NULL);
    size_t a5_start = (size_t)(a5_pos - real_read);
    size_t a3_start = (size_t)(a3_pos - real_read);

    // Case 1: one mismatch inside 5' anchor.
    char read_mm[256];
    assert(real_len + 1 < sizeof(read_mm));
    memcpy(read_mm, real_read, real_len + 1);
    read_mm[a5_start + 5] = (read_mm[a5_start + 5] == 'A') ? 'C' : 'A';

    cfg.max_error = 0;
    assert(anchor_runtime_init_range(&ar, &cfg, 20, 20) == 0);
    assert(anchor_extract_window_range(&ar, read_mm, strlen(read_mm), &start, &len) == 1);

    cfg.max_error = 1;
    assert(anchor_runtime_init_range(&ar, &cfg, 20, 20) == 0);
    assert(anchor_extract_window_range(&ar, read_mm, strlen(read_mm), &start, &len) == 0);
    assert(len == 20);
    assert(strncmp(read_mm + start, guide, len) == 0);

    // Case 2: one deletion inside 3' anchor.
    char read_del[256];
    assert(real_len + 1 < sizeof(read_del));
    size_t del_idx = a3_start + 7;
    memcpy(read_del, real_read, del_idx);
    memcpy(read_del + del_idx, real_read + del_idx + 1, real_len - del_idx);

    cfg.max_error = 1;
    assert(anchor_runtime_init_range(&ar, &cfg, 20, 20) == 0);
    assert(anchor_extract_window_range(&ar, read_del, strlen(read_del), &start, &len) == 0);
    assert(len == 20);
    assert(strncmp(read_del + start, guide, len) == 0);
}

static void test_anchor_extract_window_error_with_extra_near_hit(void) {
    AnchorConfig cfg = {0};
    AnchorRuntime ar = {0};
    size_t start = 0, len = 0;

    cfg.enabled = 1;
    cfg.anchor5 = "AAAA";
    cfg.anchor3 = "TTTT";
    cfg.max_error = 1;

    // With max_error=1 there are multiple 5' near-hits ("AAAT", "AATA", ...),
    // but only one valid anchor pair that yields a 3nt payload.
    const char *read = "AAATAAAAGGGTTTT";
    assert(anchor_runtime_init_range(&ar, &cfg, 3, 3) == 0);
    assert(anchor_extract_window_range(&ar, read, strlen(read), &start, &len) == 0);
    assert(start == 8);
    assert(len == 3);
    assert(strncmp(read + start, "GGG", len) == 0);
}

int main(void) {
    test_trie_insert_statuses();
    test_trie_search_exact_mm();
    test_find_kmer_bitset_seed_mm();
    test_find_matches_sam_line();
    test_find_matches_sam_line_soft_clip();
    test_find_matches_multihit_qname_tag();
    test_find_matches_md_tag();
    test_sam_header_strip_and_dedupe();
    test_anchor_extract_window_variants();
    test_anchor_extract_window_realworld_mm_indel();
    test_anchor_extract_window_error_with_extra_near_hit();

    fprintf(stderr, "All tests passed.\n");
    return 0;
}
