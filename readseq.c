#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>

#include "kseq.h"         
KSEQ_INIT(gzFile, gzread)  

#include "trie.h" 
#include "utils.h"
#include "kmer.h"

static char *revcomp_new_n(const char *s, size_t n) {
    char *rc = (char*)malloc(n + 1);
    if (!rc) return NULL;
    for (size_t i = 0; i < n; ++i) {
        char c = s[n - 1 - i], r = 'N';
        if (c=='A') r='T'; else if (c=='C') r='G';
        else if (c=='G') r='C'; else if (c=='T') r='A';
        rc[i] = r;
    }
    rc[n] = '\0';
    return rc;
}

const size_t STEP = 1000000;

static size_t qname_len_no_rc(const char *name) {
    if (!name) return 0;
    size_t n = strlen(name);
    if (n >= 3 && name[n-3] == '/' && name[n-2] == 'r' && name[n-1] == 'c') return n - 3;
    return n;
}

static void sam_write_unmapped(FILE *sam_fp,
                               const char *read_name,
                               const char *read_seq,
                               const char *read_qual,
                               size_t read_len) {
    if (!sam_fp || !read_name || !read_seq) return;
    const char *qual = (read_qual && strlen(read_qual) == read_len) ? read_qual : "*";
    size_t qname_len = qname_len_no_rc(read_name);
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


int load_fastq(const char *path, TrieNode *root, kh_counter_t *counts,
               int kmerlen, size_t *min_len, size_t *max_len, int k_mm, FILE *sam_fp) {

    gzFile fp = gzopen(path, "rb");
    if (fp == 0) { perror("gzopen"); return 1; }

    kseq_t *ks = kseq_init(fp);
    if (!ks) { gzclose(fp); return 1; }

    size_t nreads = 0;
    int l;
    khash_t(kset64) *kmer_hit = kmer_set_from_trie(root, (size_t)kmerlen);
    khash_t(kset64) *searched = kh_init(kset64);
    if (!kmer_hit) kmer_hit = kh_init(kset64);
    if (!kmer_hit || !searched) {
        perror("kh_init");
        if (kmer_hit) kh_destroy(kset64, kmer_hit);
        if (searched) kh_destroy(kset64, searched);
        kseq_destroy(ks);
        gzclose(fp);
        return 1;
    }

    while ((l = kseq_read(ks)) >= 0) {
        if (++nreads % STEP == 0) {
            fprintf(stderr, "processing %zu reads...\r", nreads);
            fflush(stderr); 
        }

        kvec_t(uint32_t) hits; kv_init(hits);
        find_kmer(ks->seq.s, 
                  ks->seq.l, 
                  root,
                  kmer_hit,
                  searched,
                  kmerlen,
                  &hits);

        if (hits.n == 0) {
            sam_write_unmapped(sam_fp, ks->name.s, ks->seq.s, ks->qual.s, ks->seq.l);
            kv_destroy(hits);
            continue;
        } 

        khash_t(strset) *matches = kh_init(strset);
        if (!matches) {
            perror("kh_init");
            kv_destroy(hits);
            kseq_destroy(ks);
            kh_destroy(kset64, kmer_hit);
            kh_destroy(kset64, searched);
            gzclose(fp);
            return 1;
        }

        // Find Matches Based on Hits From Kmer Search
        find_matches(
            ks->seq.s, 
            ks->seq.l, 
            &hits,
            *min_len, 
            *max_len,
            root, 
            matches,
            k_mm,
            ks->name.s,
            ks->qual.s,
            sam_fp
        );

        if (sam_fp && kh_size(matches) == 0) {
            sam_write_unmapped(sam_fp, ks->name.s, ks->seq.s, ks->qual.s, ks->seq.l);
        }

        add_to_counter(matches, counts);

        kv_destroy(hits);
        kh_destroy(strset, matches);
    }

    total_hits (counts, nreads);

    kseq_destroy(ks);
    gzclose(fp);
    kh_destroy(kset64, kmer_hit);
    kh_destroy(kset64, searched);
    return 0;
}



/* Loads reference sequences into trie.
   Returns 0 on success, 1 on I/O/memory error, 2 on duplicate when dup-policy=error.
   Works for FASTA or FASTQ, plain or .gz. */
int load_reference(const char *path, TrieNode *root, kh_counter_t *map, 
                   size_t *min_out, size_t *max_out, size_t *n_out, 
                   int add_revcomp, size_t *kmer_len, TrieDupPolicy dup_policy) 
{
    gzFile fp = gzopen(path, "rb");
    if (fp == 0) { perror("gzopen"); return 1; }

    kseq_t *ks = kseq_init(fp);
    if (!ks) { gzclose(fp); return 1; }

    size_t inserted = 0, min_len = (size_t)-1, max_len = 0;
    size_t k = kmer_len ? *kmer_len : 0;
    int l;
    while ((l = kseq_read(ks)) >= 0) {
        /* ks->name.s    : sequence name (kstring)
           ks->comment.s : optional comment (may be empty)
           ks->seq.s     : sequence string
           ks->qual.s    : FASTQ qualities (empty for FASTA)
           ks->seq.l     : sequence length
        */

        if (k > 0 && ks->seq.l < k) continue;

        char *seq = malloc(ks->seq.l + 1);
        if (!seq) {
            perror("malloc");
            kseq_destroy(ks);
            gzclose(fp);
            return 1;
        }
        memcpy(seq, ks->seq.s, ks->seq.l);
        seq[ks->seq.l] = '\0';

        normalize_acgt(seq);

        size_t L = ks->seq.l; // length of this record’s sequence
        if (L < min_len) min_len = L;
        if (L > max_len) max_len = L;

        /* Use name up to first whitespace (common FASTA behavior) */
        char name_buf[256];
        size_t i = 0;
        while (ks->name.s[i] && !isspace((unsigned char)ks->name.s[i]) && i < sizeof(name_buf)-1) {
            name_buf[i] = ks->name.s[i];
            ++i;
        }
        name_buf[i] = '\0';

        const TrieNode *existing = trie_find_exact(root, seq);
        if (existing) {
            if (dup_policy != TRIE_DUP_IGNORE) {
                fprintf(stderr,
                        "Duplicate sequence detected: new seqID '%s' already exists as '%s'%s. "
                        "Possible reverse-complement or accidental duplicate.\n",
                        name_buf,
                        existing->name ? existing->name : "(unknown)",
                        existing->rev ? " [rev]" : "");
            }

            free(seq);
            if (dup_policy == TRIE_DUP_ERROR) {
                kseq_destroy(ks);
                gzclose(fp);
                return 2;
            }
            continue;
        }

        TrieInsertStatus st = trie_insert(root, seq, name_buf, 0);
        if (st == TRIE_INSERT_OOM) {
            free(seq);
            kseq_destroy(ks);
            gzclose(fp);
            return 1;
        }
        if (st == TRIE_INSERT_INVALID_BASE) {
            fprintf(stderr, "Skipping seqID '%s': invalid base found (expected A/C/G/T only).\n", name_buf);
            free(seq);
            continue;
        }
        if (st == TRIE_INSERT_DUP) {
            free(seq);
            if (dup_policy == TRIE_DUP_ERROR) {
                kseq_destroy(ks);
                gzclose(fp);
                return 2;
            }
            continue;
        }

        counter_inc(map, name_buf);
        inserted++;

        if (add_revcomp) {
            char *rc = revcomp_new_n(seq, ks->seq.l);
            if (rc) {
                char rcname[256];
                snprintf(rcname, sizeof(rcname), "%s/rc", name_buf);
                const TrieNode *existing_rc = trie_find_exact(root, rc);
                if (existing_rc) {
                    if (dup_policy != TRIE_DUP_IGNORE) {
                        fprintf(stderr,
                                "Duplicate sequence detected: new seqID '%s' already exists as '%s'%s. "
                                "Possible reverse-complement or accidental duplicate.\n",
                                rcname,
                                existing_rc->name ? existing_rc->name : "(unknown)",
                                existing_rc->rev ? " [rev]" : "");
                    }
                    if (dup_policy == TRIE_DUP_ERROR) {
                        free(rc);
                        free(seq);
                        kseq_destroy(ks);
                        gzclose(fp);
                        return 2;
                    }
                } else {
                    TrieInsertStatus rc_st = trie_insert(root, rc, rcname, 1);
                    if (rc_st == TRIE_INSERT_OOM) {
                        free(rc);
                        free(seq);
                        kseq_destroy(ks);
                        gzclose(fp);
                        return 1;
                    }
                    if (rc_st == TRIE_INSERT_INVALID_BASE) {
                        fprintf(stderr, "Skipping seqID '%s': invalid base found (expected A/C/G/T only).\n", rcname);
                    } else if (rc_st == TRIE_INSERT_DUP) {
                        if (dup_policy == TRIE_DUP_ERROR) {
                            free(rc);
                            free(seq);
                            kseq_destroy(ks);
                            gzclose(fp);
                            return 2;
                        }
                    } else {
                        counter_inc(map, rcname);
                        inserted++;
                    }
                }
                free(rc);
            }
        }
        free(seq);
    }

    if (min_out) *min_out = min_len;
    if (max_out) *max_out = max_len;
    if (n_out)   *n_out  = inserted;

    kseq_destroy(ks);
    gzclose(fp);
    return 0;
}
