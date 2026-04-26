#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "readseq.h"

static char *revcomp_new_n(const char *s, size_t n) {
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

        size_t L = ks->seq.l;
        if (L < min_len) min_len = L;
        if (L > max_len) max_len = L;

        char name_buf[256];
        size_t i = 0;
        while (ks->name.s[i] && !isspace((unsigned char)ks->name.s[i]) && i < sizeof(name_buf) - 1) {
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
