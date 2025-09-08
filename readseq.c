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

static void normalize_acgt(char *s) {
    for (char *p = s; *p; ++p) {
        if (*p == 'u' || *p == 'U') *p = 'T';
        else *p = (char)toupper((unsigned char)*p);
    }
}

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


void strset_dump(const khash_t(strset) *S, FILE *out)
{
    fprintf(out, "# %zu items\n", (size_t)kh_size(S));
    for (khiter_t it = kh_begin(S); it != kh_end(S); ++it) {
        if (!kh_exist(S, it)) continue;
        fprintf(out, "%s\n", kh_key(S, it));
    }
}


void load_fastq(const char *path, TrieNode *root, kh_counter_t *counts,
                int kmerlen, size_t *min_len, size_t *max_len) {

    gzFile fp = gzopen(path, "rb");
    if (fp == 0) { perror("gzopen"); return 0; }

    kseq_t *ks = kseq_init(fp);
    if (!ks) { gzclose(fp); return 0; }

    int l;
    khash_t(strset) *kmer_hit = kh_init(strset);
    khash_t(strset) *searched = kh_init(strset);
    while ((l = kseq_read(ks)) >= 0) {
        normalize_acgt(ks->seq.s);
        
        // Check against previous record OR Trie (For new Kmers) if Kmer exists 
        kvec_t(uint32_t) hits; kv_init(hits);
        find_kmer(ks->seq.s, 
                  ks->seq.l, 
                  root,
                  kmer_hit,
                  searched,
                  kmerlen,
                  &hits);

        khash_t(strset) *matches = kh_init(strset);
        if (!matches) { perror("kh_init"); return 1; }

        // Find Matches Based on Hits From Kmer Search
        find_matches(
            ks->seq.s, 
            ks->seq.l, 
            &hits,
            *min_len, 
            *max_len,
            root, 
            matches
        );

        add_to_counter(matches, counts);

        kv_destroy(hits);
        kh_destroy(strset, matches);
    }

    fprintf(stderr, "Completed!\n");
    kseq_destroy(ks);
    kh_destroy(strset, kmer_hit);
    kh_destroy(strset, searched);
    return 0;
}



/* Returns number of forward sequences inserted (rev-comp not counted). 
   Works for FASTA or FASTQ, plain or .gz. */
void load_reference(const char *path, TrieNode *root, kh_counter_t *map, 
                    size_t *min_out, size_t *max_out, size_t *n_out, 
                    int add_revcomp, size_t *kmer_len) 
{
    gzFile fp = gzopen(path, "rb");
    if (fp == 0) { perror("gzopen"); return 0; }

    kseq_t *ks = kseq_init(fp);
    if (!ks) { gzclose(fp); return 0; }

    size_t inserted = 0, min_len = (size_t)-1, max_len = 0;
    int l;
    while ((l = kseq_read(ks)) >= 0) {
        /* ks->name.s    : sequence name (kstring)
           ks->comment.s : optional comment (may be empty)
           ks->seq.s     : sequence string
           ks->qual.s    : FASTQ qualities (empty for FASTA)
           ks->seq.l     : sequence length
        */

        if (ks->seq.l < kmer_len) continue;

        char *seq = malloc(ks->seq.l + 1);
        if (!seq) { perror("malloc"); /* handle */ }
        memcpy(seq, ks->seq.s, ks->seq.l);
        seq[ks->seq.l] = '\0';

        normalize_acgt(seq);

        size_t L = ks->seq.l;        // length of this record’s sequence
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

        trie_insert(root, seq, name_buf);
        counter_inc(map, name_buf);
        inserted++;

        if (add_revcomp) {
            char *rc = revcomp_new_n(seq, ks->seq.l);
            if (rc) {
                char rcname[256];
                snprintf(rcname, sizeof(rcname), "%s/rc", name_buf);
                trie_insert(root, rc, rcname);
                counter_inc(map, rcname);
                inserted++;
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
