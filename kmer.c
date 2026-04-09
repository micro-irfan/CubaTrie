#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "kmer.h"
#include "utils.h"

// Trie API
typedef struct TrieNode TrieNode;
bool trie_prefix_search(const TrieNode *root,
                        const uint64_t *b2, 
                        size_t n);

// comparator for qsort (ascending)
static int cmp_u32_asc(const void *a, const void *b) {
    uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

static inline void clear_matches_keep_capacity(kFoundVec *mv) {
    if (!mv) return;
    for (size_t i = 0; i < mv->n; ++i) {
        free(mv->a[i].name);
        free(mv->a[i].seq);
    }
    mv->n = 0;
}

static int name_is_rev(const char *name) {
    if (!name) return 0;
    size_t n = strlen(name);
    return (n >= 3 && name[n-3] == '/' && name[n-2] == 'r' && name[n-1] == 'c');
}

static size_t qname_len_no_rc(const char *name) {
    if (!name) return 0;
    size_t n = strlen(name);
    if (n >= 3 && name[n-3] == '/' && name[n-2] == 'r' && name[n-1] == 'c') return n - 3;
    return n;
}

typedef struct {
    const char *ref_name;
    size_t match_start;
    size_t match_len;
    int nm;
} SamAlignmentRecord;

static void sam_write_alignment(FILE *sam_fp,
                                const char *read_name,
                                const char *read_seq,
                                const char *read_qual,
                                size_t read_len,
                                const char *ref_name,
                                size_t match_start,
                                size_t match_len,
                                int nm,
                                int nh) {
    if (!sam_fp || !read_name || !read_seq || !ref_name) return;
    if (match_start >= read_len) return;
    if (match_len > read_len - match_start) match_len = read_len - match_start;
    if (match_len == 0) return;

    size_t hard_start = match_start;
    size_t hard_end = read_len - match_start - match_len;
    char cigar[96];

    if (hard_start > 0 && hard_end > 0) {
        snprintf(cigar, sizeof(cigar), "%zuH%zuM%zuH", hard_start, match_len, hard_end);
    } else if (hard_start > 0) {
        snprintf(cigar, sizeof(cigar), "%zuH%zuM", hard_start, match_len);
    } else if (hard_end > 0) {
        snprintf(cigar, sizeof(cigar), "%zuM%zuH", match_len, hard_end);
    } else {
        snprintf(cigar, sizeof(cigar), "%zuM", match_len);
    }

    int flag = name_is_rev(ref_name) ? 16 : 0;
    int has_full_qual = (read_qual && strlen(read_qual) == read_len);
    size_t qname_len = qname_len_no_rc(read_name);
    size_t rname_len = qname_len_no_rc(ref_name);

    // For hard clipping, output only the aligned segment in SEQ/QUAL.
    if (has_full_qual) {
        if (nh > 1) {
            fprintf(sam_fp, "%.*s\t%d\t%.*s\t1\t255\t%s\t*\t0\t0\t%.*s\t%.*s\tNM:i:%d\tNH:i:%d\n",
                    (int)qname_len, read_name, flag, (int)rname_len, ref_name, cigar,
                    (int)match_len, read_seq + match_start,
                    (int)match_len, read_qual + match_start,
                    nm, nh);
        } else {
            fprintf(sam_fp, "%.*s\t%d\t%.*s\t1\t255\t%s\t*\t0\t0\t%.*s\t%.*s\tNM:i:%d\n",
                    (int)qname_len, read_name, flag, (int)rname_len, ref_name, cigar,
                    (int)match_len, read_seq + match_start,
                    (int)match_len, read_qual + match_start,
                    nm);
        }
    } else {
        if (nh > 1) {
            fprintf(sam_fp, "%.*s\t%d\t%.*s\t1\t255\t%s\t*\t0\t0\t%.*s\t*\tNM:i:%d\tNH:i:%d\n",
                    (int)qname_len, read_name, flag, (int)rname_len, ref_name, cigar,
                    (int)match_len, read_seq + match_start,
                    nm, nh);
        } else {
            fprintf(sam_fp, "%.*s\t%d\t%.*s\t1\t255\t%s\t*\t0\t0\t%.*s\t*\tNM:i:%d\n",
                    (int)qname_len, read_name, flag, (int)rname_len, ref_name, cigar,
                    (int)match_len, read_seq + match_start,
                    nm);
        }
    }
}

// Pack text[0..n) into 2-bit words (LSB-first). Also fill ambig-mask (1 bit per base) if provided.
static inline int encode_kmer(const char *s, size_t n, uint64_t *out) {
    uint64_t code = 0;
    for (size_t i = 0; i < n; ++i) {
        int b = nt2bits(s[i]);
        if (b < 0) return -1;
        code = (code << 2) | (uint64_t)b;

    }
    *out = code;
    return 0;
}

static inline void kset64_insert(khash_t(kset64) *set, uint64_t code) {
    int ret = 0;
    khiter_t it = kh_put(kset64, set, (khint64_t)code, &ret);
    if (ret > 0) kh_key(set, it) = (khint64_t)code;
}

static void kmer_set_add_sequence(khash_t(kset64) *out, const char *seq, size_t seq_len, size_t k) {
    if (!out || !seq || k == 0 || seq_len < k) return;
    uint64_t code = 0;
    if (encode_kmer(seq, k, &code) == 0) {
        kset64_insert(out, code);
    }
}

static void trie_collect_terminal_kmers(const TrieNode *node, size_t k, khash_t(kset64) *out) {
    if (!node) return;

    if (node->end && node->seq && node->seq_len >= k) {
        kmer_set_add_sequence(out, node->seq, node->seq_len, k);
    }

    for (int i = 0; i < ALPHABET_SIZE; ++i) {
        trie_collect_terminal_kmers(node->child[i], k, out);
    }
}

khash_t(kset64) *kmer_set_from_trie(const TrieNode *root, size_t k) {
    khash_t(kset64) *out = kh_init(kset64);
    if (!out || !root || k == 0) return out;

    trie_collect_terminal_kmers(root, k, out);
    return out;
}


void find_kmer(const char *s, const size_t s_len, 
               const TrieNode *trie,
               khash_t(kset64) *kmer_hits,
               khash_t(kset64) *kmer_search,
               size_t k,
               u32vec_t *hits) {

    if (!s || k <= 0 || (size_t)k > s_len) return;                    
    int has_precomputed_kmers = (kmer_hits && kh_size(kmer_hits) > 0);

    for (size_t i = 0; i + (size_t)k <= s_len; ++i) {
        uint64_t code;
        if (encode_kmer(s + i, k, &code) < 0) continue;

        // if Kmer in kmer_hit; add pos to check; else continue
        if (has_precomputed_kmers) {
            if (kh_get(kset64, kmer_hits, (khint64_t)code) != kh_end(kmer_hits)) {
                kv_push(uint32_t, 0, *hits, i);
            } 
            continue;
        }
        
        // fallback path: no precomputed set, use kmer_search as searched cache
        // if kmer in searched cache; continue
        if (kh_get(kset64, kmer_search, (khint64_t)code) != kh_end(kmer_search)) {
            continue;
        }

        // if Kmer not in kmer_hit; add pos to check; add kmer to kmer_hit
        if (trie && trie_prefix_search(trie, &code, k)) {
            kv_push(uint32_t, 0, *hits, i);
            kset64_insert(kmer_hits, code);
        }

        if (kmer_search) kset64_insert(kmer_search, code);
    }

    if (hits->n > 1) qsort(hits->a, hits->n, sizeof(uint32_t), cmp_u32_asc);
    
    return;
}


// ---------- core routine ----------
/* sequence/seq_len : the read bases
   hit              : vector of candidate start offsets
   min_len, max_len : same semantics as in Python
   trie             : trie root
   found_sequences  : output set of names (non-owning; don't free keys here)
*/

void add_to_counter(khash_t(strset) *found_sequences,
                    kh_counter_t *map)
{
    for (khint_t i = kh_begin(found_sequences);
         i != kh_end(found_sequences); ++i)
    {
        if (!kh_exist(found_sequences, i)) continue;
        const char *name = kh_key(found_sequences, i);  // non-owning pointer
        counter_inc(map, name);
    }
}


void find_matches(const char *sequence, size_t seq_len,
                  const u32vec_t *hit,
                  uint32_t min_len, uint32_t max_len,
                  const TrieNode *trie,
                  khash_t(strset) *found_sequences,
                  int k_mm,
                  const char *read_name,
                  const char *read_qual,
                  FILE *sam_fp)
{
    khash_t(posset) *start_pos_cache = kh_init(posset);   // dedupe absolute start positions
    kFoundVec matches; kv_init(matches);
    kvec_t(SamAlignmentRecord) sam_records; kv_init(sam_records);

    for (size_t idx = 0; idx < hit->n; ++idx) {
        uint32_t h = hit->a[idx];

        if ((size_t)h + (size_t)min_len > seq_len) continue;

        size_t qlen = (size_t)max_len;
        if ((size_t)h + qlen > seq_len) 
            continue;

        // call trie, gather matches for this slice
        clear_matches_keep_capacity(&matches);  // reuse buffer without leaking old strings
        
        char *slice = malloc(max_len + 1);
        if (!slice) /* handle OOM */;
        memcpy(slice, sequence + h, max_len);
        slice[max_len] = '\0';

        normalize_acgt(slice);
        trie_search_exact(trie, slice, min_len, k_mm, &matches);

        for (size_t m = 0; m < matches.n; ++m) {
            uint32_t abs_start = matches.a[m].pos + h;

            // if abs_start in start_pos_cache: continue
            khiter_t it = kh_get(posset, start_pos_cache, (khint_t)abs_start);
            if (it != kh_end(start_pos_cache)) continue;

            // mark visited start position
            int ret;
            it = kh_put(posset, start_pos_cache, (khint_t)abs_start, &ret);
            (void)ret;

            // found_sequences.add(m.sequence_name) (non-owning)
            khiter_t jt = kh_put(strset, found_sequences,
                                 (char*)matches.a[m].name, &ret);

            if (ret > 0) {
                // insert; reuse pointer owned by trie (don’t free in this set)
                kh_key(found_sequences, jt) = strdup(matches.a[m].name);
            }

            if (sam_fp) {
                SamAlignmentRecord rec = {
                    kh_key(found_sequences, jt),
                    abs_start,
                    strlen(matches.a[m].seq),
                    matches.a[m].mm
                };
                kv_push(SamAlignmentRecord, 0, sam_records, rec);
            }
        }

        free(slice);
    }

    if (sam_fp && sam_records.n > 0) {
        int nh = (int)sam_records.n;
        for (size_t i = 0; i < sam_records.n; ++i) {
            sam_write_alignment(sam_fp,
                                read_name,
                                sequence,
                                read_qual,
                                seq_len,
                                sam_records.a[i].ref_name,
                                sam_records.a[i].match_start,
                                sam_records.a[i].match_len,
                                sam_records.a[i].nm,
                                nh);
        }
    }

    kv_destroy(sam_records);
    mv_free(&matches);
    kh_destroy(posset, start_pos_cache);
}
