#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "kmer.h"
#include "utils.h"

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

typedef struct {
    const char *ref_name;
    size_t match_start;
    size_t match_len;
    int nm;
} SamAlignmentRecord;

struct TrieCursorState {
    const TrieNode *node;
    uint8_t edge_idx;
    uint16_t edge_off;
};

#define CURSOR_AT_NODE 255u

typedef struct {
    const char *name;
    const char *seq;
    size_t seq_len;
    int mm;
} CursorMatch;

typedef kvec_t(CursorMatch) CursorMatchVec;

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
    size_t qname_len = name_len_no_rc_suffix(read_name);
    size_t rname_len = name_len_no_rc_suffix(ref_name);

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

static inline void bitset_set(KmerBitset *index, uint64_t code) {
    if (!index || !index->words) return;
    if (code >= index->nbits) return;
    index->words[code >> 6] |= (uint64_t)1u << (code & 63u);
}

int kmer_bitset_test(const KmerBitset *index, uint64_t code) {
    if (!index || !index->words) return 0;
    if (code >= index->nbits) return 0;
    return (index->words[code >> 6] & ((uint64_t)1u << (code & 63u))) != 0;
}

static int kmer_index_add_state(KmerBitset *index,
                                const TrieNode *node,
                                uint8_t edge_idx,
                                uint16_t edge_off) {
    if (!index || !node) return -1;
    if (index->nstates == index->mstates) {
        size_t nm = index->mstates ? (index->mstates << 1) : 1024;
        TrieCursorState *p = (TrieCursorState*)realloc(index->states, nm * sizeof(*p));
        if (!p) return -1;
        index->states = p;
        index->mstates = nm;
    }
    index->states[index->nstates].node = node;
    index->states[index->nstates].edge_idx = edge_idx;
    index->states[index->nstates].edge_off = edge_off;
    return (int)index->nstates++;
}

static int trie_collect_seed_states_rec(const TrieNode *node,
                                        size_t depth,
                                        uint64_t code,
                                        size_t k,
                                        KmerBitset *index) {
    if (!node) return 0;

    if (depth == k) {
        int sid = kmer_index_add_state(index, node, CURSOR_AT_NODE, 0);
        if (sid < 0) return -1;
        index->state_idx[code] = sid;
        bitset_set(index, code);
        return 0;
    }

    for (int idx = 0; idx < ALPHABET_SIZE; ++idx) {
        const char *label = node->edge_label[idx];
        const TrieNode *child = node->child[idx];
        if (!label || !child) continue;

        size_t lablen = strlen(label);
        size_t d = depth;
        uint64_t code2 = code;

        for (size_t j = 0; j < lablen && d < k; ++j) {
            int b = nt2bits(label[j]);
            if (b < 0) return -1;
            code2 = (code2 << 2) | (uint64_t)b;
            ++d;
            if (d == k) {
                int sid = -1;
                if (j + 1 == lablen) {
                    sid = kmer_index_add_state(index, child, CURSOR_AT_NODE, 0);
                } else {
                    sid = kmer_index_add_state(index, node, (uint8_t)idx, (uint16_t)(j + 1));
                }
                if (sid < 0) return -1;
                index->state_idx[code2] = sid;
                bitset_set(index, code2);
            }
        }

        if (depth + lablen < k) {
            if (trie_collect_seed_states_rec(child, depth + lablen, code2, k, index) != 0) return -1;
        }
    }
    return 0;
}

KmerBitset *kmer_bitset_from_trie(const TrieNode *root, size_t k) {
    if (!root || k == 0 || k > 31) return NULL;

    KmerBitset *index = (KmerBitset*)calloc(1, sizeof(*index));
    if (!index) return NULL;
    index->k = k;
    index->nbits = (size_t)1ull << (2u * k);
    index->nwords = (index->nbits + 63u) >> 6;
    index->words = (uint64_t*)calloc(index->nwords, sizeof(uint64_t));
    if (!index->words) {
        free(index);
        return NULL;
    }
    index->state_idx = (int32_t*)malloc(index->nbits * sizeof(int32_t));
    if (!index->state_idx) {
        free(index->words);
        free(index);
        return NULL;
    }
    for (size_t i = 0; i < index->nbits; ++i) index->state_idx[i] = -1;

    if (trie_collect_seed_states_rec(root, 0, 0, k, index) != 0) {
        free(index->state_idx);
        free(index->states);
        free(index->words);
        free(index);
        return NULL;
    }
    return index;
}

void kmer_bitset_destroy(KmerBitset *index) {
    if (!index) return;
    free(index->states);
    free(index->state_idx);
    free(index->words);
    free(index);
}

static int seed_match_with_mm(const KmerBitset *index, uint64_t code, int seed_mm) {
    size_t k = index ? index->k : 0;
    if (!index || k == 0) return 0;
    if (kmer_bitset_test(index, code)) return 1;
    if (seed_mm <= 0) return 0;

    for (size_t pos = 0; pos < k; ++pos) {
        size_t shift = 2u * (k - 1u - pos);
        uint64_t orig = (code >> shift) & 3u;
        uint64_t cleared = code & ~((uint64_t)3u << shift);
        for (uint64_t alt = 0; alt < 4u; ++alt) {
            if (alt == orig) continue;
            uint64_t neighbor = cleared | (alt << shift);
            if (kmer_bitset_test(index, neighbor)) return 1;
        }
    }
    return 0;
}

void find_kmer_bitset(const char *s,
                      size_t s_len,
                      const KmerBitset *index,
                      int seed_mm,
                      u32vec_t *hits) {
    if (!s || !index || index->k == 0 || index->k > s_len) return;
    if (seed_mm < 0) seed_mm = 0;
    if (seed_mm > 1) seed_mm = 1;

    for (size_t i = 0; i + index->k <= s_len; ++i) {
        uint64_t code;
        if (encode_kmer(s + i, index->k, &code) < 0) continue;
        if (seed_match_with_mm(index, code, seed_mm)) {
            kv_push(uint32_t, 0, *hits, i);
        }
    }

    if (hits->n > 1) qsort(hits->a, hits->n, sizeof(uint32_t), cmp_u32_asc);
}

static void cursor_collect_matches_from_state(const KmerBitset *index,
                                              int state_id,
                                              const char *sequence,
                                              size_t seq_len,
                                              size_t start,
                                              uint32_t min_len,
                                              uint32_t max_len,
                                              int mm_initial,
                                              int k_mm,
                                              CursorMatchVec *out) {
    if (!index || !sequence || !out) return;
    if (state_id < 0 || (size_t)state_id >= index->nstates) return;
    if (start + index->k > seq_len) return;
    if (mm_initial > k_mm) return;

    const TrieCursorState *st = &index->states[state_id];
    const TrieNode *node = st->node;
    uint8_t edge_idx = st->edge_idx;
    uint16_t edge_off = st->edge_off;
    size_t p = start + index->k;
    size_t consumed = index->k;
    int mm_used = mm_initial;

    for (;;) {
        if (edge_idx == CURSOR_AT_NODE) {
            if (node->end && consumed >= min_len && consumed <= max_len) {
                CursorMatch cm = { node->name, node->seq, node->seq_len, mm_used };
                kv_push(CursorMatch, 0, *out, cm);
            }

            if (consumed >= max_len || p >= seq_len) break;
            int idx = nt2bits(sequence[p]);
            if (idx < 0) break;
            if (!node->edge_label[idx] || !node->child[idx]) break;
            edge_idx = (uint8_t)idx;
            edge_off = 0;
        }

        const char *label = node->edge_label[edge_idx];
        const TrieNode *child = node->child[edge_idx];
        if (!label || !child) break;

        size_t lablen = strlen(label);
        size_t j = edge_off;
        while (j < lablen && consumed < max_len && p < seq_len) {
            if (sequence[p] != label[j]) {
                if (++mm_used > k_mm) return;
            }
            ++p;
            ++consumed;
            ++j;
        }

        if (j < lablen) break; // ran out of room in this read slice
        node = child;
        edge_idx = CURSOR_AT_NODE;
        edge_off = 0;
    }
}

void find_matches_seeded(const char *sequence, size_t seq_len,
                         const u32vec_t *hit,
                         uint32_t min_len, uint32_t max_len,
                         const KmerBitset *seed_index,
                         int seed_mm,
                         khash_t(strset) *found_sequences,
                         int k_mm,
                         const char *read_name,
                         const char *read_qual,
                         FILE *sam_fp) {
    if (!sequence || !hit || !seed_index || !found_sequences) return;
    if (seed_mm < 0) seed_mm = 0;
    if (seed_mm > 1) seed_mm = 1;

    khash_t(posset) *start_pos_cache = kh_init(posset); // dedupe absolute start positions
    CursorMatchVec matches; kv_init(matches);
    kvec_t(SamAlignmentRecord) sam_records; kv_init(sam_records);

    for (size_t idx = 0; idx < hit->n; ++idx) {
        uint32_t h = hit->a[idx];

        if ((size_t)h + (size_t)min_len > seq_len) continue;
        if ((size_t)h + (size_t)max_len > seq_len) continue;

        uint64_t code = 0;
        if (encode_kmer(sequence + h, seed_index->k, &code) < 0) continue;

        int matched_this_start = 0;
        int exact_sid = seed_index->state_idx[code];

        // exact seed first
        if (exact_sid >= 0) {
            matches.n = 0;
            cursor_collect_matches_from_state(seed_index, exact_sid, sequence, seq_len, h,
                                              min_len, max_len, 0, k_mm, &matches);
            for (size_t m = 0; m < matches.n; ++m) {
                khiter_t it = kh_get(posset, start_pos_cache, (khint_t)h);
                if (it != kh_end(start_pos_cache)) { matched_this_start = 1; break; }

                int ret = 0;
                it = kh_put(posset, start_pos_cache, (khint_t)h, &ret);
                (void)ret;

                khiter_t jt = kh_put(strset, found_sequences, (char*)matches.a[m].name, &ret);
                if (ret > 0) kh_key(found_sequences, jt) = strdup(matches.a[m].name);

                if (sam_fp) {
                    SamAlignmentRecord rec = {
                        kh_key(found_sequences, jt),
                        h,
                        matches.a[m].seq_len,
                        matches.a[m].mm
                    };
                    kv_push(SamAlignmentRecord, 0, sam_records, rec);
                }
                matched_this_start = 1;
                break; // preserve existing one-hit-per-start behavior
            }
        }

        if (matched_this_start || seed_mm == 0) continue;

        // Hamming-1 seed neighbors
        for (size_t pos = 0; pos < seed_index->k && !matched_this_start; ++pos) {
            size_t shift = 2u * (seed_index->k - 1u - pos);
            uint64_t orig = (code >> shift) & 3u;
            uint64_t cleared = code & ~((uint64_t)3u << shift);
            for (uint64_t alt = 0; alt < 4u && !matched_this_start; ++alt) {
                if (alt == orig) continue;
                uint64_t neighbor = cleared | (alt << shift);
                int sid = seed_index->state_idx[neighbor];
                if (sid < 0) continue;

                matches.n = 0;
                cursor_collect_matches_from_state(seed_index, sid, sequence, seq_len, h,
                                                  min_len, max_len, 1, k_mm, &matches);
                for (size_t m = 0; m < matches.n; ++m) {
                    khiter_t it = kh_get(posset, start_pos_cache, (khint_t)h);
                    if (it != kh_end(start_pos_cache)) { matched_this_start = 1; break; }

                    int ret = 0;
                    it = kh_put(posset, start_pos_cache, (khint_t)h, &ret);
                    (void)ret;

                    khiter_t jt = kh_put(strset, found_sequences, (char*)matches.a[m].name, &ret);
                    if (ret > 0) kh_key(found_sequences, jt) = strdup(matches.a[m].name);

                    if (sam_fp) {
                        SamAlignmentRecord rec = {
                            kh_key(found_sequences, jt),
                            h,
                            matches.a[m].seq_len,
                            matches.a[m].mm
                        };
                        kv_push(SamAlignmentRecord, 0, sam_records, rec);
                    }
                    matched_this_start = 1;
                    break; // preserve existing one-hit-per-start behavior
                }
            }
        }
    }

    if (sam_fp && sam_records.n > 0) {
        int nh = (int)sam_records.n;
        flockfile(sam_fp);
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
        funlockfile(sam_fp);
    }

    kv_destroy(matches);
    kv_destroy(sam_records);
    kh_destroy(posset, start_pos_cache);
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
        if (!slice) continue;
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
                // insert; reuse pointer owned by trie (don't free in this set)
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
        flockfile(sam_fp);
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
        funlockfile(sam_fp);
    }

    kv_destroy(sam_records);
    mv_free(&matches);
    kh_destroy(posset, start_pos_cache);
}
