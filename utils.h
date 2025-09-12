#ifndef UTILS_H
#define UTILS_H

#include <stdio.h> // fprintf, stderr
#include <stdlib.h> // exit
#include <string.h>
#include <stdint.h>
#include "khash.h"
#include "kvec.h"

// Define a string->size_t map family named `counter`
KHASH_MAP_INIT_STR(counter, size_t)

// (optional) alias to shorten the type name
typedef khash_t(counter) kh_counter_t;

// Functions that operate on the map
void counter_inc(kh_counter_t *m, const char *key);
// size_t counter_get(const kh_counter_t *m, const char *key);
// void counter_print(const kh_counter_t *m);
void counter_free(kh_counter_t *m);   // frees keys + destroys map

typedef struct {
    size_t pos; 
    char   *seq;  
    char   *name; 
    int    mm;
} Match ;

// typedef struct {
//     Match *a;
//     size_t n, m;
// } MatchVec;


// ---------- Found type returned by trie ----------
typedef struct {
    uint32_t start_pos;           // start pos relative to the query’s start
    const char *sequence_name;    // pointer owned by the trie (non-owning here)
} Found;

typedef kvec_t(Match) kFoundVec;

/* -------- MATCH RESULTS CONTAINER -------- */

void mv_init(kFoundVec *mv);
void mv_push(kFoundVec *mv, size_t pos, const char *name, const char *seq) ;
void mv_free(kFoundVec *mv);


static inline void normalize_acgt(char *s) {
    for (char *p = s; *p; ++p) {
        unsigned char c = (unsigned char)*p;
        if (c == 'u' || c == 'U') { *p = 'T'; continue; }
        if (c >= 'a' && c <= 'z') c -= 32;        // uppercase letters only
        // keep only A/C/G/T/N; collapse others to 'N'
        switch (c) {
            case 'A': case 'C': case 'G': case 'T': case 'N': *p = (char)c; break;
            default: *p = 'N'; break;
        }
    }
}


static inline char bits2nt(uint32_t b) {
    static const char LUT[4] = {'A','C','G','T'};
    return LUT[b & 3u];
}

// Return base at index i (0..k-1), where i=0 is the first char encoded
static inline char kmer2bit(uint64_t code, int k, size_t i) {
    if ((int)i >= k) return '?';  // out of range
    uint32_t bits = (uint32_t)((code >> (2 * (k - 1 - (int)i))) & 3u);
    return bits2nt(bits);
}


static inline uint32_t get2(const uint64_t *text, size_t i);
static inline int isN(const uint64_t *mask, size_t i);



// Map ASCII to 2-bit: A=0,C=1,G=2,T/U=3; return -1 for N/ambiguous
static inline int nt2bits(char c) {
    switch (c & ~32) { // uppercase
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': case 'U': return 3;
        default:  return -1;
    }
}

// static inline char bits2nt(int b) {
//     static const char map[4] = {'A','C','G','T'};
//     return (b>=0 && b<4) ? map[b] : 'N';
// }

// static size_t strnlen(const char *s, size_t n) {
//     size_t i = 0;
//     while (i < n && s[i] != '\0') {
//         i++;
//     }
//     return i;
// }

// char *strndup(const char *s, size_t n) {
//     size_t len = strnlen(s, n);
//     char *copy = malloc(len + 1);
//     if (copy) {
//         memcpy(copy, s, len);
//         copy[len] = '\0';
//     }
//     return copy;
// }

typedef struct {
    const char *key;
    size_t      count;
} counter_entry_t;

counter_write_csv_sorted_collapse_rc(const kh_counter_t *m, const char *path, size_t topk);

#endif /* UTILS_H */
