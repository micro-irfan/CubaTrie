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


/* Map nucleotide char to 2-bit code: A=0, C=1, G=2, T=3, else -1 */
static inline int nt2bits(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return -1;
    }
}


static inline char bits2nt(int b) {
    static const char map[4] = {'A','C','G','T'};
    return (b>=0 && b<4) ? map[b] : 'N';
}

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
