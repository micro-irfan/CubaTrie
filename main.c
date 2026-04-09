#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "trie.h"
#include "utils.h"
#include "ketopt.h"

/*
TODO: Allow For Mismatches
TODO: Switch to Bits
TODO: Check For Anchors - if reverse complimentary is found
TODO: Allow for threading option
*/

typedef struct {
    const char *in;          // -i/--input
    const char *ref;   // -r/--reference
    const char *out;         // -o/--output  ("-" means stdout)
    const char *sam;         // --sam FILE (or "-")
    int k;                   // -k/--kmer
    int rc; 
    int verbose;
    int mm;
    unsigned threads;
    TrieDupPolicy dup_policy;
} Options;

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] <ref.fa>\n"
        "Options:\n"
        "  -i, --input FILE        input FASTQ/FA \n"
        "  -r, --reference FILE    input FASTQ/FA \n"
        "  -o, --output FILE       output CSV (default: counts.csv)\n"
        "  -s, --sam FILE          output SAM alignments (use \"-\" for stdout)\n"
        "  -k, --kmer INT          k-mer length (default: 10)\n"
        "      --no-rc             Disable Reverse Complement (default off)\n"
        "      --indels            Include Indels (default off) (Currently Not Supported)\n"
        "  -m, --mismatch INT      Number of mismatches allowed [0] (Max 5 MMs allowed)\n"
        "  -d, --dup-policy MODE   duplicate handling: error|warn|ignore [error]\n"
        "  -t, --threads UINT      Number of threads [4] (Currently Not Supported)\n"
        "  -v                      Print Debugging Log Messages\n"
        "  -h, --help              show this help\n",
        prog);
}

static int parse_dup_policy(const char *s, TrieDupPolicy *out) {
    if (!s || !out) return -1;
    if (strcmp(s, "error") == 0)  { *out = TRIE_DUP_ERROR; return 0; }
    if (strcmp(s, "warn") == 0)   { *out = TRIE_DUP_WARN; return 0; }
    if (strcmp(s, "ignore") == 0) { *out = TRIE_DUP_IGNORE; return 0; }
    return -1;
}

int parse_args(int argc, char **argv, Options *opt, int *pos_start) {
    static ko_longopt_t longopts[] = {
        {"input",     ko_required_argument, 'i'},
        {"reference", ko_required_argument, 'r'},
        {"output",    ko_required_argument, 'o'},
        {"sam",       ko_required_argument, 's'},
        {"kmer",      ko_required_argument, 'k'},
        {"dup-policy",ko_required_argument, 'd'},
        {"no-rc",     ko_no_argument      , 301},
        {"threads",   ko_required_argument, 't'},
        {"mismatch",  ko_required_argument, 'm'},
        {"help",      ko_no_argument,       'h'},
        {NULL, 0, 0}
    };

    // Defaults
    opt->in = "-";
    opt->ref = "";
    opt->out = "counts.csv";
    opt->sam = NULL;
    opt->k = 10;
    opt->rc = 1;
    opt->verbose = 0;
    opt->threads = 1;
    opt->mm = 0;
    opt->dup_policy = TRIE_DUP_ERROR;

    ketopt_t o = KETOPT_INIT;
    int c;
    while ((c = ketopt(&o, argc, argv, 1, "i:r:o:s:k:t:m:d:vh", longopts)) >= 0) {
        switch (c) {
        case 'i': opt->in = o.arg; break;
        case 'r': opt->ref = o.arg; break;
        case 'o': opt->out = o.arg; break;
        case 's': opt->sam = o.arg; break;
        case 'k': opt->k = atoi(o.arg); break;
        case 'd':
            if (parse_dup_policy(o.arg, &opt->dup_policy) != 0) {
                fprintf(stderr, "Invalid --dup-policy '%s'. Use: error|warn|ignore\n", o.arg);
                return -5;
            }
            break;
        case 't': opt->threads = atoi(o.arg); break;
        case 'm': 
            opt->mm = atoi(o.arg); 
            if (opt->mm < 0) opt->mm = 0;
            if (opt->mm > 5) opt->mm = 5;
            break;
        case 'v': opt->verbose++; break;
        case 301: opt->rc = 0; break;
        case 'h': usage(argv[0]); return -1; // signal “showed help”
        case '?': // unknown opt
            fprintf(stderr, "Unknown option: -%c\n", o.opt ? o.opt : '?');
            usage(argv[0]); return -2;
        case ':': // missing argument
            fprintf(stderr, "Missing argument for -%c\n", o.opt ? o.opt : '?');
            usage(argv[0]); return -3;
        }
    }

    *pos_start = o.ind; // first positional argument (if any)
    if (opt->k <= 0 || opt->k > 20) {
        fprintf(stderr, "Invalid kmer len\n"); usage(argv[0]); return -4;
    }
    return 0;
}

int load_reference(const char *path, 
                   TrieNode *root, 
                   kh_counter_t *map, 
                   size_t *min_out, 
                   size_t *max_out, 
                   size_t *n_out, 
                   int add_revcomp, 
                   size_t *kmer_len,
                   TrieDupPolicy dup_policy);

int load_fastq(const char *path, 
               TrieNode *root, 
               kh_counter_t *counts,
               int kmerlen, 
               size_t *min_out, 
               size_t *max_out, 
               int k_mm,
               FILE *sam_fp);


int main(int argc, char **argv) {
    Options opt; int pos = 0;
    int pr = parse_args(argc, argv, &opt, &pos);
    if (pr < 0) return (pr == -1) ? 0 : 1;

    // Positional argument: eg. require a reference FASTA path
    // if (pos >= argc) { usage(argv[0]); return 1; }
    // const char *ref_fa = argv[pos];
    
    kh_counter_t *map = kh_init(counter);
    TrieNode *root = trie_create_node();

    // Load Reference into a Trie Structure
    size_t kmer_len = (size_t)opt.k;
    size_t min_len, max_len, n;
    int rc = load_reference(opt.ref, root, map, &min_len, &max_len, &n, opt.rc, &kmer_len, opt.dup_policy);
    if (rc != 0) {
        if (rc == 2) {
            fprintf(stderr, "Stopped: duplicate sequence found (dup-policy=error).\n");
        } else {
            fprintf(stderr, "Failed while loading reference.\n");
        }
        counter_free(map);
        trie_free_node(root);
        return 1;
    }
    printf("Inserted %zu sequences (forward & reverse count)\n", n);

    FILE *sam_fp = NULL;
    if (opt.sam) {
        sam_fp = (strcmp(opt.sam, "-") == 0) ? stdout : fopen(opt.sam, "w");
        if (!sam_fp) {
            perror("fopen");
            counter_free(map);
            trie_free_node(root);
            return 1;
        }
        trie_write_sam_header(sam_fp, root);
    }

    // Find Reference in Queries
    if (load_fastq(opt.in, root, map, opt.k, &min_len, &max_len, opt.mm, sam_fp) != 0) {
        fprintf(stderr, "Failed while loading input reads.\n");
        if (sam_fp && sam_fp != stdout) fclose(sam_fp);
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    if (sam_fp && sam_fp != stdout) fclose(sam_fp);
    counter_write_csv_sorted_collapse_rc(map, opt.out, 0);

    // Free Counter and Trie Struct
    counter_free(map);
    trie_free_node(root);
    return 0;
}

