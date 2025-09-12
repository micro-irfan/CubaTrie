#include <stdio.h>
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
    int k;                   // -k/--kmer
    int rc; 
    int verbose;
    int mm;
    unsigned threads;
} Options;

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] <ref.fa>\n"
        "Options:\n"
        "  -i, --input FILE        input FASTQ/FA \n"
        "  -r, --reference FILE    input FASTQ/FA \n"
        "  -o, --output FILE       output CSV (default: counts.csv)\n"
        "  -k, --kmer INT          k-mer length (default: 10)\n"
        "      --no-rc             Disable Reverse Complement (default off)\n"
        "      --indels            Include Indels (default off) (Currently Not Supported)\n"
        "  -m, --mismatch INT      Number of mismatches allowed [0] (Max 5 MMs allowed)\n"
        "  -t, --threads UINT      Number of threads [4] (Currently Not Supported)\n"
        "  -v                      Print Debugging Log Messages\n"
        "  -h, --help              show this help\n",
        prog);
}

int parse_args(int argc, char **argv, Options *opt, int *pos_start) {
    static ko_longopt_t longopts[] = {
        {"input",     ko_required_argument, 'i'},
        {"reference", ko_required_argument, 'r'},
        {"output",    ko_required_argument, 'o'},
        {"kmer",      ko_required_argument, 'k'},
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
    opt->k = 10;
    opt->rc = 1;
    opt->verbose = 0;
    opt->threads = 1;
    opt->mm = 0;

    ketopt_t o = KETOPT_INIT;
    int c;
    while ((c = ketopt(&o, argc, argv, 1, "i:r:o:k:t:m:h", longopts)) >= 0) {
        switch (c) {
        case 'i': opt->in = o.arg; break;
        case 'r': opt->ref = o.arg; break;
        case 'o': opt->out = o.arg; break;
        case 'k': opt->k = atoi(o.arg); break;
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

void load_reference(const char *path, TrieNode *root, kh_counter_t *map, 
                      size_t *min_out, size_t *max_out, size_t *n_out, int add_revcomp, size_t *kmer_len);

void load_fastq(const char *path, TrieNode *root, kh_counter_t *counts,
                int kmerlen, size_t *min_out, size_t *max_out, int k_mm);


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
    size_t min_len, max_len, n;
    load_reference(opt.ref, root, map, &min_len, &max_len, &n, opt.rc, opt.k);
    printf("Inserted %zu sequences (forward & reverse count)\n", n);

    // Find Reference in Queries
    load_fastq(opt.in, root, map, opt.k, &min_len, &max_len, opt.mm);

    counter_write_csv_sorted_collapse_rc(map, opt.out, 0);

    // Free Counter and Trie Struct
    counter_free(map);
    trie_free_node(root);
    return 0;
}

