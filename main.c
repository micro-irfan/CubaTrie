#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "utils.h"
#include "ketopt.h"
#include "readseq.h"

typedef enum {
    COMMAND_COUNT = 0,
    COMMAND_CUT = 1
} CommandMode;

typedef struct {
    const char *in;          // -i/--input
    const char *ref;         // -r/--reference
    const char *out;         // -o/--output  ("-" means stdout)
    const char *sam;         // --sam FILE (or "-")
    int sam_soft_clip;       // --soft-clip
    int sam_emit_unmapped;   // --no-sam-unmapped toggles this off
    int k;                   // -k/--kmer (seed k)
    int seed_mm;             // --seed-mm
    int rc;
    int verbose;
    int mm;
    int exclude_multihit;
    int anchor_enabled;
    char *anchor5;
    char *anchor3;
    int anchor_error;
    unsigned threads;
    TrieDupPolicy dup_policy;
} CountOptions;

typedef struct {
    const char *in;          // -i/--input
    const char *out;         // -o/--output (required)
    int out_explicit;
    int verbose;
    int anchor_enabled;
    char *anchor5;
    char *anchor3;
    int anchor_error;
    size_t min_len;          // -m/--min
    size_t max_len;          // -M/--max
    int min_set;
    int max_set;
    int check_revcomp;       // --check-rc
    unsigned threads;        // -t/--threads
} CutOptions;

static void free_anchor_pair(char **anchor5, char **anchor3) {
    if (!anchor5 || !anchor3) return;
    free(*anchor5);
    free(*anchor3);
    *anchor5 = NULL;
    *anchor3 = NULL;
}

static int has_only_acgt(const char *s) {
    if (!s || !*s) return 0;
    for (const char *p = s; *p; ++p) {
        if (*p != 'A' && *p != 'C' && *p != 'G' && *p != 'T') return 0;
    }
    return 1;
}

static char *dup_range(const char *s, size_t n) {
    char *out = (char*)malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

static int parse_anchor_pair(const char *arg, char **anchor5_out, char **anchor3_out) {
    if (!arg || !anchor5_out || !anchor3_out) return -1;
    const char *sep = strstr(arg, "...");
    if (!sep) return -1;
    if (strstr(sep + 3, "...")) return -1; // allow exactly one separator

    size_t left_len = (size_t)(sep - arg);
    size_t right_len = strlen(sep + 3);
    if (left_len == 0 && right_len == 0) return -1;

    char *left = NULL;
    char *right = NULL;

    if (left_len > 0) {
        left = dup_range(arg, left_len);
        if (!left) return -1;
        normalize_acgt(left);
        if (!has_only_acgt(left)) {
            free(left);
            return -1;
        }
    }
    if (right_len > 0) {
        right = dup_range(sep + 3, right_len);
        if (!right) {
            free(left);
            return -1;
        }
        normalize_acgt(right);
        if (!has_only_acgt(right)) {
            free(left);
            free(right);
            return -1;
        }
    }

    *anchor5_out = left;
    *anchor3_out = right;
    return 0;
}

static void usage_top(const char *prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s count [options]\n"
        "  %s cut [options]\n"
        "  %s [options]          (legacy alias for count)\n"
        "\n"
        "Run '%s count --help' or '%s cut --help' for command-specific options.\n",
        prog, prog, prog, prog, prog);
}

static void usage_count(const char *prog) {
    fprintf(stderr,
        "Usage: %s count [options]\n"
        "       %s [options]   (legacy alias)\n"
        "Options:\n"
        "  -i, --input FILE        input FASTQ/FA \n"
        "  -r, --reference FILE    reference FASTA/FASTQ (required)\n"
        "  -o, --output FILE       output CSV (default: counts.csv)\n"
        "  -s, --sam FILE          output SAM alignments (use \"-\" for stdout)\n"
        "      --soft-clip         use soft clipping (S) in SAM CIGAR (default: hard clip H)\n"
        "      --no-sam-unmapped   do not emit unmapped records (FLAG 4) into SAM\n"
        "  -k, --kmer INT          seed k-mer length for prefilter [4..12] (default: 8)\n"
        "      --seed-mm INT       allowed seed mismatches [0|1] (default: 0)\n"
        "      --exclude-multihit  do not count reads with >1 reference hit\n"
        "      --no-rc             Disable Reverse Complement (default off)\n"
        "  -a, --anchors STR       Anchors: 5p_adapter...3p_adapter, 5p_adapter..., or ...3p_adapter\n"
        "      --anchor-error INT  Allowed adapter edit distance [0..5] (mismatch+indel)\n"
        "  -m, --mismatch INT      Number of mismatches allowed [0..5] (default: --seed-mm)\n"
        "  -d, --dup-policy MODE   duplicate handling: error|warn|ignore [error]\n"
        "  -t, --threads UINT      Number of worker threads [1]\n"
        "  -v                      Print Debugging Log Messages\n"
        "  -h, --help              show this help\n",
        prog, prog);
}

static void usage_cut(const char *prog) {
    fprintf(stderr,
        "Usage: %s cut [options]\n"
        "Options:\n"
        "  -i, --input FILE        input FASTQ/FA (gz supported) (required)\n"
        "  -o, --output FILE       output FASTQ path (required)\n"
        "                          if FILE ends with .gz, output is gzip-compressed\n"
        "  -a, --anchors STR       Anchors: 5p_adapter...3p_adapter, 5p_adapter..., or ...3p_adapter\n"
        "      --anchor-error INT  Allowed adapter edit distance [0..5] (mismatch+indel)\n"
        "  -m, --min INT           minimum trimmed insert length (required)\n"
        "  -M, --max INT           maximum trimmed insert length (required)\n"
        "  -t, --threads UINT      number of worker threads [1]\n"
        "      --check-rc          retry anchor search on reverse-complement read\n"
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

static int parse_size_t_arg(const char *arg, size_t *out) {
    if (!arg || !out || !*arg) return -1;
    errno = 0;
    char *end = NULL;
    unsigned long long v = strtoull(arg, &end, 10);
    if (errno != 0 || end == arg || *end != '\0') return -1;
    if (v == 0ULL) return -1;
    *out = (size_t)v;
    return 0;
}

static int parse_count_args(int argc, char **argv, const char *prog, CountOptions *opt, int *pos_start) {
    static ko_longopt_t longopts[] = {
        {"input",     ko_required_argument, 'i'},
        {"reference", ko_required_argument, 'r'},
        {"output",    ko_required_argument, 'o'},
        {"sam",       ko_required_argument, 's'},
        {"kmer",      ko_required_argument, 'k'},
        {"soft-clip", ko_no_argument,       304},
        {"no-sam-unmapped", ko_no_argument, 305},
        {"seed-mm",   ko_required_argument, 302},
        {"exclude-multihit", ko_no_argument, 303},
        {"anchors",   ko_required_argument, 'a'},
        {"anchor-error", ko_required_argument, 306},
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
    opt->sam_soft_clip = 0;
    opt->sam_emit_unmapped = 1;
    opt->k = 8;
    opt->seed_mm = 0;
    opt->rc = 1;
    opt->verbose = 0;
    opt->threads = 1;
    opt->mm = 0;
    opt->exclude_multihit = 0;
    opt->anchor_enabled = 0;
    opt->anchor5 = NULL;
    opt->anchor3 = NULL;
    opt->anchor_error = 0;
    opt->dup_policy = TRIE_DUP_ERROR;

    ketopt_t o = KETOPT_INIT;
    int mm_explicit = 0;
    int c;
    while ((c = ketopt(&o, argc, argv, 1, "i:r:o:s:k:t:m:d:a:vh", longopts)) >= 0) {
        switch (c) {
        case 'i': opt->in = o.arg; break;
        case 'r': opt->ref = o.arg; break;
        case 'o': opt->out = o.arg; break;
        case 's': opt->sam = o.arg; break;
        case 'a': {
            char *a5 = NULL, *a3 = NULL;
            if (parse_anchor_pair(o.arg, &a5, &a3) != 0) {
                fprintf(stderr,
                        "Invalid --anchors format. Use: 5p_adapter...3p_adapter, 5p_adapter..., or ...3p_adapter (A/C/G/T only)\n");
                return -8;
            }
            free_anchor_pair(&opt->anchor5, &opt->anchor3);
            opt->anchor5 = a5;
            opt->anchor3 = a3;
            opt->anchor_enabled = 1;
            break;
        }
        case 'k': opt->k = atoi(o.arg); break;
        case 304:
            opt->sam_soft_clip = 1;
            break;
        case 305:
            opt->sam_emit_unmapped = 0;
            break;
        case 302:
            opt->seed_mm = atoi(o.arg);
            if (opt->seed_mm < 0 || opt->seed_mm > 1) {
                fprintf(stderr, "Invalid --seed-mm '%s'. Use: 0 or 1\n", o.arg);
                return -6;
            }
            break;
        case 303:
            opt->exclude_multihit = 1;
            break;
        case 306:
            opt->anchor_error = atoi(o.arg);
            if (opt->anchor_error < 0 || opt->anchor_error > 5) {
                fprintf(stderr, "Invalid --anchor-error '%s'. Use 0..5\n", o.arg);
                return -9;
            }
            break;
        case 'd':
            if (parse_dup_policy(o.arg, &opt->dup_policy) != 0) {
                fprintf(stderr, "Invalid --dup-policy '%s'. Use: error|warn|ignore\n", o.arg);
                return -5;
            }
            break;
        case 't': {
            int t = atoi(o.arg);
            opt->threads = (t > 0) ? (unsigned)t : 1u;
            break;
        }
        case 'm':
            mm_explicit = 1;
            opt->mm = atoi(o.arg);
            if (opt->mm < 0) opt->mm = 0;
            if (opt->mm > 5) opt->mm = 5;
            break;
        case 'v': opt->verbose++; break;
        case 301: opt->rc = 0; break;
        case 'h': usage_count(prog); return -1; // signal "showed help"
        case '?': // unknown opt
            fprintf(stderr, "Unknown option: -%c\n", o.opt ? o.opt : '?');
            usage_count(prog); return -2;
        case ':': // missing argument
            fprintf(stderr, "Missing argument for -%c\n", o.opt ? o.opt : '?');
            usage_count(prog); return -3;
        }
    }

    *pos_start = o.ind; // first positional argument (if any)
    if (opt->k <= 3 || opt->k > 12) {
        fprintf(stderr, "Invalid k-mer seed length. Allowed range: 4..12\n");
        usage_count(prog);
        return -4;
    }
    if (!opt->ref || opt->ref[0] == '\0') {
        fprintf(stderr, "Missing required option: -r/--reference\n");
        usage_count(prog);
        return -11;
    }
    if (!mm_explicit) {
        // If -m/--mismatch is omitted, inherit seed mismatch allowance.
        opt->mm = opt->seed_mm;
    }
    if (opt->seed_mm > opt->mm) {
        fprintf(stderr,
                "Invalid mismatch settings: --seed-mm=%d requires -m/--mismatch >= %d.\n",
                opt->seed_mm, opt->seed_mm);
        return -7;
    }
    if (opt->threads == 0) opt->threads = 1;
    if (opt->anchor_enabled && (!opt->anchor5 && !opt->anchor3)) {
        fprintf(stderr, "Anchor mode enabled but anchors were not parsed.\n");
        return -10;
    }
    return 0;
}

static int parse_cut_args(int argc, char **argv, const char *prog, CutOptions *opt, int *pos_start) {
    static ko_longopt_t longopts[] = {
        {"input",        ko_required_argument, 'i'},
        {"output",       ko_required_argument, 'o'},
        {"anchors",      ko_required_argument, 'a'},
        {"anchor-error", ko_required_argument, 306},
        {"min",          ko_required_argument, 'm'},
        {"max",          ko_required_argument, 'M'},
        {"threads",      ko_required_argument, 't'},
        {"check-rc",     ko_no_argument,       307},
        {"help",         ko_no_argument,       'h'},
        {NULL, 0, 0}
    };

    opt->in = NULL;
    opt->out = NULL;
    opt->out_explicit = 0;
    opt->verbose = 0;
    opt->anchor_enabled = 0;
    opt->anchor5 = NULL;
    opt->anchor3 = NULL;
    opt->anchor_error = 0;
    opt->min_len = 0;
    opt->max_len = 0;
    opt->min_set = 0;
    opt->max_set = 0;
    opt->check_revcomp = 0;
    opt->threads = 1;

    ketopt_t o = KETOPT_INIT;
    int c;
    while ((c = ketopt(&o, argc, argv, 1, "i:o:a:m:M:t:vh", longopts)) >= 0) {
        switch (c) {
        case 'i':
            opt->in = o.arg;
            break;
        case 'o':
            opt->out = o.arg;
            opt->out_explicit = 1;
            break;
        case 'a': {
            char *a5 = NULL, *a3 = NULL;
            if (parse_anchor_pair(o.arg, &a5, &a3) != 0) {
                fprintf(stderr,
                        "Invalid --anchors format. Use: 5p_adapter...3p_adapter, 5p_adapter..., or ...3p_adapter (A/C/G/T only)\n");
                return -8;
            }
            free_anchor_pair(&opt->anchor5, &opt->anchor3);
            opt->anchor5 = a5;
            opt->anchor3 = a3;
            opt->anchor_enabled = 1;
            break;
        }
        case 306:
            opt->anchor_error = atoi(o.arg);
            if (opt->anchor_error < 0 || opt->anchor_error > 5) {
                fprintf(stderr, "Invalid --anchor-error '%s'. Use 0..5\n", o.arg);
                return -9;
            }
            break;
        case 'm':
            if (parse_size_t_arg(o.arg, &opt->min_len) != 0) {
                fprintf(stderr, "Invalid --min value '%s'. Use a positive integer.\n", o.arg);
                return -12;
            }
            opt->min_set = 1;
            break;
        case 'M':
            if (parse_size_t_arg(o.arg, &opt->max_len) != 0) {
                fprintf(stderr, "Invalid --max value '%s'. Use a positive integer.\n", o.arg);
                return -13;
            }
            opt->max_set = 1;
            break;
        case 't': {
            int t = atoi(o.arg);
            opt->threads = (t > 0) ? (unsigned)t : 1u;
            break;
        }
        case 307:
            opt->check_revcomp = 1;
            break;
        case 'v':
            opt->verbose++;
            break;
        case 'h':
            usage_cut(prog);
            return -1;
        case '?':
            fprintf(stderr, "Unknown option: -%c\n", o.opt ? o.opt : '?');
            usage_cut(prog);
            return -2;
        case ':':
            fprintf(stderr, "Missing argument for -%c\n", o.opt ? o.opt : '?');
            usage_cut(prog);
            return -3;
        }
    }

    *pos_start = o.ind;

    if (!opt->in || opt->in[0] == '\0') {
        fprintf(stderr, "Missing required option: -i/--input\n");
        usage_cut(prog);
        return -14;
    }
    if (strcmp(opt->in, "-") == 0) {
        fprintf(stderr, "cut mode does not support stdin input; provide a FASTQ/FASTQ.GZ file path.\n");
        return -20;
    }
    if (!opt->out_explicit || !opt->out || opt->out[0] == '\0') {
        fprintf(stderr, "Missing required option: -o/--output\n");
        usage_cut(prog);
        return -15;
    }
    if (!opt->anchor_enabled || (!opt->anchor5 && !opt->anchor3)) {
        fprintf(stderr, "Missing required option: -a/--anchors\n");
        usage_cut(prog);
        return -16;
    }
    if (!opt->min_set) {
        fprintf(stderr, "Missing required option: -m/--min\n");
        usage_cut(prog);
        return -17;
    }
    if (!opt->max_set) {
        fprintf(stderr, "Missing required option: -M/--max\n");
        usage_cut(prog);
        return -18;
    }
    if (opt->min_len > opt->max_len) {
        fprintf(stderr, "Invalid length range: --min (%zu) must be <= --max (%zu)\n",
                opt->min_len, opt->max_len);
        return -19;
    }

    return 0;
}

static int detect_command_mode(int argc, char **argv, CommandMode *mode, int *arg_shift) {
    if (!mode || !arg_shift) return -1;
    *mode = COMMAND_COUNT;
    *arg_shift = 0;
    if (argc <= 1) return 0;

    if (strcmp(argv[1], "count") == 0) {
        *mode = COMMAND_COUNT;
        *arg_shift = 1;
        return 0;
    }
    if (strcmp(argv[1], "cut") == 0) {
        *mode = COMMAND_CUT;
        *arg_shift = 1;
        return 0;
    }
    if (strcmp(argv[1], "help") == 0) {
        usage_top(argv[0]);
        return 1;
    }
    if (argv[1][0] == '-') {
        // Legacy mode: treat as count.
        *mode = COMMAND_COUNT;
        *arg_shift = 0;
        return 0;
    }

    fprintf(stderr, "Unknown command: %s\n", argv[1]);
    usage_top(argv[0]);
    return -1;
}

int main(int argc, char **argv) {
    CommandMode mode = COMMAND_COUNT;
    int arg_shift = 0;
    int mode_rc = detect_command_mode(argc, argv, &mode, &arg_shift);
    if (mode_rc == 1) return 0;
    if (mode_rc < 0) return 1;

    char **sub_argv = argv + arg_shift;
    int sub_argc = argc - arg_shift;
    int pos = 0;

    if (mode == COMMAND_CUT) {
        CutOptions cut = {0};
        int pr = parse_cut_args(sub_argc, sub_argv, argv[0], &cut, &pos);
        (void)pos;
        if (pr < 0) {
            free_anchor_pair(&cut.anchor5, &cut.anchor3);
            return (pr == -1) ? 0 : 1;
        }

        AnchorConfig anchor_cfg = {0};
        anchor_cfg.enabled = 1;
        anchor_cfg.anchor5 = cut.anchor5;
        anchor_cfg.anchor3 = cut.anchor3;
        anchor_cfg.max_error = cut.anchor_error;

        int rc = cut_fastq_by_anchors(cut.in,
                                      cut.out,
                                      &anchor_cfg,
                                      cut.min_len,
                                      cut.max_len,
                                      cut.check_revcomp,
                                      cut.threads);
        free_anchor_pair(&cut.anchor5, &cut.anchor3);
        return rc == 0 ? 0 : 1;
    }

    CountOptions opt;
    int pr = parse_count_args(sub_argc, sub_argv, argv[0], &opt, &pos);
    (void)pos;
    if (pr < 0) {
        free_anchor_pair(&opt.anchor5, &opt.anchor3);
        return (pr == -1) ? 0 : 1;
    }

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
        free_anchor_pair(&opt.anchor5, &opt.anchor3);
        counter_free(map);
        trie_free_node(root);
        return 1;
    }
    if (opt.anchor_enabled && min_len != max_len) {
        fprintf(stderr,
                "Anchor mode requires all references to have the same length; got min=%zu max=%zu.\n",
                min_len, max_len);
        free_anchor_pair(&opt.anchor5, &opt.anchor3);
        counter_free(map);
        trie_free_node(root);
        return 1;
    }
    fprintf(stderr, "Inserted %zu sequences (forward & reverse count)\n", n);

    FILE *sam_fp = NULL;
    if (opt.sam) {
        sam_fp = (strcmp(opt.sam, "-") == 0) ? stdout : fopen(opt.sam, "w");
        if (!sam_fp) {
            perror("fopen");
            free_anchor_pair(&opt.anchor5, &opt.anchor3);
            counter_free(map);
            trie_free_node(root);
            return 1;
        }
        setvbuf(sam_fp, NULL, _IOFBF, 16 * 1024 * 1024);
        trie_write_sam_header(sam_fp, root);
    }

    AnchorConfig anchor_cfg = {0};
    if (opt.anchor_enabled) {
        anchor_cfg.enabled = 1;
        anchor_cfg.anchor5 = opt.anchor5;
        anchor_cfg.anchor3 = opt.anchor3;
        anchor_cfg.max_error = opt.anchor_error;
    }

    // Find Reference in Queries
    if (load_fastq(opt.in, root, map, opt.k, opt.seed_mm, &min_len, &max_len, opt.mm,
                   opt.exclude_multihit, sam_fp, opt.sam_soft_clip, opt.sam_emit_unmapped,
                   opt.threads, opt.anchor_enabled ? &anchor_cfg : NULL) != 0) {
        fprintf(stderr, "Failed while loading input reads.\n");
        if (sam_fp && sam_fp != stdout) fclose(sam_fp);
        free_anchor_pair(&opt.anchor5, &opt.anchor3);
        counter_free(map);
        trie_free_node(root);
        return 1;
    }

    if (sam_fp && sam_fp != stdout) fclose(sam_fp);
    counter_write_csv_sorted_collapse_rc(map, opt.out, 0);

    // Free Counter and Trie Struct
    free_anchor_pair(&opt.anchor5, &opt.anchor3);
    counter_free(map);
    trie_free_node(root);
    return 0;
}
