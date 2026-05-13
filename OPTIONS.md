# Command Line Interface Options

```text
./cubaTrie [options] 
./cubaTrie cut [options]
```

## `count` Options

| Option | Description | Default |
|---|---|---|
| `-i`, `--input FILE` | Input FASTQ/FA (gz supported). | `-` |
| `-r`, `--reference FILE` | Reference FASTA/FASTQ used to build the trie. | none |
| `-o`, `--output FILE` | Output CSV path. | `counts.csv` |
| `-s`, `--sam FILE` | Write SAM alignments. Use `-` for stdout. | off |
| `--soft-clip` | Use soft clipping (`S`) for SAM CIGAR and keep full read in `SEQ/QUAL` for debugging; if omitted, hard clipping (`H`) is used. | off (hard clip) |
| `--no-sam-unmapped` | Do not write unmapped reads (FLAG `4`) into SAM output. | off (unmapped reads included) |
| `-k`, `--kmer INT` | Seed k-mer length for prefilter. Allowed: `4..12`. | `8` |
| `--seed-mm INT` | Seed mismatches allowed. Allowed: `0` or `1`. | `0` |
| `--indel` | Enable indels during trie extension. If omitted, extension allows substitutions only. | off |
| `-m`, `--mismatch INT` | Total edit distance allowed during trie extension. With `--indel`, this is mismatch+indel; without `--indel`, substitutions-only. Clamped to `0..5`; if omitted, defaults to `--seed-mm`. | `0` or `1` (inherits `--seed-mm`) |
| `--exclude-multihit` | Exclude reads with more than one reference hit from final counting (SAM output unaffected). | off |
| `-a`, `--anchors STR` | Enable anchor-gated mapping. Format: `5p_adapter...3p_adapter`, `5p_adapter...`, or `...3p_adapter` (A/C/G/T only). | off |
| `--anchor-error INT` | Allowed anchor edit distance (mismatch + indel together). Range: `0..5`. | `0` |
| `ZA:Z` (SAM optional tag) | Added in anchor mode; reports insert position, orientation (`FWD`/`RC`), per-anchor status (`ast`), and anchor-level edit details including MD-like strings. Present on mapped records and on unmapped (`FLAG 4`) records in anchor mode. | emitted only when `--sam` and `--anchors` are both used |
| `-t`, `--threads UINT` | Number of worker threads. | `1` |
| `-d`, `--dup-policy MODE` | Duplicate reference handling: `error`, `warn`, `ignore`. | `error` |
| `--no-rc` | Do not add reverse complements of references. | off |
| `-v` | Verbose/debug output (repeatable). | off |
| `-h`, `--help` | Show help. | off |


```sh
# Basic run
./cubaTrie -r reference.fasta -i reads.fastq.gz -o counts.csv

# With SAM output and multithreading
./cubaTrie -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam -t 8

# With SAM soft clipping (default is hard clipping if omitted)
./cubaTrie -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam --soft-clip

# Skip unmapped records in SAM output
./cubaTrie -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam --no-sam-unmapped

# More sensitive seeding (allow 1 mismatch in seed and total edit distance 1 in extension)
./cubaTrie -r reference.fasta -i reads.fastq.gz -k 8 --seed-mm 1 -m 1 -o counts.csv

# Enable indels in extension (optional; default is substitutions-only)
./cubaTrie -r reference.fasta -i reads.fastq.gz -k 8 --seed-mm 1 --indel -m 1 -o counts.csv

# Exclude multi-hit reads from count table
./cubaTrie -r reference.fasta -i reads.fastq.gz --exclude-multihit -o counts.csv

# Anchor-gated mapping (all references must have same length, e.g. 20nt)
./cubaTrie -r reference.fasta -i reads.fastq.gz \
  -a GGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 \
  -o counts.csv

# One-sided 5' anchor: take exactly reference-length sequence immediately after anchor
./cubaTrie -r reference.fasta -i reads.fastq.gz \
  -a GGAAAGGACGAAACACCG... \
  --anchor-error 1 \
  -o counts.csv

# One-sided 3' anchor: take exactly reference-length sequence immediately before anchor
./cubaTrie -r reference.fasta -i reads.fastq.gz \
  -a ...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 \
  -o counts.csv
```


## `cut` Options

```text
./cubaTrie cut [options]
```

| Option | Description | Default |
|---|---|---|
| `-i`, `--input FILE` | Input FASTQ/FA (gz supported). | none (required) |
| `-o`, `--output FILE` | Output FASTQ path. If path ends with `.gz`, output is gzip-compressed. | optional |
| `-c`, `--count FILE` | Output CSV of kept insert counts with header `sequence,count` sorted by descending count. Use `-` for stdout. | optional |
| `-a`, `--anchors STR` | Trimming anchors: `5p_adapter...3p_adapter`, `5p_adapter...`, or `...3p_adapter` (A/C/G/T only). | none (required) |
| `--anchor-error INT` | Allowed anchor edit distance (mismatch + indel together). Range: `0..5`. | `0` |
| `-m`, `--min INT` | Minimum insert length to keep after anchor trimming. | `1` if omitted |
| `-M`, `--max INT` | Maximum insert length to keep after anchor trimming. | unbounded if omitted |
| `-t`, `--threads UINT` | Number of worker threads for trimming. | `1` |
| `--check-rc` | Retry anchor detection on reverse-complement read and reverse quality accordingly. | off |
| `-v` | Verbose/debug output (repeatable). | off |
| `-h`, `--help` | Show help. | off |

At least one of `-o/--output` or `-c/--count` must be provided.
Omitting `-m/-M` keeps any payload from length `1` onward for both one-sided and two-sided anchors.

## Example Commands

```sh
# Trim reads based on anchors; keep only inserts with length in [18, 24]
./cubaTrie cut -i reads.fastq.gz -o trimmed.fastq.gz \
  -a GGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 -m 18 -M 24 -t 8

# Trim and emit insert count table sorted by most frequent first
./cubaTrie cut -i reads.fastq.gz -o trimmed.fastq.gz -c trimmed.counts.csv \
  -a GGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 -m 18 -M 24 -t 8

# Count-only mode (skip FASTQ output)
./cubaTrie cut -i reads.fastq.gz -c trimmed.counts.csv \
  -a GGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 -m 18 -M 24 -t 8

# One-sided anchor mode with implicit length defaults:
# --min defaults to 1, --max is unbounded when omitted
./cubaTrie cut -i reads.fastq.gz -o trimmed.fastq.gz \
  -a GGAAAGGACGAAACACCG... \
  --anchor-error 1 -t 8

# Trim with reverse-complement retry
./cubaTrie cut -i reads.fastq.gz -o trimmed.fastq.gz \
  -a GGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 -m 18 -M 24 -t 8 --check-rc
```
