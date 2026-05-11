# cubaTrie

cubaTrie is a Short Reference (20-30bp) Mapper employing Compressed Radix Trie Structure For Querying generic Longer Reads (>50bp) or Larger Sequences like viral sequences to find conserved target sequences. Cuba (pronounced as Chew-bah and not the country Q-ba) means to try in the Malay Language. Essentially, the tool means Try Try*.

## Getting Started

```sh
git clone https://github.com/micro_irfan/CubaTrie
cd CubaTrie && make 

# Short references against merged paired data (using fastp)
./cubaTrie count -r reference.fasta -i merged.fastq.gz -o counts.csv
```

## Suggested Workflows For Paired-End Data

Users can run dedup and merge paired-end data before running cubaTrie

```sh
# Run Deduplication Step (If Required) with Fastp
fastp \
    --in1 ${read1} \
    --in2 ${read2} \
    --out1 ${sample_id}_R1.trimmed.fastq.gz \
    --out2 ${sample_id}_R2.trimmed.fastq.gz \
    --json ${sample_id}_fastp.dedup.json \
    --html ${sample_id}_fastp.dedup.html \
    --dedup \
    --thread 8

# Run Merge Step with Fastp
fastp \
    --in1 ${sample_id}_R1.trimmed.fastq.gz \
    --in2 ${sample_id}_R2.trimmed.fastq.gz \
    --merge \
    --merged_out ${sample_id}.merged.fastq.gz \
    --out1 ${sample_id}_R1.merged.fastq.gz \
    --out2 ${sample_id}_R2.merged.fastq.gz \
    --json ${sample_id}_fastp.merge.json \
    --html ${sample_id}_fastp.merge.html \
    --thread 8

./cubaTrie count -r ref.fa -i ${sample_id}.merged.fastq.gz -o ${sample_id}.counts.csv -k 8 -m 1 -t 8

# SAM File Output
./cubaTrie count -r ref.fa -i ${sample_id}.merged.fastq.gz -o ${sample_id}.counts.csv -k 8 -m 1 -t 8 --sam ${sample_id}.sam
```

## CLI

```text
./cubaTrie count [options]
./cubaTrie cut [options]
./cubaTrie [options]    # legacy alias for count
```

## `count` Options

```text
./cubaTrie count [options]
```

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

## `ZA:Z` SAM Tag Format

`ZA:Z` is emitted when both `--sam` and `--anchors` are enabled:
- on mapped SAM records in anchor mode
- on unmapped (`FLAG 4`) SAM records if anchors were detected but short-reference mapping failed
- on unmapped (`FLAG 4`) SAM records in two-sided mode when full pairing fails but a single start anchor (`a5` or `a3rc`) is confidently detected (`partial=1`)

Expected format (full paired-anchor detection):

```text
ZA:Z:ori=<FWD|RC>;ins=<insert_start_1based>,<insert_len>;ast=5:<pass|fail|na>,3:<pass|fail|na>[;a5=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>][;a3=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>]
```

Expected format (two-sided mode partial start-anchor fallback on unmapped reads):

```text
ZA:Z:ori=<FWD|RC>;partial=1;ast=5:<pass|fail|na>,3:<pass|fail|na>;a5=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>
ZA:Z:ori=<FWD|RC>;partial=1;ast=5:<pass|fail|na>,3:<pass|fail|na>;a3rc=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>
ZA:Z:ori=<FWD|RC>;partial=1;ast=5:<pass|fail|na>,3:<pass|fail|na>;a5=<...>[;a3f=<start_1based>-<end_1based>,len=<int>,ed=<int>,md=<MD_like>]
ZA:Z:ori=<FWD|RC>;partial=1;ast=5:<pass|fail|na>,3:<pass|fail|na>;a3rc=<...>[;a5f=<start_1based>-<end_1based>,len=<int>,ed=<int>,md=<MD_like>]
```

Expected format (two-sided mode unmapped diagnostics):

```text
ZA:Z:ori=<FWD|RC>;ins=<insert_start_1based>,<observed_insert_len>;exp=<reference_len>;reason=<insert_len_lt|insert_len_gt|anchor_ambiguous|anchor_window_rejected>;ast=5:<pass|fail|na>,3:<pass|fail|na>
ZA:Z:ori=<FWD|RC>;ins=<insert_start_1based>,<observed_insert_len>;exp=<reference_len>;reason=<...>;ast=5:<pass|fail|na>,3:<pass|fail|na>;a5=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>;a3=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>
ZA:Z:reason=anchor_not_found;ast=5:<pass|fail|na>,3:<pass|fail|na>;a5f=<start_1based>-<end_1based>,len=<int>,ed=<int>,md=<MD_like>[;a3f=<start_1based>-<end_1based>,len=<int>,ed=<int>,md=<MD_like>]
```

Field meaning:

- `ori`: orientation of anchor detection in the read (`FWD` or `RC`).
- `ins`: extracted insert window used for mapping.
- `a5`: 5' anchor match window in read coordinates.
- `a3`: 3' anchor match window in read coordinates.
- `a3rc`: reverse-complement 3' anchor start hit (reported only in two-sided fallback mode with `ori=RC`).
- `partial=1`: indicates two-sided anchor mode fallback where only a start anchor (`a5` or `a3rc`) was confidently detected.
- `exp`: expected insert/reference length used by `count` mode.
- `reason`: diagnostic reason for unmapped anchor attempts.
- `ast`: anchor status summary where `5` and `3` are each `pass`, `fail`, or `na`.
- `na` means that anchor was not searched in that diagnostic path (or not configured).
- `a5f` / `a3f`: best-effort failed-anchor diagnostics with position range, matched segment length, edit distance, and MD-like string.
- `ed`: anchor edit distance used by anchor matching.
- `md`: MD-like string comparing read anchor segment against expected anchor sequence (`A/C/G/T`, digits for match runs, `^` for deletions from read relative to anchor).

Notes:

- Coordinates in `ins`, `a5`, and `a3` are 1-based.
- Coordinates in `a3rc` are also 1-based.
- `a5` and/or `a3` are present depending on anchor mode (both-sided or one-sided).
- In `count` mode, two-sided anchors (`a5...a3`) must bracket an insert whose length matches the reference length.
- In `count` mode, one-sided anchors keep fixed-length behavior (insert length follows reference length).
- In two-sided mode, if full pairing fails but a single start anchor is confidently found, unmapped SAM may include `partial=1` with either `a5` or `a3rc`.
- Unmapped SAM in anchor mode includes a diagnostic `reason` when extraction fails.
- `reason=anchor_not_found` indicates no confident anchor window was recoverable; when available, failed anchors include `a5f`/`a3f` metadata for debugging.

Examples:

```text
ZA:Z:ori=FWD;ins=31,20;ast=5:pass,3:pass;a5=13-30,ed=0,md=18;a3=51-68,ed=1,md=7A10
ZA:Z:ori=RC;ins=31,20;ast=5:pass,3:pass;a5=55-72,ed=1,md=7A10;a3=9-26,ed=0,md=18
ZA:Z:ori=FWD;ins=19,20;ast=5:pass,3:na;a5=1-18,ed=0,md=18
ZA:Z:ori=FWD;ins=31,22;exp=20;reason=insert_len_gt;ast=5:pass,3:pass;a5=13-30,ed=0,md=18;a3=53-70,ed=1,md=7A10
ZA:Z:ori=FWD;partial=1;ast=5:pass,3:fail;a5=1-18,ed=0,md=18;a3f=45-62,len=18,ed=3,md=4T5A7
ZA:Z:ori=RC;partial=1;ast=5:fail,3:pass;a3rc=3-20,ed=1,md=7A10;a5f=51-68,len=18,ed=2,md=5C12
ZA:Z:reason=anchor_not_found;ast=5:fail,3:na;a5f=2-5,len=4,ed=2,md=0GG0
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
| `-m`, `--min INT` | Minimum insert length to keep after anchor trimming. | none (required) |
| `-M`, `--max INT` | Maximum insert length to keep after anchor trimming. | none (required) |
| `-t`, `--threads UINT` | Number of worker threads for trimming. | `1` |
| `--check-rc` | Retry anchor detection on reverse-complement read and reverse quality accordingly. | off |
| `-v` | Verbose/debug output (repeatable). | off |
| `-h`, `--help` | Show help. | off |

At least one of `-o/--output` or `-c/--count` must be provided.

## Example Commands

```sh
# Basic run
./cubaTrie count -r reference.fasta -i reads.fastq.gz -o counts.csv

# With SAM output and multithreading
./cubaTrie count -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam -t 8

# With SAM soft clipping (default is hard clipping if omitted)
./cubaTrie count -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam --soft-clip

# Skip unmapped records in SAM output
./cubaTrie count -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam --no-sam-unmapped

# More sensitive seeding (allow 1 mismatch in seed and total edit distance 1 in extension)
./cubaTrie count -r reference.fasta -i reads.fastq.gz -k 8 --seed-mm 1 -m 1 -o counts.csv

# Enable indels in extension (optional; default is substitutions-only)
./cubaTrie count -r reference.fasta -i reads.fastq.gz -k 8 --seed-mm 1 --indel -m 1 -o counts.csv

# Exclude multi-hit reads from count table
./cubaTrie count -r reference.fasta -i reads.fastq.gz --exclude-multihit -o counts.csv

# Anchor-gated mapping (all references must have same length, e.g. 20nt)
./cubaTrie count -r reference.fasta -i reads.fastq.gz \
  -a GGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 \
  -o counts.csv

# One-sided 5' anchor: take exactly reference-length sequence immediately after anchor
./cubaTrie count -r reference.fasta -i reads.fastq.gz \
  -a GGAAAGGACGAAACACCG... \
  --anchor-error 1 \
  -o counts.csv

# One-sided 3' anchor: take exactly reference-length sequence immediately before anchor
./cubaTrie count -r reference.fasta -i reads.fastq.gz \
  -a ...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 \
  -o counts.csv

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

# Trim with reverse-complement retry
./cubaTrie cut -i reads.fastq.gz -o trimmed.fastq.gz \
  -a GGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 -m 18 -M 24 -t 8 --check-rc

```
