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
| `-m`, `--mismatch INT` | Total mismatches allowed during extension. Clamped to `0..5`; if omitted, defaults to `--seed-mm`. | `0` or `1` (inherits `--seed-mm`) |
| `--exclude-multihit` | Exclude reads with more than one reference hit from final counting (SAM output unaffected). | off |
| `-a`, `--anchors STR` | Enable anchor-gated mapping. Format: `5p_adapter...3p_adapter`, `5p_adapter...`, or `...3p_adapter` (A/C/G/T only). | off |
| `--anchor-error INT` | Allowed anchor edit distance (mismatch + indel together). Range: `0..5`. | `0` |
| `-t`, `--threads UINT` | Number of worker threads. | `1` |
| `-d`, `--dup-policy MODE` | Duplicate reference handling: `error`, `warn`, `ignore`. | `error` |
| `--no-rc` | Do not add reverse complements of references. | off |
| `-v` | Verbose/debug output (repeatable). | off |
| `-h`, `--help` | Show help. | off |

## `cut` Options

```text
./cubaTrie cut [options]
```

| Option | Description | Default |
|---|---|---|
| `-i`, `--input FILE` | Input FASTQ/FA (gz supported). | none (required) |
| `-o`, `--output FILE` | Output FASTQ path. If path ends with `.gz`, output is gzip-compressed. | none (required) |
| `-a`, `--anchors STR` | Trimming anchors: `5p_adapter...3p_adapter`, `5p_adapter...`, or `...3p_adapter` (A/C/G/T only). | none (required) |
| `--anchor-error INT` | Allowed anchor edit distance (mismatch + indel together). Range: `0..5`. | `0` |
| `-m`, `--min INT` | Minimum insert length to keep after anchor trimming. | none (required) |
| `-M`, `--max INT` | Maximum insert length to keep after anchor trimming. | none (required) |
| `-t`, `--threads UINT` | Number of worker threads for trimming. | `1` |
| `--check-rc` | Retry anchor detection on reverse-complement read and reverse quality accordingly. | off |
| `-v` | Verbose/debug output (repeatable). | off |
| `-h`, `--help` | Show help. | off |

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

# More sensitive seeding (allow 1 mismatch in seed and query)
./cubaTrie count -r reference.fasta -i reads.fastq.gz -k 8 --seed-mm 1 -m 1 -o counts.csv

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

# Trim with reverse-complement retry
./cubaTrie cut -i reads.fastq.gz -o trimmed.fastq.gz \
  -a GGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATA \
  --anchor-error 1 -m 18 -M 24 -t 8 --check-rc

```
