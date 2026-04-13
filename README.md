# cubaTrie

cubaTrie is a Short Reference (20-30bp) Mapper employing Compressed Radix Trie Structure For Querying generic Longer Reads (>50bp) or Larger Sequences like viral sequences to find conserved target sequences. Cuba (pronounced as Chew-bah and not the country Q-ba) means to try in the Malay Language. Essentially, the tool means Try Try*.

## Getting Started

```sh
git clone https://github.com/micro_irfan/CubaTrie
cd CubaTrie && make 

# Short references against merged paired data (using fastp)
./cubaTrie -r reference.fasta -i merged.fastq.gz -o counts.csv
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

./cubaTrie -r ref.fa -i ${sample_id}.merged.fastq.gz -o ${sample_id}.counts.csv -k 8 -m 1 -t 8

# SAM File Output
./cubaTrie -r ref.fa -i ${sample_id}.merged.fastq.gz -o ${sample_id}.counts.csv -k 8 -m 1 -t 8 --sam ${sample_id}.sam
```

## CLI Options

```text
./cubaTrie [options]
```

| Option | Description | Default |
|---|---|---|
| `-i`, `--input FILE` | Input FASTQ/FA (gz supported). | `-` |
| `-r`, `--reference FILE` | Reference FASTA/FASTQ used to build the trie. | none |
| `-o`, `--output FILE` | Output CSV path. | `counts.csv` |
| `-s`, `--sam FILE` | Write SAM alignments. Use `-` for stdout. | off |
| `--soft-clip` | Use soft clipping (`S`) for SAM CIGAR and keep full read in `SEQ/QUAL` for debugging; if omitted, hard clipping (`H`) is used. | off (hard clip) |
| `--no-sam-unmapped` | Do not write unmapped reads (FLAG `4`) into SAM output. | off (unmapped reads included) |
| `-k`, `--kmer INT` | Seed k-mer length for prefilter. Allowed: `1..12`. | `8` |
| `--seed-mm INT` | Seed mismatches allowed. Allowed: `0` or `1`. | `0` |
| `-m`, `--mismatch INT` | Total mismatches allowed during extension. Clamped to `0..5`; if omitted, defaults to `--seed-mm`. | `0` or `1` (inherits `--seed-mm`) |
| `--exclude-multihit` | Exclude reads with more than one reference hit from final counting (SAM output unaffected). | off |
| `-t`, `--threads UINT` | Number of worker threads. | `1` |
| `-d`, `--dup-policy MODE` | Duplicate reference handling: `error`, `warn`, `ignore`. | `error` |
| `--no-rc` | Do not add reverse complements of references. | off |
| `-v` | Verbose/debug output (repeatable). | off |
| `-h`, `--help` | Show help. | off |

## Example Commands

```sh
# Basic run
./cubaTrie -r reference.fasta -i reads.fastq.gz -o counts.csv

# With SAM output and multithreading
./cubaTrie -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam -t 8

# With SAM soft clipping (default is hard clipping if omitted)
./cubaTrie -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam --soft-clip

# Skip unmapped records in SAM output
./cubaTrie -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam --no-sam-unmapped

# More sensitive seeding (allow 1 mismatch in seed and query)
./cubaTrie -r reference.fasta -i reads.fastq.gz -k 8 --seed-mm 1 -m 1 -o counts.csv

# Exclude multi-hit reads from count table
./cubaTrie -r reference.fasta -i reads.fastq.gz --exclude-multihit -o counts.csv

```
