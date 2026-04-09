# CubaTrie

Short Reference Mapper employing Radix Trie Structure For Querying Reads or Larger Sequences. Cuba means try in the Malay Language. Essentially its Try Try* - a personal stab at writing a tool in C.

## Getting Started

```sh
git clone https://github.com/micro_irfan/cuba_trie
cd cuba_trie && make

# Short references against merged paired data (using fastp)
./cuba_trie -r reference.fasta -i merged.fastq.gz -o counts.csv
```

## CLI Options

```text
./cuba_trie [options]
```

| Option | Description | Default |
|---|---|---|
| `-i`, `--input FILE` | Input FASTQ/FA (gz supported). | `-` |
| `-r`, `--reference FILE` | Reference FASTA/FASTQ used to build the trie. | none |
| `-o`, `--output FILE` | Output CSV path. | `counts.csv` |
| `-s`, `--sam FILE` | Write SAM alignments. Use `-` for stdout. | off |
| `-k`, `--kmer INT` | Seed k-mer length for prefilter. Allowed: `1..12`. | `8` |
| `--seed-mm INT` | Seed mismatches allowed. Allowed: `0` or `1`. | `0` |
| `-m`, `--mismatch INT` | Total mismatches allowed during extension. Clamped to `0..5`. | `0` |
| `-t`, `--threads UINT` | Number of worker threads. | `1` |
| `-d`, `--dup-policy MODE` | Duplicate reference handling: `error`, `warn`, `ignore`. | `error` |
| `--no-rc` | Do not add reverse complements of references. | off |
| `-v` | Verbose/debug output (repeatable). | off |
| `-h`, `--help` | Show help. | off |

## Example Commands

```sh
# Basic run
./cuba_trie -r reference.fasta -i reads.fastq.gz -o counts.csv

# With SAM output and multithreading
./cuba_trie -r reference.fasta -i reads.fastq.gz -o counts.csv -s output.sam -t 8

# More sensitive seeding (allow 1 mismatch in seed and query)
./cuba_trie -r reference.fasta -i reads.fastq.gz -k 8 --seed-mm 1 -m 1 -o counts.csv

```

