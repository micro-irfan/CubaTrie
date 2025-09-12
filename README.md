# CubaTrie

Short Reference Mapper employing Radix Trie Structure For Querying Reads or Larger Sequences. Cuba means try in the Malay Language. Essentially its Try Try* - a personal stab at writing a tool in C.

## Getting Started

```sh
git clone https://github.com/micro_irfan/cuba_trie
cd cuba_trie && make

# Short References against Merged Paired Data (Using Fastp)
./cuba_trie -r reference.fasta -i merged.fastq.gz -o test.csv -k 10

```

