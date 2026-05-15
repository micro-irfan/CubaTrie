# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.0.0]
* Added a cut mode which trims flanking region with `cubaTrie cut`.
* Added an inbuilt cut for count mode with `--anchor`. 
* Added a custom tag SAM tag `ZA:Z` to preserve flanking region mapping information.
* Allow For Indels during Trie Traversal using DFS with `--indel`.
* Updated Read Scanning to a Rolling Hash Style via shifting the k-mer bits to the next base
* Added a Golden Unit Test to test for CRISPR CasRfx gRNA targeting EV-A71 Genome. 
* Added Scripts used for Publications to `benchmarking`

## [v0.0.3]
* Replaced hash-set seed prefilter with a bitset seed index.
* Added seed mismatch option `--seed-mm` with allowed values `0|1`.
* Updated default seed k-mer length to a min `8` with max `12`.
* Added direct-access seed-to-trie cursor table for continuation from seed depth.
* Added multi-threaded read processing in `load_fastq` via `-t/--threads`.

## [v0.0.2]
* Create a Kmer HashSet to check if prefix exist in the Compressed Radix Trie.
* Allow Mismatches (Substituition) during Trie Traversal.
* Option to Output SAM File. NM and NH tags are included.
* Added an option to ignore reverse complementary sequences in short sequence reference.
* Added unit test to /tests directory.

## [v0.0.1]
* Initial release of cuba trie - beta version, only exact match, unoptimized.
