#!/usr/bin/env python3
"""
FASTQ Sequence Search using Trie (Prefix Tree) Algorithm

This script efficiently searches for multiple short sequences (from FASTA file)
in large FASTQ files using a trie data structure and reports read counts.
"""

import gzip
import argparse
import sys
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict, Counter


class TrieNode:
    """Node in the trie data structure."""

    def __init__(self):
        self.children = {}  # Dictionary to store child nodes
        self.is_end_of_sequence = False
        self.sequence_name = None  # Store sequence name from FASTA
        self.seq = None

class Found:

    def __init__(self, start_pos, seq_id, matched_seq):
        self.start_pos = start_pos
        self.sequence_name = seq_id
        self.matched_seq = matched_seq

class SequenceTrie:
    """Trie implementation optimized for DNA/RNA sequence searching."""

    def __init__(self):
        self.root = TrieNode()
        self.sequences = {}  # Map sequence_name to original sequence

    def insert(self, seq: str, sequence_name: str):
        node = self.root
        seq = seq.upper().strip()  # Normalize to uppercase
        self.sequences[sequence_name] = seq

        skip = False
        for pos, nucleotide in enumerate(seq):
            if nucleotide not in node.children:
                newNode = TrieNode()
                newNode.seq = seq[pos:]
                newNode.is_end_of_sequence = True
                newNode.sequence_name = sequence_name
                node.children[nucleotide] = newNode
                skip = True
                break
            
            ## Create New Nodes
            if node.children[nucleotide].seq:
                seq_in_node = node.children[nucleotide].seq
                node.children[nucleotide].seq = None
                node.children[nucleotide].is_end_of_sequence = False
                tmp_seq_name = node.children[nucleotide].sequence_name
                node.children[nucleotide].sequence_name = None

                seq_to_compare = seq[pos:]
                min_len = min(len(seq_to_compare), len(seq_in_node))
                for i in range(min_len):
                    if seq_to_compare[i] != seq_in_node[i] and i != min_len:
                        nt1 = seq_to_compare[i]
                        nt2 = seq_in_node[i]
                        seq1 = seq_to_compare[i:]
                        seq2 = seq_in_node[i:]

                        newNode1 = TrieNode()   
                        newNode1.is_end_of_sequence = True
                        newNode1.sequence_name = sequence_name
                        newNode1.seq = seq1
                        node.children[nt1] = newNode1

                        newNode2 = TrieNode()
                        newNode2.is_end_of_sequence = True
                        newNode2.sequence_name = tmp_seq_name
                        newNode2.seq = seq2
                        node.children[nt2] = newNode2
                        break

                    nucleotide = seq_in_node[i]
                    node.children[nucleotide] = TrieNode()
                    node = node.children[nucleotide]
                skip = True
                break
            
            ## Continue to the next nucleotide
            node = node.children[nucleotide]
        
        if not skip:
            node.is_end_of_sequence = True
            node.sequence_name = sequence_name


    def search_exact(self, text: str, min_len: int):
        """
        Search for exact matches in the text.
        Returns list of (position, sequence_name, matched_sequence) tuples.
        """
        matches = []
        text = text.upper()

        start_pos = 0
        for start_pos in range(len(text) - min_len + 1):
            node = self.root
            current_pos = start_pos
            while (current_pos < len(text) and text[current_pos] in node.children):
                node = node.children[text[current_pos]]
                if node.is_end_of_sequence:
                    seq = node.seq
                    seq_len = len(seq)
                    seq_to_check = text[current_pos:current_pos+seq_len]
                    if seq == seq_to_check:
                        matched_seq = self.sequences[node.sequence_name]
                        matches.append(Found(start_pos, node.sequence_name, matched_seq))

                current_pos += 1

        return matches

    def find_prefix(self, text: str):
        node = self.root

        for c, nt in enumerate(text):
            if nt not in node.children:
                ## Check Sequence
                if node.is_end_of_sequence:
                    seq_to_check = node.seq[1:len(text)-c+1]
                    seq_remainder = text[c:]
                    return seq_to_check == seq_remainder
                else:
                    return False
            node = node.children[nt]

        return True




class FASTQSearcher:
    """Main class for searching sequences in FASTQ files."""

    def __init__(self, fasta_file: str):
        """
        Initialize with a FASTA file containing target sequences.

        Args:
            fasta_file: Path to FASTA file with target sequences
        """
        self.trie = SequenceTrie()
        self.target_sequences = {}
        self.sequence_names = []

        # Read FASTA file and build the trie
        self._load_fasta(fasta_file)

    def _load_fasta(self, fasta_file: str):
        """Load sequences from FASTA file and build trie."""
        opener = gzip.open if fasta_file.endswith('.gz') else open
        mode = 'rt' if fasta_file.endswith('.gz') else 'r'

        current_name = None
        current_sequence = []

        print(f"Loading target sequences from {fasta_file}...")

        self.min_len = float("inf")
        self.max_len = 0
        try:
            with opener(fasta_file, mode) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_name and current_sequence:
                            seq = ''.join(current_sequence)
                            self.target_sequences[current_name] = seq
                            self.sequence_names.append(current_name)
                            self.trie.insert(seq, current_name)

                            current_name = f'{current_name}-rev'
                            seq = self.reverse_complement(seq)
                            self.target_sequences[current_name] = seq
                            self.sequence_names.append(current_name)
                            self.trie.insert(seq, current_name)

                            seq_len = len(seq)
                            if seq_len < self.min_len:
                                self.min_len = seq_len
                            if seq_len > self.max_len:
                                self.max_len = seq_len

                        # Start new sequence
                        current_name = line[1:]  # Remove '>'
                        current_sequence = []
                    elif line and current_name:
                        current_sequence.append(line)

                # Don't forget the last sequence
                if current_name and current_sequence:
                    seq = ''.join(current_sequence)
                    self.target_sequences[current_name] = seq
                    self.sequence_names.append(current_name)
                    self.trie.insert(seq, current_name)

                    current_name = f'{current_name}-rev'
                    seq = self.reverse_complement(seq)
                    self.target_sequences[current_name] = seq
                    self.sequence_names.append(current_name)
                    self.trie.insert(seq, current_name)

                    seq_len = len(seq)
                    if seq_len < self.min_len:
                        self.min_len = seq_len
                    if seq_len > self.max_len:
                        self.max_len = seq_len

            print(f"Loaded {len(self.target_sequences)} target sequences.")

        except FileNotFoundError:
            print(f"Error: Could not find FASTA file '{fasta_file}'")
            sys.exit(1)
        except Exception as e:
            print(f"Error reading FASTA file: {e}")
            sys.exit(1)

    def parse_fastq(self, filename: str):
        """
        Generator to parse FASTQ files (handles both regular and gzipped files).

        Yields:
            Tuple of (header, sequence, plus_line, quality)
        """
        opener = gzip.open if filename.endswith('.gz') else open
        mode = 'rt' if filename.endswith('.gz') else 'r'

        try:
            with opener(filename, mode) as f:
                while True:
                    header = f.readline().strip()
                    if not header:
                        break
                    sequence = f.readline().strip()
                    plus_line = f.readline().strip()
                    quality = f.readline().strip()

                    yield header, sequence, plus_line, quality
        except FileNotFoundError:
            print(f"Error: Could not find FASTQ file '{filename}'")
            sys.exit(1)
        except Exception as e:
            print(f"Error reading FASTQ file: {e}")
            sys.exit(1)

    def break_into_kmers(self, text, window=10):
        kmer = {}
        for i in range(len(text)-window+1):
            seq = text[i:i+window]
            if seq not in kmer.keys():
                kmer[seq] = [i]
            else:
                kmer[seq].append(i)

        return kmer
    
    def find_kmer(self, kmer_dict):
        hit = []
        for kmer, position in kmer_dict.items():
            if kmer in self.kmer_hit:
                hit += position
                continue

            if kmer in self.kmer_search:
                continue

            if self.trie.find_prefix(kmer):
                hit += position
                self.kmer_hit.add(kmer)
            
            ## Add that its been checked
            self.kmer_search.add(kmer)
        
        hit.sort()
        return hit

    def search_fastq(self, 
                     fastq_file: str, 
                     allow_mismatches: bool = False,
                     max_mismatches: int = 1) -> Dict[str, int]:
        """
        Search for target sequences in a FASTQ file and count reads.

        Args:
            fastq_file: Path to FASTQ file
            allow_mismatches: Whether to allow mismatches in search
            max_mismatches: Maximum number of mismatches allowed
            search_reverse: Whether to search reverse complement

        Returns:
            Dictionary mapping sequence names to read counts
        """
        read_counts = Counter()
        read_id_cache = dict()
        total_reads = 0
        reads_with_matches = 0

        # Initialize counts for all target sequences
        for seq_name in self.sequence_names:
            if '-rev' in seq_name: continue
            read_counts[seq_name] = 0
            read_id_cache[seq_name] = set()

        print(f"Searching in {fastq_file}...")

        self.kmer_hit = set()
        self.kmer_search = set()
        for header, sequence, plus_line, quality in self.parse_fastq(fastq_file):
            total_reads += 1
            read_has_match = False

            if total_reads % 1000000 == 0:
                print(f"Processed {total_reads} reads...")

            # Track which sequences were found in this read (avoid double counting)
            found_sequences = set()

            # Create Kmer Dict
            kmer_dict = self.break_into_kmers(sequence)

            # find kmer prefix - similar to seeding
            hit = self.find_kmer(kmer_dict)
            start_pos_cache = set()
            for h in hit: 
                if h+self.min_len > len(sequence): continue
                    
                query = sequence[h:] if h+self.max_len > len(sequence) else sequence[h:h+self.max_len]
                # if len(query) < self.min_len: continue

                match = self.trie.search_exact(query, self.min_len)
                if not match: continue 
                
                for m in match:
                    m.start_pos += h
                    if m.start_pos in start_pos_cache: continue
                    start_pos_cache.add(m.start_pos)
                    found_sequences.add(m.sequence_name)

            # Count each sequence found in this read once
            for seq_name in found_sequences:
                seq_name = seq_name.replace('-rev','') if 'rev' in seq_name else seq_name
                read_counts[seq_name] += 1
                read_id_cache[seq_name].add(header)
                read_has_match = True

            if read_has_match:
                reads_with_matches += 1

        print(f"Finished processing {total_reads} reads.")
        print(f"Reads with matches: {reads_with_matches} ({reads_with_matches/total_reads*100:.2f}%)")

        return dict(read_counts), read_id_cache

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """Generate reverse complement of a DNA sequence."""
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement_map.get(base, base) for base in reversed(sequence.upper()))

    def write_results(self, results: Dict[str, int], output_file: str = None):
        """Write search results to file or stdout."""
        output = []
        output.append("Sequence_Name\tSequence\tRead_Count")

        # Sort by read count (descending)
        sorted_results = sorted(results.items(), key=lambda x: x[1], reverse=True)

        for seq_name, count in sorted_results:
            sequence = self.target_sequences[seq_name]
            output.append(f"{seq_name}\t{sequence}\t{count}")

        output_text = '\n'.join(output)

        if output_file:
            try:
                with open(output_file, 'w') as f:
                    f.write(output_text)
                print(f"Results written to {output_file}")
            except Exception as e:
                print(f"Error writing to output file: {e}")
                print("Results:")
                print(output_text)
        else:
            print("\n=== RESULTS ===")
            print(output_text)


def write_cache(cache, output):
    with open(output, 'w') as write_file:
        write_file.write("Seq_Id,Ref_Id\n")
        for seq_id, read_id in cache.items():
            read_id = list(read_id)
            to_write = f'{seq_id},{";".join(read_id)}\n'
            write_file.write(to_write)

def main():
    """Main function with command-line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Search for target sequences (from FASTA) in FASTQ files using trie algorithm",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic exact search
  python fastq_search.py -f targets.fasta -q reads.fastq

  # Search with mismatches and save to file
  python fastq_search.py -f targets.fasta -q reads.fastq.gz -m 1 -o results.txt

  # Search forward strand only
  python fastq_search.py -f targets.fasta -q reads.fastq --no-reverse
        """
    )

    # Required arguments
    parser.add_argument('-f', '--fasta', required=True,
                       help='FASTA file containing target sequences')
    parser.add_argument('-q', '--fastq', required=True,
                       help='FASTQ file to search in (can be gzipped)')

    # Optional arguments
    parser.add_argument('-o', '--output',
                       help='Output file for results (default: stdout)')
    parser.add_argument('-m', '--mismatches', type=int, default=0,
                       help='Maximum number of mismatches allowed (default: 0 for exact match)')
    parser.add_argument('--no-reverse', action='store_true',
                       help='Do not search reverse complement strand')

    args = parser.parse_args()


    try:
        # Initialize searcher
        searcher = FASTQSearcher(args.fasta)

        # Perform search
        allow_mismatches = args.mismatches > 0
        search_reverse = not args.no_reverse

        if allow_mismatches:
            print(f"Performing search with up to {args.mismatches} mismatches...")
        else:
            print("Performing exact sequence search...")

        if search_reverse:
            print("Searching both forward and reverse complement strands...")
        else:
            print("Searching forward strand only...")

        results, cache = searcher.search_fastq(
            args.fastq,
            allow_mismatches=allow_mismatches,
            max_mismatches=args.mismatches,
        )

        # Write results
        searcher.write_results(results, args.output)
        write_cache(cache, f'readID-{args.output}')

        total_hits = sum(results.values())
        sequences_with_hits = sum(1 for count in results.values() if count > 0)
        print(f"\nSummary:")
        print(f"Total sequences searched: {len(results)}")
        print(f"Sequences with hits: {sequences_with_hits}")
        print(f"Total read hits: {total_hits}")

    except KeyboardInterrupt:
        original_print("\nSearch interrupted by user.")
        sys.exit(1)
    except Exception as e:
        original_print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

# time python3 reference_matching_adjust.py -f rbp_sgrna_short.fasta -q RHH9441.merged.fastq.gz -o test.adjusted.csv 