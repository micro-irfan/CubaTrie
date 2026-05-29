#!/usr/bin/env python3

'''
Summarize Alignment Statistics Heuristically

Paired End Mapped to unique sgRNA (19-20bp) using BWA
if one pair is unmapped - ignore the other pair
if one pair is non-exact - skip the other pair
'''


import pysam
import re
import argparse

class Reference():

    ## Missing in Control File
    missing = ["BPF_1", "tagRFP_1" , "tagRFP_2"]

    def __init__(self, reference):
        genes = set()
        self.sequence = {}
        with open(reference, 'r') as f:
            for line in f:
                if line.startswith(">"): 
                    reference_id = line.strip('\n').replace('>','')
                    genes.add(reference_id)
                seq = line.strip('\n') 
                self.sequence[reference_id] = seq

        print (f"Number of Genes: {len(genes)}")
        genes_sorted = list(genes)
        genes_sorted.sort()
        self.genes = genes_sorted

class Bam():
    
    cigar_dict = {
            0: "M", # match/mismatch
            1: "I", 
            2: "D", 
            3: "N", 
            4: "S", 
            5: "H", 
            6: "P", 
            7: "=", # match
            8: "X", # mismatch
            9: "B"
        }

    def __init__(self, bamFile, reference, sample=""):
        self.bam = bamFile
        self.sample = bamFile.split('.')[0] if bamFile else sample
        
        self.stats = { 
            'unmapped' : 0,
            'unpaired' : 0,
            'read' : 0,
            'mapped%' : 0
        }

        self.reference = reference
        genes = self.reference.genes

        self.gene_read_id = {k:set() for k in genes}
        self.gene_count = {k:0 for k in genes}


    def open_bam(self, software = 'bwa'):
        cache = set()
        
        ## Counter For Reads that are multimapped
        self.stats['multimapped'] = 0

        ## Counter For Not Proper Pair Reads - ie. different alignments
        self.stats['not_proper_pair'] = 0

        ## Counter for when mate is unmapped
        self.stats['unmapped_mate'] = 0

        ## Counter to count exact mapping
        self.stats['exact_match'] = 0
        self.stats['alignment_count'] = 0

        with pysam.AlignmentFile(self.bam, "rb") as bamfile:
            for read in bamfile.fetch(until_eof=True): 
                self.stats['alignment_count'] += 1
                read_id = read.query_name

                if read.is_secondary or read.is_supplementary: 
                    self.stats['multimapped'] += 1
                    continue

                if read_id in cache:                 
                    cache.remove(read_id)
                    continue
                
                cache.add(read_id)
                self.stats['read'] += 1

                # Unmapped Case
                if read.mate_is_unmapped: 
                    self.stats['unmapped_mate'] += 1 
                    continue

                if read.is_unmapped: 
                    if read.mate_is_mapped: 
                        self.stats['unmapped_mate'] += 1 
                    else:
                        self.stats['unmapped'] += 1
                    continue
                
                ## Check sgRNA sequence
                cigar_tuples = read.cigartuples
                reference_id = bamfile.get_reference_name(read.reference_id)
                sgrna_length = len(self.reference.sequence[reference_id])
                for op, length in cigar_tuples:
                    if self.cigar_dict[op] == 'M' and length == sgrna_length:
                        md_tag = read.get_tag('MD') # Ensure theres no mismatch
                        if md_tag == str(sgrna_length):
                            self.stats['exact_match'] += 1
                
                if not read.is_proper_pair: 
                    self.stats['not_proper_pair'] += 1 
            
                self.gene_count[reference_id] += 1

                
        self.stats['mapped'] = self.stats['read']-self.stats['unmapped']-self.stats['unmapped_mate']
        if self.stats['mapped'] == 0:
            return
        
        print (self.stats)

        mapped_percent = (self.stats['mapped']/self.stats['read']) * 100
        self.stats['mapped%'] += mapped_percent
        self.stats['unpaired'] += len(cache)
        self.stats['exact_match%'] = (self.stats['exact_match']/self.stats['mapped']) * 100
        print (self.stats)

def write_output(bam):
    sample_id = bam.sample
    with open(f'Mapping-Stats-short.{sample_id}.csv', 'w') as write_file:
        write_file.write(f'Statistics,{sample_id}\n')
        for header, value in bam.stats.items():
            write_file.write(f'{header},{value}\n')

def analyze_bam(ref, bam_files, tool):
    for b in bam_files:
        bam = Bam(b, ref)
        bam.open_bam(software=tool)
        write_output(bam)

def get_args():
    parser = argparse.ArgumentParser(
        prog = 'Gene Count For sgRNA Expression',
        description = "Summary Statistics and Gene Count Post Alignment - BWA (Preferred) or Bowtie2"
    )

    parser.add_argument('-b', '--bam', dest = "bam", required = True,
                        help="String of Bam Files, comma seperated")
    
    parser.add_argument('-r', '--reference', dest = 'ref', required = True,
                        help="Reference Fasta File")
    
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    bam = args.bam.split(',')
    ref = Reference(args.ref)
    analyze_bam(ref, bam, 'bwa')

if __name__ == "__main__":
    main()