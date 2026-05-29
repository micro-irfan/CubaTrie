#!/usr/bin/env python3

'''
Summarize Alignment Statistics 

Paired End Mapped to sgRNA w anchors (~110bp) using BWA / Bowtie2 (sensitive, very)
if one pair is unmapped - ignore the other pair
if one pair is non-exact - check the other pair
if one pair falls out of region - check the other pair 
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

class ReadWrapper:
    def __init__(self, read):
        self.read = read
        self.custom_attr = None
        self.max_gap = -1
        self.total_operation = -1

class Bam():
    
    START_SGRNA = 22 
    anchor5_len = 22
    anchon3_len = 120
    sgrna_length_max = 20

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


    def __init__(self, bam_file, reference, sample=""):
        self.bam = bam_file
        self.sample = bam_file.split('.')[0] if bam_file else sample
        
        
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
        self.gene_rpm = {}
        self.gene_log1p_rpm = {}


    def open_bam(self, software = 'bwa', allow_recheck=True, only_exact=True):
        cache = {}
        
        ## Counter For Reads that are multimapped
        self.stats['multimapped'] = 0

        ## Counter For Reads that are invalid - ie. Reference Based Lenght is not correct
        self.stats['invalid_ref_start'] = 0
        
        ## Counter For Not Proper Pair Reads - ie. different alignments
        self.stats['not_proper_pair'] = 0

        ## Counter for when mate is unmapped
        self.stats['unmapped_mate'] = 0

        ## Counter to count exact mapping
        self.stats['exact_match'] = 0
        self.stats['rescued_exact_match'] = 0
        self.stats['alignment_count'] = 0

        if software == 'bowtie2': 
            self.stats['CP'] = 0

        self.operations = {k:0 for k in range(20)}
        self.max_gap = {k:0 for k in range(20)}

        write_file = open(f"{self.sample}-error.maxgap.csv", 'w')

        with pysam.AlignmentFile(self.bam, "rb") as bamfile:
            for read in bamfile.fetch(until_eof=True): 
                self.stats['alignment_count'] += 1
                r = ReadWrapper(read)
                read_id = r.read.query_name
                to_remove = False

                if r.read.is_secondary or r.read.is_supplementary: 
                    self.stats['multimapped'] += 1
                    continue
                
                ## Only One Read Has to be exact match
                r.to_remove = False
                if read_id in cache:
                    ## Check if exact match - if False - check the paired read
                    to_remove = True
                    if cache[read_id].to_remove:
                        del cache[read_id]                 
                        continue
                else:
                    ### Add Only One Count Accordingly
                    cache[read_id] = r
                    self.stats['read'] += 1

                    if not r.read.is_proper_pair: 
                        self.stats['not_proper_pair'] += 1 
            
                    # Unmapped Case
                    if r.read.mate_is_unmapped: 
                        self.stats['unmapped_mate'] += 1 
                        cache[read_id].to_remove = True
                        continue

                    if r.read.is_unmapped: 
                        if r.read.mate_is_mapped: 
                            self.stats['unmapped_mate'] += 1 
                        else:
                            self.stats['unmapped'] += 1
                        cache[read_id].to_remove = True
                        continue
                    
                    if software == 'bowtie2':
                        tag = r.read.get_tag('YT')
                        try:
                            self.stats[tag] += 1
                        except KeyError:
                            self.stats[tag] = 1

                        if tag != 'CP' : 
                            cache[read_id].to_remove = True
                            continue

                    ## Check For Secondary Reads Alignment
                    if r.read.has_tag('SA'): 
                        sa_tag = r.read.get_tag('SA')
                        to_continue = False
                        for c, aln in enumerate(sa_tag.strip(';').split(';')):
                            if not c % 2: continue
                            aln = aln.split(',')
                            ref_start_pos = int(aln[1])
                            cigar_tuples = parse_cigar_string(aln[3])
                            sec_read_len = sum([l for l, op in cigar_tuples if op not in 'SH'])
                            ref_end_pos = ref_start_pos + sec_read_len

                            ## Potential Multimapped
                            if ref_start_pos < self.START_SGRNA and ref_end_pos > self.START_SGRNA + self.sgrna_length_max: 
                                cache[read_id].to_remove = True
                                to_continue = True
                                break

                        if to_continue: 
                            continue

                ## Check sgRNA sequence
                ref_start_pos = r.read.reference_start
                if ref_start_pos >= self.START_SGRNA: 
                    if to_remove:
                        del cache[read_id]
                    else:
                        self.stats['invalid_ref_start'] += 1
                        if not allow_recheck:
                            cache[read_id].to_remove = True
                    continue
            
                cigar_tuples = r.read.cigartuples
                reference_id = bamfile.get_reference_name(r.read.reference_id)
                md_tag = r.read.get_tag('MD')
                sgrna = self.reference.sequence[reference_id]
                adjusted_cigar = self.adjust_for_insertion(cigar_tuples, md_tag, ref_start_pos, sgrna)
                max_gap, total_operation = self._count_opened_gap(adjusted_cigar)
                
                if max_gap > 8:
                    ## For Debugging Purposes
                    to_write = f'{read_id},{reference_id},{max_gap},{adjusted_cigar},{total_operation}\n'
                    self.write_error(to_write, write_file)

                added = False
                if total_operation == 0:
                    self.stats['exact_match'] += 1
                    self.operations[total_operation] += 1
                    self.max_gap[max_gap] += 1
                    cache[read_id].to_remove = True
                    if to_remove:
                        self.stats['rescued_exact_match'] += 1
                    added = True

                if to_remove and not added: 
                    if cache[read_id].total_operation > 0:
                        max_gap = min(cache[read_id].max_gap, max_gap)
                        total_operation = min(cache[read_id].total_operation, total_operation)

                    try:
                        self.max_gap[max_gap] += 1
                    except KeyError:
                        self.max_gap[max_gap] = 1
                    
                    try:
                        self.operations[total_operation] += 1
                    except KeyError:
                        self.operations[total_operation] = 1
                    
                if to_remove:
                    del cache[read_id]
                    continue

                cache[read_id].max_gap = max_gap
                cache[read_id].total_operation = total_operation
                
                if only_exact and total_operation != 0: continue

                self.gene_count[reference_id] += 1

        self.stats['mapped'] = self.stats['read']-self.stats['unmapped']-self.stats['unmapped_mate']-self.stats['invalid_ref_start']
        # if software == 'bwa':
        #     self.stats['mapped'] = self.stats['read']-self.stats['unmapped']-self.stats['unmapped_mate']-self.stats['invalid_ref_start']
        # if software == 'bowtie2':
        #     self.stats['mapped'] = self.stats['CP'] - self.stats['invalid_ref_start']

        if self.stats['mapped'] == 0:
            return

        print (self.stats)
        
        mapped_percent = (self.stats['mapped']/self.stats['read']) * 100
        self.stats['mapped%'] += mapped_percent
        self.stats['unpaired'] += len(cache)
        self.stats['exact_match%'] = (self.stats['exact_match']/self.stats['mapped']) * 100
        self.stats['rescued_exact_match%'] = (self.stats['rescued_exact_match']/self.stats['mapped']) * 100
        
        print (self.stats)
        print (self.max_gap)
        print (self.operations)

        write_file.close()
    
    def write_error(self, to_write, write_file):
        write_file.write(to_write)

    # Reference http://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files
    def _expand_cigar(self, cigartuples):
        """Expand CIGAR string into per-reference-position operation string."""
        ## Remove SoftClipiing
        return ''.join(self.cigar_dict[op] * length for op, length in cigartuples if self.cigar_dict[op] != 'S')

    def _parse_md(self, md_string):
        """Parse MD string into per-reference-position alignment string."""
        result = []
        md_items = re.findall(r'(\d+)|(\^[A-Z]+)|([A-Z])', md_string)
        for num, deletion, mismatch in md_items:
            if num:
                result.extend(['M'] * int(num))
            elif deletion:
                result.extend(['D'] * (len(deletion)-1))  # minus 1 for the ^ char
            elif mismatch:
                result.append('X')
        return ''.join(result)

    def get_relative_start_end(self, reference_start, sgrna):

        sgrna_length = len(sgrna) - self.anchor5_len - self.anchon3_len
        offset = 0
        sgrna_start = self.START_SGRNA + offset
        sgrna_end = sgrna_start + sgrna_length
    
        relative_sgrna_start = sgrna_start - reference_start
        assert relative_sgrna_start >= 0
        relative_sgrna_end = sgrna_end - reference_start

        return relative_sgrna_start, relative_sgrna_end

    def adjust_for_insertion(self, cigartuples, md_tag, reference_start, sgrna):
        cigar_str = self._expand_cigar(cigartuples)
        md_string = self._parse_md(md_tag)
        adjusted_cigar = ''
        insert_count = 0
        relative_sgrna_start, relative_sgrna_end = self.get_relative_start_end(reference_start, sgrna)
        for c, operation in enumerate(cigar_str): 
            if operation == 'S': break
            if c > relative_sgrna_end + 1: break
            # print (c, insert_count, c-insert_count, len(md_string), cigar_str, adjusted_cigar, relative_sgrna_end)
            if operation == 'I': 
                insert_count += 1  
                adjusted_cigar += operation
                relative_sgrna_end += 1
                continue

            if operation == md_string[c-insert_count]: 
                adjusted_cigar += operation
            else:
                md_operation = md_string[c-insert_count]
                adjusted_cigar += md_operation

        adjusted_cigar = adjusted_cigar[relative_sgrna_start:relative_sgrna_end]
        return adjusted_cigar
    
    def _count_opened_gap(self, cigar):
        max_gap = 1
        gap = 0
        total_operation = 0
        for operation in cigar:
            if operation != 'M':
                total_operation += 1
                gap += 1

            if operation == 'M' and gap>max_gap:
                max_gap = gap
                gap = 0

        if total_operation == 0: 
            max_gap -= 1 # Counts should match
        return max_gap,total_operation
    
def parse_cigar_string(cigarstring):
    """Convert CIGAR string to list of (length, operation) tuples."""
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigarstring)]

color_dict = {
    'bwa' : 'navy',
    'bowtie2' : 'black'
}

tool = ['bwa', 'bowtie2']


def write_output(bam):
    sample_id = bam.sample
    with open(f'Mapping-Stats.{sample_id}.csv', 'w') as write_file:
        write_file.write(f'Statistics,{sample_id}\n')
        for header, value in bam.stats.items():
            write_file.write(f'{header},{value}\n')

    max_value = max(list(bam.max_gap.keys()) + list(bam.operations.keys()))
    with open(f'Operations-log.{sample_id}.csv', 'w') as write_file:
        write_file.write('Count,Operations,MaxGap\n')
        for i in range(max_value):
            write_file.write(f'{i},{bam.operations.get(i,0)},{bam.max_gap.get(i,0)}\n')
            

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

    parser.add_argument('--tool', default = "bwa",
                        help="Support Alignment Tools - Bowtie2 and BWA (default)")

    args = parser.parse_args()
    return args

def main():
    args = get_args()
    bam = args.bam.split(',')
    ref = Reference(args.ref)
    tool = args.tool.split('_')[0]
    analyze_bam(ref, bam, tool)

if __name__ == "__main__":
    main()
