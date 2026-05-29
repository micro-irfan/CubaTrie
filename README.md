# cubaTrie

cubaTrie is a Short Reference (20-30bp) Mapper employing Compressed Radix Trie Structure For Querying generic Longer Reads (>50bp) or Larger Sequences like viral sequences to find conserved target sequences. CubaTrie has been tested for quantification of gRNA for CRISPR Cas Pooled Screening and CRISPR Cas Therapeutics/Diagnostics, and quantification of intramolecular homologous recombination between tandem repetitive elements using Long Read Sequencing Technology [workflow](README.md#Suggested-Workflows).  

CubaTrie has two modes - cut and count - which can be performed in a single run. CubaTrie cut trims 5' and 3' flanking region. CubaTrie count maps super short references to reads or sequences of similar or longer length. All CLI options are included on [OPTIONS](OPTIONS.md). We also provide the mapping information for 5' and 3' flanking region under the custom tag `ZA:Z` for futher debugging. More information on Anchors can be found on [ANCHORS](ANCHOR.md).   

Cuba (pronounced as Chew-bah and not the country Q-ba) means to try in the Malay Language. Essentially, the tool means Try Trie*, since Trie Data Structures are not commonly found in Bioinformatics. Since this started out as a side project, I documented the process, the DSA implemented, current limitations, and ad-hoc findings, on my personal pages (Links to be added later). 

## Getting Started

```sh
git clone https://github.com/micro_irfan/CubaTrie
cd CubaTrie && make 

./cubaTrie -r ref.fa -i ${sample_id}.R1.fastq.gz -o ${sample_id}.counts.csv -k 8 -t 8

# SAM File Output
./cubaTrie -r ref.fa -i ${sample_id}.merged.fastq.gz -o ${sample_id}.counts.csv -k 8 -m 1 -t 8 --sam ${sample_id}.sam
```

## Suggested Workflows

```sh
# CRISPR Cas Pooled
./cubaTrie -r sgRNA.fasta -i ${sample_id}.R1.fastq.gz -o ${sample_id}.count.csv -k 10 -m 0 -d warn -t 8 --exclude-multihit --sam - -a ATTTTCAATTTAACGTCG...GTTTTAGAGCTAGAAATA --anchor-error 2 --soft-clip | samtools sort -o ${sample_id}.bam

# CRISPR Cas Therapeutics Or Diagnostics - gRNA on-target
./cubaTrie -i EnteroEV71.2026.fasta -r sgRNA.10.fasta -o enterovirus2026.3mm.count.csv -k 4 -m 3 --seed-mm 1 -t 4 -s - | samtools sort -o enterovirus2026.3mm.bam

# CRISPR Cas Therapeutics Or Diagnostics - gRNA off-target (Human Host)
./cubaTrie -i gencode.v33.transcripts.fa -r sgRNA.10.fasta -o offtarget.hs.count.csv -k 4 -m 3 -t 4 --seed-mm 1 --sam - --no-sam-unmapped 

# sgRNA mu6 cassette Quantification for LRS
./cubaTrie -i minion_mu6.fastq.gz -r mu6.reference.fasta -o mu6.count.csv -k 4 -m 3 --indel --seed-mm 1 -t 4 --sam  - | samtools sort -o mu6_alignment.bam 

# CubaTrie Cut Function
./cubaTrie cut -i ${sample_id}.R1.fastq.gz -o ${sample_id}.R1.trimmed.fastq.gz -a ATTTTCAATTTAACGTCG...GTTTTAGAGCTAGAAATA --anchor-error 3 -m 20 -M 20 -t 4
```

This project is licensed under the GNU General Public License v3.0.  
See [LICENSE](LICENSE) for the full text.
