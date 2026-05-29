#!/bin/bash

## gRNA mapping to conserved region of EV-A71 viral genome ##
time ./bowtie-build EnteroEV71.2026.fasta EnteroEV71
time ./bowtie -a -S -f EnteroEV71 sgRNA.10.fasta bowtie.Enterovirus2026.3mm.sam -v 3 -p 4
time ./bowtie -a -S -f EnteroEV71 sgRNA.10.fasta bowtie.Enterovirus2026.1mm.sam -v 1 -p 4

time ./cubaTrie -i EnteroEV71.2026.fasta \
                -r sgRNA.10.fasta \
                -o enterovirus2026.3mm.count.csv \
                -k 4 -m 3 --seed-mm 1 -t 4 \
                -s - | samtools sort -o enterovirus2026.3mm.bam

time ./cubaTrie -i EnteroEV71.2026.fasta \
                -r sgRNA.10.fasta \
                -o enterovirus2026.1mm.count.csv \
                -k 4 -m 1 --seed-mm 1 -t 4 \
                -s - | samtools sort -o enterovirus2026.1mm.bam

time ./bowtie-build gencode.v33.transcripts.fa HumanTranscript
time ./bowtie -a -S -f HumanTranscript sgRNA.10.fasta bowtie.OffTarget.3.sam -v 3 -p 4
time ./bowtie -a -S -f HumanTranscript sgRNA.10.fasta bowtie.OffTarget.1.sam -v 1 -p 4

time ./cubaTrie -i gencode.v33.transcripts.fa \
                -r sgRNA.10.fasta \
                -o offtarget.hs.3mm.count.csv \
                -k 4 -m 3 -t 4 --seed-mm 1 --sam - --no-sam-unmapped 

time ./cubaTrie -i gencode.v33.transcripts.fa \
                -r sgRNA.10.fasta \
                -o offtarget.hs.1mm.count.csv \
                -k 4 -m 1 -t 4 --seed-mm 1 --sam - --no-sam-unmapped 

time ./bowtie-build minion_mu6.fasta minion_mu6
time ./bowtie -a -S -f minion_mu6 mu6.reference.fasta bowtie.minion.1.sam -v 0 -p 4
time ./bowtie -a -S -f minion_mu6 mu6.reference.fasta bowtie.minion.2.sam -v 1 -p 4
time ./bowtie -a -S -f minion_mu6 mu6.reference.fasta bowtie.minion.3.sam -v 2 -p 4

time ./cubaTrie -i minion_mu6.fastq.gz \
                -r mu6.reference.fasta \
                -o mu6.indel.2mm.count.csv \
                -k 4 -m 2 --indel --seed-mm 1 -t 4 \
                --sam  - | samtools sort -o mu6_alignment.indel.2mm.bam 

time ./cubaTrie -i minion_mu6.fastq.gz \
                -r mu6.reference.fasta \
                -o mu6.indel.1mm.count.csv \
                -k 4 -m 1 --indel --seed-mm 1 -t 4 \
                --sam  - | samtools sort -o mu6_alignment.indel.1mm.bam 

time ./cubaTrie -i minion_mu6.fastq.gz \
                -r mu6.reference.fasta \
                -o mu6.0mm.count.csv \
                -k 4 -m 0 --seed-mm 0 -t 4 \
                --sam  - | samtools sort -o mu6_alignment.0mm.bam 

time ./cubaTrie -i minion_mu6.fastq.gz \
                -r mu6.reference.fasta \
                -o mu6.1mm.count.csv \
                -k 4 -m 1 --seed-mm 1 -t 4 \
                --sam  - | samtools sort -o mu6_alignment.1mm.bam 

time ./cubaTrie -i minion_mu6.fastq.gz \
                -r mu6.reference.fasta \
                -o mu6.2mm.count.csv \
                -k 4 -m 2 --seed-mm 1 -t 4 \
                --sam  - | samtools sort -o mu6_alignment.2mm.bam 