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


