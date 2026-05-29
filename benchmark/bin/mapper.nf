#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BOWTIE2_BUILD {
    label "bowtie2"
    tag "$sample_id"

    input:
    path(reference)

    output:
    path("*.bt2")

    script:
    """
    bowtie2-build ${reference} rbp_sgrna
    """
}


process BOWTIE2_ALIGN_S {
    label "bowtie2"
    tag "$sample_id"
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(reference)
    val(tool)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    int sorting_threads = Math.min((task.cpus / 4) as int, 3)
    int mapping_threads = task.cpus - sorting_threads

    """
    bowtie2 --sensitive-local \
                -p ${mapping_threads} \
                -k 1 \
                -x rbp_sgrna \
                -1 ${reads[0]} -2 ${reads[1]} | \
    samtools sort -@ ${sorting_threads} -o ${sample_id}.bam
    samtools index ${sample_id}.bam 
    """
}

process BOWTIE2_ALIGN_E2E {
    label "bowtie2"
    tag "$sample_id"
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(reference)
    val(tool)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    int sorting_threads = Math.min((task.cpus / 4) as int, 3)
    int mapping_threads = task.cpus - sorting_threads
    """
    bowtie2 --end-to-end \
            -p ${mapping_threads} \
            -k 1 \
            -x rbp_sgrna \
            -1 ${reads[0]} -2 ${reads[1]} | \
    samtools sort -@ ${sorting_threads} -o ${sample_id}.bam
    samtools index ${sample_id}.bam 
    """
}


process BOWTIE2_ALIGN_SINGLE_READ_S {
    label "bowtie2"
    tag "$sample_id"
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(reference)
    val(tool)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    int sorting_threads = Math.min((task.cpus / 4) as int, 3)
    int mapping_threads = task.cpus - sorting_threads
    """
    bowtie2 -x rbp_sgrna \
            -U ${reads} \
            --sensitive-local \
            -k 1 \
            -p ${mapping_threads} | \
    samtools sort -@ ${sorting_threads} -o ${sample_id}.bam
    samtools index ${sample_id}.bam 
    """
}

process BOWTIE2_ALIGN_SINGLE_READ {
    label "bowtie2"
    tag "$sample_id"
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(reference)
    val(tool)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    int sorting_threads = Math.min((task.cpus / 4) as int, 3)
    int mapping_threads = task.cpus - sorting_threads
    """
    bowtie2 -x rbp_sgrna \
            -U ${reads} \
            --end-to-end \
            -N 0 \
            -p ${mapping_threads} | \
    samtools sort -@ ${sorting_threads} -o ${sample_id}.bam
    samtools index ${sample_id}.bam 
    """
}

process BWA_INDEX {
    label 'bwa'
    tag "$sample_id"

    input:
    path(reference)

    output:
    path "${reference.getName()}*"

    script:
    """
    bwa index ${reference}
    """
}

process BWA_SINGLE {
    label 'bwa'
    tag "$sample_id"
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(reference)
    val(tool)
    maxRetries 1 

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    int sorting_threads = Math.min((task.cpus / 4) as int, 3)
    int mapping_threads = task.cpus - sorting_threads
    """
    bwa mem -t ${mapping_threads} ${reference.baseName[0]} ${reads} | \
    samtools sort -@ ${sorting_threads} -o ${sample_id}.bam
    samtools index ${sample_id}.bam 
    """
}

process BWA {
    label 'bwa'
    tag "$sample_id"
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(reference)
    val(tool)
    maxRetries 1 

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    int sorting_threads = Math.min((task.cpus / 4) as int, 2)
    int mapping_threads = task.cpus - sorting_threads
    """
    bwa mem -t ${mapping_threads} ${reference.baseName[0]} ${reads[0]} ${reads[1]} | \
    samtools sort -@ ${sorting_threads} -o ${sample_id}.bam
    samtools index ${sample_id}.bam 
    """
}

process SAM_TO_BAM {
    maxForks 2
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'copy'
    label "samtools_8"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(sam)
    val(tool)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    """
    samtools view -Sb ${sam} | samtools sort -@ ${task.cpus} -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """
}





process CUBATRIE_CLEAN {
    tag "$sample_id"
    publishDir "${params.outdir}/cubatrie_${method}_${n_cpus}/${sample_id}", mode: 'copy'
    cpus {n_cpus}
    label "cubaTrie_anchors"

    input:
        tuple val(sample_id), path(reads)
        path(reference)
        val(n_cpus)
        val(method)

    output:
        path "${sample_id}.counts.csv"

    script:
    """
    cubaTrie -r ${reference} -i ${reads} -o ${sample_id}.counts.csv -d warn -t ${n_cpus} --exclude-multihit
    """
}

process CUBATRIE_CLEAN_ANCHORS {
    tag "$sample_id"
    publishDir "${params.outdir}/cubatrie_${method}_${n_cpus}/${sample_id}", mode: 'copy'
    cpus {n_cpus}
    label "cubaTrie_anchors"

    input:
        tuple val(sample_id), path(reads)
        path(reference)
        val(n_cpus)
        val(method)
        val(flank5)
        val(flank3)

    output:
        path "${sample_id}.counts.csv"

    script:
    """
    cubaTrie -r ${reference} \
             -i ${reads} \
             -o ${sample_id}.counts.csv \
             -d warn \
             -t ${n_cpus} \
             --exclude-multihit \
             -a ${flank5}...${flank3} \
             --anchor-error 3
    """
}

process CUBATRIE_BAM_ANCHORS {
    publishDir "${params.outdir}/cubatrie_bam_${method}_${n_cpus}/${sample_id}", mode: 'copy'
    cpus {n_cpus}
    tag "$sample_id"
    label "cubaTrie_anchors"

    input:
        tuple val(sample_id), path(reads)
        path(reference)
        val(n_cpus)
        val(method)
        val(flank5)
        val(flank3)

    output:
        tuple path("${sample_id}.counts.csv"), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    int sorting_threads = Math.min((task.cpus / 4) as int, 2)
    int mapping_threads = task.cpus - sorting_threads
    """
    cubaTrie -r ${reference} \
             -i ${reads} \
             -o ${sample_id}.counts.csv \
             -d warn \
             -t ${mapping_threads} \
             --exclude-multihit \
             -a ${flank5}...${flank3} \
             --anchor-error 3 \
             --soft-clip \
             --sam - | \
    samtools sort -@ ${sorting_threads} -o ${sample_id}.bam
    samtools index ${sample_id}.bam 
    """
}

process CUBATRIE_CLEAN_EXTRA_ARGS {
    tag "$sample_id"
    publishDir {
        def extra_tag = extra_args
            .trim()
            .replaceAll(/\s+/, '_')      // spaces -> _
            .replaceAll(/[^A-Za-z0-9._-]/, '')  // remove unsafe chars
        "${params.outdir}/cubatrie_${n_cpus}_${extra_tag}/${sample_id}"
    }, mode: 'copy'
    cpus {n_cpus}

    input:
        tuple val(sample_id), path(reads)
        path(reference)
        val(n_cpus)
        val(extra_args)

    output:
        path "${sample_id}.counts.csv"

    script:
    """
    ${params.cubatrie} -r ${reference} -i ${reads} -o ${sample_id}.counts.csv -d warn -t ${n_cpus} ${extra_args}
    """
}

process CUBATRIE_BAM {
    publishDir "${params.outdir}/cubatrie_bam_${method}_${n_cpus}/${sample_id}", mode: 'copy'
    cpus {n_cpus}
    tag "$sample_id"
    label "cubaTrie_anchors"

    input:
        tuple val(sample_id), path(reads)
        path(reference)
        val(n_cpus)
        val(method)

    output:
        tuple path("${sample_id}.counts.csv"), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    int sorting_threads = Math.min((task.cpus / 4) as int, 2)
    int mapping_threads = task.cpus - sorting_threads
    """
    cubaTrie -r ${reference} -i ${reads} -o ${sample_id}.counts.csv -d warn -t ${mapping_threads} --exclude-multihit --soft-clip --sam - | \
    samtools sort -@ ${sorting_threads} -o ${sample_id}.bam
    samtools index ${sample_id}.bam 
    """
}

process CUBATRIE_SAM {
    publishDir "${params.outdir}/cubatrie_sam_${method}_${n_cpus}/${sample_id}", mode: 'copy'
    cpus {n_cpus}
    tag "$sample_id"
    label "cubaTrie"

    input:
        tuple val(sample_id), path(reads)
        path(reference)
        val(n_cpus)
        val(method)

    output:
        tuple path("${sample_id}.counts.csv"), path("${sample_id}.sam")

    script:
    int sorting_threads = Math.min((task.cpus / 4) as int, 2)
    int mapping_threads = task.cpus - sorting_threads
    """
    cubaTrie -r ${reference} -i ${reads} -o ${sample_id}.counts.csv -d warn -t ${mapping_threads} --exclude-multihit --sam ${sample_id}.sam
    """
}

process FASTP_DEDUP {
    publishDir "${params.outdir}/dedup/${sample_id}"
    label "fastp"
    tag "$sample_id"

    input:
    tuple val(sample_id), file(read1), file(read2)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.dedup.json", emit: json_report
    path "${sample_id}_fastp.dedup.html", emit: html_report

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}_R1.trimmed.fastq.gz \
        --out2 ${sample_id}_R2.trimmed.fastq.gz \
        --json ${sample_id}_fastp.dedup.json \
        --html ${sample_id}_fastp.dedup.html \
        --dedup \
        --thread ${task.cpus} 
    """
}


process FASTP_TRIMMED {
    publishDir "${params.outdir}/trimmed/${sample_id}", mode: 'copy'
    label "fastp"
    tag "$sample_id"

    input:
    tuple val(sample_id), file(read1), file(read2)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.trimmed.json", emit: json_report
    path "${sample_id}_fastp.trimmed.html", emit: html_report

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}_R1.trimmed.fastq.gz \
        --out2 ${sample_id}_R2.trimmed.fastq.gz \
        --json ${sample_id}_fastp.trimmed.json \
        --html ${sample_id}_fastp.trimmed.html \
        --thread ${task.cpus} 
    """
}

process FASTP_TRIMMED_SINGLE {
    publishDir "${params.outdir}/trimmed/${sample_id}", mode: 'copy'
    label "fastp"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.trimmed.json", emit: json_report
    path "${sample_id}_fastp.trimmed.html", emit: html_report

    script:
    """
    fastp \
        --in1 ${read1} \
        --out1 ${sample_id}.trimmed.fastq.gz \
        --json ${sample_id}_fastp.trimmed.json \
        --html ${sample_id}_fastp.trimmed.html \
        --thread ${task.cpus}
    """
}

process FASTP_MERGE {
    publishDir "${params.outdir}/${method}/${sample_id}"
    label "fastp"
    tag "$sample_id"

    input:
    tuple val(sample_id), file(read1), file(read2)
    val(method)

    output:
    tuple val(sample_id), path("${sample_id}.merged.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.merge.json", emit: json_report
    path "${sample_id}_fastp.merge.html", emit: html_report

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --merge \
        --merged_out ${sample_id}.merged.fastq.gz \
        --out1 ${sample_id}_R1.merged.fastq.gz \
        --out2 ${sample_id}_R2.merged.fastq.gz \
        --json ${sample_id}_fastp.merge.json \
        --html ${sample_id}_fastp.merge.html \
        --thread ${task.cpus} 
    """
}


process CUTADAPT {
    publishDir "${params.outdir}/cutadapt/${sample_id}", mode: 'copy'
    label "cutadapt"
    tag "$sample_id"

    input:
    tuple val(sample_id), file(read1), file(read2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), emit: trimmed_reads

    script:
    """
    cutadapt \
        -a GGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATA \
        -O 18 -e 0.3 \
        --discard-untrimmed \
        --revcomp \
        --max-n 0 \
        -m 20 -M 20 \
        -j ${task.cpus} \
        -o "${sample_id}_R1.trimmed.fastq.gz" \
        ${read1}
    """
}

process CUTADAPT_FLEX {
    publishDir "${params.outdir}/cutadapt/${sample_id}", mode: 'copy'
    label "cutadapt"
    tag "$sample_id"

    input:
    tuple val(sample_id), file(read)
    val(flank5)
    val(flank3)

    output:
    tuple val(sample_id), path("${sample_id}.cutadapt.fastq.gz"), emit: trimmed_reads

    script:
    """
    cutadapt \
        -a ${flank5}...${flank3} \
        -O 18 -e 0.3 \
        --discard-untrimmed \
        --revcomp \
        --max-n 0 \
        -m 20 -M 20 \
        -j ${task.cpus} \
        -o "${sample_id}.cutadapt.fastq.gz" \
        ${read}
    """
}

process COUNT_READS {
    tag "$sample_id"
    publishDir "${params.outdir}/${method}/${sample_id}", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
        val(method)

    output:
        tuple val(sample_id), path("${sample_id}.counts.txt")

    script:
    """
    zcat ${reads} | wc -l > ${sample_id}.counts.txt
    """

}



process CUBATRIE_CUT {
    tag "$sample_id"
    publishDir "${params.outdir}/cubatrie_cut/${sample_id}", mode: 'copy'
    cpus 8
    label "cubaTrie_anchors"

    input:
        tuple val(sample_id), path(reads)
        val(flank5)
        val(flank3)

    output:
        tuple val(sample_id), path("${sample_id}.cut.fastq.gz"), emit: trimmed_reads

    script:
    """
    cubaTrie cut -i ${reads} -o ${sample_id}.cut.fastq.gz \
                 -a ${flank5}...${flank3} \
                 --anchor-error 3 \
                 -m 20 -M 20 -t ${task.cpus}
    """
}




process CUBATRIEPY {
    publishDir "${params.outdir}/cubatriepy/${sample_id}", mode: 'copy'
    label 'pysam'
    tag "$sample_id"

    input:
    tuple val(sample_id), file(reads)
    path(reference)

    output:
    tuple path("${sample_id}.multimapper.txt"), emit: mapping_stats

    script:
    """
    python ${baseDir}/bin/short_multi_mapper.kmer.py -f ${reference} -q ${reads} -o ${sample_id}.multimapper.txt
    """
}


process SUMMARIZE_STATS_SHORT {
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'copy'
    label 'pysam'
    tag "$sample_id"

    input:
    tuple val(sample_id), file(bam), file(bamIdx)
    path(reference)
    val(tool)

    output:
    tuple path("Mapping-Stats-short.${sample_id}.csv"), emit: mapping_stats

    script:
    """
    python ${baseDir}/bin/summarize_short_reference.py --bam ${bam} --reference ${reference}
    """
}

process SUMMARIZE_STATS_LONG {
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'copy'
    label 'pysam'
    tag "$sample_id"

    input:
    tuple val(sample_id), file(bam), file(bamIdx)
    path(reference)
    val(tool)

    output:
    tuple val(sample_id), path("Mapping-Stats.${sample_id}.csv"), emit: mapping_stats
    tuple val(sample_id), path("Operations-log.${sample_id}.csv"), emit: operation_logs
    tuple val(sample_id), path("${sample_id}-error.maxgap.csv"), emit: max_gap_debug

    script:
    """
    python ${baseDir}/bin/summarize_mapped_reads.py --bam ${bam} --reference ${reference} --tool ${tool}
    echo "Rerun2"
    """
}

process COUNT_FILTERED_BARCODES {
    publishDir "${params.outdir}/${method}/${sample_id}", mode: 'copy'
    label 'samtools_8'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(method)

    output:
    tuple path("${sample_id}.counts.csv")

    script:
    """
    echo "gene,count" > ${sample_id}.counts.csv
    samtools view -h -F 2304 ${bam} | \
    awk 'BEGIN{OFS="\t"}
        /^@/ {next}
        {
            nm = ""
            for (i = 12; i <= NF; i++) {
                if (\$i ~ /^NM:i:/) {
                    split(\$i, a, ":")
                    nm = a[3]
                    break
                }
            }
            if (nm == 0) {
                print \$3
            }
        }
    ' | \
    LC_ALL=C sort | \
    uniq -c | \
    awk '{OFS=","; print \$2, \$1}' | \
    sort -k1,1nr \
    >> ${sample_id}.counts.csv
    """
}

