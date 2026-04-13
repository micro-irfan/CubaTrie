#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BOWTIE2_ALIGN_S {
    label "bowtie2"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_mapped.sam")

    script:
    """
    bowtie2-build ${reference} rbp_sgrna
    bowtie2 --sensitive-local \
                -p ${task.cpus} \
                -k 1 \
                -x rbp_sgrna \
                -1 ${reads[0]} -2 ${reads[1]} > ${sample_id}_mapped.sam
    """
}

process BOWTIE2_ALIGN_VS {
    label "bowtie2"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_mapped.sam")

    script:
    """
    bowtie2-build ${reference} rbp_sgrna
    bowtie2 --very-sensitive-local \
            -p ${task.cpus} \
            -k 1 \
            -x rbp_sgrna \
            -1 ${reads[0]} -2 ${reads[1]} > ${sample_id}_mapped.sam
    """
}


process BWA {
    label 'bwa'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path(reference)
    maxRetries 1 

    output:
    tuple val(sample_id), path("${sample_id}_mapped.sam")

    script:
    """
    bwa index ${reference}
    bwa mem -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} > ${sample_id}_mapped.sam
    """
}

process SAM_TO_BAM {
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
    publishDir "${params.outdir}/cubatrie_${n_cpus}/${sample_id}", mode: 'copy'
    cpus {n_cpus}

    input:
        tuple val(sample_id), path(reads)
        path(reference)
        val(n_cpus)

    output:
        path "${sample_id}.counts.csv"

    script:
    """
    ${params.cubatrie} -r ${reference} -i ${reads} -o ${sample_id}.counts.csv -d warn -t ${n_cpus}
    """
}

process CUBATRIE_SAM {
    publishDir "${params.outdir}/cubatrie_sam_${n_cpus}/${sample_id}", mode: 'copy'
    cpus {n_cpus}
    tag "$sample_id"

    input:
        tuple val(sample_id), path(reads)
        path(reference)
        val(n_cpus)

    output:
        tuple path("${sample_id}.counts.csv"), path("${sample_id}.sam")

    script:
    """
    ${params.cubatrie} -r ${reference} -i ${reads} -o ${sample_id}.counts.csv -d warn -t ${n_cpus} --sam ${sample_id}.sam
    """
}

process FASTP_DEDUP {
    publishDir "${params.outdir}/trimmed/${sample_id}"
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

process FASTP_MERGE {
    publishDir "${params.outdir}/trimmed/${sample_id}"
    label "fastp"
    tag "$sample_id"

    input:
    tuple val(sample_id), file(reads)

    output:
    tuple val(sample_id), path("${sample_id}.merged.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.merge.json", emit: json_report
    path "${sample_id}_fastp.merge.html", emit: html_report

    script:
    """
    fastp \
        --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --merge \
        --merged_out ${sample_id}.merged.fastq.gz \
        --out1 ${sample_id}_R1.merged.fastq.gz \
        --out2 ${sample_id}_R2.merged.fastq.gz \
        --json ${sample_id}_fastp.merge.json \
        --html ${sample_id}_fastp.merge.html \
        --thread ${task.cpus} 
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
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'move'
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
    publishDir "${params.outdir}/${tool}/${sample_id}", mode: 'move'
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
    """
}