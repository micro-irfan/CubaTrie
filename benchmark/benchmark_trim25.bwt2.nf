#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Default parameters
params.samplesheet = "$projectDir/samplesheet.csv"  // Add this line
params.outdir = "$projectDir/results"
params.threads = Runtime.runtime.availableProcessors()
params.help = false

// Print help message
if (params.help) {
    log.info"""
    ===================================================
    Benchmark for CubaTrie
    ===================================================
    
    Usage:
    nextflow run main.nf [options]
    
    Options:
      --samplesheet     CSV file containing sample information (default: $params.samplesheet)
      --outdir          Output directory (default: $params.outdir)
      --batchName       Batch name for processing (default: $params.batchName)
      --help            Show this message
    """
    exit 0
}

include { BOWTIE2_BUILD } from './bin/mapper'
include { BOWTIE2_ALIGN_E2E } from './bin/mapper'
include { BOWTIE2_ALIGN_S  } from './bin/mapper'
include { BOWTIE2_ALIGN_SINGLE_READ } from './bin/mapper'
include { CUTADAPT } from './bin/mapper'
include { FASTP_TRIMMED } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_2 } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_3 } from './bin/mapper'
include { SUMMARIZE_STATS_SHORT as SUMMARIZE_SHORT_2 } from './bin/mapper'

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
        def sample_id = row.sample_id
        def dir = params.workDir
        tuple(
            sample_id,
            file("${dir}/data/n1m_mrna/${sample_id}_1.fastq.gz"),
            file("${dir}/data/n1m_mrna/${sample_id}_2.fastq.gz")
        )
    }
    .set { input_reads }

workflow {
    ref_short = "${params.workDir}/ref/brunello_library.fa"
    ref_long = "${params.workDir}/ref/brunello_library_backbone.fa"

    FASTP_TRIMMED(input_reads)
    
    long_index_bwt = BOWTIE2_BUILD(ref_long)

    mapped_bam_bwt1 = BOWTIE2_ALIGN_S(FASTP_TRIMMED.out.trimmed_reads, long_index_bwt, 'bowtie2_scaffold_s')
    SUMMARIZE_LONG_2(mapped_bam_bwt1, ref_long, 'bowtie2_scaffold_s')

    mapped_bam_bwt2 = BOWTIE2_ALIGN_E2E(FASTP_TRIMMED.out.trimmed_reads, long_index_bwt, 'bowtie2_scaffold_e2e')
    SUMMARIZE_LONG_3(mapped_bam_bwt2, ref_long, 'bowtie2_scaffold_e2e')

    // CUTADAPT
    CUTADAPT(input_reads)
    mapped_bam_bwt3 = BOWTIE2_ALIGN_SINGLE_READ(CUTADAPT.out.trimmed_reads, long_index_bwt, 'bowtie2_trim25')
    SUMMARIZE_SHORT_2(mapped_bam_bwt3, ref_short, 'bowtie2_trim25')
}