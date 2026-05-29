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
include { BWA_INDEX as BWA_INDEX_LONG } from './bin/mapper'
include { BWA_INDEX as BWA_INDEX_SHORT } from './bin/mapper'
include { BWA as BWA_SCAFFOLD } from './bin/mapper'
include { BWA as BWA_SHORT_FASTP } from './bin/mapper'
include { BOWTIE2_ALIGN_E2E } from './bin/mapper'
include { BOWTIE2_ALIGN_S  } from './bin/mapper'
include { BOWTIE2_ALIGN_SINGLE_READ } from './bin/mapper'
include { CUBATRIE_CLEAN as CUBATRIE_CLEAN_8_R1 } from './bin/mapper'
include { CUBATRIE_SAM as CUBATRIE_SAM_8_R1 } from './bin/mapper'
include { CUBATRIE_CLEAN as CUBATRIE_CLEAN_4_R1 } from './bin/mapper'
include { CUBATRIE_SAM as CUBATRIE_SAM_4_R1 } from './bin/mapper'
include { CUBATRIE_CLEAN as CUBATRIE_CLEAN_8_CA } from './bin/mapper'
include { CUBATRIE_SAM as CUBATRIE_SAM_8_CA } from './bin/mapper'
include { CUTADAPT } from './bin/mapper'
include { FASTP_DEDUP } from './bin/mapper'
include { FASTP_TRIMMED } from './bin/mapper'
include { FASTP_MERGE as FASTP_MERGE_DEDUP } from './bin/mapper'
include { FASTP_MERGE as FASTP_MERGE_TRIMMED } from './bin/mapper'
include { CUBATRIEPY } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_1 } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_2 } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_3 } from './bin/mapper'
include { SUMMARIZE_STATS_SHORT as SUMMARIZE_SHORT_1 } from './bin/mapper'
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
    long_index_bwa = BWA_INDEX_LONG(ref_long)
    short_index_bwa = BWA_INDEX_SHORT(ref_short)

    mapped_bam_bwa1 = BWA_SCAFFOLD(FASTP_TRIMMED.out.trimmed_reads, long_index_bwa, 'bwa_scaffold')
    SUMMARIZE_LONG_1(mapped_bam_bwa1, ref_long, 'bwa_scaffold')

    mapped_bam_bwa2 = BWA_SHORT_FASTP(FASTP_TRIMMED.out.trimmed_reads, short_index_bwa, 'bwa_short_fastp')
    SUMMARIZE_SHORT_1(mapped_bam_bwa2, ref_short, 'bwa_short_fastp')

    mapped_bam_bwt1 = BOWTIE2_ALIGN_S(FASTP_TRIMMED.out.trimmed_reads, long_index_bwt)
    SUMMARIZE_LONG_2(mapped_bam_bwt1, ref_long, 'bowtie2_scaffold')

    mapped_bam_bwt2 = BOWTIE2_ALIGN_E2E(FASTP_TRIMMED.out.trimmed_reads, long_index_bwt)
    SUMMARIZE_LONG_3(mapped_bam_bwt2, ref_long, 'bowtie2_scaffold_e2e')

    r1_reads = FASTP_TRIMMED.out.trimmed_reads.map { sample_id, reads ->
        tuple(sample_id, reads[0])
    }

    CUBATRIE_CLEAN_4_R1(r1_reads, ref_short, 4, 'trimmed_R1')
    CUBATRIE_SAM_4_R1(r1_reads, ref_short, 4, 'trimmed_R1')

    CUBATRIE_CLEAN_8_R1(r1_reads, ref_short, 8, 'trimmed_R1')
    CUBATRIE_SAM_8_R1(r1_reads, ref_short, 8, 'trimmed_R1')

    // Python Version
    CUBATRIEPY(r1_reads, ref_short)

    // CUTADAPT
    CUTADAPT(input_reads)
    mapped_bam_bwt3 = BOWTIE2_ALIGN_SINGLE_READ(CUTADAPT.out.trimmed_reads, long_index_bwt)
    SUMMARIZE_SHORT_2(mapped_bam_bwt3, ref_short, 'bowtie2_trim25')

    CUBATRIE_CLEAN_8_CA(CUTADAPT.out.trimmed_reads, ref_short, 8, 'trimmed_CA')
    CUBATRIE_SAM_8_CA(CUTADAPT.out.trimmed_reads, ref_short, 8, 'trimmed_CA')
}