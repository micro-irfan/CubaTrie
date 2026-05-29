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

include { BOWTIE2_BUILD as BWT_INDEX_SHORT } from './bin/mapper'
include { BOWTIE2_BUILD as BWT_INDEX_LONG } from './bin/mapper'
include { BWA_INDEX as BWA_INDEX_LONG } from './bin/mapper'
include { BWA_INDEX as BWA_INDEX_SHORT } from './bin/mapper'
include { BWA_SINGLE as BWA_SCAFFOLD } from './bin/mapper'
include { BWA_SINGLE as BWA_SHORT_FASTP } from './bin/mapper'

include { BOWTIE2_ALIGN_SINGLE_READ as BWT2_SINGLE_READ_E2E } from './bin/mapper'
include { BOWTIE2_ALIGN_SINGLE_READ as BWT2_SINGLE_READ_SHORT } from './bin/mapper'
include { BOWTIE2_ALIGN_SINGLE_READ_S } from './bin/mapper'
include { BOWTIE2_ALIGN_SINGLE_READ as BWT2_SINGLE_READ_E2E_CA } from './bin/mapper'

include { CUBATRIE_CLEAN as CUBATRIE_CLEAN_8_R1 } from './bin/mapper'
include { CUBATRIE_BAM as CUBATRIE_SAM_8_R1 } from './bin/mapper'
include { CUBATRIE_CLEAN as CUBATRIE_CLEAN_8_CA } from './bin/mapper'
include { CUBATRIE_BAM as CUBATRIE_SAM_8_CA } from './bin/mapper'

include { CUTADAPT_FLEX } from './bin/mapper'
include { FASTP_TRIMMED_SINGLE } from './bin/mapper'

include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_1 } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_2 } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_3 } from './bin/mapper'
include { SUMMARIZE_STATS_SHORT as SUMMARIZE_SHORT_1 } from './bin/mapper'
include { SUMMARIZE_STATS_SHORT as SUMMARIZE_SHORT_2 } from './bin/mapper'
include { SUMMARIZE_STATS_SHORT as SUMMARIZE_SHORT_3 } from './bin/mapper'

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
        def sample_id = row.sample_id
        def dir = params.workDir
        tuple(
            sample_id,
            file("${dir}/data/${sample_id}.fastq.gz")
        )
    }
    .set { input_reads }

workflow {
    ref_short = "${params.workDir}/ref/rbp_sgrna_short.fasta"
    ref_long = "${params.workDir}/ref/rbp_sgrna_scaffold.fasta"

    FASTP_TRIMMED_SINGLE(input_reads)
    
    long_index_bwt = BWT_INDEX_LONG(ref_long)
    short_index_bwt = BWT_INDEX_SHORT(ref_short)
    long_index_bwa = BWA_INDEX_LONG(ref_long)
    short_index_bwa = BWA_INDEX_SHORT(ref_short)

    mapped_bam_bwa1 = BWA_SCAFFOLD(FASTP_TRIMMED_SINGLE.out.trimmed_reads, long_index_bwa, 'bwa_scaffold')
    SUMMARIZE_LONG_1(mapped_bam_bwa1, ref_long, 'bwa_scaffold')

    mapped_bam_bwa2 = BWA_SHORT_FASTP(FASTP_TRIMMED_SINGLE.out.trimmed_reads, short_index_bwa, 'bwa_short_fastp')
    SUMMARIZE_SHORT_1(mapped_bam_bwa2, ref_short, 'bwa_short_fastp')

    mapped_bam_bwt1 = BOWTIE2_ALIGN_SINGLE_READ_S(FASTP_TRIMMED_SINGLE.out.trimmed_reads, long_index_bwt, 'bowtie2_scaffold_s')
    SUMMARIZE_LONG_2(mapped_bam_bwt1, ref_long, 'bowtie2_scaffold_s')

    mapped_bam_bwt2 = BWT2_SINGLE_READ_E2E(FASTP_TRIMMED_SINGLE.out.trimmed_reads, long_index_bwt, 'bowtie2_scaffold_e2e')
    SUMMARIZE_LONG_3(mapped_bam_bwt2, ref_long, 'bowtie2_scaffold_e2e')

    mapped_bam_bwt4 = BWT2_SINGLE_READ_SHORT(FASTP_TRIMMED_SINGLE.out.trimmed_reads, short_index_bwt, 'bowtie2_short_e2e')
    SUMMARIZE_SHORT_3(mapped_bam_bwt4, ref_short, 'bowtie2_short_e2e')

    CUBATRIE_CLEAN_8_R1(FASTP_TRIMMED_SINGLE.out.trimmed_reads, ref_short, 8, 'single')
    CUBATRIE_SAM_8_R1(FASTP_TRIMMED_SINGLE.out.trimmed_reads, ref_short, 8, 'single')

    // CUTADAPT
    CUTADAPT_FLEX(input_reads, "ATTTTCAATTTAACGTCG", "GTTTTAGAGCTAGAAATA")
    mapped_bam_bwt3 = BWT2_SINGLE_READ_E2E_CA(CUTADAPT_FLEX.out.trimmed_reads, long_index_bwt, 'bowtie2_CA')
    SUMMARIZE_SHORT_2(mapped_bam_bwt3, ref_short, 'bowtie2_CA')

    CUBATRIE_CLEAN_8_CA(CUTADAPT_FLEX.out.trimmed_reads, ref_short, 8, 'cutadapt')
    CUBATRIE_SAM_8_CA(CUTADAPT_FLEX.out.trimmed_reads, ref_short, 8, 'cutadapt')
}