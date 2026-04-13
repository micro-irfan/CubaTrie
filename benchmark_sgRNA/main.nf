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
      --workDir         Path to Working Directory with Data and Reference
      --cubatrie        Path to cubaTrie exe
      --samplesheet     CSV file containing sample information (default: $params.samplesheet)
      --outdir          Output directory (default: $params.outdir)
      --batchName       Batch name for processing (default: $params.batchName)
      --help            Show this message
    """
    exit 0
}

include { BWA as BWA_SCAFFOLD } from './bin/mapper'
include { BWA as BWA_SHORT } from './bin/mapper'
include { BOWTIE2_ALIGN_VS } from './bin/mapper'
include { BOWTIE2_ALIGN_S  } from './bin/mapper'
include { CUBATRIE_CLEAN  } from './bin/mapper'
include { CUBATRIE_SAM  } from './bin/mapper'
include { SAM_TO_BAM as SAM_TO_BAM_1  } from './bin/mapper'
include { SAM_TO_BAM as SAM_TO_BAM_2  } from './bin/mapper'
include { SAM_TO_BAM as SAM_TO_BAM_3  } from './bin/mapper'
include { SAM_TO_BAM as SAM_TO_BAM_4  } from './bin/mapper'
include { FASTP_DEDUP  } from './bin/mapper'
include { FASTP_MERGE  } from './bin/mapper'
include { CUBATRIEPY  } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_1 } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_2 } from './bin/mapper'
include { SUMMARIZE_STATS_LONG as SUMMARIZE_LONG_3 } from './bin/mapper'
include { SUMMARIZE_STATS_SHORT as SUMMARIZE_SHORT } from './bin/mapper'

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
        def sample_id = row.sample_id
        def dir = params.workDir
        tuple(
            sample_id,
            file("${dir}/data/${sample_id}/${sample_id}_External-MUX12864_DKDL250003494-1A_232KHCLT3_L3_1.fq.gz"),
            file("${dir}/data/${sample_id}/${sample_id}_External-MUX12864_DKDL250003494-1A_232KHCLT3_L3_2.fq.gz")
        )
    }
    .set { input_reads }


workflow {
    ref_short = "${params.workDir}/ref/rbp_sgrna_short.fasta"
    ref_long = "${params.workDir}/ref/rbp_sgrna_scaffold.fasta"

    FASTP_DEDUP(input_reads)

    bwa1_sam_ch = BWA_SCAFFOLD(FASTP_DEDUP.out.trimmed_reads, ref_long)
    mapped_bam_bwa1 = SAM_TO_BAM_1(bwa1_sam_ch, 'bwa_scaffold')
    SUMMARIZE_LONG_1(mapped_bam_bwa1, ref_long, 'bwa_scaffold')

    bwa2_sam_ch = BWA_SHORT(FASTP_DEDUP.out.trimmed_reads, ref_short)
    mapped_bam_bwa2 = SAM_TO_BAM_2(bwa2_sam_ch, 'bwa_short')
    SUMMARIZE_SHORT(mapped_bam_bwa2, ref_short, 'bwa_short')

    bwt1_sam_ch = BOWTIE2_ALIGN_S(FASTP_DEDUP.out.trimmed_reads, ref_long)
    mapped_bam_bwt1 = SAM_TO_BAM_3(bwt1_sam_ch, 'bowtie2_scaffold')
    SUMMARIZE_LONG_2(mapped_bam_bwt1, ref_long, 'bowtie2_scaffold')

    bwt2_sam_ch = BOWTIE2_ALIGN_VS(FASTP_DEDUP.out.trimmed_reads, ref_long)
    mapped_bam_bwt2 = SAM_TO_BAM_4(bwt2_sam_ch, 'bowtie2_scaffold_vs')
    SUMMARIZE_LONG_3(mapped_bam_bwt2, ref_long, 'bowtie2_scaffold_vs')

    // Merge Reads
    FASTP_MERGE(FASTP_DEDUP.out.trimmed_reads)

    // Python Version
    CUBATRIEPY(FASTP_MERGE.out.trimmed_reads, ref_short)

    // C Version
    CUBATRIE_CLEAN(FASTP_MERGE.out.trimmed_reads, ref_short, 4)
    CUBATRIE_SAM(FASTP_MERGE.out.trimmed_reads, ref_short, 4)
}