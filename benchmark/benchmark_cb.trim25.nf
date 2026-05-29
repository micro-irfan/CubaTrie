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

include { FASTP_TRIMMED  } from './bin/mapper'
include { CUTADAPT_FLEX  } from './bin/mapper'
include { CUBATRIE_CLEAN as CUBATRIE_CLEAN_R1 } from './bin/mapper'
include { CUBATRIE_CLEAN as CUBATRIE_CLEAN_CA } from './bin/mapper'
include { CUBATRIE_BAM as CUBATRIE_SAM_R1 } from './bin/mapper'
include { CUBATRIE_BAM as CUBATRIE_SAM_CA } from './bin/mapper'
include { CUBATRIE_CLEAN_ANCHORS  } from './bin/mapper'
include { CUBATRIE_BAM_ANCHORS  } from './bin/mapper'
include { BOWTIE2_ALIGN_SINGLE_READ as BWT2_SINGLE_READ_E2E_CA } from './bin/mapper'
include { COUNT_FILTERED_BARCODES } from './bin/mapper'
include { BOWTIE2_BUILD as BWT_INDEX_LONG } from './bin/mapper'

workflow {
    ref_short = "${params.workDir}/ref/brunello_library.fa"
    ref_long = "${params.workDir}/ref/brunello_library_backbone.fa"

    FASTP_TRIMMED(input_reads)

    flank_5 = "GGAAAGGACGAAACACCG"
    flank_3 = "GTTTTAGAGCTAGAAATA"

    r1_reads = FASTP_TRIMMED.out.trimmed_reads.map { sample_id, reads ->
        tuple(sample_id, reads[0])
    }

    CUBATRIE_CLEAN_ANCHORS(r1_reads, ref_short, 8, 'trimmed_CB', flank_5, flank_3)
    CUBATRIE_BAM_ANCHORS(r1_reads, ref_short, 8, 'trimmed_CB', flank_5, flank_3)

    CUBATRIE_CLEAN_R1(r1_reads, ref_short, 8, 'trimmed_R1')
    CUBATRIE_SAM_R1(r1_reads, ref_short, 8, 'trimmed_R1')
    
    long_index_bwt = BWT_INDEX_LONG(ref_long)

    CUTADAPT_FLEX(r1_reads, flank_5, flank_3)
    mapped_bam_bwt2 = BWT2_SINGLE_READ_E2E_CA(CUTADAPT_FLEX.out.trimmed_reads, long_index_bwt, 'trimmed_bowtie2_CA')
    COUNT_FILTERED_BARCODES(mapped_bam_bwt2, 'trimmed_bowtie2_CA')
    
    CUBATRIE_CLEAN_CA(CUTADAPT_FLEX.out.trimmed_reads, ref_short, 8, 'trimmed_CA')
    CUBATRIE_SAM_CA(CUTADAPT_FLEX.out.trimmed_reads, ref_short, 8, 'trimmed_CA')
}
