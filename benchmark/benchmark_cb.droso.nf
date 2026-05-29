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
            file("${dir}/data/${sample_id}.fastq.gz")
        )
    }
    .set { input_reads }

include { FASTP_TRIMMED_SINGLE  } from './bin/mapper'
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
    ref_short = "${params.workDir}/ref/rbp_sgrna_short.fasta"
    ref_long = "${params.workDir}/ref/rbp_sgrna_scaffold.fasta"

    flank_5 = "ATTTTCAATTTAACGTCG"
    flank_3 = "GTTTTAGAGCTAGAAATA"

    FASTP_TRIMMED_SINGLE(input_reads)

    CUBATRIE_CLEAN_ANCHORS(FASTP_TRIMMED_SINGLE.out.trimmed_reads, ref_short, 8, 'trimmed_CB_anchors', flank_5, flank_3)
    CUBATRIE_BAM_ANCHORS(FASTP_TRIMMED_SINGLE.out.trimmed_reads, ref_short, 8, 'trimmed_CB_anchors', flank_5, flank_3)

    // CUBATRIE_CLEAN_R1(FASTP_TRIMMED_SINGLE.out.trimmed_reads, ref_short, 8, 'trimmed_R1_anchors')
    // CUBATRIE_SAM_R1(FASTP_TRIMMED_SINGLE.out.trimmed_reads, ref_short, 8, 'trimmed_R1_anchors')

    // CUTADAPT_FLEX(input_reads, flank_5, flank_3)
    // CUBATRIE_CLEAN_CA(CUTADAPT_FLEX.out.trimmed_reads, ref_short, 8, 'trimmed_CA_anchors')
    // CUBATRIE_SAM_CA(CUTADAPT_FLEX.out.trimmed_reads, ref_short, 8, 'trimmed_CA_anchors')

    // CUBATRIE_CLEAN_ANCHORS(input_reads, ref_short, 8, 'CB_anchors', flank_5, flank_3)
    // CUBATRIE_CLEAN_R1(input_reads, ref_short, 8, 'R1_no_anchors')

    // CUTADAPT_FLEX(FASTP_TRIMMED_SINGLE.out.trimmed_reads, flank_5, flank_3)
    // CUBATRIE_CLEAN_CA(CUTADAPT_FLEX.out.trimmed_reads, ref_short, 8, 'trimmed_CA_anchors')
    
    // long_index_bwt = BWT_INDEX_LONG(ref_long)
    // mapped_bam_bwt2 = BWT2_SINGLE_READ_E2E_CA(CUTADAPT_FLEX.out.trimmed_reads, long_index_bwt, 'trimmed_bowtie2_CA')
    // COUNT_FILTERED_BARCODES(mapped_bam_bwt2, 'trimmed_bowtie2_CA')

}
