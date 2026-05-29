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
    .fromPath(params.samplesheet_trim25)
    .splitCsv(header: true)
    .map { row ->
        def sample_id = row.sample_id
        def dir = params.workDir
        tuple(
            sample_id,
            file("${dir}/trim25/data/n1m_mrna/${sample_id}_1.fastq.gz"),
            file("${dir}/trim25/data/n1m_mrna/${sample_id}_2.fastq.gz")
        )
    }
    .set { input_reads_trim25 }

Channel
    .fromPath(params.samplesheet_droso)
    .splitCsv(header: true)
    .map { row ->
        def sample_id = row.sample_id
        def dir = params.workDir
        tuple(
            sample_id,
            file("${dir}/drosophila_paper/data/${sample_id}.fastq.gz")
        )
    }
    .set { input_reads_droso }

include { FASTP_TRIMMED_SINGLE as FASTP_DROSO  } from './bin/mapper'
include { FASTP_TRIMMED as FASTP_TRIM25  } from './bin/mapper'

include { CUTADAPT_FLEX as CA_DROSO  } from './bin/mapper'
include { CUTADAPT_FLEX as CA_TRIM25  } from './bin/mapper'

include { CUBATRIE_CUT as CB_DROSO  } from './bin/mapper'
include { CUBATRIE_CUT as CB_TRIM25  } from './bin/mapper'

include { COUNT_READS as COUNT_CA_DROSO  } from './bin/mapper'
include { COUNT_READS as COUNT_CB_DROSO  } from './bin/mapper'
include { COUNT_READS as COUNT_CA_TRIM25  } from './bin/mapper'
include { COUNT_READS as COUNT_CB_TRIM25  } from './bin/mapper'

workflow {
    droso_flank_5 = "ATTTTCAATTTAACGTCG"
    droso_flank_3 = "GTTTTAGAGCTAGAAATA"

    FASTP_DROSO(input_reads_droso)

    CA_DROSO(FASTP_DROSO.out.trimmed_reads, droso_flank_5, droso_flank_3)
    CB_DROSO(FASTP_DROSO.out.trimmed_reads, droso_flank_5, droso_flank_3)

    COUNT_CA_DROSO(CA_DROSO.out.trimmed_reads, 'cutadapt')
    COUNT_CB_DROSO(CB_DROSO.out.trimmed_reads, 'cubatrie_cut')

    trim25_flank_5 = "GGAAAGGACGAAACACCG"
    trim25_flank_3 = "GTTTTAGAGCTAGAAATA"

    FASTP_TRIM25(input_reads_trim25)

    r1_reads = FASTP_TRIM25.out.trimmed_reads.map { sample_id, reads ->
        tuple(sample_id, reads[0])
    }

    CA_TRIM25(r1_reads, trim25_flank_5, trim25_flank_3)
    CB_TRIM25(r1_reads, trim25_flank_5, trim25_flank_3)
    
    COUNT_CA_TRIM25(CA_TRIM25.out.trimmed_reads, 'cutadapt')
    COUNT_CB_TRIM25(CB_TRIM25.out.trimmed_reads, 'cubatrie_cut')
}