// main.nf
// I-BLESS Pipeline Month 2: FASTQ → QC → BAM → Break Calling
// ---------------------------------------------------------

nextflow.enable.dsl = 2

// Include modules
include { FASTQC }      from '../../modules/qc/fastqc.nf'
include { ALIGN_BWA }   from '../../modules/mapping/bwa_mem.nf'
include { MARKDUP }     from '../../modules/mapping/markdup.nf'
include { ALIGN_STATS } from '../../modules/mapping/aln_stats.nf'


// ---------------------------------------------------------
// Input channel — SINGLE-END FASTQ FILES
// ---------------------------------------------------------
Channel
    .fromPath("${params.fastq_dir}/${params.sample_pattern}")
    .map { file ->
        def sample = file.getName().replace("_R1.fastq.gz", "")
        tuple(sample, file)
    }
    .set { reads_ch }


// ---------------------------------------------------------
// Break Calling Process
// ---------------------------------------------------------
process BREAK_CALLING {
    tag "$sample_id"
    publishDir "${params.outdir}/breaks", mode: 'copy'

    input:
    tuple val(sample_id), file(bam)

    output:
    tuple val(sample_id), file("${sample_id}.breaks.bed")

    script:
    // Build extra options for SgrDI filtering
    def extra_opts = ''
    if( params.break_calling.remove_sgrdi ) {
        extra_opts = "--remove-sgrdi --reference-fasta ${params.genome_fasta}"
    }

    """
    python ${baseDir}/modules/break_calling/break_caller.py \
        --bam ${bam} \
        --out ${sample_id}.breaks.bed \
        --mode ${params.break_calling.mode} \
        --bin-size ${params.break_calling.bin_size} \
        ${extra_opts}
    """
}

// ---------------------------------------------------------
// Workflow
// ---------------------------------------------------------
workflow {

    // 1. QC on raw reads
    fastqc_out   = FASTQC(reads_ch)

    // 2. Mapping (sorted BAM)
    mapped_bam_ch = ALIGN_BWA(reads_ch)

    // 3. Deduplication
    dedup_bam_ch  = MARKDUP(mapped_bam_ch)

    // 4. Alignment stats on deduplicated BAM
    stats_out     = ALIGN_STATS(dedup_bam_ch)

    // 5. Break calling from deduplicated BAM
    breaks_out    = BREAK_CALLING(dedup_bam_ch)
}