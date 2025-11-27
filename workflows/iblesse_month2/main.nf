// main.nf for Month 2: I-BLESS pipeline (FASTQ → QC → BAM → Breaks)
// ---------------------------------------------------------------

nextflow.enable.dsl = 2

include { FASTQC }      from '../../modules/qc/fastqc.nf'
include { ALIGN_BWA }   from '../../modules/mapping/bwa_mem.nf'

// BREAK_CALLING will be inline for now
// Later we can extract it to modules/break_calling/break_call.nf

// ---------------------------------------------------------------
// Input channel: single-end FASTQs (I-BLESS is SE)
// ---------------------------------------------------------------

Channel
    .fromFilePairs("${params.fastq_dir}/${params.sample_pattern}", flat: true)
    .set { reads_ch }

// ---------------------------------------------------------------
// Process: Break Calling (inline)
// ---------------------------------------------------------------

process BREAK_CALLING {
    tag "${sample_id}"

    input:
    tuple sample_id, file(bam) from bam_ch

    output:
    tuple sample_id, file("${sample_id}.breaks.bed") into breaks_bed_ch

    script:
    """
    python ${projectDir}/modules/break_calling/break_caller.py \
        --bam ${bam} \
        --out ${sample_id}.breaks.bed \
        --mode ${params.break_calling.mode} \
        --bin-size ${params.break_calling.bin_size} \
        ${params.break_calling.remove_sgrdi ? "--remove-sgrdi" : ""}
    """
}

// ---------------------------------------------------------------
// Workflow definition
// ---------------------------------------------------------------

workflow {
    take: reads_ch

    main:

        // 1) QC
        fastqc_results = FASTQC(reads_ch)

        // 2) Mapping
        bam_ch = ALIGN_BWA(reads_ch)

        // 3) Break calling
        break_results = BREAK_CALLING(bam_ch)

    emit:
        fastqc_results
        bam_ch
        break_results
}
