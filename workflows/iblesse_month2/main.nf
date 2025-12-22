// main.nf
// I-BLESS Pipeline Month 2 + Month 3 Steps 1–2 (Visualization + Normalization)
// FASTQ → QC → BAM → Break Calling → bedGraph/bigWig → normalized tracks
// ---------------------------------------------------------

nextflow.enable.dsl = 2

// Include modules
include { FASTQC }      from '../../modules/qc/fastqc.nf'
include { ALIGN_BWA }   from '../../modules/mapping/bwa_mem.nf'
include { MARKDUP }     from '../../modules/mapping/markdup.nf'
include { ALIGN_STATS } from '../../modules/mapping/aln_stats.nf'

// Month 3 Step 1
include { BREAKS_TO_TRACKS; MAKE_IGV_SESSION } from '../../modules/viz/breaks_to_tracks.nf'
// Month 3 Step 2
include { NORMALIZE_TRACKS } from '../../modules/viz/normalize_tracks.nf'
include { DIFF_BREAKS } from '../../modules/stats/diff_breaks.nf'


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
    def extra_opts = ''
    if( params.break_calling.remove_sgrdi ) {
        extra_opts = "--remove-sgrdi --reference-fasta ${params.genome_fasta}"
    }

    """
    python ${projectDir}/../../modules/break_calling/break_caller.py \
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

    // Month 3 Step 1 + Step 2
    if ( params.viz?.enabled ) {

        // Step 1: create raw binned bedGraphs + bigWigs + total_dsb.txt
        tracks_out = BREAKS_TO_TRACKS(breaks_out, file(params.genome_fasta))

        /*
          tracks_out tuple:
          (sample_id,
           plus.bedGraph, minus.bedGraph,
           plus.bw, minus.bw,
           total_dsb.txt)
        */

        // IGV sessions (raw tracks)
        igv_in = tracks_out.map { sample_id, plus_bg, minus_bg, plus_bw, minus_bw, total_txt ->
            tuple(sample_id, plus_bw, minus_bw)
        }
        sessions_out = MAKE_IGV_SESSION(igv_in)

        // Step 2: normalize bedGraphs and create normalized bigWigs
        norm_in = tracks_out.map { sample_id, plus_bg, minus_bg, plus_bw, minus_bw, total_txt ->
            tuple(sample_id, plus_bg, minus_bg, total_txt)
        }
        norm_tracks_out = NORMALIZE_TRACKS(norm_in, file(params.genome_fasta))
    }
}
