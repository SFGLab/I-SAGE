// main.nf
// I-BLESS Pipeline Month 2 + Month 3 + Month 4 updates
// FASTQ → QC → BAM → Break Calling → bedGraph/bigWig → normalized tracks → differential stats → validation → bin-size sweep summary
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
// Month 3 Step 3
include { DIFF_BREAKS } from '../../modules/stats/diff_breaks.nf'
// Month 4 (P1 Task 4): bin-size sweep summary
include { BIN_SWEEP_SUMMARY } from '../../modules/stats/bin_sweep_summary.nf'
// Month 3 Step 4
include { VALIDATE_DIFF } from '../../modules/validation/validate_diff.nf'


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
    python ${projectDir}/modules/break_calling/break_caller.py \
        --bam ${bam} \
        --out ${sample_id}.breaks.bed \
        --mode ${params.break_calling.mode} \
        --bin-size ${params.break_calling.bin_size} \
        ${extra_opts}
    """
}


// ---------------------------------------------------------
// Helper: build per-bin-size samples TSV for stats module
// Input:  (bin_size, [ \"sample\\tplus_bg\\tminus_bg\", ... ]) 
// Output: (bin_size, samples.tsv)
// ---------------------------------------------------------
process MAKE_SAMPLES_TSV {
    tag "bin_${bin_size}"
    publishDir "${params.outdir}/stats/bin_${bin_size}", mode: 'copy'

    input:
    tuple val(bin_size), val(lines)

    output:
    tuple val(bin_size), file("samples_for_stats_bin_${bin_size}.tsv")

    script:
    // lines is a Groovy List<String>
    def content = (lines as List).join("\n") + "\n"
    """
    cat > samples_for_stats_bin_${bin_size}.tsv <<'EOT'
    ${content}
EOT
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

    // ---------------------------------------------------------
    // Month 3 Step 1 + Step 2 (Visualization + Normalization)
    // ---------------------------------------------------------
    if ( params.viz?.enabled ) {

        // Bin-size sweep support:
        // - if viz.bin_sizes provided: run viz→stats/validation per bin size
        // - else fall back to viz.bin_size (or 500 if unset)
        def bin_sizes = (params.viz?.bin_sizes ?: [ params.viz?.bin_size ?: 500 ]) as List
        Channel
            .fromList(bin_sizes.collect { it as int })
            .set { bin_sizes_ch }

        // Fan-out breaks across bin sizes WITHOUT re-invoking the module in a loop
        def breaks_with_bin = breaks_out
            .cross(bin_sizes_ch)
            .map { br, bs ->
                def (sample_id, breaks_bed) = br
                tuple(bs as int, sample_id, breaks_bed)
            }

        // Step 1: create raw binned bedGraphs + bigWigs + total_dsb.txt
        def tracks_out = BREAKS_TO_TRACKS(breaks_with_bin, file(params.genome_fasta))
        /*
          tracks_out tuple:
          (bin_size,
           sample_id,
           plus.bedGraph, minus.bedGraph,
           plus.bw, minus.bw,
           total_dsb.txt)
        */

        // IGV sessions (raw tracks)
        def igv_in = tracks_out.map { bin_size, sample_id, plus_bg, minus_bg, plus_bw, minus_bw, total_txt ->
            tuple(bin_size, sample_id, plus_bw, minus_bw)
        }
        MAKE_IGV_SESSION(igv_in)

        // Step 2: normalize bedGraphs and create normalized bigWigs
        def norm_in = tracks_out.map { bin_size, sample_id, plus_bg, minus_bg, plus_bw, minus_bw, total_txt ->
            tuple(bin_size, sample_id, plus_bg, minus_bg, total_txt)
        }
        NORMALIZE_TRACKS(norm_in, file(params.genome_fasta))

        // ---------------------------------------------------------
        // Month 3 Step 3 + Step 4 (Stats + Validation)
        // ---------------------------------------------------------
        if ( params.stats?.enabled ) {

            // Build per-bin-size samples TSVs from raw bedGraphs
            def samples_lines_grouped = tracks_out
                .map { bin_size, sample_id, plus_bg, minus_bg, plus_bw, minus_bw, total_txt ->
                    tuple(bin_size, "${sample_id}\t${plus_bg}\t${minus_bg}")
                }
                .groupTuple(by: 0)  // (bin_size, [line1, line2, ...])

            def samples_tsv_by_bin = MAKE_SAMPLES_TSV(samples_lines_grouped) // (bin_size, samples.tsv)

            // Build a channel of tuples: (bin_size, samples_tsv, contrast_name, contrast_spec)
            def contrasts  = params.stats.contrasts ?: []
            def conditions = params.stats.conditions ?: [:]

            def resolve_group = { obj ->
                if (obj == null) return []
                if (obj instanceof List) return obj as List
                def s = obj as String
                if (conditions.containsKey(s)) return (conditions[s] as List)
                return [s]
            }

            def stats_in_ch = samples_tsv_by_bin.flatMap { bs, samples_tsv ->
                contrasts.collect { c ->
                    def cname = c.name as String
                    def case_list = []
                    def ctrl_list = []

                    if (c.containsKey('case_condition') || c.containsKey('control_condition')) {
                        def cc = c.case_condition as String
                        def ct = c.control_condition as String
                        if (!conditions.containsKey(cc) || !conditions.containsKey(ct)) {
                            throw new IllegalArgumentException(
                                "Unknown condition key in contrast '${cname}': case_condition='${cc}', control_condition='${ct}'. Available: ${conditions.keySet()}"
                            )
                        }
                        case_list = (conditions[cc] as List)
                        ctrl_list = (conditions[ct] as List)
                    } else {
                        case_list = resolve_group(c.case)
                        ctrl_list = resolve_group(c.control)
                    }

                    if (!case_list || !ctrl_list) {
                        throw new IllegalArgumentException("Contrast '${cname}' must define non-empty case and control groups")
                    }

                    def case_ids = (case_list as List).join(',')
                    def ctrl_ids = (ctrl_list as List).join(',')
                    def spec = "${cname}:${case_ids}:${ctrl_ids}"
                    tuple(bs as int, samples_tsv, cname, spec)
                }
            }

            // Step 3: differential testing (once; runs per (bin_size, contrast) tuple)
            def diff_out = DIFF_BREAKS(stats_in_ch)

            // Step 4: validation (optional; runs per (bin_size, contrast) tuple)
            if ( params.validation?.enabled ) {
                VALIDATE_DIFF(stats_in_ch)
            }

            // After stats finish, write a single summary table across all bin sizes/contrasts
            def summary_in = diff_out.map { bin_size, contrast_name, tsv, sig, sig_up, sig_down, volcano_png, volcano_pdf, ma_png, ma_pdf, summary_txt ->
                tuple(bin_size, contrast_name, summary_txt)
            }

            def manifest_ch = summary_in.collectFile(
                name: "bin_sweep_inputs.tsv",
                storeDir: "${params.outdir}/stats",
                newLine: true
            ) { row ->
                "${row[0]}\t${row[1]}\t${row[2]}"
            }

            BIN_SWEEP_SUMMARY(manifest_ch)
        }
    }
}