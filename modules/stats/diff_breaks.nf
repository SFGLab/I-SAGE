// modules/stats/diff_breaks.nf
// Month 3 Step 3: differential testing + volcano + tables
nextflow.enable.dsl = 2

process DIFF_BREAKS {
    tag "${contrast_name}"
    publishDir "${params.outdir}/stats", mode: 'copy'

    /*
      Accept a SINGLE tuple input per contrast:
      (samples_tsv, contrast_name, contrast_spec)
      This allows calling DIFF_BREAKS(stats_in_ch) once for many contrasts.
    */
    input:
    tuple path(samples_tsv), val(contrast_name), val(contrast_spec)   // "name:CASE1,CASE2:CTRL1,CTRL2"

    output:
    path("${contrast_name}.tsv")
    path("${contrast_name}.volcano.png")
    path("${contrast_name}.summary.txt")

    script:
    def eps = params.stats?.eps ?: 0.5

    """
    set -euo pipefail

    python ${projectDir}/../../modules/stats/differential_breaks.py \
      --samples-tsv ${samples_tsv} \
      --contrast "${contrast_spec}" \
      --out-prefix ${contrast_name} \
      --eps ${eps}
    """
}