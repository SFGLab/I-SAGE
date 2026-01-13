// modules/stats/diff_breaks.nf
// Differential testing + volcano + tables + significant-bin split outputs
nextflow.enable.dsl = 2

process DIFF_BREAKS {
    tag "${contrast_name} (bin=${bin_size})"
    publishDir "${params.outdir}/stats/bin_${bin_size}", mode: 'copy'

    input:
    tuple val(bin_size), path(samples_tsv), val(contrast_name), val(contrast_spec)   // "name:CASE1,CASE2:CTRL1,CTRL2"

    output:
    tuple val(bin_size), val(contrast_name),
          path("${contrast_name}.tsv"),
          path("${contrast_name}.sig.tsv"),
          path("${contrast_name}.sig_up.tsv"),
          path("${contrast_name}.sig_down.tsv"),
          path("${contrast_name}.volcano.png"),
          path("${contrast_name}.volcano.pdf"),
          path("${contrast_name}.ma.png"),
          path("${contrast_name}.ma.pdf"),
          path("${contrast_name}.summary.txt")

    script:
    def eps = params.stats?.eps ?: 0.5
    def rep_method = params.stats?.replicate_method ?: 'pooled'
    def ebv_regex = params.stats?.ebv_regex ?: ''
    def ebv_arg = ebv_regex ? "--ebv-regex '${ebv_regex}'" : ''
    def fdr = params.stats?.fdr ?: 0.05

    def regions_bed = params.stats?.regions_bed
    def regions_arg = regions_bed ? "--regions-bed ${regions_bed}" : ""

    """
    set -euo pipefail

    python ${projectDir}/modules/stats/differential_breaks.py \\
      --samples-tsv ${samples_tsv} \\
      --contrast \"${contrast_spec}\" \\
      --out-prefix ${contrast_name} \\
      --eps ${eps} \\
      --fdr ${fdr} \\
      --replicate-method ${rep_method} \\
      ${ebv_arg} \
      ${regions_arg}
    """
}