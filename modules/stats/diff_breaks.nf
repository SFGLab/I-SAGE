// modules/stats/diff_breaks.nf
// Month 3 Step 3: differential testing + volcano + tables
nextflow.enable.dsl = 2

process DIFF_BREAKS {
    tag "${contrast_name}"
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    path samples_tsv
    val contrast_name
    val contrast_spec   // "name:CASE1,CASE2:CTRL1,CTRL2"

    output:
    path("${contrast_name}.tsv")
    path("${contrast_name}.volcano.png")
    path("${contrast_name}.summary.txt")

    script:
    """
    set -euo pipefail

    python ${projectDir}/../../modules/stats/differential_breaks.py \
      --samples-tsv ${samples_tsv} \
      --contrast "${contrast_spec}" \
      --out-prefix ${contrast_name}
    """
}
