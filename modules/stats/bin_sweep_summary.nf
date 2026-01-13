nextflow.enable.dsl = 2

process BIN_SWEEP_SUMMARY {
    tag "bin_sweep_summary"
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    path(manifest_tsv)

    output:
    path("bin_sweep_summary.tsv")

    script:
    """
    set -euo pipefail

    cat "${manifest_tsv}" | \
      python ${projectDir}/modules/stats/bin_sweep_summary.py --out bin_sweep_summary.tsv
    """
}
