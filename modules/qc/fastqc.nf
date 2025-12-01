process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), file(reads)

    output:
    tuple val(sample_id), file("*.html"), file("*.zip")

    script:
    """
    fastqc -o . ${reads}
    """
}
