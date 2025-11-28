process FASTQC {
    tag "$sample_id"

    input:
    tuple sample_id, file(reads)

    output:
    tuple sample_id, file("*.html"), file("*.zip")

    script:
    """
    fastqc -o . ${reads}
    """
}
