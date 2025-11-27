process FASTQC {
    tag "${sample_id}"

    input:
    set sample_id, file(reads) from reads_ch

    output:
    set sample_id, file("*.html"), file("*.zip") into fastqc_results

    script:
    """
    fastqc -o . ${reads}
    """
}
