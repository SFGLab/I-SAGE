process MARKDUP {
    tag "$sample_id"
    publishDir "${params.outdir}/mapping/dedup", mode: 'copy'

    input:
    tuple val(sample_id), file(bam)

    output:
    tuple val(sample_id), file("${sample_id}.dedup.bam")

    script:
    """
    samtools markdup -r ${bam} ${sample_id}.dedup.bam
    samtools index ${sample_id}.dedup.bam
    """
}
