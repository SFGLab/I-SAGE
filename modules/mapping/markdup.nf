process MARKDUP {
    tag "$sample_id"

    input:
    tuple sample_id, file(bam)

    output:
    tuple sample_id, file("${sample_id}.dedup.bam")

    script:
    """
    # Remove duplicates with samtools markdup
    samtools markdup -r ${bam} ${sample_id}.dedup.bam
    samtools index ${sample_id}.dedup.bam
    """
}
