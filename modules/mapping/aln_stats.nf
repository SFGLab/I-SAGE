process ALIGN_STATS {
    tag "$sample_id"

    input:
    tuple sample_id, file(bam)

    output:
    tuple sample_id, 
         file("${sample_id}.flagstat.txt"),
         file("${sample_id}.stats.txt")

    script:
    """
    samtools flagstat ${bam} > ${sample_id}.flagstat.txt
    samtools stats ${bam} > ${sample_id}.stats.txt
    """
}
