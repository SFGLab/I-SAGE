process ALIGN_BWA {
    tag "$sample_id"

    input:
    tuple sample_id, file(reads)

    output:
    tuple sample_id, file("${sample_id}.sorted.bam")

    script:
    """
    bwa mem -t ${task.cpus} ${params.genome_index} ${reads} \
      | samtools sort -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}
