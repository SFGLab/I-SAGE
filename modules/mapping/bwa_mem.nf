process ALIGN_BWA {
    tag "${sample_id}"

    input:
    set sample_id, file(reads) from reads_ch

    output:
    tuple sample_id, file("${sample_id}.sorted.bam") into bam_ch

    script:
    """
    bwa mem -t ${task.cpus} ${params.genome_index} ${reads} | \
        samtools sort -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}
