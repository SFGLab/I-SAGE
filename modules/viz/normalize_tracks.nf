// modules/viz/normalize_tracks.nf
// Month 3 Step 2: DSB-CPM normalization of binned bedGraphs + normalized bigWig
nextflow.enable.dsl = 2

process NORMALIZE_TRACKS {
    tag "${sample_id}"
    publishDir "${params.outdir}/viz/normalized", mode: 'copy'

    input:
    tuple val(sample_id),
          path(plus_bg),
          path(minus_bg),
          path(total_txt)
    path genome_fasta

    output:
    tuple val(sample_id),
          path("${sample_id}.plus.norm.bedGraph"),
          path("${sample_id}.minus.norm.bedGraph"),
          path("${sample_id}.plus.norm.bw"),
          path("${sample_id}.minus.norm.bw")

    script:
    def scale_to = params.normalization?.scale ?: 1000000

    """
    set -euo pipefail

    command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found on PATH"; exit 127; }
    command -v bedGraphToBigWig >/dev/null 2>&1 || { echo "ERROR: bedGraphToBigWig not found on PATH"; exit 127; }

    # Build chrom.sizes
    if [ ! -f "${genome_fasta}.fai" ]; then
      samtools faidx "${genome_fasta}"
    fi
    cut -f1,2 "${genome_fasta}.fai" > chrom.sizes

    TOTAL=\$(cat "${total_txt}")
    if [ "\$TOTAL" -le 0 ]; then
      echo "ERROR: total_dsb is <= 0 for ${sample_id}"
      exit 1
    fi

    # DSB-CPM scaling factor
    SCALE=\$(awk -v t="\$TOTAL" -v s=${scale_to} 'BEGIN { printf "%.12f", s/t }')

    # Apply scaling to bedGraph value column (col4)
    awk -v OFS='\\t' -v S="\$SCALE" '{ \$4=\$4*S; print }' "${plus_bg}"  > "${sample_id}.plus.norm.bedGraph"
    awk -v OFS='\\t' -v S="\$SCALE" '{ \$4=\$4*S; print }' "${minus_bg}" > "${sample_id}.minus.norm.bedGraph"

    # Ensure sorted (safety)
    sort -k1,1 -k2,2n "${sample_id}.plus.norm.bedGraph"  -o "${sample_id}.plus.norm.bedGraph"
    sort -k1,1 -k2,2n "${sample_id}.minus.norm.bedGraph" -o "${sample_id}.minus.norm.bedGraph"

    # Convert normalized bedGraph to bigWig
    bedGraphToBigWig "${sample_id}.plus.norm.bedGraph"  chrom.sizes "${sample_id}.plus.norm.bw"
    bedGraphToBigWig "${sample_id}.minus.norm.bedGraph" chrom.sizes "${sample_id}.minus.norm.bw"
    """
}
