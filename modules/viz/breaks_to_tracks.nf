// modules/viz/breaks_to_tracks.nf
// Month 3 Step 1: breaks.bed -> strand bedGraph + bigWig + IGV session
// (Updated for Step 2: also emits total_dsb.txt)
nextflow.enable.dsl = 2

process BREAKS_TO_TRACKS {
    tag "${sample_id}"
    publishDir "${params.outdir}/viz/tracks", mode: 'copy'

    input:
    tuple val(sample_id), path(breaks_bed)
    path genome_fasta

    output:
    tuple val(sample_id),
          path("${sample_id}.plus.bedGraph"),
          path("${sample_id}.minus.bedGraph"),
          path("${sample_id}.plus.bw"),
          path("${sample_id}.minus.bw"),
          path("${sample_id}.total_dsb.txt")

    script:
    def bin_size = params.viz?.bin_size ?: 500

    """
    set -euo pipefail

    command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found on PATH"; exit 127; }
    command -v bedGraphToBigWig >/dev/null 2>&1 || { echo "ERROR: bedGraphToBigWig not found on PATH"; exit 127; }

    # Build chrom.sizes from FASTA index
    if [ ! -f "${genome_fasta}.fai" ]; then
      samtools faidx "${genome_fasta}"
    fi
    cut -f1,2 "${genome_fasta}.fai" > chrom.sizes

    # Total DSB count per sample (sum of column 5)
    awk 'BEGIN{sum=0} !/^#/ {sum += \$5} END{print sum}' "${breaks_bed}" > "${sample_id}.total_dsb.txt"

    # Re-bin counts into fixed windows of size ${bin_size}
    # Input BED columns: chrom start end name score strand
    awk -v OFS='\\t' -v B=${bin_size} '
      {
        chrom=$1; start=$2; score=$5; strand=$6;
        binStart = int(start / B) * B;
        key = chrom OFS binStart OFS strand;
        sum[key] += score;
      }
      END {
        for (k in sum) {
          split(k, a, OFS);
          chrom=a[1]; binStart=a[2]; strand=a[3];
          binEnd = binStart + B;
          if (strand == "+") print chrom, binStart, binEnd, sum[k] > "'"${sample_id}"'.plus.bedGraph";
          else if (strand == "-") print chrom, binStart, binEnd, sum[k] > "'"${sample_id}"'.minus.bedGraph";
        }
      }
    ' "${breaks_bed}"

    # Sort bedGraphs (required by bedGraphToBigWig)
    sort -k1,1 -k2,2n "${sample_id}.plus.bedGraph"  -o "${sample_id}.plus.bedGraph"
    sort -k1,1 -k2,2n "${sample_id}.minus.bedGraph" -o "${sample_id}.minus.bedGraph"

    # Convert to bigWig
    bedGraphToBigWig "${sample_id}.plus.bedGraph"  chrom.sizes "${sample_id}.plus.bw"
    bedGraphToBigWig "${sample_id}.minus.bedGraph" chrom.sizes "${sample_id}.minus.bw"
    """
}

process MAKE_IGV_SESSION {
    tag "${sample_id}"
    publishDir "${params.outdir}/viz/sessions", mode: 'copy'

    input:
    tuple val(sample_id), path(plus_bw), path(minus_bw)

    output:
    path("${sample_id}.igv_session.xml")

    script:
    """
    cat > ${sample_id}.igv_session.xml << 'EOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg38" hasGeneTrack="true" hasSequenceTrack="true" version="8">
  <Resources>
    <Resource path="${plus_bw}"/>
    <Resource path="${minus_bw}"/>
  </Resources>
</Session>
EOF
    """
}
