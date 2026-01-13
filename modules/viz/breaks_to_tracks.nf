// modules/viz/breaks_to_tracks.nf
// Month 3 Step 1: breaks.bed -> strand bedGraph + bigWig + IGV session
// (Updated for Step 2: also emits total_dsb.txt)
// Month 4 update: accept bin_size as an explicit input (supports bin-size sweeps)
nextflow.enable.dsl = 2

process BREAKS_TO_TRACKS {
    tag "${sample_id} (bin=${bin_size})"
    publishDir "${params.outdir}/viz/tracks/bin_${bin_size}", mode: 'copy'

    input:
    tuple val(bin_size), val(sample_id), path(breaks_bed)
    path genome_fasta

    output:
    tuple val(bin_size),
          val(sample_id),
          path("${sample_id}.plus.bedGraph"),
          path("${sample_id}.minus.bedGraph"),
          path("${sample_id}.plus.bw"),
          path("${sample_id}.minus.bw"),
          path("${sample_id}.total_dsb.txt")

    script:

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

    # Re-bin counts into fixed windows of size ${bin_size} and CLIP binEnd to chromosome length.
    # Input BED columns: chrom start end name score strand
    awk -v OFS='\\t' -v B=${bin_size} '
      NR==FNR { len[\$1]=\$2; next }   # first file = chrom.sizes
      {
        chrom=\$1; start=\$2; score=\$5; strand=\$6;

        # skip chromosomes not present in chrom.sizes
        if (!(chrom in len)) next;

        binStart = int(start / B) * B;
        chromLen = len[chrom];

        # skip bins that start beyond contig length
        if (binStart >= chromLen) next;

        binEnd = binStart + B;
        if (binEnd > chromLen) binEnd = chromLen;

        # must have positive width
        if (binEnd <= binStart) next;

        key = chrom OFS binStart OFS binEnd OFS strand;
        sum[key] += score;
      }
      END {
        for (k in sum) {
          split(k, a, OFS);
          chrom=a[1]; binStart=a[2]; binEnd=a[3]; strand=a[4];
          if (strand == "+") print chrom, binStart, binEnd, sum[k] > "'"${sample_id}"'.plus.bedGraph";
          else if (strand == "-") print chrom, binStart, binEnd, sum[k] > "'"${sample_id}"'.minus.bedGraph";
        }
      }
    ' chrom.sizes "${breaks_bed}"

    # Ensure bedGraphs exist even if one strand has zero bins
    [ -f "${sample_id}.plus.bedGraph" ]  || : > "${sample_id}.plus.bedGraph"
    [ -f "${sample_id}.minus.bedGraph" ] || : > "${sample_id}.minus.bedGraph"

    # Sort bedGraphs (required by bedGraphToBigWig)
    sort -k1,1 -k2,2n "${sample_id}.plus.bedGraph"  -o "${sample_id}.plus.bedGraph"
    sort -k1,1 -k2,2n "${sample_id}.minus.bedGraph" -o "${sample_id}.minus.bedGraph"

    # Convert to bigWig
    bedGraphToBigWig "${sample_id}.plus.bedGraph"  chrom.sizes "${sample_id}.plus.bw"
    bedGraphToBigWig "${sample_id}.minus.bedGraph" chrom.sizes "${sample_id}.minus.bw"
    """
}

process MAKE_IGV_SESSION {
    tag "${sample_id} (bin=${bin_size})"
    publishDir "${params.outdir}/viz/sessions/bin_${bin_size}", mode: 'copy'

    input:
    tuple val(bin_size), val(sample_id), path(plus_bw), path(minus_bw)

    output:
    tuple val(bin_size), val(sample_id), path("${sample_id}.igv_session.xml")

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
