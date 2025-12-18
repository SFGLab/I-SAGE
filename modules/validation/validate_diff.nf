// modules/validation/validate_diff.nf
// Month 3 Step 4: downsampling + spike-in validation per contrast
nextflow.enable.dsl = 2

process VALIDATE_DIFF {
    tag "${contrast_name}"
    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    path samples_tsv
    val contrast_name
    val contrast_spec

    output:
    path("${contrast_name}.downsample.tsv")
    path("${contrast_name}.downsample.png")
    path("${contrast_name}.spikein.tsv")
    path("${contrast_name}.validation.summary.txt")

    script:
    def fdr = params.validation?.fdr ?: 0.05
    def fracs = params.validation?.downsample_fracs ?: "1.0,0.5,0.25,0.1"
    def reps = params.validation?.downsample_reps ?: 1
    def spike_bins = params.validation?.spikein_bins ?: 1000
    def spike_mult = params.validation?.spikein_mult ?: 3.0
    def seed = params.validation?.seed ?: 123

    """
    set -euo pipefail

    python ${projectDir}/../../modules/validation/validate_diff.py \
      --samples-tsv ${samples_tsv} \
      --contrast "${contrast_spec}" \
      --out-prefix ${contrast_name} \
      --fdr ${fdr} \
      --downsample-fracs "${fracs}" \
      --downsample-reps ${reps} \
      --spikein-bins ${spike_bins} \
      --spikein-mult ${spike_mult} \
      --seed ${seed}
    """
}