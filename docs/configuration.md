````markdown
# Configuration Guide — `iblesse.yaml`

All behavior of the I-SAGE pipeline is controlled via a single YAML configuration file (typically `configs/iblesse.yaml`). This file defines input locations, pipeline behavior, statistical settings, and validation strategies.

This page documents **every configuration section**, what it controls, and how parameters interact.

---

## Global Parameters

### `fastq_dir`
Directory containing input FASTQ files.

```yaml
fastq_dir: "/path/to/fastq"
````

All FASTQ files matching `sample_pattern` are discovered recursively under this directory.

---

### `sample_pattern`

Glob pattern used to identify FASTQ files.

```yaml
sample_pattern: "*_R1.fastq.gz"
```

This pattern must uniquely identify one file per sample. Incorrect patterns may lead to missing or duplicated samples.

---

### `outdir`

Base directory for all pipeline outputs.

```yaml
outdir: "/path/to/output"
```

Subdirectories (`viz/`, `stats/`, `validation/`) are created automatically.

---

### `logdir`

Directory for logs, reports, and execution traces.

```yaml
logdir: "/path/to/logs"
```

---

### `genome_fasta` and `genome_index`

Reference genome FASTA and corresponding index.

```yaml
genome_fasta: "/path/to/genome.fna"
genome_index: "/path/to/genome.fna"
```

The reference may include additional contigs (e.g. EBV) if present in the alignment.

---

## Break Calling

### `break_calling`

Controls how DNA breaks are identified from aligned reads.

```yaml
break_calling:
  mode: "per_base"
  bin_size: 1000
  remove_sgrdi: true
```

#### Parameters

* `mode`

  * `per_base` (recommended): break calls at single-base resolution
  * `binned`: aggregate breaks directly into bins

* `bin_size`
  Used only when `mode: binned`. Ignored in `per_base` mode.

* `remove_sgrdi`
  Whether to remove SgrDI restriction enzyme sites from break calls.

---

## Visualization and Binning

### `viz`

Controls how break calls are aggregated into bins for visualization and statistics.

```yaml
viz:
  enabled: true
  bin_sizes: [250, 500, 1000]
```

#### Parameters

* `enabled`
  Enables generation of binned bedGraph tracks.

* `bin_size`
  Single bin size (legacy mode).

* `bin_sizes`
  List of bin sizes to evaluate simultaneously.
  When set, the pipeline performs a **bin-size sweep** and runs downstream steps independently for each bin size.

Each bin size produces separate outputs under:

```
viz/bin_<size>/
stats/bin_<size>/
validation/bin_<size>/
```

---

## Normalization

### `normalization`

Controls normalization of binned break counts.

```yaml
normalization:
  method: dsb_cpm
  scale: 1000000
```

#### Parameters

* `method`
  Currently supported:

  * `dsb_cpm`: counts-per-million normalization

* `scale`
  Scaling factor used during normalization.

---

## Differential Statistics

### `stats`

Controls differential break analysis.

```yaml
stats:
  enabled: true
  fdr: 0.05
  replicate_method: meta_fisher
```

---

### Replicate Handling

#### `replicate_method`

Controls how biological replicates are handled.

* `pooled`
  Counts are summed across replicates and tested once.

* `meta_fisher`
  Per-replicate tests are performed and p-values combined using Fisher’s method.
  Effect size is reported as the median log2 fold-change.

```yaml
replicate_method: meta_fisher
```

If only one replicate is present, the pipeline automatically falls back to pooled behavior.

---

### Conditions

Conditions define biological groups and their replicates.

```yaml
conditions:
  APH:  ["sample1_APH", "sample2_APH"]
  DMSO: ["sample1_DMSO", "sample2_DMSO"]
```

Each sample ID must correspond to a discovered FASTQ.

---

### Contrasts

Contrasts define comparisons between conditions.

```yaml
contrasts:
  - name: "APH_vs_DMSO"
    case_condition: "APH"
    control_condition: "DMSO"
```

Contrasts are evaluated independently for each bin size (if bin-size sweep is enabled).

---

### FDR Threshold

```yaml
fdr: 0.05
```

Controls significance threshold for:

* significant bins
* up/down split
* reported summary statistics

---

## Region-Restricted Testing (Optional)

Statistical testing can be restricted to specific genomic regions.

```yaml
stats:
  regions_bed: "/path/to/regions.bed"
```

Only bins overlapping these regions are tested.
Totals are computed **within the restricted region set**, not genome-wide.

---

## EBV Annotation and Enrichment (Optional)

If EBV contigs are present in the reference, EBV-specific reporting can be enabled.

```yaml
stats:
  ebv_regex: "(?i)^chrEBV$"
```

This:

* Annotates bins as EBV vs non-EBV
* Reports EBV enrichment among significant bins
* Adds EBV metrics to summary files

If `ebv_regex` is omitted or empty, EBV analysis is disabled.

---

## Validation

### `validation`

Controls robustness and sensitivity analyses.

```yaml
validation:
  enabled: true
  fdr: 0.05
  downsample_fracs: "1.0,0.5,0.25,0.1"
  downsample_reps: 1
  spikein_bins: 1000
  spikein_mult: 3.0
  seed: 123
```

#### Parameters

* `downsample_fracs`
  Fractions of data retained during downsampling.

* `downsample_reps`
  Number of replicates per downsampling fraction.

* `spikein_bins`
  Number of bins with artificial signal added.

* `spikein_mult`
  Fold-change applied to spike-in bins.

* `seed`
  Random seed for reproducibility.

Validation is performed **per bin size** if bin-size sweep is enabled.

---

## Configuration Interactions (Important)

* Bin-size sweeps multiply runtime by number of bin sizes
* Replicate-aware testing requires correctly defined conditions
* EBV analysis requires EBV contigs in the reference
* Region-restricted testing reduces multiple testing burden

---

## Recommended Workflow

1. Start with a single bin size (e.g. 500 bp)
2. Validate contrasts and outputs
3. Enable bin-size sweep
4. Enable validation and EBV analysis
5. Interpret results jointly

---

## Next Section

**Pipeline Modules** — detailed explanation of each pipeline stage and its implementation.

```

