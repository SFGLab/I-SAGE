# Pipeline Modules

This section documents the internal structure of the I-SAGE pipeline. Each stage is implemented as a modular Nextflow process, with clear inputs, outputs, and responsibilities.

Understanding these modules is useful for:

- Interpreting intermediate outputs
- Debugging failed runs
- Extending the pipeline with new functionality

## High-Level Module Flow

```
FASTQ
  ↓
Alignment & Deduplication
  ↓
Break Calling
  ↓
Visualization Tracks
  ↓
Normalization
  ↓
Differential Statistics
  ↓
Validation & Sensitivity Analysis
```

Each step is described below.

## Alignment & Deduplication

**Purpose**

Align raw iBLESS reads to the reference genome and remove PCR duplicates.

**Implementation**

- Uses standard short-read alignment tools (configured via Nextflow profile)
- Produces coordinate-sorted BAM files

**Outputs**

- Deduplicated BAM files per sample

**Notes**

- Alignment parameters are defined in `nextflow.config`
- Reference genome must match downstream chromosome naming

## Break Calling

**Module Location**

`modules/break_calling/`

**Purpose**

Identify DNA double-strand breaks from aligned reads.

**Modes**

- `per_base` (default) — breaks are identified at single-base resolution
- `binned` — breaks are directly aggregated into bins (optional)

**Key Parameters**

```yaml
break_calling:
  mode: "per_base"
  bin_size: 1000
  remove_sgrdi: true
```

**Outputs**

- Strand-specific break calls
- Per-base break coordinates

**Notes**

- In `per_base` mode, binning is deferred to the visualization stage
- Restriction enzyme site removal reduces background artifacts

## Visualization Tracks

**Module Location**

`modules/viz/`

**Purpose**

Aggregate per-base break calls into fixed-width genomic bins and generate bedGraph tracks for visualization and statistics.

**Key Parameters**

```yaml
viz:
  enabled: true
  bin_sizes: [250, 500, 1000]
```

**Behavior**

- Break counts are summed per bin
- Plus and minus strands are processed separately
- Multiple bin sizes can be evaluated in parallel

**Outputs**

- Strand-specific bedGraph files
- Combined total-signal bedGraph files

**Notes**

- Outputs are organized per bin size
- Tracks are suitable for IGV/UCSC Genome Browser

## Normalization

**Module Location**

`modules/viz/normalize_tracks.nf`

**Purpose**

Normalize binned break counts to enable comparisons across samples.

**Method**

- Counts-per-million (CPM)–style normalization

**Key Parameters**

```yaml
normalization:
  method: dsb_cpm
  scale: 1000000
```

**Outputs**

- Normalized bedGraph tracks

## Differential Statistics

**Module Location**

- `modules/stats/differential_breaks.py`
- `modules/stats/diff_breaks.nf`

**Purpose**

Identify genomic bins with statistically significant differences in break frequency between conditions.

**Statistical Framework**

- Per-bin Fisher exact test
- Benjamini–Hochberg FDR correction
- Optional replicate-aware meta-analysis

**Key Features**

- Replicate-aware testing (`meta_fisher`)
- Upregulated vs downregulated bin separation
- Optional region-restricted testing (BED)
- Optional EBV annotation and enrichment

**Outputs**

- TSV files of all bins
- Significant bins (all / up / down)
- Volcano and MA plots (PNG + PDF)
- Per-contrast summary files

## Bin-Size Sweep Summary

**Module Location**

- `modules/stats/bin_sweep_summary.py`
- `modules/stats/bin_sweep_summary.nf`

**Purpose**

Summarize differential statistics across multiple bin sizes.

**Behavior**

- Aggregates per-bin-size summaries
- Reports number of tested bins and significant bins
- Tracks up/down counts and EBV metrics

**Outputs**

- `bin_sweep_summary.tsv`

**Notes**

- Facilitates sensitivity analysis
- Enables informed bin-size selection

## Validation & Sensitivity Analysis

**Module Location**

`modules/validation/`

**Purpose**

Assess robustness and reproducibility of differential results.

**Validation Strategies**

**Downsampling**

- Re-run stats on subsets of the data
- Measure stability of significant bins

**Spike-In**

- Artificially add signal to random bins
- Evaluate recovery performance

**Key Parameters**

```yaml
validation:
  downsample_fracs: "1.0,0.5,0.25,0.1"
  downsample_reps: 1
  spikein_bins: 1000
  spikein_mult: 3.0
```

**Outputs**

- Validation reports
- Downsampling plots (PNG + PDF)

## Workflow Orchestration

**Module Location**

`workflows/iblesse_month2/main.nf`

**Purpose**

Coordinate execution of all modules, handle parameter propagation, and manage bin-size sweeps.

**Responsibilities**

- Fan-out across bin sizes
- Ensure consistent inputs to stats and validation
- Organize outputs into structured directories

## Extending the Pipeline

To add a new module:

1. Create a new process under `modules/`
2. Define explicit inputs and outputs
3. Wire it into `main.nf`
4. Document the module here

---

**Next:** See [Statistical Methods](statistical_methods.md) for detailed explanation of the statistical models, assumptions, and limitations.
