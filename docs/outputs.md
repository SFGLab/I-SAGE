````markdown
# Outputs and Interpretation

This section describes the outputs produced by the I-SAGE pipeline, their directory structure, file formats, and how to interpret them correctly.

Understanding these outputs is critical for drawing valid biological conclusions.

---

## Output Directory Structure

All outputs are written under the directory specified by:

```yaml
outdir: "/path/to/output"
````

When bin-size sweep is enabled, outputs are organized per bin size.

Typical structure:

```
outdir/
├── viz/
│   ├── bin_250/
│   ├── bin_500/
│   └── bin_1000/
├── stats/
│   ├── bin_250/
│   ├── bin_500/
│   └── bin_1000/
├── validation/
│   ├── bin_250/
│   ├── bin_500/
│   └── bin_1000/
└── logs/
```

Each bin-size directory contains results specific to that resolution.

---

## Visualization Outputs (`viz/`)

### BedGraph Tracks

Location:

```
viz/bin_<size>/
```

Files include:

* strand-specific bedGraphs
* total-signal bedGraphs
* normalized tracks (if normalization is enabled)

**Use cases**

* Load into IGV or UCSC Genome Browser
* Inspect local break patterns
* Compare signal distribution across conditions

**Important**

* BedGraphs are binned representations
* They do not represent single-base resolution after binning

---

## Differential Statistics Outputs (`stats/`)

Location:

```
stats/bin_<size>/
```

Each contrast produces a consistent set of files.

---

### Main Results Table

**File**

```
<contrast>.tsv
```

**Columns**

* `chrom`, `start`, `end`: genomic bin coordinates
* `count_case`: break counts in case condition
* `count_ctrl`: break counts in control condition
* `log2FC`: log2 fold change
* `pvalue`: raw Fisher test p-value
* `padj`: FDR-adjusted p-value
* `is_EBV` (optional): EBV annotation

**Interpretation**

* This file includes **all tested bins**
* Most bins will not be significant after FDR correction

---

### Significant Bins

**Files**

```
<contrast>.sig.tsv
<contrast>.sig_up.tsv
<contrast>.sig_down.tsv
```

* `sig.tsv`: all significant bins (FDR ≤ threshold)
* `sig_up.tsv`: bins with `log2FC > 0`
* `sig_down.tsv`: bins with `log2FC < 0`

**Interpretation**

* Directional split is essential for biological interpretation
* Absolute counts of significant bins depend strongly on bin size

---

### Summary File

**File**

```
<contrast>.summary.txt
```

**Key Fields**

* `bins_tested`
* `bins_FDR_le_threshold`
* `sig_up`
* `sig_down`
* `replicate_method`
* `fdr_threshold`
* `ebv_*` metrics (if enabled)

**Interpretation**

* This file is the authoritative numerical summary
* Use it for reporting counts and comparisons

---

### Volcano Plots

**Files**

```
<contrast>.volcano.png
<contrast>.volcano.pdf
```

**Axes**

* X-axis: log2 fold change
* Y-axis: −log10(p-value)

**Features**

* FDR threshold indicated
* Up/down significant bins highlighted
* Subtitle reports replicate method and sample counts

**Interpretation**

* Use to assess effect size vs significance
* Large fold-change without significance indicates low support
* Dense vertical bands suggest global shifts or normalization effects

---

### MA Plots

**Files**

```
<contrast>.ma.png
<contrast>.ma.pdf
```

**Axes**

* X-axis: log10 mean break count
* Y-axis: log2 fold change

**Interpretation**

* Highlights intensity-dependent biases
* Helps identify low-count artifacts

---

## Bin-Size Sweep Summary

**File**

```
stats/bin_sweep_summary.tsv
```

**Contents**

* Bin size
* Contrast
* Number of tested bins
* Number of significant bins
* Up/down counts
* EBV metrics (if enabled)

**Interpretation**

* Used to assess robustness across bin sizes
* Large variability across bin sizes indicates resolution sensitivity
* Stable signals across sizes are more reliable

---

## Validation Outputs (`validation/`)

Location:

```
validation/bin_<size>/
```

---

### Downsampling Results

**Files**

```
<contrast>.downsample.png
<contrast>.downsample.pdf
```

**Interpretation**

* Shows how significance degrades as data is reduced
* Stable signals persist across downsampling
* Rapid collapse suggests marginal significance

---

### Spike-In Validation

**Outputs**

* Reported in validation summaries
* Evaluates sensitivity and false-negative behavior

**Interpretation**

* High recovery indicates adequate power
* Poor recovery suggests insufficient depth or overly strict thresholds

---

## EBV-Specific Outputs (Optional)

When enabled:

* `is_EBV` column in TSVs
* EBV enrichment statistics in summary files

**Interpretation**

* EBV enrichment reflects biological consistency
* Should not be interpreted as mechanistic proof
* Sensitive to reference composition

---

## Common Interpretation Pitfalls

* Comparing raw counts across bin sizes
* Treating all significant bins as equally meaningful
* Ignoring replicate consistency
* Overinterpreting EBV enrichment without context
* Assuming genome-wide significance implies local biological relevance

---

## Recommended Interpretation Workflow

1. Inspect summary files
2. Review volcano and MA plots
3. Compare results across bin sizes
4. Check validation stability
5. Visualize key regions in a genome browser

---

## Next Section

**Developer Guide** — internal structure, extension points, and best practices for maintaining the pipeline.

```

