# Statistical Methods

This section documents the statistical framework used by I-SAGE to identify **differential DNA double-strand break (DSB) regions** from iBLESS data.

The emphasis is on **explicit assumptions**, **transparent tests**, and **robustness assessment**, rather than black-box modeling.

## Overview

I-SAGE performs **bin-level differential analysis** of DSB counts between experimental conditions.

The statistical workflow consists of:

1. Aggregation of break counts into genomic bins
2. Hypothesis testing per bin
3. Multiple testing correction
4. Optional replicate-aware aggregation
5. Robustness and sensitivity validation

## Data Representation

For each genomic bin, the pipeline computes:

- `count_case` — total break counts in the case condition
- `count_ctrl` — total break counts in the control condition

Counts are derived from strand-aware bedGraph tracks after normalization.

## Bin-Level Hypothesis Test

### Fisher Exact Test

For each bin, I-SAGE performs a **two-sided Fisher exact test**, comparing break enrichment in the bin versus the remaining tested genomic space.

The contingency table is:

```
          In bin     Outside bin
Case           k_c       T_c − k_c
Control        k_t       T_t − k_t
```

Where:

- `k_c`, `k_t` — bin-level counts
- `T_c`, `T_t` — total counts across all tested bins

This formulation tests whether the **relative enrichment of breaks** differs between conditions.

### Why Fisher Exact Test?

- Counts are sparse and non-normally distributed
- Bin sizes may vary
- No assumption of equal variance
- Exact test appropriate for low counts

This choice prioritizes **interpretability and correctness** over parametric efficiency.

## Effect Size: Log2 Fold Change

For each bin, effect size is reported as:

$$\text{log2FC} = \log_2 \left( \frac{(k_c + \varepsilon) / (T_c + \varepsilon)}{(k_t + \varepsilon) / (T_t + \varepsilon)} \right)$$

Where:

- $\varepsilon$ — small pseudocount (default: 0.5)

This represents the **relative enrichment of breaks** in the bin, normalized by global signal.

## Multiple Testing Correction

P-values are corrected using the **Benjamini–Hochberg procedure** to control the **false discovery rate (FDR)**.

Bins with:

```
FDR ≤ threshold
```

are considered statistically significant.

The default threshold is:

```yaml
stats:
  fdr: 0.05
```

## Replicate-Aware Analysis

### Motivation

Pooling replicates can allow a single replicate to dominate significance. I-SAGE therefore supports an optional **replicate-aware meta-analysis**.

### `pooled` Mode

- Counts are summed across replicates
- A single Fisher test is performed per bin

This mode is simple but ignores replicate variability.

### `meta_fisher` Mode

In replicate-aware mode:

1. For each bin, a Fisher test is performed **per replicate** against the pooled control
2. Per-replicate p-values are combined using **Fisher's method**
3. Effect size is reported as the **median log2 fold change** across replicates

This approach enforces **replicate agreement** without fitting a full variance model.

### Limitations

- Does not model dispersion explicitly
- Assumes independence between replicate tests
- Not equivalent to a negative binomial GLM

This design is intentional and favors **transparency over complexity**.

## Directional Classification

Significant bins are classified as:

- **Upregulated** — $\text{log2FC} > 0$
- **Downregulated** — $\text{log2FC} < 0$

Separate output files are generated for each class.

## Region-Restricted Testing

If a BED file is provided:

```yaml
stats:
  regions_bed: "/path/to/regions.bed"
```

- Only bins overlapping these regions are tested
- Total counts are computed **within the restricted region set**

This reduces multiple testing burden and enables hypothesis-driven analyses.

## EBV Annotation and Enrichment

If EBV contigs are present and enabled via:

```yaml
stats:
  ebv_regex: "(?i)^chrEBV$"
```

Bins are annotated as EBV vs non-EBV.

### Enrichment Test

A 2×2 contingency table is constructed:

```
                Significant   Not significant
EBV bins            a               c
Non-EBV bins        b               d
```

Enrichment is assessed using a two-sided Fisher exact test.

Reported metrics include:

- Percentage of EBV bins among tested bins
- Percentage of EBV bins among significant bins
- Enrichment ratio
- Enrichment p-value

## Validation and Sensitivity Analysis

Statistical significance alone is insufficient. I-SAGE therefore includes built-in validation.

### Downsampling

- Break counts are randomly downsampled to fixed fractions
- Differential testing is re-run
- Stability of significant bins is assessed

### Spike-In Analysis

- Artificial signal is added to random bins
- Recovery of spiked bins is measured
- Evaluates sensitivity and false-negative behavior

## Bin-Size Sensitivity

When multiple bin sizes are evaluated:

- All statistics are computed independently per bin size
- Results are summarized in a bin-size sweep table

This guards against bin-size–specific artifacts.

## Interpretation Guidelines

- Statistical significance does not imply causality
- Effect size and reproducibility should be considered jointly
- EBV enrichment is a biological consistency check, not proof of mechanism

## Summary

I-SAGE uses a **simple, explicit statistical framework**:

- Exact tests
- Clear assumptions
- Built-in robustness checks

This design prioritizes **scientific defensibility and reproducibility** over black-box modeling.

---

**Next:** See [Outputs & Interpretation](outputs.md) for how to read and interpret the results produced by the pipeline.
