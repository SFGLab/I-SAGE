# I-SAGE — iBLESS Analysis Pipeline

I-SAGE is a reproducible, Nextflow-based pipeline for analyzing **iBLESS** sequencing data and performing **genome-wide differential DNA double-strand break (DSB) analysis**, particularly under replication stress conditions.

The pipeline integrates alignment, break calling, normalization, visualization, differential statistics, validation, and sensitivity analyses into a single configurable workflow.

---

## Table of Contents

1. [Key Features](#key-features)
2. [Pipeline Overview](#pipeline-overview)
3. [Requirements](#requirements)
4. [Quick Start](#quick-start)
5. [Configuration](#configuration)
6. [Differential Statistics](#differential-statistics)
7. [EBV Annotation](#ebv-annotation-optional)
8. [Outputs](#outputs)
9. [Documentation](#documentation)
10. [Project Status](#project-status)
11. [License](#license)

---

## Key Features

- Per-base iBLESS break calling
- Strand-aware break aggregation
- Configurable binning for visualization and statistics
- Replicate-aware differential break analysis
- Genome-wide or region-restricted testing (BED)
- Automated validation via downsampling and spike-in
- Bin-size sensitivity analysis
- Optional EBV contig annotation and enrichment
- Publication- and genome-browser–ready outputs

---

## Documentation
[![Documentation](https://img.shields.io/badge/docs-online-blue)](https://sfglab.github.io/I-SAGE/)

The full documentation for **I-SAGE** is available as a hosted website:

👉 **https://sfglab.github.io/I-SAGE/**

The documentation includes:
- Getting started guide
- Full configuration reference
- Pipeline module descriptions
- Statistical methods and assumptions
- Output interpretation
- Developer and contribution guidelines



## Pipeline Overview

```
FASTQ
  ↓
Alignment & Deduplication
  ↓
Break Calling (per-base)
  ↓
Visualization Tracks (binned bedGraph)
  ↓
Normalization
  ↓
Differential Break Statistics
  ↓
Validation & Sensitivity Analyses
```

---

## Requirements

### Software

- Nextflow ≥ 22
- Java ≥ 11
- Python ≥ 3.9

### Python Packages

- numpy
- pandas
- scipy
- matplotlib

**Note:** Packages are typically managed via environment/profile on HPC systems.

---

## Quick Start

```bash
nextflow run workflows/iblesse_month2/main.nf \
  -profile eden_local \
  -params-file configs/iblesse.yaml
```

For detailed setup instructions, see the [Configuration](#configuration) section below.

---

## Configuration

Pipeline behavior is controlled via a YAML configuration file (e.g., `configs/iblesse.yaml`).

### Main Configuration Sections

#### 1. Break Calling (`break_calling`)
- Controls per-base or binned break calling parameters

#### 2. Visualization (`viz`)
- Visualization bin size or bin-size sweep settings

#### 3. Statistics (`stats`)
- Contrasts specification
- Replicate handling
- FDR thresholds
- EBV annotation options

#### 4. Validation (`validation`)
- Downsampling parameters
- Spike-in validation settings

A fully working example is provided in `configs/iblesse.yaml`.

---

## Differential Statistics

The pipeline performs robust statistical analysis across genomic bins:

### Method
- **Per-bin Fisher exact test** for break count differences
- **Benjamini–Hochberg FDR correction** for multiple testing

### Advanced Options
- Replicate-aware analysis via Fisher meta-analysis

### Output Files
- All tested bins (full results)
- Significant bins (filtered by FDR threshold)
- Upregulated / downregulated bins (directional results)
- Volcano and MA plots (PNG + PDF formats)

---

## EBV Annotation (Optional)

If EBV contigs are present in the reference genome, I-SAGE can:

### Capabilities
- Annotate bins as EBV vs. non-EBV
- Quantify EBV enrichment among significant bins
- Report enrichment statistics in summary files

### Enable EBV Analysis

Add the following to your `configs/iblesse.yaml`:

```yaml
stats:
  ebv_regex: "(?i)^chrEBV$"
```

---

## Outputs

Results are organized under `outdir/`:

### Output Directory Structure

```
outdir/
├── viz/              # bedGraph tracks (per bin size)
├── stats/            # Differential statistics, plots, summaries
├── validation/       # Robustness and sensitivity analyses
├── reports/          # HTML reports and execution traces
└── logs/             # Pipeline logs
```

### Key Output Files

- **Visualization tracks** – BigWig/bedGraph format for genome browsers
- **Statistical tables** – TSV files with bin-level results
- **Plots** – Volcano plots, MA plots, and heatmaps (PNG + PDF)
- **Summary reports** – HTML and text-based summaries

---

## Documentation

Full documentation is under development and will be available in the `documentation/` directory, including:

- Pipeline architecture and design
- Module-level descriptions
- Statistical methods and validation
- Configuration guide and best practices
- Output interpretation and usage

---

## Project Status

- **Phase:** Active development (Month 4)
- **Stability:** Core APIs and outputs are stabilizing
- **Development roadmap:** See `CHANGELOG.md` for recent updates

---

## Citation / Acknowledgment

If this tool supports your work, please cite the repository and acknowledge:
**“Developed by Pranjul Mishra, under the guidance of Dr. Joanna Borkowska and Prof. Dariusz Plewczyński (Structural and Functional Genomics Laboratory).”**

## License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.

---

## Support & Contribution

For issues, questions, or contributions, please refer to the project repository.