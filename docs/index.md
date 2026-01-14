```markdown
# I-SAGE Documentation

Welcome to the documentation for **I-SAGE**, a reproducible Nextflow-based pipeline for analyzing **iBLESS** sequencing data and performing **genome-wide differential DNA double-strand break (DSB) analysis**.

This documentation is intended for:
- Computational biologists running iBLESS analyses
- Developers extending or maintaining the pipeline
- Reviewers seeking clarity on statistical methods and assumptions

The goal is to make every component of the pipeline **transparent, reproducible, and interpretable**.

---

## What is I-SAGE?

I-SAGE is an end-to-end analysis framework that transforms raw iBLESS FASTQ files into:
- strand-aware break tracks
- normalized genome-wide signals
- statistically validated differential break calls
- robustness and sensitivity assessments

It is designed to support **replication stress experiments**, including treatments such as HU and aphidicolin, and to scale from single contrasts to replicate-aware, multi-parameter analyses.

---

## Design Principles

I-SAGE is built around the following principles:

- **Reproducibility first**  
  All steps are parameterized, logged, and traceable via Nextflow.

- **Explicit statistics**  
  Statistical tests, assumptions, and thresholds are documented and configurable.

- **Modularity**  
  Each pipeline stage is implemented as an independent module.

- **Scientific defensibility**  
  Validation steps (downsampling, spike-ins, bin-size sensitivity) are part of the pipeline, not an afterthought.

---

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

Each stage is described in detail in the following sections.

---

## Documentation Structure

This documentation is organized as follows:

### 🔹 Getting Started
- Installation and requirements
- Quick-start examples
- Configuration basics

### 🔹 Pipeline Modules
- Alignment and deduplication
- Break calling
- Visualization and normalization
- Differential statistics
- Validation and sensitivity analysis

### 🔹 Statistical Methods
- Bin-level testing framework
- Replicate-aware analysis
- Multiple testing correction
- EBV annotation and enrichment

### 🔹 Configuration Guide
- Full explanation of `iblesse.yaml`
- Parameter interactions
- Recommended settings

### 🔹 Outputs & Interpretation
- Output directory structure
- TSV files and plots
- Genome browser tracks
- Common pitfalls in interpretation

### 🔹 Developer Guide
- Code structure
- Adding new modules
- Testing expectations
- Contribution workflow

---

## Intended Usage

I-SAGE is suitable for:
- Genome-wide DSB profiling under replication stress
- Comparing treated vs control conditions
- Assessing reproducibility across replicates
- Sensitivity analyses across bin sizes or genomic regions

It is **not** intended to replace low-level exploratory notebooks, but to formalize analyses that must be **re-run, reviewed, and trusted**.

---

## Project Status

The pipeline is under active development.  
Core functionality is stable; documentation is being expanded.

Major changes are tracked via Git commits and documented in the repository.

---

## Getting Help

- For usage questions: consult the relevant documentation sections
- For bugs or feature requests: open a GitHub issue
- For design discussions: see `CONTRIBUTING.md`

---

## Citation

If you use I-SAGE in your work, please cite the relevant methodological references and acknowledge the pipeline.

(A formal citation entry will be added upon publication.)

---

Next: proceed to **Getting Started → Installation & Requirements**.
```

---
