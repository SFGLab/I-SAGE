# Developer Guide

## Overview

This guide is a reference for developers and maintainers of the I-SAGE pipeline. It documents design principles, code layout, statistical choices, extension points, and testing expectations.

## Table of Contents

1. [Design Philosophy](#design-philosophy)
2. [Code Organization](#code-organization)
3. [Why Nextflow](#why-nextflow)
4. [Statistical Decisions](#statistical-decisions)
5. [Replicate Handling](#replicate-handling)
6. [Bin-Size Policy](#bin-size-policy)
7. [EBV Annotation](#ebv-annotation)
8. [Region-Restricted Testing](#region-restricted-testing)
9. [Validation](#validation)
10. [Adding a Module](#adding-a-module)
11. [Testing & Release Checklist](#testing--release-checklist)
12. [When to Say No](#when-to-say-no)
13. [Contact & Maintainers](#contact--maintainers)

## Design Philosophy

Principles that should guide all development:

- **Reproducibility over convenience** — All results must be fully reproducible
- **Transparency and interpretability over black-box methods** — Users must understand what happened
- **Modularity** — Each stage is a self-contained, reviewable unit
- **Explicit assumptions** — Document assumptions alongside code

## Code Organization

Recommended repository layout:

```
workflows/
  iblesse_month2/
    main.nf
modules/
  break_calling/
  viz/
  stats/
  validation/
configs/
docs/
```

Guidelines:

- `main.nf` orchestrates processes and propagates parameters
- Each directory in `modules/` implements a single conceptual stage with explicit inputs/outputs
- Statistical and validation logic lives in Python scripts under `modules/` with unit tests

## Why Nextflow

Rationale:

- **Reproducible execution** — Works on local and HPC environments
- **Clear separation** — Workflow logic distinct from tool implementation
- **Built-in features** — Caching, provenance, and traceability

**Important:** Do not bypass Nextflow by embedding execution logic in ad-hoc scripts.

## Statistical Decisions

### Primary Choices

- **Per-bin Fisher exact test** for differential break counts
- **Benjamini–Hochberg FDR** for multiple-testing correction

### Rationale

- iBLESS counts are sparse; many bins have low or zero counts
- Fisher's test is exact and assumption-light, making it interpretable and robust for sparse data

**Note:** Document any deviation from this approach with motivation and tests demonstrating improved behavior.

## Replicate Handling

### Supported Strategies

- **`pooled`** — Sum counts across replicates
- **`meta_fisher`** — Combine per-replicate p-values via Fisher's method

### Rationale

- Avoids fitting unstable dispersion models with few replicates
- Keeps analysis conservative and interpretable

**Note:** If introducing GLMs or shrinkage estimators, include benchmarks and documented justification.

## Bin-Size Policy

- Treat bin size as an explicit parameter (biological and statistical)
- Support bin-size sweeps to assess sensitivity
- Avoid hardcoded single-bin assumptions in downstream logic

Provide scripts to aggregate and summarize results across bin sizes when adding new modules.

## EBV Annotation

- EBV support is **optional and opt-in** via configuration (e.g., `stats.ebv_regex`)
- EBV is used for **annotation and enrichment reporting**, not for filtering or mechanistic claims
- When EBV is enabled, include clear reporting on how EBV bins were identified

## Region-Restricted Testing

- Region restriction narrows the tested universe and recomputes totals within that universe
- Distinguish region-restricted testing from post-hoc enrichment or overlap analyses in documentation and UIs

## Validation

Built-in validation strategies:

- **Downsampling experiments** — Assess robustness with reduced data
- **Spike-in recovery tests** — Validate detection sensitivity
- **Bin-size sensitivity analysis** — Understand parameter sensitivity

**Important:** Validation outputs are diagnostic; they inform interpretation but should not be used to automatically filter results without explicit justification.

## Adding a Module

Follow these steps:

1. Define scientific purpose and acceptance criteria
2. Implement the process under `modules/<name>/` with clear inputs/outputs
3. Add or update unit tests and example data
4. Wire the process in `main.nf` with configuration knobs
5. Document module in [Pipeline Modules](pipeline_modules.md) and update config example
6. Run integration test (small dataset) and add to CI if applicable

**Merge condition:** Documentation and tests must be present before merging.

## Testing & Release Checklist

Before merging:

- [ ] Run pipeline on a small, representative dataset
- [ ] Confirm outputs match expectations and reproducibility (same inputs → same outputs)
- [ ] Update documentation for any user-facing or config changes
- [ ] Add/adjust unit tests and integration tests
- [ ] Record changes in CHANGELOG.md

## When to Say No

Reject changes that:

- Add complexity without biological benefit
- Introduce opaque statistical methods without documentation and tests
- Encourage overfitting or p-hacking through undocumented parameters

**Philosophy:** Prefer conservative, explainable solutions.

## Contact & Maintainers

**Primary Maintainer:** Pranjul Mishra

**Supervisors:** Dr. Joanna Borkowska, Prof. Dariusz Plewczyński

For design decisions, open an issue and tag the maintainers. Include rationale, tests, and example outputs.

---

**See Also:** [Contributing Guidelines](../CONTRIBUTING.md) | [Code of Conduct](../CODE_OF_CONDUCT.md)
