# Developer Guide

Purpose
- Reference for developers and maintainers of the I-SAGE pipeline.
- Documents design principles, code layout, statistical choices, extension points, and testing expectations.

Table of contents
1. Design philosophy
2. Code organization
3. Why Nextflow
4. Statistical decisions
5. Replicate handling
6. Bin-size policy
7. EBV annotation
8. Region-restricted testing
9. Validation
10. Adding a module
11. Testing & release checklist
12. When to say no
13. Contact

---

## 1. Design philosophy
Principles that should guide all development:
- Reproducibility over convenience
- Transparency and interpretability over black-box methods
- Modularity: each stage is a self-contained, reviewable unit
- Explicit assumptions documented alongside code

---

## 2. Code organization

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
documentation/
```

Guidelines:
- main.nf orchestrates processes and propagates parameters.
- Each directory in modules/ implements a single conceptual stage with explicit inputs/outputs.
- Statistical and validation logic lives in Python scripts under modules/ with unit tests.

---

## 3. Why Nextflow
Rationale:
- Reproducible execution on local and HPC environments
- Clear separation of workflow logic and tool implementation
- Built-in caching, provenance, and traceability

Do not bypass Nextflow by embedding execution logic in ad-hoc scripts.

---

## 4. Statistical decisions

Primary choices:
- Per-bin Fisher exact test for differential break counts
- Benjamini–Hochberg FDR for multiple-testing correction

Rationale:
- iBLESS counts are sparse; many bins have low/zero counts
- Fisher's test is exact and assumption-light, hence interpretable and robust for sparse data

Document any deviation from this approach with motivation and tests demonstrating improved behaviour.

---

## 5. Replicate handling

Supported strategies:
- pooled: sum counts across replicates
- meta_fisher: combine per-replicate p-values via Fisher's method

Rationale:
- Avoids fitting unstable dispersion models with few replicates
- Keeps analysis conservative and interpretable

If introducing GLMs or shrinkage estimators, include benchmarks and documented justification.

---

## 6. Bin-size policy

- Treat bin size as an explicit parameter (biological and statistical).
- Support bin-size sweeps to assess sensitivity.
- Avoid hardcoded single-bin assumptions in downstream logic.

Provide scripts to aggregate and summarize results across bin sizes when adding new modules.

---

## 7. EBV annotation

- EBV support is optional and opt-in via configuration (e.g., stats.ebv_regex).
- EBV is used for annotation and enrichment reporting, not for filtering or mechanistic claims.
- When EBV is enabled, include clear reporting on how EBV bins were identified.

---

## 8. Region-restricted testing

- Region restriction narrows the tested universe and recomputes totals within that universe.
- Distinguish region-restricted testing from post-hoc enrichment or overlap analyses in docs and UIs.

---

## 9. Validation

Built-in validation strategies:
- Downsampling experiments
- Spike-in recovery tests
- Bin-size sensitivity analysis

Validation outputs are diagnostic; they inform interpretation but should not be used to automatically filter results without explicit justification.

---

## 10. Adding a module

Steps:
1. Define scientific purpose and acceptance criteria.
2. Implement the process under modules/<name>/ with clear inputs/outputs.
3. Add or update unit tests and example data.
4. Wire the process in main.nf with configuration knobs.
5. Document module in documentation/pipeline_modules.md and update configs example.
6. Run the integration test (small dataset) and add to CI if applicable.

Merge only if documentation and tests are present.

---

## 11. Testing & release checklist

Before merging:
- Run pipeline on a small, representative dataset.
- Confirm outputs match expectations and reproducibility (same inputs → same outputs).
- Update documentation for any user-facing or config changes.
- Add/adjust unit tests and integration tests.
- Record changes in CHANGELOG.md.

---

## 12. When to say no

Reject changes that:
- Add complexity without biological benefit
- Introduce opaque statistical methods without documentation and tests
- Encourage overfitting or p-hacking through undocumented parameters

Prefer conservative, explainable solutions.

---

## 13. Contact / Maintainers

Primary maintainer: Pranjul Mishra  
Supervisors: Dr. Joanna Borkowska, Prof. Dariusz Plewczyński

For design decisions, open an issue and tag the maintainers. Include rationale, tests, and example outputs.

---
<!-- End of developer guide -->

