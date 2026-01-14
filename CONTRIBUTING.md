# Contributing to I-SAGE

Thank you for your interest in contributing to I-SAGE.

This document describes how to propose changes, report issues, and contribute code in a way that maintains pipeline correctness, reproducibility, and scientific validity.

---

## Contribution Principles

- Correctness over speed
- Reproducibility over convenience
- Explicit assumptions over implicit behavior
- Incremental, reviewable changes

---

## Types of Contributions

We welcome:
- Bug fixes
- Performance improvements
- New analysis features
- Documentation improvements
- Validation and robustness extensions

All contributions should be aligned with the pipeline’s scientific goals.

---

## Development Workflow

1. Create a feature branch from `main`
2. Make **atomic commits** (one logical change per commit)
3. Use descriptive commit messages
4. Open a pull request with:
   - clear description of the change
   - motivation and expected impact
   - notes on backward compatibility

---

## Coding Standards

### Nextflow
- Keep processes modular
- Do not hardcode paths
- Avoid implicit globals
- Ensure outputs are explicitly declared

### Python
- Python ≥ 3.9
- Prefer explicit imports
- Avoid side effects at import time
- Add comments for non-obvious logic
- Ensure scripts can run headless (HPC-safe plotting)

---

## Configuration Changes

Any change that modifies:
- default parameters
- output formats
- statistical behavior

**must** be reflected in:
- `configs/iblesse.yaml`
- documentation
- commit message

---

## Testing Expectations

Before submitting:
- Run the pipeline (or affected module) on a small test dataset
- Ensure outputs are generated as expected
- Confirm no silent changes to statistical logic

---

## Scientific Changes

Changes affecting:
- statistical tests
- replicate handling
- normalization
- thresholds

must include:
- clear rationale
- references where appropriate
- discussion of limitations

---

## Documentation

All new features must be documented.
Undocumented features will not be merged.

---

## Questions

For questions about design or scope, open an issue before implementing major changes.
