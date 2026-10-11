# Changelog

## 0.3.0

The manuscript now uses the serialized module clocks from the original study.
The module clocks retrained for v0.2.0 are withdrawn and are not in this tree.
Training data, serialized pipelines, gene-level coefficients and module membership
are not redistributed; `scripts/setup_v030.py --verify-only` checks the checksums
of copies you already have. Result tables are in `results/`, figures for the
revised manuscript are in `figures/v17/`, and `python verify/verify_claims.py --r2`
recomputes the numbers that could be tied to a supplementary-table cell.

## 0.2.0

Revised manuscript: `scripts/revision/`, `results/revision/`, `figures_revision/`.

## 0.1.0

Original submission.
