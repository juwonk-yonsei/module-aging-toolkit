#!/usr/bin/env python3
"""Quickstart: apply the attribution ladder to your own two-cohort data.

This runs on small synthetic data (no download required) to show the API. To
use it on real data, replace the synthetic arrays with your own per-module
effect vectors (one value per module, per cohort) and per-donor module-score
matrices.

Run:
    python examples/reuse_quickstart.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# make `module_aging` importable from a source checkout without installing
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from module_aging import (  # noqa: E402
    classify,
    concordance_permutation_p,
    effect_concordance,
    split_half_ceiling,
)

rng = np.random.default_rng(0)
N_MODULES = 23

# --------------------------------------------------------------------------- #
# Scenario A: a lymphoid-like cell type that REPLICATES across cohorts
# --------------------------------------------------------------------------- #
truth = rng.normal(0, 0.5, N_MODULES)
eff_A = truth + rng.normal(0, 0.15, N_MODULES)
eff_B = truth + rng.normal(0, 0.15, N_MODULES)  # correlated with A

conc = effect_concordance(eff_A, eff_B)
p = concordance_permutation_p(eff_A, eff_B, n_perm=2000, seed=1)
print("Scenario A (replicates):")
print(f"  cross-cohort Pearson r = {conc.pearson_r:.2f}, sign = {conc.sign_concordance:.0%}, "
      f"permutation p = {p:.4f}")
print("  verdict:", classify(conc.pearson_r, p).conclusion)

# --------------------------------------------------------------------------- #
# Scenario B: a myeloid-like cell type that FAILS across cohorts but is
# stable within-cohort -> genuine between-cohort difference (not power)
# --------------------------------------------------------------------------- #
eff_A2 = rng.normal(0, 0.5, N_MODULES)
eff_B2 = rng.normal(0, 0.5, N_MODULES)  # independent -> no cross-cohort concordance
conc2 = effect_concordance(eff_A2, eff_B2)
p2 = concordance_permutation_p(eff_A2, eff_B2, n_perm=2000, seed=2)

# within-cohort split-half: 40 donors, strong stable module signal in cohort A
n_donors = 40
labels = np.r_[np.ones(20), np.zeros(20)].astype(int)
signal = rng.normal(0, 1, N_MODULES)
donor_scores = (labels[:, None] * signal[None, :]) + rng.normal(0, 0.4, (n_donors, N_MODULES))
sh = split_half_ceiling(donor_scores, labels, n_repeats=200, seed=3)

verdict = classify(conc2.pearson_r, p2, splithalf_median=sh["median"])
print("\nScenario B (fails cross-cohort, stable within-cohort):")
print(f"  cross-cohort Pearson r = {conc2.pearson_r:.2f}, permutation p = {p2:.4f}")
print(f"  within-cohort split-half ceiling = {sh['median']:.2f} "
      f"[{sh['lo']:.2f}, {sh['hi']:.2f}]")
print(f"  verdict [{verdict.stage}]: {verdict.conclusion}")

print("\nDone. Swap in your own effect vectors / donor matrices to diagnose "
      "non-replication in your data.")
