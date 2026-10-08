#!/usr/bin/env python3
"""Quickstart: apply the diagnostic checklist (manuscript Table 2) to two cohorts.

Runs on small synthetic data (no download required) to show the API. For real data,
replace the arrays with your own sample x module score matrices and case labels; pass
``strata`` (fine cell type) and ``donors`` when samples are donor x cell-type pseudobulks.

Run:
    python examples/reuse_quickstart.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# make `module_aging` importable from a source checkout without installing
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from module_aging import checklist, cross_cohort_test, split_half_reliability  # noqa: E402

rng = np.random.default_rng(0)
N_MODULES = 23


def cohort(effect, n_case, n_ctrl, noise=1.0):
    """Donor x module scores with a case-control effect vector."""
    case = np.r_[np.ones(n_case), np.zeros(n_ctrl)]
    return case[:, None] * effect[None, :] + rng.normal(0, noise, (len(case), N_MODULES)), case


def run(name, scores_a, case_a, scores_b, case_b):
    r, p = cross_cohort_test(scores_a, case_a, scores_b, case_b, n_perm=1000, seed=1)
    rs, tau = split_half_reliability(scores_a, case_a, n_split=100, n_perm=100, seed=2)
    res = checklist(r, p, rs, tau)
    print(f"{name}\n  cross-cohort r = {r:.2f} (two-sided p = {p:.3f}); split-half r = {rs:.2f}, tau = {tau:.2f}"
          f"\n  step {res.step}: {res.decision}")


shared = rng.normal(0, 1, N_MODULES)
run("Shared effect in both cohorts", *cohort(shared, 30, 30), *cohort(shared, 12, 10))
run("Stable effect in cohort A, different effect in cohort B",
    *cohort(shared, 30, 30), *cohort(rng.normal(0, 1, N_MODULES), 12, 10))
run("No effect", *cohort(np.zeros(N_MODULES), 30, 30), *cohort(np.zeros(N_MODULES), 12, 10))
