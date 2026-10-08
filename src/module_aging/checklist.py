"""Diagnostic checklist for cross-cohort non-replication (revised manuscript, Table 2).

The steps are applied in order and stop at the first decision:

1. Do the effect vectors agree between cohorts?  Pearson r between the module effect
   vectors of the two cohorts, with a two-sided p from permuting disease labels within
   each cohort. r > 0 and p < alpha: "replicated".
2. Is the primary cohort's effect vector reliable at this sample size?  Median r between
   disease-stratified random halves of the primary cohort, against tau, the 95th
   percentile of the same statistic after label permutation. r < tau: "inconclusive (power)".
3. Does matching cell sub-states restore agreement?  Cross-cohort r within marker-defined
   sub-states (Bonferroni-adjusted p) compared with the r of random splits of equal size.
   p < alpha and r above the random-split r: "composition".
4. Otherwise: "between-cohort difference (unresolved)". The checklist cannot separate
   technical from biological sampling causes.

Inputs are sample x module score matrices with a case indicator per sample. A sample can
be a donor or a donor x fine-cell-type pseudobulk: pass ``strata`` (fine cell type) to
estimate the effect with stratum fixed effects, and ``donors`` to permute and split at the
donor level, as in the manuscript. In simulation the checklist called weak shared effects
a between-cohort difference in 41% of cohort pairs and a pure composition shift replicated
in 44% (Section 3.8), so its decisions are diagnostic, not proof.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

REPLICATED = "replicated"
INCONCLUSIVE = "inconclusive (power)"
COMPOSITION = "composition"
UNRESOLVED = "between-cohort difference (unresolved)"


def _donor_index(case: np.ndarray, donors) -> tuple[np.ndarray, np.ndarray]:
    """Map samples to donors; return (donor index per sample, case label per donor)."""
    if donors is None:
        return np.arange(len(case)), case.copy()
    _, inv = np.unique(np.asarray(donors), return_inverse=True)
    lab = np.zeros(inv.max() + 1)
    lab[inv] = case
    if not np.array_equal(lab[inv], case):
        raise ValueError("the case label must be constant within each donor")
    return inv, lab


def effect_vector(scores, case, strata=None) -> np.ndarray:
    """Case-minus-control effect per module (OLS, with stratum fixed effects if given)."""
    X = np.asarray(scores, float)
    y = np.asarray(case, float)
    if X.ndim != 2 or X.shape[0] != y.shape[0]:
        raise ValueError("scores must be (n_samples, n_modules) and match case")
    if not np.isfinite(X).all():
        raise ValueError("scores must be finite")
    if strata is None:
        D = np.ones((len(y), 1))
    else:
        _, s = np.unique(np.asarray(strata), return_inverse=True)
        D = np.eye(s.max() + 1)[s]
    return np.linalg.lstsq(np.column_stack([y, D]), X, rcond=None)[0][0]


def _r(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a, b)[0, 1])


def _permuted(case: np.ndarray, donors, rng) -> np.ndarray:
    inv, lab = _donor_index(case, donors)
    return rng.permutation(lab)[inv]


def cross_cohort_test(scores_a, case_a, scores_b, case_b, *, strata_a=None, strata_b=None,
                      donors_a=None, donors_b=None, n_perm: int = 2000, seed: int = 0) -> tuple[float, float]:
    """Step 1: Pearson r between the cohorts' effect vectors and its two-sided permutation p."""
    rng = np.random.default_rng(seed)
    ya, yb = np.asarray(case_a, float), np.asarray(case_b, float)
    r = _r(effect_vector(scores_a, ya, strata_a), effect_vector(scores_b, yb, strata_b))
    null = np.array([_r(effect_vector(scores_a, _permuted(ya, donors_a, rng), strata_a),
                        effect_vector(scores_b, _permuted(yb, donors_b, rng), strata_b))
                     for _ in range(n_perm)])
    return r, float((1 + np.sum(np.abs(null) >= abs(r))) / (n_perm + 1))


def _split_half(scores, y, strata, donors, n_split, rng) -> float:
    X = np.asarray(scores, float)
    inv, lab = _donor_index(y, donors)
    cases, ctrls = np.where(lab == 1)[0], np.where(lab == 0)[0]
    st = None if strata is None else np.asarray(strata)
    rs = []
    for _ in range(n_split):
        a, b = rng.permutation(cases), rng.permutation(ctrls)
        half = np.zeros(len(lab), bool)
        half[a[: len(a) // 2]] = True
        half[b[: len(b) // 2]] = True
        m = half[inv]
        e1 = effect_vector(X[m], y[m], None if st is None else st[m])
        e2 = effect_vector(X[~m], y[~m], None if st is None else st[~m])
        rs.append(_r(e1, e2))
    return float(np.median(rs))


def split_half_reliability(scores, case, *, strata=None, donors=None, n_split: int = 200,
                           n_perm: int = 200, n_split_null: int = 20, q: float = 95,
                           seed: int = 0) -> tuple[float, float]:
    """Step 2: median split-half r of the primary cohort and its permutation threshold tau."""
    rng = np.random.default_rng(seed)
    y = np.asarray(case, float)
    r = _split_half(scores, y, strata, donors, n_split, rng)
    null = [_split_half(scores, _permuted(y, donors, rng), strata, donors, n_split_null, rng)
            for _ in range(n_perm)]
    return r, float(np.nanpercentile(null, q))


@dataclass
class ChecklistResult:
    step: int
    decision: str
    detail: dict = field(default_factory=dict)


def checklist(r_cross: float, p_cross: float, r_splithalf: float, tau: float, *,
              r_substate: Optional[float] = None, p_substate: Optional[float] = None,
              r_random_split: Optional[float] = None, alpha: float = 0.05) -> ChecklistResult:
    """Apply Table 2 to precomputed statistics.

    ``p_substate`` is the Bonferroni-adjusted p over sub-states and ``r_substate`` the mean
    cross-cohort r over sub-states; ``r_random_split`` is the mean r of random splits of the
    same sizes. Omit the sub-state arguments if no sub-state split is available.
    """
    if r_cross > 0 and p_cross < alpha:
        return ChecklistResult(1, REPLICATED, {"r_cross": r_cross, "p_cross": p_cross})
    if r_splithalf < tau:
        return ChecklistResult(2, INCONCLUSIVE, {"r_splithalf": r_splithalf, "tau": tau})
    if (r_substate is not None and p_substate is not None and r_random_split is not None
            and p_substate < alpha and r_substate > r_random_split):
        return ChecklistResult(3, COMPOSITION, {"r_substate": r_substate, "p_substate": p_substate,
                                                "r_random_split": r_random_split})
    return ChecklistResult(4, UNRESOLVED, {"r_cross": r_cross, "r_splithalf": r_splithalf, "tau": tau,
                                           "r_substate": r_substate, "r_random_split": r_random_split})
