"""Reusable attribution-ladder primitives (manuscript Table 1).

A clock-agnostic protocol to diagnose *why* a cell-type-resolved module effect
vector fails to replicate across two cohorts. Given per-module effect vectors
(and, for the power stage, donor-level effect matrices) this module provides the
statistical primitives used in the paper and a decision helper that returns the
Table-1 verdict.

The functions here operate on plain NumPy/pandas inputs and carry **no**
dependency on the transcriptomic clocks themselves, so they can be reused to
diagnose non-replication of any decomposed single-cell signature.

Stages (applied in order, stop at the first that resolves the discrepancy):
    0. observe      cross-cohort concordance vs disease-label permutation null
    1. technical    re-score on the common detected-gene support (symmetric mask)
    2. power        within-cohort split-half "ceiling" at the replication n
    3. composition  matched marker-defined sub-states
    4. resolution   whole-donor pseudobulk / independent bulk cohort
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

try:
    from scipy.stats import pearsonr, spearmanr
except Exception:  # pragma: no cover - scipy always present in practice
    pearsonr = spearmanr = None


# --------------------------------------------------------------------------- #
# Concordance metrics
# --------------------------------------------------------------------------- #
@dataclass
class Concordance:
    pearson_r: float
    spearman_r: float
    sign_concordance: float
    n: int


def effect_concordance(a: np.ndarray, b: np.ndarray) -> Concordance:
    """Concordance between two per-module effect vectors.

    Parameters
    ----------
    a, b : array-like
        Module effect estimates (e.g. IPF-vs-control beta) for the *same* set of
        modules, from cohort A and cohort B. NaNs are dropped pairwise.

    Returns
    -------
    Concordance
        Pearson r, Spearman r, sign-concordance (fraction of modules with the
        same sign), and the number of modules compared.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if a.size < 3:
        raise ValueError("need >= 3 finite paired values")
    pr = float(pearsonr(a, b)[0]) if pearsonr else float(np.corrcoef(a, b)[0, 1])
    sr = float(spearmanr(a, b)[0]) if spearmanr else float("nan")
    sign = float(np.mean(np.sign(a) == np.sign(b)))
    return Concordance(pearson_r=pr, spearman_r=sr, sign_concordance=sign, n=int(a.size))


def concordance_permutation_p(
    a: np.ndarray, b: np.ndarray, n_perm: int = 2000, seed: int = 0
) -> float:
    """One-sided permutation p-value for positive cross-cohort concordance.

    The module labels of ``b`` are shuffled ``n_perm`` times to build the null
    distribution of Pearson r under no module-wise correspondence. Mirrors the
    disease-label permutation null used in the paper at the effect-vector level.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    obs = np.corrcoef(a, b)[0, 1]
    rng = np.random.default_rng(seed)
    ge = 1  # +1 for the observed statistic (avoids p = 0)
    for _ in range(n_perm):
        if np.corrcoef(a, rng.permutation(b))[0, 1] >= obs:
            ge += 1
    return ge / (n_perm + 1)


def split_half_ceiling(
    donor_effects_a: np.ndarray,
    labels: np.ndarray,
    n_repeats: int = 100,
    seed: int = 0,
) -> dict:
    """Within-cohort split-half reproducibility ceiling.

    Repeatedly partition donors into two disease-matched halves, recompute the
    per-module effect vector in each half (mean over cases minus mean over
    controls), and correlate the two halves. The median gives the best
    achievable cross-cohort concordance at this sample size (the "power ceiling"
    of Table-1 stage 2).

    Parameters
    ----------
    donor_effects_a : ndarray, shape (n_donors, n_modules)
        Per-donor module scores for one cohort.
    labels : array-like of {0,1}, shape (n_donors,)
        Disease label per donor (1 = case, 0 = control).
    """
    X = np.asarray(donor_effects_a, float)
    y = np.asarray(labels).astype(int)
    if X.ndim != 2 or X.shape[0] != y.shape[0]:
        raise ValueError("donor_effects_a must be (n_donors, n_modules) matching labels")
    rng = np.random.default_rng(seed)
    case, ctrl = np.where(y == 1)[0], np.where(y == 0)[0]
    if len(case) < 2 or len(ctrl) < 2:
        raise ValueError("need >= 2 cases and >= 2 controls")

    def _effect(idx_case, idx_ctrl):
        return X[idx_case].mean(0) - X[idx_ctrl].mean(0)

    rs = []
    for _ in range(n_repeats):
        c = rng.permutation(case)
        k = rng.permutation(ctrl)
        e1 = _effect(c[: len(c) // 2], k[: len(k) // 2])
        e2 = _effect(c[len(c) // 2 :], k[len(k) // 2 :])
        ok = np.isfinite(e1) & np.isfinite(e2)
        if ok.sum() >= 3:
            rs.append(np.corrcoef(e1[ok], e2[ok])[0, 1])
    rs = np.array(rs)
    return {
        "median": float(np.median(rs)),
        "lo": float(np.percentile(rs, 2.5)),
        "hi": float(np.percentile(rs, 97.5)),
        "n_repeats": int(len(rs)),
    }


# --------------------------------------------------------------------------- #
# Table-1 decision helper
# --------------------------------------------------------------------------- #
@dataclass
class LadderVerdict:
    stage: str
    conclusion: str
    detail: dict = field(default_factory=dict)


def classify(
    r_cross: float,
    p_cross: float,
    *,
    r_harmonized: Optional[float] = None,
    splithalf_median: Optional[float] = None,
    r_substate: Optional[float] = None,
    r_bulk: Optional[float] = None,
    p_bulk: Optional[float] = None,
    r_recover: float = 0.3,
    p_alpha: float = 0.05,
    splithalf_high: float = 0.5,
) -> LadderVerdict:
    """Apply the Table-1 decision tree to precomputed stage statistics.

    Only supply the stages you have run; the ladder stops at the first stage
    that resolves the discrepancy. Thresholds are exposed as keyword arguments.
    """
    # Stage 0 -- is there genuine non-replication to explain?
    if r_cross >= r_recover and p_cross < p_alpha:
        return LadderVerdict("0. observe", "replicates (no attribution needed)",
                             {"r_cross": r_cross, "p_cross": p_cross})

    # Stage 1 -- technical (gene coverage / imputation)
    if r_harmonized is not None and r_harmonized >= r_recover:
        return LadderVerdict("1. technical",
                             "recovered on common gene support -> technical/coverage artifact",
                             {"r_harmonized": r_harmonized})

    # Stage 2 -- power vs genuine between-cohort difference
    power_ruled_out = False
    if splithalf_median is not None:
        if splithalf_median < splithalf_high:
            return LadderVerdict("2. power",
                                 "low within-cohort ceiling -> power-limited (insufficient n "
                                 "to estimate the effect)",
                                 {"splithalf_median": splithalf_median})
        # high ceiling: signal is stable within-cohort -> NOT power-limited; keep going
        power_ruled_out = True

    # Stage 3 -- composition (sub-state mixing)
    if r_substate is not None and r_substate >= r_recover:
        return LadderVerdict("3. composition",
                             "recovered on matched sub-states -> composition/sub-state mixing",
                             {"r_substate": r_substate})

    # Stage 4 -- resolution (recover at bulk / whole-donor)
    if r_bulk is not None and r_bulk >= r_recover and (p_bulk is None or p_bulk < p_alpha):
        return LadderVerdict("4. resolution",
                             "recovers without cell-type resolution -> resolution-dependent "
                             "(report at bulk/protocol-matched level only)",
                             {"r_bulk": r_bulk, "p_bulk": p_bulk})

    if power_ruled_out:
        return LadderVerdict("2. power (ruled out)",
                             "stable within-cohort (high split-half) yet fails across cohorts "
                             "-> NOT power-limited; genuine between-cohort (sampling) difference. "
                             "Supply composition/resolution stages to localize further.",
                             {"splithalf_median": splithalf_median, "r_cross": r_cross})

    return LadderVerdict("unresolved",
                         "not resolved by the supplied stages; report as genuine non-replication",
                         {"r_cross": r_cross, "p_cross": p_cross})
