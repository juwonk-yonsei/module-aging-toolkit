"""Operating characteristics of the non-replication diagnostic checklist (Section 3.8, Supplementary Methods S3,
Supplementary Table S15a).

Stylized simulation at the donor-pseudobulk level for one cell type with two sub-states.
23 module scores per donor with the module correlation observed in the IPF data
(rev03_corr_ipf.csv); cohort sizes as for macrophages (primary 32 IPF / 27 Control,
replication 12 / 9); sub-state fractions as observed (monocyte-derived share 0.42 in
Control, 0.66 in IPF).

True scenarios
  shared             identical disease effect in both sub-states and cohorts
  power              as 'shared' but a small effect (0.25 SD)
  composition        opposite effects in the two sub-states and different sampling of
                     sub-states (monocyte-derived share 0.25 vs 0.80)
  composition_shift  no within-state effect; IPF raises the monocyte-derived share
                     (0.42 -> 0.66) in both cohorts
  confound           a stable disease-correlated shift present only in the primary
                     cohort (e.g. tissue source), no shared effect
  null               no effect

Checklist rule (explicit, applied in this order)
  1 cross-cohort r > 0 and significant (two-sided permutation p < 0.05)  -> replicated
  2 primary-cohort split-half r below the 95th percentile under 'null'   -> inconclusive (power)
  3 matched sub-state r significant and larger than a random-split r     -> composition
  4 otherwise                                                            -> between-cohort difference, unresolved

Outputs: rev12_checklist_confusion.csv, rev12_checklist_raw.csv
"""
import numpy as np
import pandas as pd

import rev_common as rc

N_REP = 300
B_PERM = 200
N_SPLIT = 20
NI1, NC1, NI2, NC2 = 32, 27, 12, 9


def group_diff(Y, lab_mat):
    """Y: n x M; lab_mat: K x n boolean IPF indicators -> K x M mean differences."""
    L = lab_mat.astype(float)
    ni = L.sum(1, keepdims=True)
    nc = (1 - L).sum(1, keepdims=True)
    return (L @ Y) / ni - ((1 - L) @ Y) / nc


def rowcorr(A, Bm):
    A = A - A.mean(1, keepdims=True)
    Bm = Bm - Bm.mean(1, keepdims=True)
    return (A * Bm).sum(1) / np.sqrt((A ** 2).sum(1) * (Bm ** 2).sum(1))


def perm_mats(lab, K, rng):
    return np.array([rng.permutation(lab) for _ in range(K)])


def simulate_cohort(rng, n_i, n_c, pi_c, pi_i, mu_mono, d_res, d_mono, Lc, shift=None, sd=1.0):
    lab = np.r_[np.ones(n_i, bool), np.zeros(n_c, bool)]
    n, M = len(lab), len(mu_mono)
    pi = np.clip(np.where(lab, pi_i, pi_c) + rng.normal(0, 0.08, n), 0.05, 0.95)
    eps_res = rng.standard_normal((n, M)) @ Lc.T * sd * 1.2
    eps_mono = rng.standard_normal((n, M)) @ Lc.T * sd * 1.2
    D = lab[:, None].astype(float)
    y_res = D * d_res + eps_res
    y_mono = mu_mono + D * d_mono + eps_mono
    if shift is not None:
        y_res = y_res + D * shift
        y_mono = y_mono + D * shift
    pooled = pi[:, None] * y_mono + (1 - pi[:, None]) * y_res
    # random split of the same cells: two halves with independent noise, same mixture
    half_noise = rng.standard_normal((n, M)) @ Lc.T * sd * 0.6
    rand_half = pooled + half_noise
    return lab, pooled, y_res, y_mono, rand_half


SCENARIOS = {  # effect size, opposite sub-state effects, (pi_ctrl, pi_ipf) cohort 1 / cohort 2, confound
    "shared": (0.6, False, (0.4, 0.4), (0.4, 0.4), False),
    "power": (0.25, False, (0.4, 0.4), (0.4, 0.4), False),
    "composition": (0.6, True, (0.25, 0.25), (0.8, 0.8), False),
    "composition_shift": (0.0, False, (0.42, 0.66), (0.42, 0.66), False),
    "confound": (0.0, False, (0.4, 0.4), (0.4, 0.4), True),
    "null": (0.0, False, (0.4, 0.4), (0.4, 0.4), False),
}


def run(scenario, rng, Lc, M):
    eff, opposite, pi1, pi2, confound = SCENARIOS[scenario]
    base = rng.normal(0, 1, M)
    base = base / np.sqrt(np.mean(base ** 2))
    mu_mono = rng.normal(0, 1.0, M)
    d = base * eff
    d_res, d_mono = (d, -d) if opposite else (d, d)
    shift = rng.normal(0, 1, M) * 0.6 if confound else None
    l1, P1, R1, Mo1, H1 = simulate_cohort(rng, NI1, NC1, *pi1, mu_mono, d_res, d_mono, Lc, shift)
    l2, P2, R2, Mo2, H2 = simulate_cohort(rng, NI2, NC2, *pi2, mu_mono, d_res, d_mono, Lc)
    e1, e2 = group_diff(P1, l1[None])[0], group_diff(P2, l2[None])[0]
    r = np.corrcoef(e1, e2)[0, 1]
    n1 = group_diff(P1, perm_mats(l1, B_PERM, rng))
    n2 = group_diff(P2, perm_mats(l2, B_PERM, rng))
    p_cross = (1 + np.sum(np.abs(rowcorr(n1, n2)) >= abs(r))) / (B_PERM + 1)
    # split-half in the primary cohort (disease-stratified halves)
    idx_i, idx_c = np.where(l1)[0], np.where(~l1)[0]
    rs = []
    for _ in range(N_SPLIT):
        a_i, a_c = rng.permutation(idx_i), rng.permutation(idx_c)
        h1 = np.r_[a_i[:len(a_i) // 2], a_c[:len(a_c) // 2]]
        h2 = np.r_[a_i[len(a_i) // 2:], a_c[len(a_c) // 2:]]
        rs.append(np.corrcoef(group_diff(P1[h1], l1[h1][None])[0], group_diff(P1[h2], l1[h2][None])[0])[0, 1])
    r_sh = float(np.median(rs))
    # matched sub-states
    r_sub, p_sub = [], []
    for Y1, Y2 in [(R1, R2), (Mo1, Mo2)]:
        a, b = group_diff(Y1, l1[None])[0], group_diff(Y2, l2[None])[0]
        rr = np.corrcoef(a, b)[0, 1]
        nn = rowcorr(group_diff(Y1, perm_mats(l1, B_PERM, rng)), group_diff(Y2, perm_mats(l2, B_PERM, rng)))
        r_sub.append(rr)
        p_sub.append((1 + np.sum(np.abs(nn) >= abs(rr))) / (B_PERM + 1))
    r_rand = np.corrcoef(group_diff(H1, l1[None])[0], group_diff(H2, l2[None])[0])[0, 1]
    return dict(scenario=scenario, r_cross=r, p_cross=p_cross, r_splithalf=r_sh,
                r_substate=float(np.mean(r_sub)), p_substate=float(np.min(p_sub) * 2),
                r_random_split=r_rand)


def classify(row, tau_sh):
    if row.p_cross < 0.05 and row.r_cross > 0:
        return "replicated"
    if row.r_splithalf < tau_sh:
        return "inconclusive (power)"
    if row.p_substate < 0.05 and row.r_substate > row.r_random_split:
        return "composition"
    return "between-cohort difference (unresolved)"


def main():
    rng = np.random.default_rng(rc.SEED)
    C = pd.read_csv(rc.OUT / "rev03_corr_ipf.csv", index_col=0).values
    w, V = np.linalg.eigh(C)
    Lc = V @ np.diag(np.sqrt(np.clip(w, 1e-6, None)))
    M = C.shape[0]
    rows = [run(s, rng, Lc, M) for s in SCENARIOS for _ in range(N_REP)]
    R = pd.DataFrame(rows)
    tau_sh = float(np.percentile(R[R.scenario == "null"].r_splithalf, 95))
    R["attribution"] = [classify(r, tau_sh) for r in R.itertuples()]
    rc.save(R, "rev12_checklist_raw.csv")
    Cm = pd.crosstab(R.scenario, R.attribution, normalize="index").reindex(list(SCENARIOS))
    Cm["tau_splithalf"] = tau_sh
    rc.save(Cm, "rev12_checklist_confusion.csv", index=True)
    print(f"split-half threshold (95th pct under null) = {tau_sh:.2f}")
    print(Cm.round(2).to_string())
    print(R.groupby("scenario")[["r_cross", "r_splithalf", "r_substate", "r_random_split"]].median().round(2))


if __name__ == "__main__":
    main()
