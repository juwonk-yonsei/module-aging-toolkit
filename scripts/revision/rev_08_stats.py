"""Disease contrasts and inference (Sections 2.6, 3.4 and 3.5, Supplementary Table S9).

Per cell type x module (253 tests; plus composite reported separately):
  * IPF vs Control, score ~ IPF + age + sex: classical OLS, HC3 small-sample robust SE,
    Freedman-Lane permutation p (B = 5000); 95% CIs; BH-FDR per family.
  * Direct IPF vs COPD contrast from score ~ IPF + COPD + age + sex (Control reference),
    classical and HC3; BH-FDR. "IPF-specific" requires FDR < 0.10 for IPF vs Control AND
    for the direct IPF - COPD contrast (replaces "significant in IPF, not in COPD").
Global per module (all 11 cell types, IPF vs Control, cell-type fixed effects):
  * donor-clustered robust OLS (as in the original submission), linear mixed model with a
    donor random intercept, and wild cluster restricted bootstrap (Rademacher, B = 9999).

Outputs: rev08_celltype_stats.csv, rev08_global_stats.csv, rev08_summary.csv
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

import rev_common as rc

B_PERM = 5000
B_WILD = 9999


def ols_full(Y, X):
    """Y n x k. Returns beta (p x k), classical cov diag, HC3 cov diag, dof, (X'X)^-1, resid."""
    XtX_inv = np.linalg.pinv(X.T @ X)
    Bt = XtX_inv @ X.T @ Y
    E = Y - X @ Bt
    n, p = X.shape
    dof = n - p
    s2 = (E ** 2).sum(0) / dof
    h = np.einsum("ij,jk,ik->i", X, XtX_inv, X)
    w = (E / (1 - np.clip(h, None, 0.999))[:, None]) ** 2
    A = XtX_inv @ X.T                                   # p x n
    hc3 = np.einsum("pi,ik,pi->pk", A, w, A)            # p x k  (diag only)
    cls = np.outer(np.diag(XtX_inv), s2)
    return Bt, cls, hc3, dof, XtX_inv, E


def contrast(Bt, XtX_inv, E, X, c, dof):
    """Linear contrast c'beta with classical and HC3 SE for every column."""
    est = c @ Bt
    s2 = (E ** 2).sum(0) / dof
    se_cls = np.sqrt((c @ XtX_inv @ c) * s2)
    h = np.einsum("ij,jk,ik->i", X, XtX_inv, X)
    a = c @ XtX_inv @ X.T                               # n
    w = (E / (1 - np.clip(h, None, 0.999))[:, None]) ** 2
    se_hc3 = np.sqrt((a[:, None] ** 2 * w).sum(0))
    return est, se_cls, se_hc3


def p_t(t, dof):
    return 2 * stats.t.sf(np.abs(t), dof)


def freedman_lane(Y, X, col, rng, B):
    keep = [j for j in range(X.shape[1]) if j != col]
    Xr = X[:, keep]
    H = Xr @ np.linalg.pinv(Xr.T @ Xr) @ Xr.T
    fit, res = H @ Y, Y - H @ Y
    XtX_inv = np.linalg.pinv(X.T @ X)
    A = (XtX_inv @ X.T)[col]

    def tstat(Yy):
        b = XtX_inv @ X.T @ Yy
        e = Yy - X @ b
        s2 = (e ** 2).sum(0) / (X.shape[0] - X.shape[1])
        return b[col] / np.sqrt(XtX_inv[col, col] * s2)

    t0 = tstat(Y)
    cnt = np.zeros(Y.shape[1])
    for _ in range(B):
        cnt += np.abs(tstat(fit + res[rng.permutation(len(Y))])) >= np.abs(t0) - 1e-12
    return (cnt + 1) / (B + 1)


def wild_cluster(y, X, col, clusters, rng, B):
    """Wild cluster restricted bootstrap-t (Rademacher) for H0: beta[col] = 0."""
    uc, cid = np.unique(clusters, return_inverse=True)
    G = len(uc)

    def cr1_t(yy):
        XtX_inv = np.linalg.pinv(X.T @ X)
        b = XtX_inv @ X.T @ yy
        e = yy - X @ b
        S = np.zeros((G, X.shape[1]))
        np.add.at(S, cid, X * e[:, None])
        meat = S.T @ S
        n, k = X.shape
        V = XtX_inv @ meat @ XtX_inv * (G / (G - 1)) * ((n - 1) / (n - k))
        return b[col] / np.sqrt(V[col, col]), b[col], np.sqrt(V[col, col])

    t0, b0, se0 = cr1_t(y)
    keep = [j for j in range(X.shape[1]) if j != col]
    Xr = X[:, keep]
    br = np.linalg.lstsq(Xr, y, rcond=None)[0]
    fit, res = Xr @ br, y - Xr @ br
    ts = np.empty(B)
    for i in range(B):
        w = rng.choice([-1.0, 1.0], size=G)[cid]
        ts[i] = cr1_t(fit + res * w)[0]
    return b0, se0, t0, float(np.mean(np.abs(ts) >= abs(t0)))


def main():
    rng = np.random.default_rng(rc.SEED)
    S = rc.adams_scores()
    mods = rc.module_names(S) + ["Composite"]
    rows = []
    for ct, g in S.groupby("celltype"):
        base = g[g.Disease_Identity.isin(["IPF", "Control"])]
        X = np.column_stack([np.ones(len(base)), (base.Disease_Identity == "IPF").astype(float),
                             base.age, base.sexM])
        Y = base[mods].values.astype(float)
        Bt, cls, hc3, dof, _, _ = ols_full(Y, X)
        q = stats.t.ppf(0.975, dof)
        p_perm = freedman_lane(Y, X, 1, rng, B_PERM)
        n_copd = int((g.Disease_Identity == "COPD").sum())
        if n_copd >= 4:
            three = g[g.Disease_Identity.isin(["IPF", "Control", "COPD"])]
            X3 = np.column_stack([np.ones(len(three)), (three.Disease_Identity == "IPF").astype(float),
                                  (three.Disease_Identity == "COPD").astype(float), three.age, three.sexM])
            Y3 = three[mods].values.astype(float)
            B3, cls3, hc33, dof3, XtX3, E3 = ols_full(Y3, X3)
            est, se_c, se_h = contrast(B3, XtX3, E3, X3, np.array([0, 1, -1, 0, 0.]), dof3)
            copd_b, copd_p = B3[2], p_t(B3[2] / np.sqrt(cls3[2]), dof3)
        else:
            est = se_c = se_h = copd_b = copd_p = np.full(len(mods), np.nan)
            dof3 = np.nan
        for k, m in enumerate(mods):
            se, seh = np.sqrt(cls[1, k]), np.sqrt(hc3[1, k])
            rows.append(dict(celltype=ct, module=m, n_IPF=int((base.Disease_Identity == "IPF").sum()),
                             n_Control=int((base.Disease_Identity == "Control").sum()), n_COPD=n_copd,
                             beta_IPF=Bt[1, k], se=se, ci_lo=Bt[1, k] - q * se, ci_hi=Bt[1, k] + q * se,
                             p_IPF=p_t(Bt[1, k] / se, dof), se_hc3=seh, ci_lo_hc3=Bt[1, k] - q * seh,
                             ci_hi_hc3=Bt[1, k] + q * seh, p_IPF_hc3=p_t(Bt[1, k] / seh, dof),
                             p_IPF_perm=p_perm[k], beta_COPD=copd_b[k], p_COPD=copd_p[k],
                             IPF_minus_COPD=est[k], se_IPF_minus_COPD=se_c[k],
                             p_IPF_vs_COPD=p_t(est[k] / se_c[k], dof3) if n_copd >= 4 else np.nan,
                             p_IPF_vs_COPD_hc3=p_t(est[k] / se_h[k], dof3) if n_copd >= 4 else np.nan))
    R = pd.DataFrame(rows)
    mm = R.module != "Composite"
    for c in ["p_IPF", "p_IPF_hc3", "p_IPF_perm", "p_COPD", "p_IPF_vs_COPD", "p_IPF_vs_COPD_hc3"]:
        R.loc[mm, c.replace("p_", "fdr_")] = rc.bh(R.loc[mm, c])
    R["IPF_specific_direct"] = (R.fdr_IPF < .1) & (R.fdr_IPF_vs_COPD < .1)
    rc.save(R, "rev08_celltype_stats.csv")

    # ---- global
    D = S[S.Disease_Identity.isin(["IPF", "Control"])].copy()
    D["ipf"] = (D.Disease_Identity == "IPF").astype(float)
    cts = sorted(D.celltype.unique())
    Xg = np.column_stack([np.ones(len(D)), D.ipf, D.age, D.sexM] +
                         [(D.celltype == c).astype(float) for c in cts[1:]])
    grows = []
    for m in mods:
        y = D[m].values.astype(float)
        b0, se0, t0, p_w = wild_cluster(y, Xg, 1, D.Subject_Identity.values, rng, B_WILD)
        dd = D.rename(columns={m: "y"})
        try:
            mx = smf.mixedlm("y ~ ipf + age + sexM + C(celltype)", dd, groups=dd["Subject_Identity"]).fit(reml=True)
            b_mx, p_mx, lo_mx, hi_mx = mx.params["ipf"], mx.pvalues["ipf"], *mx.conf_int().loc["ipf"].values
        except Exception:
            b_mx = p_mx = lo_mx = hi_mx = np.nan
        G = D.Subject_Identity.nunique()
        q = stats.t.ppf(0.975, G - 1)
        grows.append(dict(module=m, beta=b0, se_cr1=se0, ci_lo=b0 - q * se0, ci_hi=b0 + q * se0,
                          p_cr1_t=p_t(t0, G - 1), p_cr1_normal=2 * stats.norm.sf(abs(t0)),
                          p_wild_cluster=p_w, beta_mixed=b_mx, ci_lo_mixed=lo_mx, ci_hi_mixed=hi_mx,
                          p_mixed=p_mx))
    Gt = pd.DataFrame(grows)
    for c in ["p_cr1_normal", "p_cr1_t", "p_wild_cluster", "p_mixed"]:
        Gt[c.replace("p_", "fdr_")] = rc.bh(Gt[c])
    rc.save(Gt, "rev08_global_stats.csv")

    hits = R[mm & (R.fdr_IPF < .1)]
    Sm = pd.DataFrame([dict(
        n_tests=int(mm.sum()), hits_classical=int((R[mm].fdr_IPF < .1).sum()),
        hits_hc3=int((R[mm].fdr_IPF_hc3 < .1).sum()), hits_perm=int((R[mm].fdr_IPF_perm < .1).sum()),
        hits_all_three=int(((R[mm].fdr_IPF < .1) & (R[mm].fdr_IPF_hc3 < .1) & (R[mm].fdr_IPF_perm < .1)).sum()),
        hits_with_copd_testable=int(hits.n_COPD.ge(4).sum()),
        hits_ipf_specific_direct=int(hits.IPF_specific_direct.sum()),
        hits_old_rule_copd_ns=int((hits.n_COPD.ge(4) & (hits.p_COPD > .05)).sum()),
        global_fdr10_cr1=int((Gt[Gt.module != "Composite"].fdr_cr1_normal < .1).sum()),
        global_fdr10_wild=int((Gt[Gt.module != "Composite"].fdr_wild_cluster < .1).sum()),
        global_fdr10_mixed=int((Gt[Gt.module != "Composite"].fdr_mixed < .1).sum()))])
    rc.save(Sm, "rev08_summary.csv")
    print(Sm.T.to_string())
    print(Gt.round(4).to_string())
    print(hits[["celltype", "module", "beta_IPF", "fdr_IPF", "fdr_IPF_hc3", "fdr_IPF_perm", "n_COPD",
                "IPF_minus_COPD", "fdr_IPF_vs_COPD", "IPF_specific_direct"]].round(3).to_string())


if __name__ == "__main__":
    main()
