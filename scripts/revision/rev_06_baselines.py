"""Baselines and negative controls for the module scores (Sections 2.7 and 3.6, Supplementary Table S10).

A. IPF (GSE136831), every cell type x module (253 tests), covariate-adjusted IPF-vs-Control
   t statistic for:
     clock            trained elastic-net module score (as in the paper)
     geneset          equal-weight mean of the standardized module genes (unsigned)
     geneset_signed   equal-weight mean with the sign of the trained coefficient
   and two competitive nulls (N_NULL draws each):
     random_genes     size-matched random genes carrying the module's coefficients
     perm_coef        the module's coefficients shuffled across its own genes
   p_competitive = P(|t_null| >= |t_clock|).

B. Rodent training data: out-of-fold accuracy of the signed gene-set score (signs learned
   within training folds) versus the trained module clock on identical dataset-grouped folds.

Outputs: rev06_ipf_baselines.csv, rev06_ipf_baselines_summary.csv, rev06_rodent_geneset_cv.csv
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import GroupKFold

import rev_common as rc

N_NULL = 1000


def t_disease(Y, md):
    """Y: n x k matrix. OLS score ~ IPF + age + sexM; returns t (k,) for IPF and df."""
    X = np.column_stack([np.ones(len(md)), (md.Disease_Identity == "IPF").astype(float),
                         md.age.values, md.sexM.values])
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ X.T @ Y
    R = Y - X @ B
    dof = len(md) - X.shape[1]
    s2 = (R ** 2).sum(0) / dof
    se = np.sqrt(XtX_inv[1, 1] * s2)
    return B[1] / np.where(se > 0, se, np.nan), B[1], dof


def part_a():
    pb, meta = rc.load_adams()
    clocks = rc.load_module_clocks()
    gstats = rc.train_gene_stats()
    rng = np.random.default_rng(rc.SEED)
    rows = []
    for ct in rc.testable_celltypes(meta):
        X, md = rc.preprocess_ct(pb, meta, ct, diseases=("IPF", "Control"))
        md = md.loc[X.index]
        ns = rc.NullScorer(X, gstats)
        for lab, d in clocks.items():
            y_clock = d["pipeline"].predict(X.reindex(columns=[str(g) for g in d["genes"]]).values)
            ys = np.column_stack([y_clock, ns.geneset(d), ns.geneset(d, signed=True)])
            t_obs, b_obs, dof = t_disease(ys, md)
            t_rand, *_ = t_disease(ns.random_genes(d, N_NULL, rng), md)
            t_perm, *_ = t_disease(ns.perm_coef(d, N_NULL, rng), md)
            p_param = 2 * stats.t.sf(np.abs(t_obs), dof)
            rows.append(dict(celltype=ct, module=lab, n_genes=len(d["genes"]), df=dof,
                             t_clock=t_obs[0], beta_clock=b_obs[0], p_clock=p_param[0],
                             t_geneset=t_obs[1], p_geneset=p_param[1],
                             t_geneset_signed=t_obs[2], p_geneset_signed=p_param[2],
                             p_competitive_random=(1 + np.nansum(np.abs(t_rand) >= abs(t_obs[0]))) / (N_NULL + 1),
                             p_competitive_perm=(1 + np.nansum(np.abs(t_perm) >= abs(t_obs[0]))) / (N_NULL + 1),
                             null_random_abs_t_median=np.nanmedian(np.abs(t_rand)),
                             null_perm_abs_t_median=np.nanmedian(np.abs(t_perm))))
    R = pd.DataFrame(rows)
    for c in ["p_clock", "p_geneset", "p_geneset_signed", "p_competitive_random", "p_competitive_perm"]:
        R[c.replace("p_", "fdr_")] = rc.bh(R[c])
    old = pd.read_csv(rc.RES / "ipf_module_stats_adjusted.csv")
    R = R.merge(old[["celltype", "module", "fdr_IPF", "perm_p_IPF"]], on=["celltype", "module"], how="left")
    rc.save(R, "rev06_ipf_baselines.csv")

    hits = R.fdr_IPF < 0.10
    S = pd.DataFrame([dict(
        n_tests=len(R),
        hits_clock_fdr10_published=int(hits.sum()),
        hits_geneset_fdr10=int((R.fdr_geneset < .1).sum()),
        hits_geneset_signed_fdr10=int((R.fdr_geneset_signed < .1).sum()),
        corr_t_clock_geneset=stats.spearmanr(R.t_clock, R.t_geneset)[0],
        corr_t_clock_geneset_signed=stats.spearmanr(R.t_clock, R.t_geneset_signed)[0],
        hits_beating_random_p05=int((hits & (R.p_competitive_random < .05)).sum()),
        hits_beating_perm_p05=int((hits & (R.p_competitive_perm < .05)).sum()),
        hits_competitive_random_fdr10=int((R.fdr_competitive_random < .1).sum()),
        hits_competitive_perm_fdr10=int((R.fdr_competitive_perm < .1).sum()),
        median_abs_t_hits=R.loc[hits, "t_clock"].abs().median(),
        median_null_abs_t_random_hits=R.loc[hits, "null_random_abs_t_median"].median(),
    )])
    rc.save(S, "rev06_ipf_baselines_summary.csv")
    print(S.T.round(3).to_string())
    return R


def part_b():
    from rev_02_species_cv import TARGETS, load_training, safe_r
    ann, expr = load_training()
    E = expr.T
    E.columns = E.columns.astype(str)
    clocks = rc.load_module_clocks()
    cv = pd.read_csv(rc.OUT / "rev02_species_cv.csv")
    groups = ann.Source.values
    rows = []
    for oc, col in TARGETS.items():
        y_all = pd.to_numeric(ann[col], errors="coerce").values
        ok = np.isfinite(y_all)
        for lab, d in clocks.items():
            X = E[[str(g) for g in d["genes"]]].values[ok]
            y, g = y_all[ok], groups[ok]
            oof = np.full(len(y), np.nan)
            for tr, te in GroupKFold(n_splits=5).split(X, y, g):
                med = np.nanmedian(X[tr], 0)
                Xtr = np.where(np.isnan(X[tr]), med, X[tr])
                mu, sd = Xtr.mean(0), Xtr.std(0)
                sd[sd == 0] = 1
                Ztr = (Xtr - mu) / sd
                w = np.sign([np.corrcoef(Ztr[:, j], y[tr])[0, 1] if Ztr[:, j].std() > 0 else 0
                             for j in range(Ztr.shape[1])])
                w = np.nan_to_num(w)
                Zte = (np.where(np.isnan(X[te]), med, X[te]) - mu) / sd
                oof[te] = Zte @ w / max(1, (w != 0).sum())
            r_cl = cv[(cv.outcome == oc) & (cv.module == lab)].r_oof_all.values
            rows.append(dict(outcome=oc, module=lab, n_genes=X.shape[1], r_oof_geneset_signed=safe_r(oof, y),
                             r_oof_clock=float(r_cl[0]) if len(r_cl) else np.nan))
    R = pd.DataFrame(rows)
    R["clock_minus_geneset"] = R.r_oof_clock - R.r_oof_geneset_signed
    rc.save(R, "rev06_rodent_geneset_cv.csv")
    print(R.groupby("outcome")[["r_oof_clock", "r_oof_geneset_signed", "clock_minus_geneset"]].median().round(3))
    print("clock > signed gene set:", R.groupby("outcome").apply(lambda d: int((d.clock_minus_geneset > 0).sum())).to_dict())


def part_c(n_null=10000):
    """Higher-resolution competitive nulls for the published FDR<0.10 hits; BH across hits."""
    R = pd.read_csv(rc.OUT / "rev06_ipf_baselines.csv")
    hits = R[R.fdr_IPF < 0.10]
    pb, meta = rc.load_adams()
    clocks = rc.load_module_clocks()
    gstats = rc.train_gene_stats()
    rng = np.random.default_rng(rc.SEED + 1)
    rows = []
    for ct, sub in hits.groupby("celltype"):
        X, md = rc.preprocess_ct(pb, meta, ct, diseases=("IPF", "Control"))
        md = md.loc[X.index]
        ns = rc.NullScorer(X, gstats)
        for lab in sub.module:
            d = clocks[lab]
            y = d["pipeline"].predict(X.reindex(columns=[str(g) for g in d["genes"]]).values)
            t_obs = t_disease(y[:, None], md)[0][0]
            t_rand = t_disease(ns.random_genes(d, n_null, rng), md)[0]
            t_perm = t_disease(ns.perm_coef(d, n_null, rng), md)[0]
            rows.append(dict(celltype=ct, module=lab, t_clock=t_obs,
                             p_comp_random=(1 + np.nansum(np.abs(t_rand) >= abs(t_obs))) / (n_null + 1),
                             p_comp_perm=(1 + np.nansum(np.abs(t_perm) >= abs(t_obs))) / (n_null + 1)))
    H = pd.DataFrame(rows)
    H["fdr_comp_random_within_hits"] = rc.bh(H.p_comp_random)
    H["fdr_comp_perm_within_hits"] = rc.bh(H.p_comp_perm)
    H["module_specific"] = (H.fdr_comp_random_within_hits < 0.10) & (H.fdr_comp_perm_within_hits < 0.10)
    rc.save(H.sort_values("p_comp_random"), "rev06_hits_competitive.csv")
    print(H.sort_values("p_comp_random").round(4).to_string())


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
