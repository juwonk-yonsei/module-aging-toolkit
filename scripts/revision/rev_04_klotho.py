"""Klotho-KO validation with full statistics and negative controls (Section 3.2, Supplementary Table S7).

Per tissue (kidney, skeletal muscle; 6 KO vs 6 WT, male):
  * per-module KO-WT difference, Welch 95% CI, Hedges' g, exact two-sided permutation p
    over all 924 genotype relabelings, BH-FDR across modules x tissues;
  * a priori expectation: progeroid KO -> positive (pro-mortality) shift; count of modules
    shifting positive, tested against the exact relabeling null (accounts for module
    correlation);
  * negative controls: (i) WT-vs-WT 3:3 splits (no biological contrast);
    (ii) coefficient-permuted clocks (same genes, shuffled weights);
    (iii) size-matched random genes carrying the module coefficients.

Outputs: rev04_klotho_modules.csv, rev04_klotho_summary.csv, rev04_klotho_wtwt.csv
"""
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

import rev_common as rc

N_NULL = 1000


def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    J = 1 - 3 / (4 * (na + nb) - 9)
    return J * (a.mean() - b.mean()) / sp


def welch_ci(a, b):
    d = a.mean() - b.mean()
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = np.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    q = stats.t.ppf(0.975, df)
    return d, d - q * se, d + q * se, 2 * stats.t.sf(abs(d / se), df)


def main():
    import tage_prep as tp
    ext = rc.TAGE_DIR / "inst" / "extdata"
    counts = pd.read_csv(ext / "Exprs_example.csv", index_col=0)
    meta = pd.read_csv(ext / "Metadata_example.csv", index_col=0)
    clocks = rc.load_module_clocks()
    comp = rc.load_composite()
    gstats = rc.train_gene_stats()
    rng = np.random.default_rng(rc.SEED)

    rows, summ, wtwt = [], [], []
    for tissue, mt in meta.groupby("Tissue"):
        pp = tp.preprocess(counts[mt.index], mt, species="mouse", gene_mapping_type="Ensembl",
                           control_group_column="Genotype", control_group_label="WT")
        X = pp["scaled_diff"]
        S = rc.score_modules(X, clocks, comp)
        ko = (mt.loc[S.index, "Genotype"] == "Klotho KO").values
        n = len(ko)
        relabels = [np.isin(np.arange(n), c) for c in combinations(range(n), int(ko.sum()))]
        D = np.array([[S[c].values[m].mean() - S[c].values[~m].mean() for c in S.columns] for m in relabels])
        obs = np.array([S[c].values[ko].mean() - S[c].values[~ko].mean() for c in S.columns])
        ns = rc.NullScorer(X, gstats)
        for j, c in enumerate(S.columns):
            a, b = S[c].values[ko], S[c].values[~ko]
            d, lo, hi, p_w = welch_ci(a, b)
            rec = dict(tissue=tissue, module=c, n_KO=int(ko.sum()), n_WT=int((~ko).sum()), delta=d,
                       ci_lo=lo, ci_hi=hi, hedges_g=hedges_g(a, b), p_welch=p_w,
                       p_perm_exact=float(np.mean(np.abs(D[:, j]) >= abs(obs[j]) - 1e-12)))
            if c != "Composite":
                for kind in ["perm_coef", "random_genes"]:
                    Z = getattr(ns, kind)(clocks[c], N_NULL, rng)
                    dn = Z[ko].mean(0) - Z[~ko].mean(0)
                    rec[f"p_vs_{kind}"] = (1 + np.sum(np.abs(dn) >= abs(d))) / (N_NULL + 1)
                    rec[f"null_{kind}_sd"] = dn.std()
            rows.append(rec)
        mods = [i for i, c in enumerate(S.columns) if c != "Composite"]
        cnt_obs = int((obs[mods] > 0).sum())
        cnt_null = (D[:, mods] > 0).sum(1)
        summ.append(dict(tissue=tissue, n_modules=len(mods), n_positive=cnt_obs,
                         p_count_exact=float(np.mean(cnt_null >= cnt_obs)),
                         null_count_median=float(np.median(cnt_null)),
                         composite_delta=obs[list(S.columns).index("Composite")]))
        wt = np.where(~ko)[0]
        for c3 in combinations(wt, 3):
            if 0 not in c3:
                continue
            m = np.isin(wt, c3)
            dd = S.iloc[wt][S.columns].values
            wtwt.append(dict(tissue=tissue, split=",".join(map(str, c3)),
                             **{c: dd[m, j].mean() - dd[~m, j].mean() for j, c in enumerate(S.columns)}))

    R = pd.DataFrame(rows)
    mod_mask = R.module != "Composite"
    R.loc[mod_mask, "fdr_perm"] = rc.bh(R.loc[mod_mask, "p_perm_exact"])
    R.loc[mod_mask, "fdr_welch"] = rc.bh(R.loc[mod_mask, "p_welch"])
    both = R[mod_mask].pivot(index="module", columns="tissue", values="delta")
    Sm = pd.DataFrame(summ)
    Sm["n_positive_both_tissues"] = int((both > 0).all(1).sum())
    W = pd.DataFrame(wtwt)
    for t in R.tissue.unique():
        cols = [c for c in W.columns if c not in ("tissue", "split", "Composite")]
        wmax = W[W.tissue == t][cols].abs().values
        Sm.loc[Sm.tissue == t, "wtwt_median_abs_delta"] = float(np.median(wmax))
        Sm.loc[Sm.tissue == t, "obs_median_abs_delta"] = float(R[(R.tissue == t) & mod_mask].delta.abs().median())
    rc.save(R, "rev04_klotho_modules.csv")
    rc.save(Sm, "rev04_klotho_summary.csv")
    rc.save(W, "rev04_klotho_wtwt.csv")
    print(Sm.round(3).to_string())
    print(R[R.fdr_perm < 0.1][["tissue", "module", "delta", "ci_lo", "ci_hi", "hedges_g", "p_perm_exact",
                                "fdr_perm", "p_vs_perm_coef", "p_vs_random_genes"]].round(3).to_string())


if __name__ == "__main__":
    main()
