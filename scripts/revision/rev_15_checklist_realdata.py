"""Diagnostic checklist applied to the observed data (Section 3.7, Table 2, Supplementary Table S15b)
and imputation denominators (Supplementary Table S15c-e).

Checklist (same rule and thresholds as the simulation in rev_12):
  1 cross-cohort r > 0 and two-sided p < 0.05                         -> replicated
  2 discovery-cohort split-half r below its label-permutation 95th pct -> inconclusive (power)
  3 matched sub-state r significant and above the random-split r       -> composition
  4 otherwise                                                          -> between-cohort difference (unresolved)
Split-half: 200 disease-stratified donor halves of GSE136831; effect = IPF-Control
difference with fine-cell-type fixed effects (as in rev_10). Null: donor labels permuted
(200 permutations x 20 splits each).

Imputation: for each of the 11 testable cell types, the number of module-clock input genes
(1,995) and composite input genes present after preprocessing (i.e. not median-imputed),
and the number present in every cell type.

RLE support: per cell type, the genes retained by the detection filter and the subset with
non-zero counts in every pseudobulk, which are the genes the RLE size factors are computed from.

Outputs: rev15_checklist_realdata.csv, rev15_detected_features.csv, rev15_detected_by_module.csv,
rev15_rle_support.csv
"""
import numpy as np
import pandas as pd

import rev_common as rc
from rev_10_replication import COARSE_A, coarse_effects

N_SPLIT = 200
B_NULL = 200
N_SPLIT_NULL = 20


def split_half(S, mods, labels, n_split, rng):
    don = labels.index.values
    ipf, ctl = don[labels.values == "IPF"], don[labels.values == "Control"]
    rs = []
    for _ in range(n_split):
        a, b = rng.permutation(ipf), rng.permutation(ctl)
        h1 = set(np.r_[a[: len(a) // 2], b[: len(b) // 2]])
        m1 = S.Subject_Identity.isin(h1).values
        e1 = coarse_effects(S[m1], mods, labels)
        e2 = coarse_effects(S[~m1], mods, labels)
        co = S.coarse.iloc[0]
        if co in e1.index and co in e2.index:
            rs.append(np.corrcoef(e1.loc[co].values, e2.loc[co].values)[0, 1])
    return float(np.median(rs)) if rs else np.nan


def checklist(rng):
    A = rc.adams_scores()
    mods = rc.module_names(A)
    A = A[A.Disease_Identity.isin(["IPF", "Control"])].copy()
    A["coarse"] = A.celltype.map(COARSE_A)
    cc = pd.read_csv(rc.OUT / "rev10_coarse_concordance.csv")
    ss = pd.read_csv(rc.OUT / "rev11_substate_concordance.csv")
    rows = []
    for r in cc.itertuples():
        S = A[A.coarse == r.coarse]
        lab = S.groupby("Subject_Identity").Disease_Identity.first()
        r_sh = split_half(S, mods, lab, N_SPLIT, rng)
        null = [split_half(S, mods, pd.Series(rng.permutation(lab.values), index=lab.index), N_SPLIT_NULL, rng)
                for _ in range(B_NULL)]
        tau = float(np.nanpercentile(null, 95))
        sub = ss[(ss.celltype == r.coarse) & (ss.split == "marker sub-state")]
        rnd = ss[(ss.celltype == r.coarse) & (ss.split == "random split")]
        r_sub = sub.r.mean() if len(sub) else np.nan
        p_sub = min(1.0, sub.p_two_sided.min() * 2) if len(sub) else np.nan
        r_rnd = rnd.r.mean() if len(rnd) else np.nan
        if r.p_two_sided < 0.05 and r.r > 0:
            att = "replicated"
        elif r_sh < tau:
            att = "inconclusive (power)"
        elif len(sub) and p_sub < 0.05 and r_sub > r_rnd:
            att = "composition"
        else:
            att = "between-cohort difference (unresolved)"
        rows.append(dict(coarse=r.coarse, n_IPF_A=r.n_IPF_A, n_Ctrl_A=r.n_Ctrl_A, n_IPF_H=r.n_IPF_H,
                         n_Ctrl_H=r.n_Ctrl_H, r_cross=r.r, p_cross=r.p_two_sided, fdr_cross=r.fdr_two_sided,
                         r_splithalf=r_sh, tau_splithalf=tau, p_splithalf=(1 + np.sum(np.array(null) >= r_sh)) / (B_NULL + 1),
                         r_substate=r_sub, p_substate=p_sub, r_random_split=r_rnd, attribution=att))
        print(rows[-1])
    return pd.DataFrame(rows)


def detected_features():
    pb, meta = rc.load_adams()
    clocks = rc.load_module_clocks()
    comp_genes = set(rc.load_composite()[1])
    mod_genes = {m: set(map(str, d["genes"])) for m, d in clocks.items()}
    union = set().union(*mod_genes.values())
    present, rows, by_mod = {}, [], []
    for ct in rc.testable_celltypes(meta):
        X, _ = rc.preprocess_ct(pb, meta, ct)
        X.columns = X.columns.map(str)
        pres = set(X.columns[X.notna().any(axis=0).values])
        present[ct] = pres
        rows.append(dict(celltype=ct, n_samples=len(X), module_genes_total=len(union),
                         module_genes_detected=len(union & pres), module_frac=len(union & pres) / len(union),
                         composite_genes_total=len(comp_genes), composite_genes_detected=len(comp_genes & pres),
                         composite_frac=len(comp_genes & pres) / max(1, len(comp_genes))))
        for m, g in mod_genes.items():
            by_mod.append(dict(celltype=ct, module=m, n_genes=len(g), detected=len(g & pres), frac=len(g & pres) / len(g)))
    common = set.intersection(*present.values())
    D = pd.DataFrame(rows)
    D.loc[len(D)] = dict(celltype="detected in all cell types", n_samples=np.nan, module_genes_total=len(union),
                         module_genes_detected=len(union & common), module_frac=len(union & common) / len(union),
                         composite_genes_total=len(comp_genes), composite_genes_detected=len(comp_genes & common),
                         composite_frac=len(comp_genes & common) / max(1, len(comp_genes)))
    return D, pd.DataFrame(by_mod)


def rle_support():
    """Genes entering the RLE size factors (non-zero in every pseudobulk) per testable cell type."""
    pb, meta = rc.load_adams()
    rows = []
    for ct in rc.testable_celltypes(meta):
        sub = meta[(meta.Manuscript_Identity == ct) & meta.Disease_Identity.isin(["IPF", "Control", "COPD"])]
        m = rc.tp.map_genes(rc.tp.filter_genes(pb[sub.group.tolist()]), "human", "Ensembl")
        rows.append(dict(celltype=ct, n_samples=m.shape[1], genes_retained=m.shape[0],
                         genes_nonzero_all=int((m.values > 0).all(axis=1).sum())))
    return pd.DataFrame(rows)


def main():
    rng = np.random.default_rng(rc.SEED)
    D, Dm = detected_features()
    rc.save(D, "rev15_detected_features.csv")
    rc.save(Dm, "rev15_detected_by_module.csv")
    rc.save(rle_support(), "rev15_rle_support.csv")
    print(D.round(3).to_string())
    print(Dm.groupby("module").frac.agg(["min", "median", "max"]).round(2).to_string())
    C = checklist(rng)
    rc.save(C, "rev15_checklist_realdata.csv")
    print(C.round(3).to_string())


if __name__ == "__main__":
    main()
