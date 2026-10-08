"""Human chronological-age transfer shown explicitly (Section 3.3, Supplementary Table S8).

(1) HLCA core (CELLxGENE; 42 healthy donors, 20-81 y): within each coarse cell type
    (>=12 donors), Spearman rho between module / composite score and donor age, raw and
    partial (sex, log cell number, log library size), BH-FDR per outcome.
(2) Adams control donors (GSE136831; 28 donors, 20-80 y): the same test within each of
    the 11 analysed cell types.
(3) GTEx v8 lung bulk (578 donors): per-module confounder-adjusted rho (existing output of
    scripts/gtex_validate.py). Superseded by rev_05b_gtex.py, which scores all 23 module
    clocks; the manuscript reports rev_05b.

Outputs: rev05_hlca_age.csv, rev05_adams_control_age.csv, rev05_human_age_summary.csv,
         figs/rev05_human_age.png
"""
import numpy as np
import pandas as pd
from scipy import stats

import rev_common as rc

ADAMS_TO_HLCA = {"B": "B", "Ciliated": "Ciliated", "Macrophage": "Macrophage",
                 "Macrophage_Alveolar": "Macrophage", "NK": "NK", "T": "T", "T_Cytotoxic": "T",
                 "cDC2": "DC", "cMonocyte": "Monocyte", "ncMonocyte": "Monocyte",
                 "Fibroblast": "Stromal"}


def partial_spearman(y, x, Z):
    """Spearman partial correlation of y and x given covariates Z (rank-residual method)."""
    ry, rx = stats.rankdata(y), stats.rankdata(x)
    Zr = np.column_stack([np.ones(len(y))] + [stats.rankdata(z) for z in Z.T])
    ey = ry - Zr @ np.linalg.lstsq(Zr, ry, rcond=None)[0]
    ex = rx - Zr @ np.linalg.lstsq(Zr, rx, rcond=None)[0]
    r = np.corrcoef(ey, ex)[0, 1]
    dof = len(y) - 2 - Z.shape[1]
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), dof)


def age_tests(S, age, Z, cols, extra):
    rows = []
    for c in cols:
        y = S[c].values
        rho, p = stats.spearmanr(y, age)
        prho, pp = partial_spearman(y, age, Z)
        rows.append(dict(**extra, module=c, n=len(y), rho=rho, p=p, rho_adj=prho, p_adj=pp))
    return rows


def score_both(X, clocks_m, clocks_c, comp_m, comp_c):
    Sm = rc.score_modules(X, clocks_m, comp_m)
    Sc = rc.score_modules(X, clocks_c, comp_c)
    return {"Mortality": Sm, "Chrono": Sc}


def main():
    import tage_prep as tp
    clocks = {"Mortality": rc.load_module_clocks("Mortality"), "Chrono": rc.load_module_clocks("Chrono")}
    comps = {"Mortality": rc.load_composite("Mortality"), "Chrono": rc.load_composite("Chrono")}

    # ---- (1) HLCA
    H = rc.DATA / "hlca"
    pb = pd.read_pickle(H / "pseudobulk_counts.pkl")
    meta = pd.read_csv(H / "pseudobulk_meta.csv")
    meta["key"] = meta["key"].astype(str)
    meta["lib_size"] = pb[meta.key].sum(0).values
    rows = []
    for ct, sub in meta.groupby("celltype"):
        if len(sub) < 12:
            continue
        md = sub.set_index("key")
        X = tp.preprocess(pb[sub.key.tolist()], md, species="human", gene_mapping_type="Ensembl")["scaled_diff"]
        md = md.loc[X.index]
        age = md.age.values.astype(float)
        Z = np.column_stack([(md.sex == "male").astype(float), np.log(md.n_cells), np.log(md.lib_size)])
        for oc in ["Mortality", "Chrono"]:
            S = rc.score_modules(X, clocks[oc], comps[oc])
            rows += age_tests(S, age, Z, list(S.columns), dict(dataset="HLCA", celltype=ct, outcome=oc))
    Hh = pd.DataFrame(rows)

    # ---- (2) Adams control donors
    pbA, metaA = rc.load_adams()
    rows = []
    for ct in rc.testable_celltypes(metaA):
        X, md = rc.preprocess_ct(pbA, metaA, ct)
        ctrl = md.loc[X.index].Disease_Identity.eq("Control").values
        Xc, mc = X[ctrl], md.loc[X.index][ctrl]
        if len(mc) < 6:
            continue
        age = mc.age.values.astype(float)
        Z = np.column_stack([mc.sexM.values, np.log(mc.n_cells), np.log(mc.lib_size)])
        for oc in ["Mortality", "Chrono"]:
            S = rc.score_modules(Xc, clocks[oc], comps[oc])
            rows += age_tests(S, age, Z, list(S.columns), dict(dataset="Adams controls", celltype=ct, outcome=oc))
    Ha = pd.DataFrame(rows)

    for T in (Hh, Ha):
        for oc in ["Mortality", "Chrono"]:
            m = (T.outcome == oc) & (T.module != "Composite")
            T.loc[m, "fdr"] = rc.bh(T.loc[m, "p"])
            T.loc[m, "fdr_adj"] = rc.bh(T.loc[m, "p_adj"])
            mc_ = (T.outcome == oc) & (T.module == "Composite")
            T.loc[mc_, "fdr"] = rc.bh(T.loc[mc_, "p"])
            T.loc[mc_, "fdr_adj"] = rc.bh(T.loc[mc_, "p_adj"])
    rc.save(Hh, "rev05_hlca_age.csv")
    rc.save(Ha, "rev05_adams_control_age.csv")

    # ---- summary
    rows = []
    for name, T in [("HLCA", Hh), ("Adams controls", Ha)]:
        for oc in ["Mortality", "Chrono"]:
            mods = T[(T.outcome == oc) & (T.module != "Composite")]
            comp = T[(T.outcome == oc) & (T.module == "Composite")]
            rows.append(dict(dataset=name, outcome=oc, n_celltypes=mods.celltype.nunique(),
                             n_tests=len(mods), donors_range=f"{T.n.min()}–{T.n.max()}",
                             module_median_rho=mods.rho.median(),
                             module_pos_fdr10=int(((mods.rho > 0) & (mods.fdr < .1)).sum()),
                             module_neg_fdr10=int(((mods.rho < 0) & (mods.fdr < .1)).sum()),
                             module_pos_fdr10_adj=int(((mods.rho_adj > 0) & (mods.fdr_adj < .1)).sum()),
                             module_neg_fdr10_adj=int(((mods.rho_adj < 0) & (mods.fdr_adj < .1)).sum()),
                             composite_median_rho=comp.rho.median(),
                             composite_pos_fdr10=int(((comp.rho > 0) & (comp.fdr < .1)).sum()),
                             composite_neg_fdr10=int(((comp.rho < 0) & (comp.fdr < .1)).sum())))
    # concordance of HLCA and Adams-control age effects (mortality modules, matched cell types)
    a = Ha[(Ha.outcome == "Mortality") & (Ha.module != "Composite")].copy()
    a["hlca_ct"] = a.celltype.map(ADAMS_TO_HLCA)
    h = Hh[(Hh.outcome == "Mortality") & (Hh.module != "Composite")]
    j = a.merge(h, left_on=["hlca_ct", "module"], right_on=["celltype", "module"], suffixes=("_adams", "_hlca"))
    r, p = stats.spearmanr(j.rho_adams, j.rho_hlca)
    Sm = pd.DataFrame(rows)
    Sm["hlca_vs_adams_rho_concordance"] = r
    Sm["hlca_vs_adams_concordance_p"] = p
    Sm["hlca_vs_adams_n_pairs"] = len(j)
    rc.save(Sm, "rev05_human_age_summary.csv")
    print(Sm.round(3).T.to_string())

    # ---- figure
    plt = rc.setup_mpl()
    g = pd.read_csv(rc.RES / "gtex_validation.csv")
    g = g[g.outcome == "Mortality"]
    fig = plt.figure(figsize=(7.2, 8.4), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.8])
    for k, (T, title) in enumerate([(Hh, "a  HLCA healthy lung, within cell type (42 donors, 20–81 y)"),
                                    (Ha, "b  GSE136831 control donors, within cell type (28 donors, 20–80 y)")]):
        ax = fig.add_subplot(gs[k])
        M = T[T.outcome == "Mortality"].pivot_table(index="celltype", columns="module", values="rho")
        F = T[T.outcome == "Mortality"].pivot_table(index="celltype", columns="module", values="fdr")
        cols = ["Composite"] + [c for c in clocks["Mortality"] if c in M.columns]
        M, F = M[cols], F.reindex(index=M.index, columns=cols)
        im = ax.imshow(M.values, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
        for i in range(M.shape[0]):
            for jj in range(M.shape[1]):
                if F.values[i, jj] < 0.10:
                    ax.text(jj, i, "*", ha="center", va="center", fontsize=8)
        ax.set_xticks(range(M.shape[1]))
        ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=6)
        ax.set_yticks(range(M.shape[0]))
        ax.set_yticklabels(M.index, fontsize=6)
        ax.axvline(0.5, color="k", lw=0.8)
        ax.set_title(title, loc="left")
        fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman ρ with age")
    ax = fig.add_subplot(gs[2])
    g = g.assign(label=g.module.replace({"COMPOSITE": "Composite"})).sort_values("adj_spearman")
    colors = ["#D55E00" if l == "Composite" else ("#0072B2" if f < 0.05 else "#999999")
              for l, f in zip(g.label, g.adj_fdr)]
    ax.barh(range(len(g)), g.adj_spearman, color=colors)
    ax.set_yticks(range(len(g)))
    ax.set_yticklabels(g.label, fontsize=5.5)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("Confounder-adjusted Spearman ρ with age (GTEx lung bulk, n = 578)")
    ax.set_title("c  GTEx v8 lung, bulk tissue (mortality clocks)", loc="left")
    fig.savefig(rc.FIG / "rev05_human_age.png")


if __name__ == "__main__":
    main()
