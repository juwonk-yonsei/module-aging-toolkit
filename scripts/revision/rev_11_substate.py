"""Myeloid sub-state analysis with a random-split control (Section 3.7, Supplementary Table S14).

Sub-states: cells of each myeloid cell type split at the cohort-specific median of a
monocyte-derived-minus-resident marker signature (monoHi / resHi; markers in
scripts/stageB_substate.py), aggregated to donor x cell type x sub-state pseudobulks.

  * module scores per sub-state; IPF effect within each sub-state and the
    sub-state x disease interaction (linear mixed model, donor random intercept, age, sex);
  * composition: per-donor monoHi fraction by disease, and the main macrophage effects
    re-estimated with the donor's monoHi fraction as a covariate;
  * cross-cohort concordance of matched sub-states versus random splits of the same
    cells (two-sided donor-permutation p).

Outputs: rev11_substate_effects.csv, rev11_interaction.csv, rev11_composition_adjusted.csv,
         rev11_substate_concordance.csv
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

import rev_common as rc

B = 2000


def score_stageB(tag, gmap, cov=None):
    import tage_prep as tp
    pb = pd.read_pickle(rc.RES / "mil" / f"stageB_pb_{tag}.pkl")
    gm = pd.read_csv(rc.RES / "mil" / f"stageB_gm_{tag}.csv")
    gm = gm[gm.Disease_Identity.isin(["IPF", "Control"]) & (gm.n_cells >= 40)]
    clocks, comp = rc.load_module_clocks(), rc.load_composite()
    recs = []
    for ct, sub in gm.groupby("celltype"):
        md = sub.set_index("group")
        X = tp.preprocess(pb[sub.group.tolist()], md, species="human", gene_mapping_type=gmap,
                          control_group_column="Disease_Identity", control_group_label="Control")["scaled_diff"]
        S = rc.score_modules(X, clocks, comp).join(md[["donor", "celltype", "substate", "Disease_Identity", "n_cells"]])
        recs.append(S)
    S = pd.concat(recs)
    if cov is not None:
        S = S.merge(cov[["age", "sexM"]], left_on="donor", right_index=True, how="left").set_index(S.index)
    return S


def eff_vectors(S, mods, labels=None):
    out = {}
    for (ct, ss), g in S.groupby(["celltype", "substate"]):
        d = g.donor.map(labels).values if labels is not None else g.Disease_Identity.values
        if (d == "IPF").sum() < 3 or (d == "Control").sum() < 3:
            continue
        out[(ct, ss)] = g[mods][d == "IPF"].mean().values - g[mods][d == "Control"].mean().values
    return out


def main():
    rng = np.random.default_rng(rc.SEED)
    cov = pd.read_csv(rc.DATA / "ipf" / "subject_covariates.csv").set_index("Subject_Identity")
    cov["sexM"] = (cov.Sex == "M").astype(float)
    SA = score_stageB("adams", "Ensembl", cov)
    mods = [c for c in SA.columns if c in rc.load_module_clocks()]

    # ---- within-sub-state effects and interaction (GSE136831)
    rows, irows = [], []
    for ct, g in SA.groupby("celltype"):
        g = g.dropna(subset=["age", "sexM"]).copy()
        g["ipf"] = (g.Disease_Identity == "IPF").astype(float)
        g["mono"] = (g.substate == "monoHi").astype(float)
        for m in mods + ["Composite"]:
            d = g.rename(columns={m: "y"})
            try:
                fit = smf.mixedlm("y ~ ipf * mono + age + sexM", d, groups=d["donor"]).fit(reml=True)
                pr, pv, cv = fit.params, fit.pvalues, fit.cov_params()
                b_res, p_res = pr["ipf"], pv["ipf"]
                b_mono = pr["ipf"] + pr["ipf:mono"]
                se_mono = np.sqrt(cv.loc["ipf", "ipf"] + cv.loc["ipf:mono", "ipf:mono"] + 2 * cv.loc["ipf", "ipf:mono"])
                p_mono = 2 * stats.norm.sf(abs(b_mono / se_mono))
                irows.append(dict(celltype=ct, module=m, beta_IPF_resHi=b_res, p_resHi=p_res,
                                  beta_IPF_monoHi=b_mono, p_monoHi=p_mono,
                                  interaction=pr["ipf:mono"], p_interaction=pv["ipf:mono"],
                                  substate_main=pr["mono"], p_substate=pv["mono"]))
            except Exception as e:
                irows.append(dict(celltype=ct, module=m, error=str(e)[:80]))
    I = pd.DataFrame(irows)
    mm = I.module != "Composite"
    for c in ["p_interaction", "p_resHi", "p_monoHi", "p_substate"]:
        I.loc[mm, c.replace("p_", "fdr_")] = rc.bh(I.loc[mm, c])
    rc.save(I, "rev11_interaction.csv")
    print(I[mm].groupby("celltype")[["fdr_interaction", "fdr_resHi", "fdr_monoHi", "fdr_substate"]]
          .apply(lambda d: (d < 0.1).sum()).to_string())

    # ---- composition: monoHi fraction per donor
    gm = pd.read_csv(rc.RES / "mil" / "stageB_gm_adams.csv")
    gm = gm[gm.Disease_Identity.isin(["IPF", "Control"])]
    frac = (gm.pivot_table(index=["donor", "celltype"], columns="substate", values="n_cells", aggfunc="sum")
            .fillna(0))
    frac["mono_frac"] = frac["monoHi"] / frac.sum(1)
    frac = frac.reset_index().merge(gm.groupby("donor").Disease_Identity.first(), left_on="donor", right_index=True)
    crow = []
    for ct, g in frac.groupby("celltype"):
        a, b = g[g.Disease_Identity == "IPF"].mono_frac, g[g.Disease_Identity == "Control"].mono_frac
        crow.append(dict(celltype=ct, mono_frac_IPF=a.median(), mono_frac_Control=b.median(),
                         p_mannwhitney=stats.mannwhitneyu(a, b).pvalue))
    # main macrophage hits adjusted for monoHi fraction
    S = rc.adams_scores()
    S = S[S.Disease_Identity.isin(["IPF", "Control"])]
    mf = frac[frac.celltype == "Macrophage"].set_index("donor").mono_frac
    st = pd.read_csv(rc.OUT / "rev08_celltype_stats.csv")
    hits = st[(st.module != "Composite") & (st.fdr_IPF < .1) & st.celltype.str.startswith("Macrophage")]
    arow = []
    for _, h in hits.iterrows():
        g = S[S.celltype == h.celltype].copy()
        g["mono_frac"] = g.Subject_Identity.map(mf)
        g = g.dropna(subset=["mono_frac"])
        e0 = rc.disease_effect(g, h.module)
        e1 = rc.disease_effect(g, h.module, covars=("age", "sexM", "mono_frac"))
        arow.append(dict(celltype=h.celltype, module=h.module, beta=e0["beta"], p=e0["p"],
                         beta_adj_monofrac=e1["beta"], p_adj_monofrac=e1["p"], ratio=e1["beta"] / e0["beta"]))
    Cadj = pd.DataFrame(arow)
    rc.save(pd.concat([pd.DataFrame(crow).assign(kind="monoHi fraction"), Cadj.assign(kind="adjusted effect")]),
            "rev11_composition_adjusted.csv")
    print(pd.DataFrame(crow).round(3).to_string())
    print(Cadj.round(3).to_string())

    # ---- cross-cohort concordance: matched sub-states vs random split (two-sided)
    rows = []
    SH = score_stageB("hab", "Gene.Symbol")
    for kind, ta, th in [("marker sub-state", "adams", "hab"), ("random split", "adams_rand1", "hab_rand1")]:
        A_ = SA if ta == "adams" else score_stageB(ta, "Ensembl")
        H_ = SH if th == "hab" else score_stageB(th, "Gene.Symbol")
        eA, eH = eff_vectors(A_, mods), eff_vectors(H_, mods)
        labA = A_.groupby("donor").Disease_Identity.first()
        labH = H_.groupby("donor").Disease_Identity.first()
        nulls = []
        for _ in range(B):
            nulls.append((eff_vectors(A_, mods, pd.Series(rng.permutation(labA.values), index=labA.index)),
                          eff_vectors(H_, mods, pd.Series(rng.permutation(labH.values), index=labH.index))))
        for key in sorted(set(eA) & set(eH)):
            r = np.corrcoef(eA[key], eH[key])[0, 1]
            nr = np.array([np.corrcoef(a[key], h[key])[0, 1] for a, h in nulls if key in a and key in h])
            rows.append(dict(split=kind, celltype=key[0], substate=key[1], r=r,
                             p_two_sided=(1 + np.sum(np.abs(nr) >= abs(r))) / (1 + len(nr))))
    Cc = pd.DataFrame(rows)
    rc.save(Cc, "rev11_substate_concordance.csv")
    print(Cc.round(3).to_string())


if __name__ == "__main__":
    main()
