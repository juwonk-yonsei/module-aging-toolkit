"""Cross-cohort replication, re-analysed (Section 3.7, Supplementary Table S13).

A. Habermann (GSE135893) module scores per fine cell type (>=50 cells, >=4 donors/group),
   harmonized to a shared coarse vocabulary; site recorded (only Vanderbilt contributes
   both diagnoses: 7 Control / 7 IPF; DNA = controls only, NTI = IPF only).
B. Replication of the 30 specific GSE136831 effects (FDR<0.10) in the matched coarse
   cell type: sign agreement and two-sided donor-permutation p; count of sign-concordant
   effects against its permutation null. Coarse-level 23-module concordance: Pearson r
   (invariant to a common shift), sign concordance of centred vectors, partial r adjusted
   for each module's response to a uniform expression shift, two-sided permutation p.
C. Whole-tissue comparison (GSE134692 bulk; whole-donor pseudobulks of both single-cell
   cohorts) with two-sided p, direction of the globally significant modules, and a
   decomposition of the whole-donor effect into a composition-only and an expression-only
   component.
D. Cell-type composition per donor in both cohorts.
E. Overlap of macrophage identity markers with module genes.

Outputs: rev10_*.csv
"""
import numpy as np
import pandas as pd
from scipy import stats

import rev_common as rc

B = 2000
COARSE_H = {"Macrophages": "Macrophage", "Proliferating Macrophages": "Macrophage",
            "Monocytes": "Monocyte", "NK Cells": "NK", "T Cells": "T cell",
            "Proliferating T Cells": "T cell", "B Cells": "B cell", "Plasma Cells": "B cell",
            "cDCs": "DC", "pDCs": "DC", "Ciliated": "Ciliated", "Differentiating Ciliated": "Ciliated",
            "AT2": "AT2", "Transitional AT2": "AT2", "AT1": "AT1", "Fibroblasts": "Fibroblast",
            "HAS1 High Fibroblasts": "Fibroblast", "PLIN2+ Fibroblasts": "Fibroblast",
            "Myofibroblasts": "Fibroblast", "Mast Cells": "Mast", "Basal": "Basal"}
COARSE_A = {"Macrophage": "Macrophage", "Macrophage_Alveolar": "Macrophage", "cMonocyte": "Monocyte",
            "ncMonocyte": "Monocyte", "NK": "NK", "T": "T cell", "T_Cytotoxic": "T cell",
            "T_Regulatory": "T cell", "B": "B cell", "B_Plasma": "B cell", "cDC1": "DC", "cDC2": "DC",
            "DC_Langerhans": "DC", "DC_Mature": "DC", "pDC": "DC", "Ciliated": "Ciliated",
            "ATII": "AT2", "ATI": "AT1", "Fibroblast": "Fibroblast", "Myofibroblast": "Fibroblast",
            "Mast": "Mast", "Basal": "Basal"}
MARKERS = {"monocyte-derived": ["SPP1", "FCN1", "VCAN", "CTHRC1", "INHBA", "MERTK", "SPARC", "CD14", "S100A8"],
           "resident alveolar": ["FABP4", "MARCO", "MRC1", "PPARG", "INHBA", "CD36", "APOC1", "C1QA"]}


# ----------------------------------------------------------------------------- helpers
def hab_site():
    m = pd.read_csv(rc.DATA / "ipf2" / "GSE135893_IPF_metadata.csv.gz", index_col=0,
                    usecols=["Unnamed: 0", "Sample_Name", "Sample_Source"] if False else None)
    return m.groupby("Sample_Name").Sample_Source.first()


def hab_scores():
    pb = pd.read_pickle(rc.RES / "hab_pseudobulk_counts.pkl")
    gm = pd.read_csv(rc.RES / "hab_pseudobulk_meta.csv")
    gm["group"] = gm["group"].astype(str)
    gm = gm[gm.n_cells >= 50]
    import tage_prep as tp
    clocks, comp = rc.load_module_clocks(), rc.load_composite()
    site = hab_site()
    recs = []
    for ct, sub in gm.groupby("celltype"):
        n = sub.Disease_Identity.value_counts()
        if n.get("IPF", 0) < 4 or n.get("Control", 0) < 4:
            continue
        md = sub.set_index("group")
        X = tp.preprocess(pb[sub.group.tolist()], md, species="human", gene_mapping_type="Gene.Symbol",
                          control_group_column="Disease_Identity", control_group_label="Control")["scaled_diff"]
        S = rc.score_modules(X, clocks, comp).join(md[["Subject_Identity", "Disease_Identity", "n_cells"]])
        S["celltype"] = ct
        recs.append(S)
    S = pd.concat(recs)
    S["coarse"] = S.celltype.map(COARSE_H)
    S["site"] = S.Subject_Identity.map(site)
    return S


def coarse_effects(S, mods, labels=None):
    """IPF - Control effect per module at coarse level with fine-cell-type fixed effects.
    labels: optional Series donor -> disease (for permutation)."""
    dis = S.Subject_Identity.map(labels) if labels is not None else S.Disease_Identity
    out = {}
    for co, g in S.groupby("coarse"):
        d = dis.loc[g.index].values if labels is not None else g.Disease_Identity.values
        ipf = (d == "IPF").astype(float)
        if ipf.sum() < 3 or (1 - ipf).sum() < 3:
            continue
        fine = pd.get_dummies(g.celltype).astype(float).values
        X = np.column_stack([ipf, fine])
        Bt = np.linalg.lstsq(X, g[mods].values.astype(float), rcond=None)[0]
        out[co] = Bt[0]
    return pd.DataFrame(out, index=mods).T


def donor_labels(S):
    return S.groupby("Subject_Identity").Disease_Identity.first()


def perm_labels(lab, rng):
    return pd.Series(rng.permutation(lab.values), index=lab.index)


def uniform_shift_response(clocks, mods):
    """Score change per unit uniform increase of every standardized input gene."""
    u = {}
    for m in mods:
        med, mu, sc, coef, b = rc.linear_parts(clocks[m]["pipeline"])
        u[m] = float(np.sum(coef / sc))
    return pd.Series(u)


def partial_r(a, b, z):
    Z = np.column_stack([np.ones(len(z)), z])
    ra = a - Z @ np.linalg.lstsq(Z, a, rcond=None)[0]
    rb = b - Z @ np.linalg.lstsq(Z, b, rcond=None)[0]
    return np.corrcoef(ra, rb)[0, 1]


# ----------------------------------------------------------------------------- main
def main():
    rng = np.random.default_rng(rc.SEED)
    clocks = rc.load_module_clocks()
    A = rc.adams_scores()
    mods = rc.module_names(A)
    A = A[A.Disease_Identity.isin(["IPF", "Control"])].copy()
    A["coarse"] = A.celltype.map(COARSE_A)
    H = hab_scores()
    H = H[H.Disease_Identity.isin(["IPF", "Control"])]
    rc.save(H.reset_index().rename(columns={"index": "group"}), "rev10_hab_scores.csv")
    u = uniform_shift_response(clocks, mods)
    labA, labH = donor_labels(A), donor_labels(H)
    Hv = H[H.site == "Vanderbilt"]
    labHv = donor_labels(Hv)

    # ---- B1: coarse-level concordance
    eA, eH, eHv = coarse_effects(A, mods), coarse_effects(H, mods), coarse_effects(Hv, mods)
    common = [c for c in eA.index if c in eH.index]
    nullA = [coarse_effects(A, mods, perm_labels(labA, rng)) for _ in range(B)]
    nullH = [coarse_effects(H, mods, perm_labels(labH, rng)) for _ in range(B)]
    rows = []
    for co in common:
        a, h = eA.loc[co].values, eH.loc[co].values
        r = np.corrcoef(a, h)[0, 1]
        nr = np.array([np.corrcoef(nA.loc[co].values, nH.loc[co].values)[0, 1]
                       for nA, nH in zip(nullA, nullH) if co in nA.index and co in nH.index])
        ac, hc = a - np.median(a), h - np.median(h)
        rec = dict(coarse=co, n_IPF_A=int(A[A.coarse == co].query("Disease_Identity=='IPF'").Subject_Identity.nunique()),
                   n_Ctrl_A=int(A[A.coarse == co].query("Disease_Identity=='Control'").Subject_Identity.nunique()),
                   n_IPF_H=int(H[H.coarse == co].query("Disease_Identity=='IPF'").Subject_Identity.nunique()),
                   n_Ctrl_H=int(H[H.coarse == co].query("Disease_Identity=='Control'").Subject_Identity.nunique()),
                   r=r, p_two_sided=(1 + np.sum(np.abs(nr) >= abs(r))) / (1 + len(nr)),
                   p_one_sided=(1 + np.sum(nr >= r)) / (1 + len(nr)),
                   sign_conc_raw=float(np.mean(np.sign(a) == np.sign(h))),
                   sign_conc_centred=float(np.mean(np.sign(ac) == np.sign(hc))),
                   r_partial_uniform_shift=partial_r(a, h, u.values),
                   r_A_vs_shift=np.corrcoef(a, u.values)[0, 1], r_H_vs_shift=np.corrcoef(h, u.values)[0, 1])
        if co in eHv.index:
            rec["r_vanderbilt_only"] = np.corrcoef(a, eHv.loc[co].values)[0, 1]
        rows.append(rec)
    C = pd.DataFrame(rows)
    C["fdr_two_sided"] = rc.bh(C.p_two_sided)
    rc.save(C, "rev10_coarse_concordance.csv")
    print(C.round(3).to_string())

    # ---- B2: the 30 specific effects
    st = pd.read_csv(rc.OUT / "rev08_celltype_stats.csv")
    hits = st[(st.module != "Composite") & (st.fdr_IPF < 0.10)].copy()
    hits["coarse"] = hits.celltype.map(COARSE_A)
    rows = []
    for _, h in hits.iterrows():
        if h.coarse not in eH.index:
            rows.append(dict(celltype=h.celltype, module=h.module, coarse=h.coarse, beta_A=h.beta_IPF))
            continue
        bH = eH.loc[h.coarse, h.module]
        nb = np.array([n.loc[h.coarse, h.module] for n in nullH if h.coarse in n.index])
        bHv = eHv.loc[h.coarse, h.module] if h.coarse in eHv.index else np.nan
        rows.append(dict(celltype=h.celltype, module=h.module, coarse=h.coarse, beta_A=h.beta_IPF,
                         fdr_A=h.fdr_IPF, beta_H=bH, same_sign=bool(np.sign(bH) == np.sign(h.beta_IPF)),
                         p_H_two_sided=(1 + np.sum(np.abs(nb) >= abs(bH))) / (1 + len(nb)),
                         beta_H_vanderbilt=bHv,
                         same_sign_vanderbilt=bool(np.sign(bHv) == np.sign(h.beta_IPF)) if np.isfinite(bHv) else np.nan))
    R = pd.DataFrame(rows)
    ok = R.beta_H.notna()
    R.loc[ok, "fdr_H"] = rc.bh(R.loc[ok, "p_H_two_sided"])
    n_same = int(R.loc[ok, "same_sign"].sum())
    T = R[ok]
    null_same = np.array([sum(np.sign(n.loc[c, m]) == np.sign(b) for c, m, b in zip(T.coarse, T.module, T.beta_A))
                          for n in nullH if set(T.coarse) <= set(n.index)])
    R.attrs["summary"] = dict(testable=int(ok.sum()), same_sign=n_same,
                              p_count=float((1 + np.sum(null_same >= n_same)) / (len(null_same) + 1)),
                              null_same_median=float(np.median(null_same)), n_perm_valid=len(null_same))
    rc.save(R, "rev10_specific_effects.csv")
    print(f"specific effects testable in Habermann: {ok.sum()}, same sign: {n_same}, "
          f"null median {np.median(null_same):.0f}, permutation p(count) = {R.attrs['summary']['p_count']:.3f}; "
          f"nominal p<0.05 in Habermann: {(R.p_H_two_sided < .05).sum()}")
    pd.DataFrame([R.attrs["summary"] | dict(n_p05_H=int((R.p_H_two_sided < .05).sum()),
                                            same_sign_vanderbilt=int(R.same_sign_vanderbilt.fillna(False).sum()))]
                 ).to_csv(rc.OUT / "rev10_specific_effects_summary.csv", index=False)

    # ---- C: whole-tissue comparison
    bulk_rows, dir_rows, decomp_rows = whole_tissue(A, mods, clocks, rng)
    rc.save(pd.DataFrame(bulk_rows), "rev10_bulk_concordance.csv")
    rc.save(pd.DataFrame(dir_rows), "rev10_bulk_direction.csv")
    rc.save(pd.DataFrame(decomp_rows), "rev10_bulk_decomposition.csv")

    # ---- D: composition
    rc.save(composition(), "rev10_composition.csv")

    # ---- E: marker overlap
    mem = pd.read_csv(rc.DATA / "supp" / "module_membership_rodent.csv")
    lab = rc.make_labels([f"{m}|{a}" for m, a in mem[["module", "annotation"]].drop_duplicates().values])
    mem["label"] = [lab[f"{m}|{a}"] for m, a in mem[["module", "annotation"]].values]
    sym = mem.assign(SYM=mem.symbol.str.upper()).set_index("SYM")
    rows = []
    for kind, genes in MARKERS.items():
        for g in genes:
            rows.append(dict(marker_set=kind, gene=g, module=sym.label.get(g, None) if g in sym.index else None))
    rc.save(pd.DataFrame(rows), "rev10_marker_overlap.csv")


def build_bulk_from_celltypes(pb, meta, ct_col, mode, ref_dis="Control"):
    """Whole-donor pseudobulk from donor x cell-type pseudobulks.
    mode: observed | composition_only | expression_only."""
    meta = meta[meta.Disease_Identity.isin(["IPF", "Control"])]
    cnt = pb[meta.group]
    cts = meta[ct_col].unique()
    ref = meta[meta.Disease_Identity == ref_dis]
    prof = {c: cnt[ref[ref[ct_col] == c].group].sum(1) / max(1, ref[ref[ct_col] == c].n_cells.sum()) for c in cts}
    ref_frac = ref.groupby(ct_col).n_cells.sum() / ref.n_cells.sum()
    out = {}
    for d, g in meta.groupby("Subject_Identity"):
        N = g.n_cells.sum()
        if mode == "observed":
            v = cnt[g.group].sum(1)
        elif mode == "composition_only":
            v = sum(prof[c] * n for c, n in zip(g[ct_col], g.n_cells))
        else:
            own = {c: cnt[gg].values[:, 0] / n for c, gg, n in zip(g[ct_col], g.group.map(lambda x: [x]), g.n_cells)}
            v = sum((own.get(c, prof[c].values) if c in own else prof[c].values) * ref_frac.get(c, 0) * N
                    for c in cts)
            v = pd.Series(v, index=cnt.index)
        out[d] = v
    dis = meta.groupby("Subject_Identity").Disease_Identity.first()
    return pd.DataFrame(out), dis


def score_bulk(counts, dis, gmap, clocks, comp):
    import tage_prep as tp
    md = pd.DataFrame({"Disease_Identity": dis.reindex(counts.columns).values}, index=counts.columns.astype(str))
    counts = counts.copy()
    counts.columns = counts.columns.astype(str)
    X = tp.preprocess(counts, md, species="human", gene_mapping_type=gmap,
                      control_group_column="Disease_Identity", control_group_label="Control")["scaled_diff"]
    S = rc.score_modules(X, clocks, comp)
    return S, md.loc[S.index, "Disease_Identity"].values


def whole_tissue(A, mods, clocks, rng):
    comp = rc.load_composite()
    bundles = {}
    sd = pd.read_pickle(rc.RES / "mil" / "stageD_pb_adams.pkl")
    bundles["GSE136831 whole-donor"] = score_bulk(sd["counts"], sd["dis"], "Ensembl", clocks, comp)
    sh = pd.read_pickle(rc.RES / "mil" / "stageD_pb_hab.pkl")
    bundles["GSE135893 whole-donor"] = score_bulk(sh["counts"], sh["dis"], "Gene.Symbol", clocks, comp)
    cnt = pd.read_csv(rc.DATA / "ipf_bulk" / "GSE134692_raw_counts.txt.gz", sep="\t", index_col=0)
    des = pd.read_csv(rc.DATA / "ipf_bulk" / "GSE134692_design.txt.gz", sep="\t", index_col="sample_id").reindex(cnt.columns)
    keep = des.DiseaseStatus.isin(["IPF", "Normal"]).values
    cnt = cnt.loc[:, keep]
    cnt.index = [str(g).split(".")[0] for g in cnt.index]
    dis = des.DiseaseStatus[keep].map({"IPF": "IPF", "Normal": "Control"})
    bundles["GSE134692 bulk"] = score_bulk(cnt, dis, "Ensembl", clocks, comp)

    def eff(S, lab):
        return S[mods][lab == "IPF"].mean().values - S[mods][lab == "Control"].mean().values

    keys = list(bundles)
    rows = []
    for i in range(3):
        for j in range(i + 1, 3):
            (Sa, la), (Sb, lb) = bundles[keys[i]], bundles[keys[j]]
            a, b = eff(Sa, la), eff(Sb, lb)
            r = np.corrcoef(a, b)[0, 1]
            nr = np.array([np.corrcoef(eff(Sa, rng.permutation(la)), eff(Sb, rng.permutation(lb)))[0, 1]
                           for _ in range(B)])
            ac, bc = a - np.median(a), b - np.median(b)
            rows.append(dict(pair=f"{keys[i]} vs {keys[j]}", r=r,
                             p_two_sided=(1 + np.sum(np.abs(nr) >= abs(r))) / (B + 1),
                             p_one_sided=(1 + np.sum(nr >= r)) / (B + 1),
                             sign_conc_centred=float(np.mean(np.sign(ac) == np.sign(bc)))))
    G = pd.read_csv(rc.OUT / "rev08_global_stats.csv")
    sig = G[(G.module != "Composite") & (G.fdr_cr1_normal < 0.10)]
    dir_rows = []
    for k, (S, lab) in bundles.items():
        e = pd.Series(eff(S, lab), index=mods)
        for _, g in sig.iterrows():
            y = S[g.module].values
            p = stats.mannwhitneyu(y[lab == "IPF"], y[lab == "Control"]).pvalue
            dir_rows.append(dict(cohort=k, module=g.module, beta_global_A=g.beta, delta=e[g.module],
                                 same_direction=bool(np.sign(e[g.module]) == np.sign(g.beta)), p_mannwhitney=p,
                                 n_IPF=int((lab == "IPF").sum()), n_Control=int((lab == "Control").sum())))
    # decomposition of whole-donor effects (GSE136831; GSE135893)
    decomp = []
    pbA = pd.read_pickle(rc.RES / "pseudobulk_counts.pkl")
    mA = pd.read_csv(rc.RES / "pseudobulk_meta.csv")
    mA["group"] = mA.group.astype(str)
    pbH = pd.read_pickle(rc.RES / "hab_pseudobulk_counts.pkl")
    mH = pd.read_csv(rc.RES / "hab_pseudobulk_meta.csv")
    mH["group"] = mH.group.astype(str)
    for name, pb, mm, ctc, gmap in [("GSE136831", pbA, mA, "Manuscript_Identity", "Ensembl"),
                                    ("GSE135893", pbH, mH, "celltype", "Gene.Symbol")]:
        effs = {}
        for mode in ["observed", "composition_only", "expression_only"]:
            cnt, dis = build_bulk_from_celltypes(pb, mm, ctc, mode)
            S, lab = score_bulk(cnt, dis, gmap, clocks, comp)
            effs[mode] = eff(S, lab)
        for mode in ["composition_only", "expression_only"]:
            decomp.append(dict(cohort=name, component=mode,
                               r_with_observed=np.corrcoef(effs["observed"], effs[mode])[0, 1],
                               slope_on_observed=np.polyfit(effs["observed"], effs[mode], 1)[0]))
    return rows, dir_rows, decomp


def composition():
    mA = pd.read_csv(rc.RES / "pseudobulk_meta.csv")
    mH = pd.read_csv(rc.RES / "hab_pseudobulk_meta.csv")
    rows = []
    for name, m, ctc, cmap in [("GSE136831", mA, "Manuscript_Identity", COARSE_A),
                               ("GSE135893", mH, "celltype", COARSE_H)]:
        m = m[m.Disease_Identity.isin(["IPF", "Control"])].copy()
        m["coarse"] = m[ctc].map(cmap).fillna("Other")
        tot = m.groupby("Subject_Identity").n_cells.sum()
        frac = m.groupby(["Subject_Identity", "coarse"]).n_cells.sum().unstack(fill_value=0).div(tot, axis=0)
        if name == "GSE136831":
            am = m[m[ctc] == "Macrophage_Alveolar"].groupby("Subject_Identity").n_cells.sum()
            mac = m[m.coarse == "Macrophage"].groupby("Subject_Identity").n_cells.sum()
            frac["Alveolar share of macrophages"] = (am / mac).reindex(frac.index)
        dis = m.groupby("Subject_Identity").Disease_Identity.first()
        for c in frac.columns:
            a, b = frac.loc[dis == "IPF", c].dropna(), frac.loc[dis == "Control", c].dropna()
            rows.append(dict(cohort=name, cell_type=c, median_IPF=a.median(), median_Control=b.median(),
                             p_mannwhitney=stats.mannwhitneyu(a, b).pvalue if len(a) and len(b) else np.nan,
                             n_IPF=len(a), n_Control=len(b)))
    D = pd.DataFrame(rows)
    D["fdr"] = rc.bh(D.p_mannwhitney)
    return D


if __name__ == "__main__":
    main()
