"""Retrained versus published rodent module clocks (Section 3.1, Supplementary Methods S2,
Supplementary Table S5).

The published rodent module clocks (Tyshkovskiy et al., Supp. Table 5B; 23 modules) are
distributed as coefficient lists only (no imputation medians, feature scaling or fitted
objects). Checks:
  A. reproduction: published coefficients applied to the released training matrix
     (Expression_data_relative_rodents_Scaled) under several feature conventions; the composite
     coefficients from Supp. Table 5A are applied the same way as a positive control;
  B. coefficient agreement with the retrained clocks;
  C. IPF effects from published-coefficient scores vs retrained scores (same pseudobulks).

Outputs: rev14_published_reproduction.csv, rev14_coef_agreement.csv, rev14_ipf_agreement.csv,
         rev14_summary.csv
"""
import numpy as np
import pandas as pd
from scipy import stats

import rev_common as rc
from rev_02_species_cv import TARGETS, load_training

OUTCOME_PUB = {"Mortality": "Mortality", "Chrono": "Chronological"}
PERF_OUTCOME = {"Mortality": "Expected mortality", "Chrono": "Chronological age"}
SUPP5 = rc.DATA / "supp" / "SuppTable5_coefficients.xlsx"


def conventions(X):
    Xf = np.where(np.isnan(X), 0.0, X)
    gz = (Xf - np.nanmean(X, 0)) / (np.nanstd(X, 0) + 1e-12)
    wz = (Xf - Xf.mean(1, keepdims=True)) / (Xf.std(1, keepdims=True) + 1e-12)
    return {"as_released": Xf, "gene_standardized": gz, "within_module_sample_z": wz}


def reproduction(ann, E):
    rows = []
    A = pd.read_excel(SUPP5, "(A) Composite clocks")
    A = A[A["Entrez ID"] != "Intercept"]
    A["Entrez ID"] = A["Entrez ID"].astype(str)
    for col, oc in [("Mortality, Rodents\nMulti-tissue, Scaling", "Mortality"),
                    ("Chronological Age, Rodents\nMulti-tissue, Scaling", "Chrono")]:
        cp = A.set_index("Entrez ID")[col].astype(float)
        cp = cp[cp != 0]
        g = [x for x in cp.index if x in E.columns]
        y = pd.to_numeric(ann[TARGETS[oc]], errors="coerce").values
        ok = np.isfinite(y)
        for conv, Z in conventions(E[g].values).items():
            rows.append(dict(table="5A composite", outcome=oc, module="composite", n_genes=len(g),
                             convention=conv, r_insample=stats.pearsonr((Z @ cp.loc[g].values)[ok], y[ok])[0]))
    pub = pd.read_csv(rc.DATA / "supp" / "module_coef_rodent.csv")
    pub["entrez"] = pub.entrez.astype(str)
    perf = pd.read_excel(rc.DATA / "supp" / "SuppTable4_performance.xlsx", "(B) Module rodent clocks")
    perf["Outcome"] = perf["Outcome"].ffill()
    for oc in ["Mortality", "Chrono"]:
        y = pd.to_numeric(ann[TARGETS[oc]], errors="coerce").values
        ok = np.isfinite(y)
        for m, p in pub[pub.outcome == OUTCOME_PUB[oc]].groupby("module"):
            cp = p.set_index("entrez").coef
            rep = perf[(perf.Module == m) & (perf.Outcome == PERF_OUTCOME[oc])]
            for conv, Z in conventions(E[list(cp.index)].values).items():
                rows.append(dict(table="5B module", outcome=oc, module=m, n_genes=len(cp), convention=conv,
                                 r_insample=stats.pearsonr((Z @ cp.values)[ok], y[ok])[0],
                                 r_reported=float(rep["Pearson r"].iloc[0]) if len(rep) else np.nan))
    return pd.DataFrame(rows)


def main():
    ann, expr = load_training()
    E = expr.T
    E.columns = E.columns.astype(str)
    Rp = reproduction(ann, E)
    rc.save(Rp, "rev14_published_reproduction.csv")
    print(Rp[Rp.convention == "as_released"].groupby(["table", "outcome"])[["r_insample", "r_reported"]]
          .median().round(3))

    pub = pd.read_csv(rc.DATA / "supp" / "module_coef_rodent.csv")
    pub["entrez"] = pub.entrez.astype(str)
    cv = pd.read_csv(rc.OUT / "rev02_species_cv.csv")
    rows = []
    for oc in ["Mortality", "Chrono"]:
        y = pd.to_numeric(ann[TARGETS[oc]], errors="coerce").values
        ok = np.isfinite(y)
        for lab, d in rc.load_module_clocks(oc).items():
            genes = [str(g) for g in d["genes"]]
            med, mu, sc, coef, b = rc.linear_parts(d["pipeline"])
            cp = pub[(pub.module == d["module"]) & (pub.outcome == OUTCOME_PUB[oc])].set_index("entrez").coef
            cp = cp.reindex(genes).fillna(0.0).values
            both = (cp != 0) & (coef != 0)
            Xm = E.reindex(columns=genes).values
            Z = (np.where(np.isnan(Xm), med[None, :], Xm) - mu[None, :]) / sc[None, :]
            rcv = cv[(cv.outcome == oc) & (cv.kind == "module") & (cv.module == lab)].r_oof_all
            rows.append(dict(outcome=oc, module=lab, color=d["module"], n_genes=len(genes),
                             n_nonzero_published=int((cp != 0).sum()), n_nonzero_retrained=int((coef != 0).sum()),
                             coef_spearman_rawscale=stats.spearmanr(cp, coef / sc)[0],
                             sign_agreement_both_nonzero=float(np.mean(np.sign(cp[both]) == np.sign(coef[both])))
                             if both.any() else np.nan,
                             r_insample_retrained=stats.pearsonr((Z @ coef + b)[ok], y[ok])[0],
                             r_cv_retrained=float(rcv.iloc[0]) if len(rcv) else np.nan))
    Cf = pd.DataFrame(rows)
    rc.save(Cf, "rev14_coef_agreement.csv")

    pb, meta = rc.load_adams()
    clocks = rc.load_module_clocks("Mortality")
    p_m = pub[pub.outcome == "Mortality"]
    recs = []
    for ct in rc.testable_celltypes(meta):
        X, md = rc.preprocess_ct(pb, meta, ct)
        md = md.loc[X.index]
        X.columns = X.columns.map(str)
        for lab, d in clocks.items():
            genes = [str(g) for g in d["genes"]]
            med, mu, sc, coef, b = rc.linear_parts(d["pipeline"])
            Xm = X.reindex(columns=genes).values
            Xm = np.where(np.isnan(Xm), med[None, :], Xm)
            cp = p_m[p_m.module == d["module"]].set_index("entrez").coef.reindex(genes).fillna(0.0).values
            df = md[["Disease_Identity", "age", "sexM"]].copy()
            df["pub"] = Xm @ cp
            df["new"] = ((Xm - mu[None, :]) / sc[None, :]) @ coef + b
            e_pub, e_new = rc.disease_effect(df, "pub"), rc.disease_effect(df, "new")
            sub = df[df.Disease_Identity.isin(["IPF", "Control"])]
            recs.append(dict(celltype=ct, module=lab, beta_retrained=e_new["beta"], p_retrained=e_new["p"],
                             t_retrained=e_new["beta"] / e_new["se"], p_published=e_pub["p"],
                             t_published=e_pub["beta"] / e_pub["se"],
                             score_corr=np.corrcoef(sub.pub, sub.new)[0, 1]))
    R = pd.DataFrame(recs)
    R["fdr_retrained"], R["fdr_published"] = rc.bh(R.p_retrained), rc.bh(R.p_published)
    rc.save(R, "rev14_ipf_agreement.csv")
    hits = R[R.fdr_retrained < 0.10]
    rep = Rp[(Rp.convention == "as_released")]
    S = pd.DataFrame([dict(
        composite_r_insample_mortality=rep[(rep.table == "5A composite") & (rep.outcome == "Mortality")].r_insample.iloc[0],
        allmodule_r_insample_mortality=rep[(rep.table == "5B module") & (rep.outcome == "Mortality")
                                           & (rep.module == "All module genes")].r_insample.iloc[0],
        allmodule_r_reported_mortality=rep[(rep.table == "5B module") & (rep.outcome == "Mortality")
                                           & (rep.module == "All module genes")].r_reported.iloc[0],
        module_r_reported_median_mortality=rep[(rep.table == "5B module") & (rep.outcome == "Mortality")
                                               & (rep.module != "All module genes")].r_reported.median(),
        module_r_insample_median_mortality=rep[(rep.table == "5B module") & (rep.outcome == "Mortality")
                                               & (rep.module != "All module genes")].r_insample.median(),
        module_r_insample_max_any_convention=Rp[(Rp.table == "5B module")].r_insample.max(),
        retrained_r_cv_median_mortality=Cf[Cf.outcome == "Mortality"].r_cv_retrained.median(),
        coef_spearman_median=Cf.coef_spearman_rawscale.median(),
        ipf_t_spearman=stats.spearmanr(R.t_retrained, R.t_published)[0],
        ipf_tests=len(R), ipf_hits_retrained=len(hits),
        ipf_hits_same_sign_published=int((np.sign(hits.t_published) == np.sign(hits.t_retrained)).sum()),
        ipf_hits_published_fdr10=int((R.fdr_published < .1).sum()),
        ipf_hits_overlap=int(((R.fdr_published < .1) & (R.fdr_retrained < .1)).sum()))])
    rc.save(S, "rev14_summary.csv")
    print(S.T.round(3).to_string())


if __name__ == "__main__":
    main()
