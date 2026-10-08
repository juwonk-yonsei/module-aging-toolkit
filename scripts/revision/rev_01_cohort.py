"""Cohort characteristics and cell-type bookkeeping: 13 scored versus 11 tested cell types
(Section 2.4, Supplementary Tables S1-S2).

Outputs: rev01_cohort_donors.csv, rev01_celltype_table.csv, rev01_depth_by_disease.csv,
         rev01_hab_celltype_table.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import rev_common as rc


def summarize_age(x):
    x = x.dropna()
    return f"{x.mean():.1f} ± {x.std():.1f} ({x.min():.0f}–{x.max():.0f})"


def main():
    cov = pd.read_csv(rc.DATA / "ipf" / "subject_covariates.csv")
    rows = []
    for dis, g in cov.groupby("disease"):
        rows.append(dict(cohort="GSE136831 (Adams)", group=dis, donors=len(g),
                         age=summarize_age(g.age), male=int((g.Sex == "M").sum()),
                         female=int((g.Sex == "F").sum()), tissue=", ".join(sorted(g.tissue.unique()))))
    hab = pd.read_csv(rc.RES / "hab_pseudobulk_meta.csv")
    site = (pd.read_csv(rc.DATA / "ipf2" / "GSE135893_IPF_metadata.csv.gz", index_col=0)
            .groupby("Sample_Name").Sample_Source.first())
    for dis, g in hab.groupby("Disease_Identity"):
        donors = pd.Series(g.Subject_Identity.unique())
        by_site = (donors.map(site).replace({"DNA": "Donor Network of Arizona", "NTI": "Norton Thoracic Institute"})
                   .fillna("unknown").value_counts())
        rows.append(dict(cohort="GSE135893 (Habermann)", group=dis,
                         donors=len(donors), age="not available",
                         male=np.nan, female=np.nan,
                         tissue="lung tissue; " + ", ".join(f"{s} {n}" for s, n in by_site.items())))
    rc.save(pd.DataFrame(rows), "rev01_cohort_donors.csv")

    # all pseudobulks (no threshold) to reconcile cell-type counts
    meta_all = pd.read_csv(rc.RES / "pseudobulk_meta.csv")
    pb_all = pd.read_pickle(rc.RES / "pseudobulk_counts.pkl")
    meta_all["lib_size"] = pb_all[meta_all.group.astype(str)].sum(0).values
    m50 = meta_all[meta_all.n_cells >= 50]
    ct3, ct4 = set(), set()
    for ct, sub in m50.groupby("Manuscript_Identity"):
        n = sub.Disease_Identity.value_counts()
        if n.get("IPF", 0) >= 3 and n.get("Control", 0) >= 3:
            ct3.add(ct)
        if n.get("IPF", 0) >= 4 and n.get("Control", 0) >= 4:
            ct4.add(ct)
    rows = []
    for ct in sorted(ct3):
        sub = m50[m50.Manuscript_Identity == ct]
        rec = dict(celltype=ct, in_coverage_set=True, in_statistics_set=ct in ct4)
        for dis in ["Control", "IPF", "COPD"]:
            s = sub[sub.Disease_Identity == dis]
            rec[f"donors_{dis}"] = len(s)
            rec[f"cells_median_{dis}"] = float(s.n_cells.median()) if len(s) else np.nan
            rec[f"cells_total_{dis}"] = int(s.n_cells.sum())
        rows.append(rec)
    T = pd.DataFrame(rows).sort_values(["in_statistics_set", "celltype"], ascending=[False, True])
    rc.save(T, "rev01_celltype_table.csv")
    print(f"coverage set (>=3 donors/group): {len(ct3)}; statistics set (>=4): {len(ct4)}")
    print("only in coverage set:", sorted(ct3 - ct4))

    # depth / cell number by disease within the statistics cell types
    rows = []
    for ct in sorted(ct4):
        sub = m50[(m50.Manuscript_Identity == ct) & m50.Disease_Identity.isin(["IPF", "Control"])]
        a, b = sub[sub.Disease_Identity == "IPF"], sub[sub.Disease_Identity == "Control"]
        for var in ["n_cells", "lib_size"]:
            p = mannwhitneyu(a[var], b[var]).pvalue
            rows.append(dict(celltype=ct, variable=var, median_IPF=a[var].median(),
                             median_Control=b[var].median(), ratio=a[var].median() / b[var].median(),
                             p_mannwhitney=p))
    D = pd.DataFrame(rows)
    D["fdr"] = rc.bh(D.p_mannwhitney)
    rc.save(D, "rev01_depth_by_disease.csv")

    # Habermann
    h50 = hab[hab.n_cells >= 50]
    rows = []
    for ct, sub in h50.groupby("celltype"):
        n = sub.Disease_Identity.value_counts()
        rows.append(dict(celltype=ct, donors_Control=n.get("Control", 0), donors_IPF=n.get("IPF", 0),
                         cells_median=float(sub.n_cells.median())))
    rc.save(pd.DataFrame(rows), "rev01_hab_celltype_table.csv")


if __name__ == "__main__":
    main()
