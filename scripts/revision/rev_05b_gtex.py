"""GTEx v8 lung bulk (578 donors): composite and all 23 module clocks vs chronological age,
raw and adjusted for ischaemic time, RIN, Hardy death scale and sex.

Same design as scripts/gtex_validate.py (profiles centred on donors aged <= 35 y), but the
module clocks are keyed by module colour so that the two cholesterol modules
(darkred: cholesterol/mTOR; salmon: cholesterol/platelet) are both scored.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

import rev_common as rc

GTEX = rc.DATA / "gtex"


def load():
    counts = pd.read_pickle(GTEX / "lung_counts.pkl")
    meta = pd.read_csv(GTEX / "lung_meta.csv", index_col=0).dropna(subset=["age_mid"])
    counts = counts[meta.index]
    meta["ref"] = np.where(meta["age_mid"] <= 35, "young", "old")
    pp = rc.tp.preprocess(counts, meta, species="human", gene_mapping_type="Ensembl",
                          control_group_column="ref", control_group_label="young")
    X = pp["scaled_diff"]
    X.columns = X.columns.map(str)
    meta = meta.loc[X.index]
    tech = pd.read_csv(GTEX / "lung_tech.tsv", sep="\t", header=None,
                       names=["SAMPID", "isch", "rin"]).set_index("SAMPID")
    ph = pd.read_csv(GTEX / "SubjectPhenotypes.txt", sep="\t").set_index("SUBJID")
    cov = pd.DataFrame(index=X.index)
    cov["isch"] = tech.reindex(X.index)["isch"].astype(float)
    cov["rin"] = tech.reindex(X.index)["rin"].astype(float)
    subj = pd.Index(X.index).str.split("-").str[:2].str.join("-")
    cov["hardy"] = pd.Series(ph.reindex(subj)["DTHHRDY"].values, index=X.index)
    cov["sexM"] = meta["sexM"].values
    C = pd.get_dummies(cov, columns=["hardy"], dummy_na=True).astype(float)
    C = C.fillna(C.median(numeric_only=True))
    C.insert(0, "intercept", 1.0)
    return X, meta["age_mid"].values, C.values


def main():
    X, age, Cmat = load()
    print(f"GTEx lung: {X.shape[0]} samples")

    def adj(pred):
        beta, *_ = np.linalg.lstsq(Cmat, pred, rcond=None)
        return spearmanr(pred - Cmat @ beta, age)

    rows = []
    for outcome in ["Mortality", "Chrono"]:
        preds = {name: d["pipeline"].predict(X.reindex(columns=[str(g) for g in d["genes"]]).values)
                 for name, d in rc.load_module_clocks(outcome).items()}
        m, genes = rc.load_composite(outcome)
        preds["Composite"] = m.predict(X.reindex(columns=genes).values)
        for name, pred in preds.items():
            r, p = spearmanr(pred, age)
            ra, pa = adj(pred)
            rows.append(dict(outcome=outcome, module=name, n=len(age), spearman=r, p=p, adj_spearman=ra, adj_p=pa))
    T = pd.DataFrame(rows)
    for oc in T.outcome.unique():
        k = (T.outcome == oc) & (T.module != "Composite")
        T.loc[k, "fdr"] = multipletests(T.loc[k, "p"], method="fdr_bh")[1]
        T.loc[k, "adj_fdr"] = multipletests(T.loc[k, "adj_p"], method="fdr_bh")[1]
    rc.save(T, "rev05b_gtex_validation.csv")
    for oc in T.outcome.unique():
        s = T[(T.outcome == oc) & (T.module != "Composite")]
        c = T[(T.outcome == oc) & (T.module == "Composite")].iloc[0]
        print(f"[{oc}] composite adj rho {c.adj_spearman:.3f} (p {c.adj_p:.1e}); modules: median adj rho "
              f"{s.adj_spearman.median():.3f}, adj FDR<0.05 positive {((s.adj_spearman > 0) & (s.adj_fdr < .05)).sum()}, "
              f"negative {((s.adj_spearman < 0) & (s.adj_fdr < .05)).sum()} of {len(s)}")
    print(T[T.outcome == "Mortality"].sort_values("adj_spearman", ascending=False).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
