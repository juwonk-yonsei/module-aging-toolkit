"""Step 3 (bulk compartment): does the WHOLE-TISSUE IPF module signature replicate
across three INDEPENDENT cohorts when cell composition is NOT resolved?

  - GSE134692 : true bulk RNA-seq, 46 IPF / 26 Normal lung homogenates (+ ALI, dropped).
  - Adams      : whole-donor pseudobulk (ALL cells summed per donor), IPF/Control.
  - Habermann  : whole-donor pseudobulk, IPF/Control.

Each -> 23-module IPF-vs-Control effect vector; pairwise Pearson concordance with a
donor/sample disease-permutation null. Contrast with the cell-type-resolved myeloid
result (Stage 2 / stageC): the manuscript claims bulk replicates where single-cell
myeloid does not.

Outputs: results/mil/stageD_bulk.csv, results/figures/stageD_bulk.png
"""
import os
import sys, gzip, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import pearsonr
REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
from mil_score import score_matrix, MODULES

IPFD = REPO / "data" / "ipf"
D2 = REPO / "data" / "ipf2"
BULK = REPO / "data" / "ipf_bulk"
MIL = REPO / "results" / "mil"
FIG = REPO / "results" / "figures"


def stream_donor_pb(mtx, group_of_cell, ngroups, nrows):
    pb = np.zeros(nrows * ngroups, dtype=np.float64)
    f = gzip.open(mtx, "rt")
    for line in f:
        if line.startswith("%"):
            continue
        break
    for ch in pd.read_csv(f, sep=r"\s+", header=None, names=["g", "c", "v"],
                          dtype={"g": np.int32, "c": np.int32, "v": np.float64}, chunksize=40_000_000):
        gi = ch["g"].values - 1; ci = ch["c"].values - 1; v = ch["v"].values
        grp = group_of_cell[ci]; keep = grp >= 0
        pb += np.bincount(gi[keep].astype(np.int64) * ngroups + grp[keep],
                          weights=v[keep], minlength=nrows * ngroups)
    return pb.reshape(nrows, ngroups)


def adams_bulk():
    cache = MIL / "stageD_pb_adams.pkl"
    if cache.exists():
        return pd.read_pickle(cache)
    genes = pd.read_csv(IPFD / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    genes.columns = [c.strip('"') for c in genes.columns]
    gid = genes["Ensembl_GeneID"].str.strip('"').values
    barc = pd.read_csv(IPFD / "GSE136831_AllCells.cellBarcodes.txt.gz", header=None)[0].values
    meta = pd.read_csv(IPFD / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
    meta.columns = [c.strip('"') for c in meta.columns]
    for c in meta.select_dtypes("object"):
        meta[c] = meta[c].str.strip('"')
    meta = meta.set_index("CellBarcode_Identity").reindex(barc)
    ok = meta.Disease_Identity.isin(["IPF", "Control"]).values
    donor = meta.Subject_Identity.astype(str).values
    uniq = pd.Index(pd.unique(donor[ok]))
    gm = {d: i for i, d in enumerate(uniq)}
    goc = np.array([gm[donor[i]] if ok[i] else -1 for i in range(len(barc))], dtype=np.int64)
    pb = stream_donor_pb(IPFD / "GSE136831_RawCounts_Sparse.mtx.gz", goc, len(uniq), len(gid))
    dis = pd.Series(meta.Disease_Identity.values, index=donor).groupby(level=0).first().reindex(uniq)
    df = pd.DataFrame(pb, index=gid, columns=uniq)
    out = {"counts": df, "dis": dis, "gmap": "Ensembl"}
    pd.to_pickle(out, cache); return out


def hab_bulk():
    cache = MIL / "stageD_pb_hab.pkl"
    if cache.exists():
        return pd.read_pickle(cache)
    gid = pd.read_csv(D2 / "GSE135893_genes.tsv.gz", header=None)[0].values
    barc = pd.read_csv(D2 / "GSE135893_barcodes.tsv.gz", header=None)[0].values
    meta = pd.read_csv(D2 / "GSE135893_IPF_metadata.csv.gz", index_col=0).reindex(barc)
    dcol = "Diagnosis"; scol = "Sample_Name"
    ok = meta[dcol].isin(["IPF", "Control"]).values
    donor = meta[scol].astype(str).values
    uniq = pd.Index(pd.unique(donor[ok]))
    gm = {d: i for i, d in enumerate(uniq)}
    goc = np.array([gm[donor[i]] if ok[i] else -1 for i in range(len(barc))], dtype=np.int64)
    pb = stream_donor_pb(D2 / "GSE135893_matrix.mtx.gz", goc, len(uniq), len(gid))
    dis = pd.Series(meta[dcol].values, index=donor).groupby(level=0).first().reindex(uniq)
    df = pd.DataFrame(pb, index=gid, columns=uniq)
    out = {"counts": df, "dis": dis, "gmap": "Gene.Symbol"}
    pd.to_pickle(out, cache); return out


def gse134692():
    cnt = pd.read_csv(BULK / "GSE134692_raw_counts.txt.gz", sep="\t", index_col=0)
    des = pd.read_csv(BULK / "GSE134692_design.txt.gz", sep="\t", index_col="sample_id")
    des = des.reindex(cnt.columns)
    keep = des.DiseaseStatus.isin(["IPF", "Normal"])
    cnt = cnt.loc[:, keep.values]
    dis = des.DiseaseStatus[keep.values].map({"IPF": "IPF", "Normal": "Control"})
    cnt.index = [str(g).split(".")[0] for g in cnt.index]
    return {"counts": cnt, "dis": dis, "gmap": "Ensembl"}


def effect_vec(bundle):
    dis = bundle["dis"]
    md = pd.DataFrame({"Disease_Identity": dis.values}, index=[str(x) for x in dis.index])
    cnt = bundle["counts"].copy(); cnt.columns = [str(x) for x in cnt.columns]
    pp = tp.preprocess(cnt, md, species="human", gene_mapping_type=bundle["gmap"],
                       control_group_column="Disease_Identity", control_group_label="Control")
    S = score_matrix(pp["scaled_diff"])
    lab = md.loc[S.index, "Disease_Identity"].values
    return S[MODULES].copy(), lab


def eff(S, lab):
    return S[lab == "IPF"].mean().values - S[lab == "Control"].mean().values


def conc(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return pearsonr(a[ok], b[ok])[0] if ok.sum() >= 3 else np.nan


def perm_p(Sa, la, Sb, lb, r_obs, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    null = np.empty(n)
    for i in range(n):
        null[i] = conc(eff(Sa, rng.permutation(la)), eff(Sb, rng.permutation(lb)))
    null = null[np.isfinite(null)]
    return (1 + np.sum(null >= r_obs)) / (1 + len(null)), null.std()


def main():
    print("scoring GSE134692 bulk ..."); B = effect_vec(gse134692())
    print("building/scoring Adams whole-donor pseudobulk ..."); A = effect_vec(adams_bulk())
    print("building/scoring Habermann whole-donor pseudobulk ..."); H = effect_vec(hab_bulk())
    named = {"GSE134692(bulk)": B, "Adams(sc-bulk)": A, "Habermann(sc-bulk)": H}
    keys = list(named)
    rows = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            (Sa, la), (Sb, lb) = named[keys[i]], named[keys[j]]
            r = conc(eff(Sa, la), eff(Sb, lb))
            p, sd = perm_p(Sa, la, Sb, lb, r)
            rows.append({"pair": f"{keys[i]} vs {keys[j]}", "r": r, "p_perm": p, "null_sd": sd})
            print(f"  {keys[i]:18s} vs {keys[j]:18s}  r={r:+.2f}  p={p:.3f}  (null sd {sd:.2f})")
    R = pd.DataFrame(rows); R.to_csv(MIL / "stageD_bulk.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(R))
    ax.bar(x, R.r, 0.6, color=["#c0392b" if p < 0.05 else "#95a5a6" for p in R.p_perm])
    for xi, r in zip(x, R.itertuples()):
        ax.annotate(f"p={r.p_perm:.3f}", (xi, r.r), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=8)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(R.pair, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("whole-tissue module-effect concordance (Pearson r)")
    ax.set_title("Step 3: bulk-compartment IPF module signature across 3 independent cohorts")
    plt.tight_layout(); fig.savefig(FIG / "stageD_bulk.png", dpi=130); plt.close(fig)
    print(f"\nsaved figure: {FIG/'stageD_bulk.png'}")
    print(f"saved table:  {MIL/'stageD_bulk.csv'}")


if __name__ == "__main__":
    main()
