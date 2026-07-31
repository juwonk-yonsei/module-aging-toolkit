"""IPF (GSE136831) module-resolved transcriptomic aging pipeline.

Phase 1: stream the sparse MatrixMarket counts -> donor x cell-type pseudobulk
         (memory-safe chunked bincount aggregation). Cached to disk.
Phase 2: per cell type, preprocess (human->mouse Entrez, relative to Control),
         apply composite multispecies chrono + mortality clocks, decompose the
         mortality clock into per-module contributions, and test IPF vs Control.

Outputs (results/):
  pseudobulk_counts.parquet, pseudobulk_meta.csv
  ipf_tage_percell.csv          (per donor x celltype tAge + module contributions)
  ipf_celltype_summary.csv      (IPF vs Control effect per cell type)
"""
import os
import sys, gzip, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")
REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "python"))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
from tage_predict import predict_tAge

IPF = REPO / "data" / "ipf"
MODELS = REPO / "data" / "clock_models"
SUPP = REPO / "data" / "supp"
RES = REPO / "results"
RES.mkdir(exist_ok=True)

MIN_CELLS = 50            # min cells per (donor x celltype) pseudobulk
CELLTYPE_COL = "Manuscript_Identity"
MTX = IPF / "GSE136831_RawCounts_Sparse.mtx.gz"


def load_meta():
    m = pd.read_csv(IPF / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
    m.columns = [c.strip('"') for c in m.columns]
    for c in m.select_dtypes("object"):
        m[c] = m[c].str.strip('"')
    m = m.set_index("CellBarcode_Identity")
    return m


def build_pseudobulk():
    genes = pd.read_csv(IPF / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    genes.columns = [c.strip('"') for c in genes.columns]
    ens = genes["Ensembl_GeneID"].str.strip('"').values
    barcodes = pd.read_csv(IPF / "GSE136831_AllCells.cellBarcodes.txt.gz", header=None)[0].values
    meta = load_meta().reindex(barcodes)

    # define groups = (Subject x celltype), exclude Multiplet / missing
    ok = (meta["CellType_Category"] != "Multiplet") & meta[CELLTYPE_COL].notna()
    grp_key = meta["Subject_Identity"].astype(str) + "||" + meta[CELLTYPE_COL].astype(str)
    grp_key = grp_key.where(ok, other=np.nan)
    uniq = pd.Index(grp_key.dropna().unique())
    gid_map = {k: i for i, k in enumerate(uniq)}
    group_of_cell = grp_key.map(gid_map).fillna(-1).astype(np.int64).values
    ngroups = len(uniq)
    print(f"cells={len(barcodes)} genes={len(ens)} groups={ngroups}")

    # stream mtx
    f = gzip.open(MTX, "rt")
    for line in f:
        if line.startswith("%"):
            continue
        nrows, ncols, nnz = map(int, line.split())
        break
    print(f"mtx dims: {nrows} x {ncols}, nnz={nnz}")
    assert nrows == len(ens) and ncols == len(barcodes), "orientation mismatch"

    pb = np.zeros(nrows * ngroups, dtype=np.float64)
    reader = pd.read_csv(f, sep=r"\s+", header=None, names=["g", "c", "v"],
                         dtype={"g": np.int32, "c": np.int32, "v": np.float64},
                         chunksize=40_000_000)
    seen = 0
    for ch in reader:
        gi = ch["g"].values - 1
        ci = ch["c"].values - 1
        v = ch["v"].values
        grp = group_of_cell[ci]
        keep = grp >= 0
        lin = gi[keep].astype(np.int64) * ngroups + grp[keep]
        pb += np.bincount(lin, weights=v[keep], minlength=nrows * ngroups)
        seen += len(ch)
        print(f"  processed {seen:,}/{nnz:,} nnz", end="\r")
    print()
    pb = pb.reshape(nrows, ngroups)

    ncells = pd.Series(group_of_cell[group_of_cell >= 0]).value_counts().reindex(range(ngroups)).fillna(0)
    gmeta = pd.DataFrame({"group": uniq})
    gmeta["Subject_Identity"] = [k.split("||")[0] for k in uniq]
    gmeta[CELLTYPE_COL] = [k.split("||")[1] for k in uniq]
    gmeta["n_cells"] = ncells.values.astype(int)
    # disease per subject
    subj_disease = load_meta().groupby("Subject_Identity")["Disease_Identity"].first()
    gmeta["Disease_Identity"] = gmeta["Subject_Identity"].map(subj_disease)

    pbdf = pd.DataFrame(pb, index=ens, columns=uniq)
    pbdf.columns = [str(c) for c in pbdf.columns]
    pbdf.to_pickle(RES / "pseudobulk_counts.pkl")
    gmeta.to_csv(RES / "pseudobulk_meta.csv", index=False)
    print("saved pseudobulk:", pbdf.shape)
    return pbdf, gmeta


def load_module_membership():
    mem = pd.read_csv(SUPP / "module_membership_multispecies.csv")
    mem["entrez"] = mem["entrez"].astype(str)
    return dict(zip(mem["entrez"], mem["module"])), mem


def decompose_mortality(model_path, X_df, module_map):
    """Return per-sample per-module contribution (in raw estimator units, pre species-adj)."""
    import joblib, sklearn.pipeline
    from tage_predict import _patch_simple_imputer
    m = joblib.load(model_path)
    genes = list(m.feature_names_in_)
    X = X_df.reindex(columns=[str(g) for g in genes])
    imp = m.named_steps["imputation"]; sca = m.named_steps["scaler"]
    _patch_simple_imputer(imp)
    sel = m.named_steps["featureselection"]; est = m.named_steps["estimator"]
    Xi = imp.transform(X.values)
    Xs = sca.transform(Xi)
    support = sel.get_support()
    Xsel = Xs[:, support]
    sel_genes = [str(g) for g, k in zip(genes, support) if k]
    contrib = Xsel * est.coef_[None, :]           # samples x selected
    mods = np.array([module_map.get(g, "unassigned") for g in sel_genes])
    out = {}
    for mod in np.unique(mods):
        out[f"mod_{mod}"] = contrib[:, mods == mod].sum(axis=1)
    return pd.DataFrame(out, index=X_df.index)


def run_clocks(pbdf, gmeta):
    mod_map, _ = load_module_membership()
    clocks = {
        "Chrono": MODELS / "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",
        "Mortality": MODELS / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl",
    }
    gmeta = gmeta[gmeta["n_cells"] >= MIN_CELLS].copy()
    results = []
    for ct, sub in gmeta.groupby(CELLTYPE_COL):
        disease = sub["Disease_Identity"]
        if (disease == "Control").sum() < 3 or (disease == "IPF").sum() < 3:
            continue
        cols = sub["group"].astype(str).tolist()
        counts = pbdf[cols]
        m_ct = sub.set_index(sub["group"].astype(str))
        pp = tp.preprocess(counts, m_ct, species="human", gene_mapping_type="Ensembl",
                           control_group_column="Disease_Identity", control_group_label="Control")
        X = pp["scaled_diff"]
        ann = m_ct.copy()
        for name, path in clocks.items():
            ann = predict_tAge(str(path), X, ann, species="human", prefix=f"{name}_")
        dec = decompose_mortality(clocks["Mortality"], X, mod_map) * 122.5  # species adj
        ann = pd.concat([ann, dec.reindex(ann.index)], axis=1)
        ann["celltype"] = ct
        results.append(ann)
    allr = pd.concat(results)
    allr.to_csv(RES / "ipf_tage_percell.csv")

    # summary IPF vs Control per celltype
    rows = []
    modcols = [c for c in allr.columns if c.startswith("mod_")]
    for ct, sub in allr.groupby("celltype"):
        ipf = sub[sub["Disease_Identity"] == "IPF"]
        ctl = sub[sub["Disease_Identity"] == "Control"]
        for clock in ["Chrono", "Mortality"]:
            col = f"{clock}_tAge"
            try:
                u, p = mannwhitneyu(ipf[col], ctl[col], alternative="two-sided")
            except ValueError:
                p = np.nan
            rows.append({"celltype": ct, "metric": col,
                         "IPF_mean": ipf[col].mean(), "Control_mean": ctl[col].mean(),
                         "delta": ipf[col].mean() - ctl[col].mean(),
                         "p": p, "n_IPF": len(ipf), "n_Control": len(ctl)})
    summ = pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests
    mmask = summ["metric"] == "Mortality_tAge"
    summ.loc[mmask, "padj"] = multipletests(summ.loc[mmask, "p"].fillna(1), method="fdr_bh")[1]
    summ = summ.sort_values(["metric", "delta"], ascending=[True, False])
    summ.to_csv(RES / "ipf_celltype_summary.csv", index=False)

    print("\n=== IPF vs Control: Mortality tAge by cell type (top accelerated) ===")
    mt = summ[summ["metric"] == "Mortality_tAge"].sort_values("delta", ascending=False)
    print(mt[["celltype", "Control_mean", "IPF_mean", "delta", "p", "padj", "n_IPF", "n_Control"]]
          .head(20).to_string(index=False))
    # module contribution: mean IPF-Control delta per module across cell types
    print("\n=== Mean module contribution to IPF mortality-aging (IPF - Control) ===")
    dd = []
    for ct, sub in allr.groupby("celltype"):
        ipf = sub[sub["Disease_Identity"] == "IPF"][modcols].mean()
        ctl = sub[sub["Disease_Identity"] == "Control"][modcols].mean()
        dd.append(ipf - ctl)
    mod_delta = pd.concat(dd, axis=1).mean(axis=1).sort_values(ascending=False)
    print(mod_delta.round(3).to_string())


if __name__ == "__main__":
    if (RES / "pseudobulk_counts.pkl").exists():
        print("loading cached pseudobulk")
        pbdf = pd.read_pickle(RES / "pseudobulk_counts.pkl")
        gmeta = pd.read_csv(RES / "pseudobulk_meta.csv")
    else:
        pbdf, gmeta = build_pseudobulk()
    run_clocks(pbdf, gmeta)
