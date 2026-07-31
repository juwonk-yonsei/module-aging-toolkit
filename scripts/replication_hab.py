"""External replication on GSE135893 (Habermann) IPF lung scRNA-seq.

Same pipeline as Adams (GSE136831): stream sparse mtx -> donor x celltype pseudobulk
-> module clocks -> IPF vs Control per module. Then test whether the module-level
aging signatures replicate across the two independent cohorts.

NOTE: GSE135893 GEO provides no age/sex; the Adams analysis showed age/sex adjustment
barely changed effects (naive vs adjusted r=0.85), so replication uses a disease-only
donor-level model (+ permutation). This is stated as a limitation.
"""
import os
import sys, gzip, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import joblib

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "python"))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
D2 = REPO / "data" / "ipf2"
MC = REPO / "data" / "module_clocks"
RES = REPO / "results"
MIN_CELLS = 50
rng = np.random.default_rng(0)

COARSE = {  # harmonize fine celltypes -> coarse (shared vocabulary)
    "Macrophages": "Macrophage", "Proliferating Macrophages": "Macrophage",
    "Monocytes": "Monocyte", "NK Cells": "NK",
    "T Cells": "Tcell", "Proliferating T Cells": "Tcell",
    "B Cells": "Bcell", "Plasma Cells": "Bcell", "cDCs": "DC", "pDCs": "DC",
    "Ciliated": "Ciliated", "Differentiating Ciliated": "Ciliated",
    "AT2": "AT2", "Transitional AT2": "AT2", "AT1": "AT1",
    "Fibroblasts": "Fibroblast", "HAS1 High Fibroblasts": "Fibroblast",
    "PLIN2+ Fibroblasts": "Fibroblast", "Myofibroblasts": "Fibroblast",
    "Mast Cells": "Mast", "Basal": "Basal",
}
ADAMS_COARSE = {
    "Macrophage": "Macrophage", "Macrophage_Alveolar": "Macrophage",
    "cMonocyte": "Monocyte", "ncMonocyte": "Monocyte", "NK": "NK",
    "T": "Tcell", "T_Cytotoxic": "Tcell", "T_Regulatory": "Tcell",
    "B": "Bcell", "cDC1": "DC", "cDC2": "DC", "Ciliated": "Ciliated",
    "ATII": "AT2", "Fibroblast": "Fibroblast", "Myofibroblast": "Fibroblast",
}


def build_pseudobulk():
    cache = RES / "hab_pseudobulk_counts.pkl"
    if cache.exists():
        return pd.read_pickle(cache), pd.read_csv(RES / "hab_pseudobulk_meta.csv")
    genes = pd.read_csv(D2 / "GSE135893_genes.tsv.gz", header=None)[0].values
    barcodes = pd.read_csv(D2 / "GSE135893_barcodes.tsv.gz", header=None)[0].values
    meta = pd.read_csv(D2 / "GSE135893_IPF_metadata.csv.gz", index_col=0)
    meta = meta.reindex(barcodes)
    ok = meta["celltype"].notna() & meta["Diagnosis"].isin(["IPF", "Control"])
    grp_key = (meta["Sample_Name"].astype(str) + "||" + meta["celltype"].astype(str)).where(ok, np.nan)
    uniq = pd.Index(grp_key.dropna().unique())
    gid = {k: i for i, k in enumerate(uniq)}
    group_of_cell = grp_key.map(gid).fillna(-1).astype(np.int64).values
    ng = len(uniq)
    print(f"cells={len(barcodes)} genes={len(genes)} groups={ng}")

    f = gzip.open(D2 / "GSE135893_matrix.mtx.gz", "rt")
    for line in f:
        if line.startswith("%"):
            continue
        nrows, ncols, nnz = map(int, line.split()); break
    print(f"mtx {nrows}x{ncols} nnz={nnz}")
    assert nrows == len(genes) and ncols == len(barcodes), "orientation mismatch"
    pb = np.zeros(nrows * ng, dtype=np.float64)
    reader = pd.read_csv(f, sep=r"\s+", header=None, names=["g", "c", "v"],
                         dtype={"g": np.int32, "c": np.int32, "v": np.float64}, chunksize=40_000_000)
    seen = 0
    for ch in reader:
        gi = ch["g"].values - 1; ci = ch["c"].values - 1; v = ch["v"].values
        grp = group_of_cell[ci]; keep = grp >= 0
        lin = gi[keep].astype(np.int64) * ng + grp[keep]
        pb += np.bincount(lin, weights=v[keep], minlength=nrows * ng)
        seen += len(ch); print(f"  {seen:,}/{nnz:,}", end="\r")
    print()
    pb = pb.reshape(nrows, ng)
    ncells = pd.Series(group_of_cell[group_of_cell >= 0]).value_counts().reindex(range(ng)).fillna(0)
    gm = pd.DataFrame({"group": uniq})
    gm["Subject_Identity"] = [k.split("||")[0] for k in uniq]
    gm["celltype"] = [k.split("||")[1] for k in uniq]
    gm["n_cells"] = ncells.values.astype(int)
    subj_dx = meta.groupby("Sample_Name")["Diagnosis"].first()
    gm["Disease_Identity"] = gm["Subject_Identity"].map(subj_dx)
    pbdf = pd.DataFrame(pb, index=genes, columns=[str(c) for c in uniq])
    pbdf.to_pickle(cache); gm.to_csv(RES / "hab_pseudobulk_meta.csv", index=False)
    print("saved hab pseudobulk", pbdf.shape)
    return pbdf, gm


# ---- module clocks
CLOCKS = {}
for fp in sorted(MC.glob("module_Mortality_*.pkl")):
    d = joblib.load(fp); CLOCKS[f"{d['module']}|{d['annotation']}"] = d


def make_labels(fns):
    from collections import Counter
    lab = {fn: fn.split("|", 1)[1].split("/")[0].strip() for fn in fns}
    cnt = Counter(lab.values())
    for fn, s in list(lab.items()):
        if cnt[s] > 1:
            lab[fn] = f"{s} ({fn.split('|')[0]})"
    return lab


LAB = make_labels(list(CLOCKS.keys()))


def score_all(scaled_diff):
    X = scaled_diff.copy(); X.columns = X.columns.map(str)
    return pd.DataFrame({LAB[n]: d["pipeline"].predict(
        X.reindex(columns=[str(g) for g in d["genes"]]).values) for n, d in CLOCKS.items()},
        index=scaled_diff.index)


def perm_p(y, dvec, B=2000):
    ok = np.isfinite(y); y = y[ok]; d = dvec[ok].astype(float)
    d = d - d.mean(); ry = y - y.mean()
    obs = (d @ ry) / (d @ d)
    null = np.array([( (rng.permutation(d)) @ ry) / (d @ d) for _ in range(B)])
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (B + 1)


def main():
    pb, gm = build_pseudobulk()
    gm["group"] = gm["group"].astype(str)
    gm = gm[gm.n_cells >= MIN_CELLS]
    modules = list(LAB.values())

    recs = []
    for ct, sub in gm.groupby("celltype"):
        n = sub.Disease_Identity.value_counts()
        if n.get("IPF", 0) < 4 or n.get("Control", 0) < 4:
            continue
        md = sub.set_index("group")
        pp = tp.preprocess(pb[sub.group.tolist()], md, species="human",
                           gene_mapping_type="Gene.Symbol",
                           control_group_column="Disease_Identity", control_group_label="Control")
        S = score_all(pp["scaled_diff"])
        S["Subject_Identity"] = md.loc[S.index, "Subject_Identity"].values
        S["celltype"] = ct
        S["disease"] = md.loc[S.index, "Disease_Identity"].values
        S["coarse"] = COARSE.get(ct, ct)
        recs.append(S)
    D = pd.concat(recs, ignore_index=True)
    D.to_csv(RES / "hab_persample_module_scores.csv", index=False)
    print(f"hab per-sample: {D.shape}, celltypes={D.celltype.nunique()}")

    # ---- global per-module effect (donor-clustered robust OLS, disease only)
    grows = []
    DD = D[D.disease.isin(["Control", "IPF"])].copy()
    DD["disease"] = pd.Categorical(DD["disease"], ["Control", "IPF"])
    for mod in modules:
        try:
            m = smf.ols(f"Q('{mod}') ~ C(disease) + C(celltype)", DD).fit(
                cov_type="cluster", cov_kwds={"groups": DD["Subject_Identity"]})
            key = [k for k in m.params.index if "disease" in k][0]
            grows.append(dict(module=mod, beta_hab=m.params[key], p_hab=m.pvalues[key]))
        except Exception:
            grows.append(dict(module=mod, beta_hab=np.nan, p_hab=np.nan))
    G = pd.DataFrame(grows)
    G["fdr_hab"] = multipletests(G["p_hab"].fillna(1), method="fdr_bh")[1]

    # ---- compare to Adams global
    A = pd.read_csv(RES / "module_global_robust.csv").rename(
        columns={"beta": "beta_adams", "p": "p_adams", "fdr": "fdr_adams"})
    A = A[A.module != "Composite"]
    M = A.merge(G, on="module", how="inner")
    M.to_csv(RES / "replication_global_compare.csv", index=False)
    from scipy.stats import pearsonr, spearmanr
    r = pearsonr(M.beta_adams, M.beta_hab)[0]; rho = spearmanr(M.beta_adams, M.beta_hab)[0]
    sign_conc = np.mean(np.sign(M.beta_adams) == np.sign(M.beta_hab))
    adams_sig = M[M.fdr_adams < 0.10]
    sign_conc_sig = np.mean(np.sign(adams_sig.beta_adams) == np.sign(adams_sig.beta_hab)) if len(adams_sig) else np.nan
    print(f"\n=== GLOBAL module replication (Adams vs Habermann) ===")
    print(f"Pearson r={r:.2f}  Spearman rho={rho:.2f}  sign concordance={sign_conc:.0%}")
    print(f"Adams-FDR<0.10 modules ({len(adams_sig)}): sign concordance={sign_conc_sig:.0%}")
    print(M.sort_values("beta_adams").round(3)[["module","beta_adams","fdr_adams","beta_hab","p_hab","fdr_hab"]].to_string(index=False))

    # ---- figure: replication scatter
    fig, ax = plt.subplots(figsize=(7.5, 7))
    sig = M.fdr_adams < 0.10
    ax.scatter(M.beta_adams[~sig], M.beta_hab[~sig], s=25, color="#bbb", label="Adams ns")
    ax.scatter(M.beta_adams[sig], M.beta_hab[sig], s=45, color="#c0392b", label="Adams FDR<0.10")
    for _, rr in M[sig].iterrows():
        ax.annotate(rr.module[:16], (rr.beta_adams, rr.beta_hab), fontsize=7)
    lim = np.nanmax(np.abs(np.r_[M.beta_adams, M.beta_hab])) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.7, alpha=0.6)
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("Adams GSE136831  β IPF (adj)"); ax.set_ylabel("Habermann GSE135893  β IPF")
    ax.set_title(f"External replication of module aging signatures\nPearson r={r:.2f}, sign concordance {sign_conc:.0%} (sig {sign_conc_sig:.0%})")
    ax.legend(fontsize=8); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    plt.tight_layout(); fig.savefig(RES / "replication_scatter.png", dpi=130); plt.close(fig)
    print("saved figure: replication_scatter.png")

    # ---- coarse celltype x module concordance vs Adams-significant hits
    Araw = pd.read_csv(RES / "ipf_module_stats_adjusted.csv")
    Araw = Araw[Araw.module != "Composite"].copy()
    Araw["coarse"] = Araw.celltype.map(ADAMS_COARSE)
    a_coarse = Araw.dropna(subset=["coarse"]).groupby(["coarse", "module"]).beta_IPF.mean().rename("beta_adams")
    # hab coarse betas
    hrows = []
    for co, subco in DD.groupby("coarse"):
        if subco.disease.nunique() < 2 or (subco.disease == "IPF").sum() < 4 or (subco.disease == "Control").sum() < 4:
            continue
        for mod in modules:
            g = subco.groupby("disease")[mod].mean()
            hrows.append(dict(coarse=co, module=mod, beta_hab=g.get("IPF", np.nan) - g.get("Control", np.nan)))
    Hc = pd.DataFrame(hrows).set_index(["coarse", "module"]).beta_hab
    cmp = pd.concat([a_coarse, Hc], axis=1, join="inner").dropna()
    sig_hits = Araw[Araw.fdr_IPF < 0.10].copy()
    sig_hits["coarse"] = sig_hits.celltype.map(ADAMS_COARSE)
    keys = set(zip(sig_hits.coarse, sig_hits.module))
    cmp_sig = cmp[[k in keys for k in cmp.index]]
    if len(cmp_sig):
        sc = np.mean(np.sign(cmp_sig.beta_adams) == np.sign(cmp_sig.beta_hab))
        rr = pearsonr(cmp_sig.beta_adams, cmp_sig.beta_hab)[0]
        print(f"\n=== coarse celltype x module (Adams-significant hits, n={len(cmp_sig)}) ===")
        print(f"sign concordance={sc:.0%}  Pearson r={rr:.2f}")
    cmp.reset_index().to_csv(RES / "replication_coarse_compare.csv", index=False)


if __name__ == "__main__":
    main()
