"""Option B: marker-based cross-cohort SUB-STATE alignment.

Hypothesis: myeloid non-replication is driven by cell-composition mismatch
(whole-lung dissociation vs biopsy), so comparing the coarse cell type across
cohorts mixes different populations. If instead we split each myeloid cell type
into cohort-invariant sub-states by the SAME marker rule and compare matched
sub-states, the IPF-vs-Control module effects should replicate.

Pipeline (per cohort):
  pass 1 (marker rows only): per-cell mono-vs-resident signature -> median split
          within cohort x cell type -> sub-state label {A, B}.
  pass 2 (full counts):      aggregate into (donor x cell type x sub-state)
          pseudobulks -> module scores -> IPF-vs-Control effect per module.
Then compare Adams<->Habermann concordance: coarse (pooled) vs matched sub-state,
plus a random-split negative control.

Outputs: results/mil/stageB_substate.csv, results/figures/stageB_substate.png
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
MIL = REPO / "results" / "mil"
FIG = REPO / "results" / "figures"

MONO = ["SPP1", "FCN1", "VCAN", "CTHRC1", "INHBA", "MERTK", "SPARC"]
RES = ["FABP4", "MARCO", "MRC1", "PPARG"]
ADAMS_CT = {"Macrophage": ["Macrophage", "Macrophage_Alveolar"],
            "Monocyte": ["cMonocyte", "ncMonocyte"]}
HAB_CT = {"Macrophage": ["Macrophages", "Proliferating Macrophages"],
          "Monocyte": ["Monocytes"]}
MIN_CELLS = 40          # min cells per (donor x ct x substate) pseudobulk


def sym2ens():
    g = pd.read_csv(IPFD / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    g.columns = [c.strip('"') for c in g.columns]
    for c in g.select_dtypes("object"):
        g[c] = g[c].str.strip('"')
    return dict(zip(g["HGNC_EnsemblAlt_GeneID"], g["Ensembl_GeneID"]))


def stream_marker_counts(mtx, marker_rows, ncols):
    """Return (n_markers x ncols) dense array of counts for given 0-based rows."""
    row_set = {r: i for i, r in enumerate(marker_rows)}
    out = np.zeros((len(marker_rows), ncols), dtype=np.float64)
    f = gzip.open(mtx, "rt")
    for line in f:
        if line.startswith("%"):
            continue
        nr, nc, nnz = map(int, line.split()); break
    reader = pd.read_csv(f, sep=r"\s+", header=None, names=["g", "c", "v"],
                         dtype={"g": np.int32, "c": np.int32, "v": np.float64}, chunksize=40_000_000)
    for ch in reader:
        gi = ch["g"].values - 1
        keep = np.isin(gi, marker_rows)
        if keep.any():
            gk = gi[keep]; ck = ch["c"].values[keep] - 1; vk = ch["v"].values[keep]
            for r in np.unique(gk):
                out[row_set[r], ck[gk == r]] += vk[gk == r]
    return out


def stream_group_pb(mtx, group_of_cell, ngroups, nrows):
    pb = np.zeros(nrows * ngroups, dtype=np.float64)
    f = gzip.open(mtx, "rt")
    for line in f:
        if line.startswith("%"):
            continue
        break
    reader = pd.read_csv(f, sep=r"\s+", header=None, names=["g", "c", "v"],
                         dtype={"g": np.int32, "c": np.int32, "v": np.float64}, chunksize=40_000_000)
    for ch in reader:
        gi = ch["g"].values - 1; ci = ch["c"].values - 1; v = ch["v"].values
        grp = group_of_cell[ci]; keep = grp >= 0
        lin = gi[keep].astype(np.int64) * ngroups + grp[keep]
        pb += np.bincount(lin, weights=v[keep], minlength=nrows * ngroups)
    return pb.reshape(nrows, ngroups)


def zscore(a):
    return (a - np.nanmean(a)) / (np.nanstd(a) + 1e-9)


def build_cohort(cohort, randomize=False, seed=0):
    """Return dict: pseudobulk counts (genes x groups), group meta with
    (donor, celltype, substate, disease), gene ids, gene_mapping. Cached to disk."""
    tag = f"{cohort}{'_rand'+str(seed) if randomize else ''}"
    cache_pb = MIL / f"stageB_pb_{tag}.pkl"; cache_gm = MIL / f"stageB_gm_{tag}.csv"
    gmap = "Ensembl" if cohort == "adams" else "Gene.Symbol"
    if cache_pb.exists() and cache_gm.exists():
        return {"pb": pd.read_pickle(cache_pb), "gm": pd.read_csv(cache_gm), "gmap": gmap}
    rng = np.random.default_rng(seed)
    if cohort == "adams":
        genes = pd.read_csv(IPFD / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
        genes.columns = [c.strip('"') for c in genes.columns]
        gid = genes["Ensembl_GeneID"].str.strip('"').values
        barc = pd.read_csv(IPFD / "GSE136831_AllCells.cellBarcodes.txt.gz", header=None)[0].values
        meta = pd.read_csv(IPFD / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
        meta.columns = [c.strip('"') for c in meta.columns]
        for c in meta.select_dtypes("object"):
            meta[c] = meta[c].str.strip('"')
        meta = meta.set_index("CellBarcode_Identity").reindex(barc)
        ctcol, donorcol, discol, libcol = "Manuscript_Identity", "Subject_Identity", "Disease_Identity", "nUMI"
        ctmap = ADAMS_CT; mtx = IPFD / "GSE136831_RawCounts_Sparse.mtx.gz"; gmap = "Ensembl"
        S2E = sym2ens()
        mono_ids = [S2E[m] for m in MONO if m in S2E]
        res_ids = [S2E[m] for m in RES if m in S2E]
        id_index = pd.Index(gid)
    else:
        gid = pd.read_csv(D2 / "GSE135893_genes.tsv.gz", header=None)[0].values
        barc = pd.read_csv(D2 / "GSE135893_barcodes.tsv.gz", header=None)[0].values
        meta = pd.read_csv(D2 / "GSE135893_IPF_metadata.csv.gz", index_col=0).reindex(barc)
        ctcol, donorcol, discol, libcol = "celltype", "Sample_Name", "Diagnosis", "nCount_RNA"
        ctmap = HAB_CT; mtx = D2 / "GSE135893_matrix.mtx.gz"; gmap = "Gene.Symbol"
        mono_ids = [m for m in MONO]; res_ids = [m for m in RES]
        id_index = pd.Index(gid)

    fine2coarse = {f: c for c, fs in ctmap.items() for f in fs}
    coarse = meta[ctcol].map(fine2coarse)
    lib = pd.to_numeric(meta[libcol], errors="coerce").values

    # pass 1: marker counts -> per-cell signature
    mono_rows = [id_index.get_loc(x) for x in mono_ids if x in id_index]
    res_rows = [id_index.get_loc(x) for x in res_ids if x in id_index]
    marker_rows = sorted(set(mono_rows) | set(res_rows))
    mc = stream_marker_counts(mtx, marker_rows, len(barc))
    rowpos = {r: i for i, r in enumerate(marker_rows)}
    cpm = np.log1p(mc / (lib[None, :] + 1e-9) * 1e6)
    mono_sig = np.nanmean([cpm[rowpos[r]] for r in mono_rows], axis=0)
    res_sig = np.nanmean([cpm[rowpos[r]] for r in res_rows], axis=0)

    # substate label per cell: median split within cohort x coarse celltype
    substate = np.array(["NA"] * len(barc), dtype=object)
    for ct in ctmap:
        m = (coarse == ct).values & np.isfinite(lib)
        if m.sum() < 10:
            continue
        score = zscore(mono_sig[m]) - zscore(res_sig[m])
        if randomize:
            rng.shuffle(score)
        thr = np.nanmedian(score)
        lab = np.where(score >= thr, "monoHi", "resHi")
        substate[np.where(m)[0]] = lab

    # groups = donor x coarse ct x substate, disease in IPF/Control
    donor = meta[donorcol].astype(str).values
    dis = meta[discol].astype(str).values
    ok = pd.Series(coarse).notna().values & np.isin(dis, ["IPF", "Control"]) & (substate != "NA")
    key = pd.Series([f"{donor[i]}|{coarse.values[i]}|{substate[i]}" if ok[i] else np.nan
                     for i in range(len(barc))])
    uniq = pd.Index(key.dropna().unique())
    gmapd = {k: i for i, k in enumerate(uniq)}
    group_of_cell = key.map(gmapd).fillna(-1).astype(np.int64).values
    pb = stream_group_pb(mtx, group_of_cell, len(uniq), len(gid))
    ncells = pd.Series(group_of_cell[group_of_cell >= 0]).value_counts().reindex(range(len(uniq))).fillna(0)

    gm = pd.DataFrame({"group": uniq})
    gm["donor"] = [k.split("|")[0] for k in uniq]
    gm["celltype"] = [k.split("|")[1] for k in uniq]
    gm["substate"] = [k.split("|")[2] for k in uniq]
    subj_dis = pd.Series(dis, index=donor).groupby(level=0).first()
    gm["Disease_Identity"] = gm["donor"].map(subj_dis)
    gm["n_cells"] = ncells.values.astype(int)
    pbdf = pd.DataFrame(pb, index=gid, columns=[str(u) for u in uniq])
    pbdf.to_pickle(cache_pb); gm.to_csv(cache_gm, index=False)
    return {"pb": pbdf, "gm": gm, "gmap": gmap}


def effect_by(pbdf, gm, gmap, by_substate):
    """Return dict {(celltype[,substate]): 23-module IPF-Control effect}."""
    gm = gm[gm.n_cells >= MIN_CELLS].copy()
    out = {}
    group_cols = ["celltype", "substate"] if by_substate else ["celltype"]
    for keys, sub in gm.groupby(group_cols):
        keyt = keys if isinstance(keys, tuple) else (keys,)
        n = sub.Disease_Identity.value_counts()
        if n.get("IPF", 0) < 3 or n.get("Control", 0) < 3:
            continue
        cols = sub.group.astype(str).tolist()
        md = sub.set_index(sub.group.astype(str))
        pp = tp.preprocess(pbdf[cols], md, species="human", gene_mapping_type=gmap,
                           control_group_column="Disease_Identity", control_group_label="Control")
        S = score_matrix(pp["scaled_diff"])
        S["dis"] = md.loc[S.index, "Disease_Identity"].values
        eff = S[S.dis == "IPF"][MODULES].mean() - S[S.dis == "Control"][MODULES].mean()
        out[keyt] = eff.values
    return out


def conc(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return pearsonr(a[ok], b[ok])[0] if ok.sum() >= 3 else np.nan


def main():
    print("building Adams sub-state pseudobulks ...")
    A = build_cohort("adams")
    print("building Habermann sub-state pseudobulks ...")
    H = build_cohort("hab")
    print("building Adams RANDOM-split control ...")
    Ar = build_cohort("adams", randomize=True, seed=1)
    Hr = build_cohort("hab", randomize=True, seed=1)

    # coarse (pooled) effects
    Ac = effect_by(A["pb"], A["gm"], A["gmap"], by_substate=False)
    Hc = effect_by(H["pb"], H["gm"], H["gmap"], by_substate=False)
    # matched sub-state effects
    As = effect_by(A["pb"], A["gm"], A["gmap"], by_substate=True)
    Hs = effect_by(H["pb"], H["gm"], H["gmap"], by_substate=True)
    # random-split effects
    Ars = effect_by(Ar["pb"], Ar["gm"], Ar["gmap"], by_substate=True)
    Hrs = effect_by(Hr["pb"], Hr["gm"], Hr["gmap"], by_substate=True)

    rows = []
    for ct in ["Macrophage", "Monocyte"]:
        r_coarse = conc(Ac.get((ct,), np.full(len(MODULES), np.nan)),
                        Hc.get((ct,), np.full(len(MODULES), np.nan)))
        for ss in ["monoHi", "resHi"]:
            r_ss = conc(As.get((ct, ss), np.full(len(MODULES), np.nan)),
                        Hs.get((ct, ss), np.full(len(MODULES), np.nan)))
            r_rand = conc(Ars.get((ct, ss), np.full(len(MODULES), np.nan)),
                          Hrs.get((ct, ss), np.full(len(MODULES), np.nan)))
            rows.append({"celltype": ct, "substate": ss, "r_coarse": r_coarse,
                         "r_substate": r_ss, "r_random": r_rand})
            print(f"  {ct:11s} {ss:7s}  coarse r={r_coarse:+.2f}  matched-substate r={r_ss:+.2f}  "
                  f"random-split r={r_rand:+.2f}")
    R = pd.DataFrame(rows)
    R.to_csv(MIL / "stageB_substate.csv", index=False)

    # figure
    fig, ax = plt.subplots(figsize=(10, 5.5))
    labels = [f"{r.celltype}\n{r.substate}" for _, r in R.iterrows()]
    x = np.arange(len(R)); w = 0.27
    ax.bar(x - w, R.r_coarse, w, label="coarse (pooled)", color="#2c3e50")
    ax.bar(x, R.r_substate, w, label="matched sub-state", color="#c0392b")
    ax.bar(x + w, R.r_random, w, label="random-split control", color="#bbbbbb")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("Adams <-> Habermann module-effect concordance (Pearson r)")
    ax.set_title("Option B: does marker-based sub-state alignment rescue myeloid replication?")
    ax.legend()
    plt.tight_layout(); fig.savefig(FIG / "stageB_substate.png", dpi=130); plt.close(fig)
    print(f"\nsaved figure: {FIG/'stageB_substate.png'}")
    print(f"saved table:  {MIL/'stageB_substate.csv'}")


if __name__ == "__main__":
    main()
