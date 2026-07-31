"""Attention-MIL infrastructure, step 1: build INSTANCE pseudobulks.

Instead of one pseudobulk per (donor x cell type), we split each
(donor x cell type) into mini-bags of ~BAG_SIZE cells (a "sub-pseudobulk"),
so a bag = (donor x cell type) and its instances = the mini-bags. This keeps
module-clock scoring statistically stable (each instance still has enough cells)
while giving attention something to weight over.

Streams the GSE136831 sparse MatrixMarket once, aggregates counts by
(donor x celltype x minibag) via chunked bincount (memory-safe), and caches:
  results/mil/inst_counts_{tag}.pkl   genes x n_instances (float32)
  results/mil/inst_meta_{tag}.csv     one row per instance

Cell -> mini-bag assignment is precomputed from the metadata (deterministic seed),
so the streaming only needs an int group-id per cell.
"""
import os
import sys, gzip, argparse, warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
IPF = REPO / "data" / "ipf"
OUT = REPO / "results" / "mil"
OUT.mkdir(parents=True, exist_ok=True)

CELLTYPE_COL = "Manuscript_Identity"
MTX = IPF / "GSE136831_RawCounts_Sparse.mtx.gz"

# cell types kept by default: myeloid/lymphoid/epithelial/stromal needed for
# Stage 1 (macrophage) + Stage 2 replication (coarse) + Stage 3 (mac subclass).
DEFAULT_CELLTYPES = [
    "Macrophage", "Macrophage_Alveolar", "cMonocyte", "ncMonocyte",
    "NK", "T", "T_Cytotoxic", "T_Regulatory", "B",
    "cDC1", "cDC2", "Ciliated", "ATII", "Fibroblast", "Myofibroblast",
]


def load_meta():
    m = pd.read_csv(IPF / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
    m.columns = [c.strip('"') for c in m.columns]
    for c in m.select_dtypes("object"):
        m[c] = m[c].str.strip('"')
    return m


def assign_instances(meta, celltypes, bag_size, min_cells, seed):
    """Return a Series (index=row order of meta) of instance ids (or -1),
    plus a per-instance metadata frame."""
    rng = np.random.default_rng(seed)
    inst_of_cell = np.full(len(meta), -1, dtype=np.int64)
    keep = (meta["CellType_Category"] != "Multiplet") & meta[CELLTYPE_COL].isin(celltypes) \
        & meta["Disease_Identity"].isin(["IPF", "Control", "COPD"])
    rows = []
    next_id = 0
    idx_all = np.arange(len(meta))
    grp = meta[CELLTYPE_COL].astype(str) + "||" + meta["Subject_Identity"].astype(str)
    grp = grp.where(keep, other=np.nan)
    for key, pos in pd.Series(idx_all).groupby(grp.values, sort=False):
        cells = pos.values
        if len(cells) < min_cells:
            continue
        rng.shuffle(cells)
        n_bags = max(1, int(round(len(cells) / bag_size)))
        for b, chunk in enumerate(np.array_split(cells, n_bags)):
            if len(chunk) == 0:
                continue
            inst_of_cell[chunk] = next_id
            ct, subj = key.split("||")
            sub = meta.iloc[chunk]
            subclass = sub["Subclass_Cell_Identity"].mode()
            rows.append({
                "instance": next_id, "celltype": ct, "Subject_Identity": subj,
                "Disease_Identity": sub["Disease_Identity"].iloc[0],
                "n_cells": len(chunk),
                "subclass_majority": subclass.iloc[0] if len(subclass) else "NA",
            })
            next_id += 1
    return inst_of_cell, pd.DataFrame(rows)


def build(celltypes, bag_size, min_cells, seed, tag):
    genes = pd.read_csv(IPF / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    genes.columns = [c.strip('"') for c in genes.columns]
    ens = genes["Ensembl_GeneID"].str.strip('"').values
    barcodes = pd.read_csv(IPF / "GSE136831_AllCells.cellBarcodes.txt.gz", header=None)[0].values
    meta = load_meta().set_index("CellBarcode_Identity").reindex(barcodes).reset_index()

    inst_of_cell, inst_meta = assign_instances(meta, celltypes, bag_size, min_cells, seed)
    ninst = int(inst_meta["instance"].max()) + 1 if len(inst_meta) else 0
    print(f"cells={len(barcodes)} genes={len(ens)} instances={ninst}")
    print(inst_meta.groupby(["celltype", "Disease_Identity"]).size().unstack(fill_value=0).to_string())

    f = gzip.open(MTX, "rt")
    for line in f:
        if line.startswith("%"):
            continue
        nrows, ncols, nnz = map(int, line.split())
        break
    print(f"mtx dims: {nrows} x {ncols}, nnz={nnz}")
    assert nrows == len(ens) and ncols == len(barcodes), "orientation mismatch"

    pb = np.zeros(nrows * ninst, dtype=np.float64)
    reader = pd.read_csv(f, sep=r"\s+", header=None, names=["g", "c", "v"],
                         dtype={"g": np.int32, "c": np.int32, "v": np.float64},
                         chunksize=40_000_000)
    seen = 0
    for ch in reader:
        gi = ch["g"].values - 1
        ci = ch["c"].values - 1
        v = ch["v"].values
        inst = inst_of_cell[ci]
        keep = inst >= 0
        lin = gi[keep].astype(np.int64) * ninst + inst[keep]
        pb += np.bincount(lin, weights=v[keep], minlength=nrows * ninst)
        seen += len(ch)
        print(f"  processed {seen:,}/{nnz:,} nnz", end="\r")
    print()
    pb = pb.reshape(nrows, ninst).astype(np.float32)

    cols = [str(i) for i in range(ninst)]
    pbdf = pd.DataFrame(pb, index=ens, columns=cols)
    pbdf.to_pickle(OUT / f"inst_counts_{tag}.pkl")
    inst_meta.to_csv(OUT / f"inst_meta_{tag}.csv", index=False)
    print(f"saved instance pseudobulk: {pbdf.shape} -> {OUT}/inst_counts_{tag}.pkl")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--bag_size", type=int, default=60)
    ap.add_argument("--min_cells", type=int, default=60)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="adams_bs60")
    ap.add_argument("--celltypes", nargs="*", default=DEFAULT_CELLTYPES)
    a = ap.parse_args()
    build(a.celltypes, a.bag_size, a.min_cells, a.seed, a.tag)
