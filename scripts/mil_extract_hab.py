"""Attention-MIL infrastructure: build INSTANCE pseudobulks for the external
replication cohort GSE135893 (Habermann). Mirrors mil_extract.py but for the
Habermann file layout (gene symbols; metadata columns celltype/Diagnosis/Sample_Name).

Instances = mini-bags of ~BAG_SIZE cells within each (donor x fine-celltype).
"""
import os
import sys, gzip, argparse, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
D2 = REPO / "data" / "ipf2"
OUT = REPO / "results" / "mil"
OUT.mkdir(parents=True, exist_ok=True)

DEFAULT_CELLTYPES = [
    "Macrophages", "Proliferating Macrophages", "Monocytes",
    "NK Cells", "T Cells", "Proliferating T Cells", "B Cells", "Plasma Cells",
    "cDCs", "pDCs", "Ciliated", "Differentiating Ciliated",
    "AT2", "Transitional AT2", "Fibroblasts", "Myofibroblasts",
    "PLIN2+ Fibroblasts", "HAS1 High Fibroblasts",
]


def assign_instances(meta, celltypes, bag_size, min_cells, seed):
    rng = np.random.default_rng(seed)
    inst_of_cell = np.full(len(meta), -1, dtype=np.int64)
    keep = meta["celltype"].isin(celltypes) & meta["Diagnosis"].isin(["IPF", "Control"])
    grp = (meta["celltype"].astype(str) + "||" + meta["Sample_Name"].astype(str)).where(keep, np.nan)
    rows = []; next_id = 0
    idx_all = np.arange(len(meta))
    for key, pos in pd.Series(idx_all).groupby(grp.values, sort=False):
        cells = pos.values
        if len(cells) < min_cells:
            continue
        rng.shuffle(cells)
        n_bags = max(1, int(round(len(cells) / bag_size)))
        for chunk in np.array_split(cells, n_bags):
            if len(chunk) == 0:
                continue
            inst_of_cell[chunk] = next_id
            ct, subj = key.split("||")
            sub = meta.iloc[chunk]
            rows.append({"instance": next_id, "celltype": ct, "Subject_Identity": subj,
                         "Disease_Identity": sub["Diagnosis"].iloc[0], "n_cells": len(chunk),
                         "subclass_majority": ct})
            next_id += 1
    return inst_of_cell, pd.DataFrame(rows)


def build(celltypes, bag_size, min_cells, seed, tag):
    genes = pd.read_csv(D2 / "GSE135893_genes.tsv.gz", header=None)[0].values
    barcodes = pd.read_csv(D2 / "GSE135893_barcodes.tsv.gz", header=None)[0].values
    meta = pd.read_csv(D2 / "GSE135893_IPF_metadata.csv.gz", index_col=0).reindex(barcodes).reset_index()

    inst_of_cell, inst_meta = assign_instances(meta, celltypes, bag_size, min_cells, seed)
    ninst = int(inst_meta["instance"].max()) + 1 if len(inst_meta) else 0
    print(f"cells={len(barcodes)} genes={len(genes)} instances={ninst}")
    print(inst_meta.groupby(["celltype", "Disease_Identity"]).size().unstack(fill_value=0).to_string())

    f = gzip.open(D2 / "GSE135893_matrix.mtx.gz", "rt")
    for line in f:
        if line.startswith("%"):
            continue
        nrows, ncols, nnz = map(int, line.split()); break
    print(f"mtx {nrows}x{ncols} nnz={nnz}")
    assert nrows == len(genes) and ncols == len(barcodes), "orientation mismatch"

    pb = np.zeros(nrows * ninst, dtype=np.float64)
    reader = pd.read_csv(f, sep=r"\s+", header=None, names=["g", "c", "v"],
                         dtype={"g": np.int32, "c": np.int32, "v": np.float64}, chunksize=40_000_000)
    seen = 0
    for ch in reader:
        gi = ch["g"].values - 1; ci = ch["c"].values - 1; v = ch["v"].values
        inst = inst_of_cell[ci]; keep = inst >= 0
        lin = gi[keep].astype(np.int64) * ninst + inst[keep]
        pb += np.bincount(lin, weights=v[keep], minlength=nrows * ninst)
        seen += len(ch); print(f"  {seen:,}/{nnz:,}", end="\r")
    print()
    pb = pb.reshape(nrows, ninst).astype(np.float32)
    pbdf = pd.DataFrame(pb, index=genes, columns=[str(i) for i in range(ninst)])
    pbdf.to_pickle(OUT / f"inst_counts_{tag}.pkl")
    inst_meta.to_csv(OUT / f"inst_meta_{tag}.csv", index=False)
    print(f"saved instance pseudobulk: {pbdf.shape} -> {OUT}/inst_counts_{tag}.pkl")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--bag_size", type=int, default=60)
    ap.add_argument("--min_cells", type=int, default=60)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="hab_bs60")
    ap.add_argument("--celltypes", nargs="*", default=DEFAULT_CELLTYPES)
    a = ap.parse_args()
    build(a.celltypes, a.bag_size, a.min_cells, a.seed, a.tag)
