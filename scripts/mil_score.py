"""Attention-MIL infrastructure, step 2: score INSTANCE pseudobulks with the
23 retrained mortality module clocks.

For each cell type separately, preprocess its instance pseudobulks with the
ported tAge pipeline (control-subtraction relative to Control-donor instances)
and apply every module clock -> instance x 23-module score matrix.

Output: results/mil/inst_scores_{tag}.csv  (instance meta + one column per module)
"""
import os
import sys, argparse, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import joblib

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
MC = REPO / "data" / "module_clocks"
OUT = REPO / "results" / "mil"

# ---- load mortality module clocks
CLOCKS = {}
for fp in sorted(MC.glob("module_Mortality_*.pkl")):
    d = joblib.load(fp)
    CLOCKS[f"{d['module']}|{d['annotation']}"] = d


def make_labels(fns):
    from collections import Counter
    lab = {fn: fn.split("|", 1)[1].split("/")[0].strip() for fn in fns}
    cnt = Counter(lab.values())
    for fn, s in list(lab.items()):
        if cnt[s] > 1:
            lab[fn] = f"{s} ({fn.split('|')[0]})"
    return lab


LAB = make_labels(list(CLOCKS.keys()))
MODULES = list(LAB.values())


def score_matrix(scaled_diff):
    X = scaled_diff.copy(); X.columns = X.columns.map(str)
    return pd.DataFrame(
        {LAB[n]: d["pipeline"].predict(X.reindex(columns=[str(g) for g in d["genes"]]).values)
         for n, d in CLOCKS.items()},
        index=scaled_diff.index)


def score_celltype(counts, meta, gene_mapping_type="Ensembl", shared_genes=None):
    """counts: genes x instances; meta indexed by instance-id string with
    Disease_Identity. Returns instances x 23 module scores.
    If shared_genes (set of Entrez id strings) is given, genes not in it are
    masked to NaN so both cohorts impute an identical gene support (Option A)."""
    pp = tp.preprocess(counts, meta, species="human", gene_mapping_type=gene_mapping_type,
                       control_group_column="Disease_Identity", control_group_label="Control")
    sd = pp["scaled_diff"]                                   # instances x gene_list
    if shared_genes is not None:
        drop = [c for c in sd.columns if str(c) not in shared_genes]
        sd = sd.copy(); sd[drop] = np.nan
    return score_matrix(sd)


def run(tag, gene_mapping_type="Ensembl", min_inst_per_class=8):
    counts = pd.read_pickle(OUT / f"inst_counts_{tag}.pkl")
    meta = pd.read_csv(OUT / f"inst_meta_{tag}.csv")
    meta["instance"] = meta["instance"].astype(str)
    meta = meta.set_index("instance")
    counts.columns = counts.columns.map(str)

    out = []
    for ct, sub in meta.groupby("celltype"):
        n = sub["Disease_Identity"].value_counts()
        if n.get("Control", 0) < min_inst_per_class or n.get("IPF", 0) < min_inst_per_class:
            print(f"  skip {ct}: too few instances {n.to_dict()}")
            continue
        cols = sub.index.tolist()
        S = score_celltype(counts[cols], sub, gene_mapping_type=gene_mapping_type)
        S = sub.join(S)
        out.append(S)
        print(f"  scored {ct}: {len(sub)} instances")
    allr = pd.concat(out)
    allr.to_csv(OUT / f"inst_scores_{tag}.csv")
    print(f"saved {OUT}/inst_scores_{tag}.csv  shape={allr.shape}  celltypes={allr.celltype.nunique()}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="adams_bs60")
    ap.add_argument("--gene_mapping", default="Ensembl", help="Ensembl or Gene.Symbol")
    a = ap.parse_args()
    run(a.tag, gene_mapping_type=a.gene_mapping)
