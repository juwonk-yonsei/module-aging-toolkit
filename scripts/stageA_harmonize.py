"""Option A (tech-first baseline): re-score both cohorts on a COMMON detectable
gene support, to test whether removing the gene-coverage / imputation confound
changes cross-cohort replication.

Shared support = mouse-Entrez genes that map AND pass the detection filter
(filter_genes) in BOTH cohorts' instance pseudobulks. Genes outside it are
masked to NaN identically in both cohorts (so imputation is symmetric).

Outputs: results/mil/inst_scores_adams_bs60_shared.csv,
         results/mil/inst_scores_hab_bs60_shared.csv
"""
import os
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
import mil_score as ms
MIL = REPO / "results" / "mil"


def detectable(tag, gene_mapping_type):
    counts = pd.read_pickle(MIL / f"inst_counts_{tag}.pkl")
    counts.columns = counts.columns.map(str)
    f = tp.filter_genes(counts)
    m = tp.map_genes(f, "human", gene_mapping_type)
    return set(m.index.astype(str))


def main():
    da = detectable("adams_bs60", "Ensembl")
    dh = detectable("hab_bs60", "Gene.Symbol")
    shared = da & dh
    print(f"detectable Entrez: Adams={len(da)}  Habermann={len(dh)}  shared={len(shared)}")

    for tag, gmap in [("adams_bs60", "Ensembl"), ("hab_bs60", "Gene.Symbol")]:
        counts = pd.read_pickle(MIL / f"inst_counts_{tag}.pkl"); counts.columns = counts.columns.map(str)
        meta = pd.read_csv(MIL / f"inst_meta_{tag}.csv"); meta["instance"] = meta["instance"].astype(str)
        meta = meta.set_index("instance")
        out = []
        for ct, sub in meta.groupby("celltype"):
            n = sub["Disease_Identity"].value_counts()
            if n.get("Control", 0) < 8 or n.get("IPF", 0) < 8:
                continue
            S = ms.score_celltype(counts[sub.index.tolist()], sub, gene_mapping_type=gmap,
                                  shared_genes=shared)
            out.append(sub.join(S))
        allr = pd.concat(out)
        allr.to_csv(MIL / f"inst_scores_{tag}_shared.csv")
        print(f"  saved inst_scores_{tag}_shared.csv  shape={allr.shape}  celltypes={allr.celltype.nunique()}")


if __name__ == "__main__":
    main()
