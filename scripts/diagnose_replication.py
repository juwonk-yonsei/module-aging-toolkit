import os
import sys
from pathlib import Path
import numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore")
from scipy.stats import pearsonr, spearmanr
REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
RES = REPO / "results"

# 1. which celltypes passed + sample sizes in Habermann
gm = pd.read_csv(RES / "hab_pseudobulk_meta.csv")
gm = gm[gm.n_cells >= 50]
tab = gm[gm.Disease_Identity.isin(["IPF","Control"])].groupby(["celltype","Disease_Identity"]).size().unstack(fill_value=0)
tab["total"] = tab.sum(axis=1)
print("=== Habermann pseudobulk sizes (n_cells>=50) ===")
print(tab.sort_values("total", ascending=False).to_string())

# 2. gene coverage (symbol -> mouse Entrez) for Habermann
pb = pd.read_pickle(RES / "hab_pseudobulk_counts.pkl")
gmap = tp._load_gene_mapping("human", "Gene.Symbol")
mapped = pd.Index(pb.index.astype(str)).map(gmap)
print(f"\n=== gene coverage ===")
print(f"Habermann genes total: {pb.shape[0]}; mapped to mouse Entrez: {mapped.notna().sum()} ({100*mapped.notna().mean():.1f}%)")
# Adams coverage comparison
pba = pd.read_pickle(RES / "pseudobulk_counts.pkl")
gmap_e = tp._load_gene_mapping("human", "Ensembl")
mapped_a = pd.Index(pba.index.astype(str)).map(gmap_e)
print(f"Adams genes total: {pba.shape[0]}; mapped: {mapped_a.notna().sum()} ({100*mapped_a.notna().mean():.1f}%)")

# 3. shared coarse celltype x module full concordance (all, not just sig)
cmp = pd.read_csv(RES / "replication_coarse_compare.csv")
print(f"\n=== shared coarse celltype x module (all, n={len(cmp)}) ===")
print("shared coarse celltypes:", sorted(cmp.coarse.unique()))
r = pearsonr(cmp.beta_adams, cmp.beta_hab)[0]; rho = spearmanr(cmp.beta_adams, cmp.beta_hab)[0]
sc = np.mean(np.sign(cmp.beta_adams)==np.sign(cmp.beta_hab))
print(f"Pearson r={r:.2f}  Spearman={rho:.2f}  sign concordance={sc:.0%}")
# per coarse celltype
print("\nper coarse celltype concordance:")
for co, g in cmp.groupby("coarse"):
    if len(g) >= 5:
        rr = pearsonr(g.beta_adams, g.beta_hab)[0]
        s = np.mean(np.sign(g.beta_adams)==np.sign(g.beta_hab))
        print(f"  {co:12s} n={len(g):2d}  r={rr:+.2f}  sign={s:.0%}")

# 4. restrict global to SHARED coarse celltypes for fair comparison
D_h = pd.read_csv(RES / "hab_persample_module_scores.csv")
D_a = pd.read_csv(RES / "persample_module_scores.csv")
from B_module_aging.scripts.replication_hab import COARSE, ADAMS_COARSE, LAB
mods = list(LAB.values())
D_a["coarse"] = D_a.celltype.map(ADAMS_COARSE)
shared = sorted(set(D_h.coarse.dropna()) & set(D_a.coarse.dropna()))
print(f"\n=== fair global on shared coarse celltypes {shared} ===")
def global_beta(D):
    D = D[D.disease.isin(["Control","IPF"]) & D.coarse.isin(shared)].copy()
    out={}
    for m in mods:
        ipf=D[D.disease=="IPF"][m].mean(); ctl=D[D.disease=="Control"][m].mean()
        out[m]=ipf-ctl
    return pd.Series(out)
ba=global_beta(D_a); bh=global_beta(D_h)
r=pearsonr(ba,bh)[0]; sc=np.mean(np.sign(ba)==np.sign(bh))
print(f"Pearson r={r:.2f}  sign concordance={sc:.0%}")
