import os
import sys
from pathlib import Path
import numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import pearsonr
REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
RES = REPO / "results"

cmp = pd.read_csv(RES / "replication_coarse_compare.csv")
order = cmp.groupby("coarse").apply(lambda g: pearsonr(g.beta_adams, g.beta_hab)[0]).sort_values(ascending=False)
cells = order.index.tolist()

n = len(cells)
fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 3.4), sharex=False, sharey=False)
if n == 1: axes = [axes]
for ax, co in zip(axes, cells):
    g = cmp[cmp.coarse == co]
    r = pearsonr(g.beta_adams, g.beta_hab)[0]
    sc = np.mean(np.sign(g.beta_adams) == np.sign(g.beta_hab))
    ok = r > 0.3
    ax.scatter(g.beta_adams, g.beta_hab, s=22, color="#c0392b" if ok else "#888", alpha=0.8)
    lim = np.nanmax(np.abs(np.r_[g.beta_adams, g.beta_hab])) * 1.15
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.7, alpha=0.6)
    ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_title(f"{co}\nr={r:+.2f}, sign {sc:.0%}", fontsize=10,
                 color="#c0392b" if ok else "#333")
    ax.set_xlabel("Adams β", fontsize=8)
    if co == cells[0]: ax.set_ylabel("Habermann β", fontsize=9)
    ax.tick_params(labelsize=7)
fig.suptitle("Cross-cohort replication of module-aging signatures per cell type\n"
             "(each point = 1 of 23 modules; Adams=whole-lung dissociate, Habermann=biopsy)",
             fontsize=11)
plt.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(RES / "replication_percelltype.png", dpi=300)
print("saved:", RES / "replication_percelltype.png")

# overall summary line
allr = pearsonr(cmp.beta_adams, cmp.beta_hab)[0]
print(f"overall shared-celltype r={allr:.2f}, n={len(cmp)}")
for co in cells:
    g = cmp[cmp.coarse == co]
    print(f"  {co:11s} r={pearsonr(g.beta_adams,g.beta_hab)[0]:+.2f} sign={np.mean(np.sign(g.beta_adams)==np.sign(g.beta_hab)):.0%}")
