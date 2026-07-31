"""Retrain the 23 rodent module-specific clocks from public Zenodo data.

Reproduces the paper's relative elastic-net approach: input = relative (control-
subtracted) scaled expression; targets = the *difference* hazard / chronological age.
Group-aware CV by dataset (Source) for honest out-of-fold performance.
Outputs: per-module sklearn pipelines (.pkl) + LOFO performance table.
"""
import os
import sys, warnings, json
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import ElasticNetCV
from sklearn.base import clone
from sklearn.model_selection import GroupKFold
from scipy.stats import pearsonr, spearmanr
import joblib

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
ROD = REPO / "data" / "rodent"
SUPP = REPO / "data" / "supp"
OUT = REPO / "data" / "module_clocks"
OUT.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "Mortality": "Difference.Hazard.log10",
    "Chrono": "Difference.Chronological_age.Normalized_by_species_max_lifespan",
}

print("loading annotation + membership ...")
ann = pd.read_excel(ROD / "Data_annotation_relative_rodents.xlsx").set_index("Sample")
mem = pd.read_csv(SUPP / "module_membership_rodent.csv")
mem["entrez"] = mem["entrez"].astype(str)
mod2genes = {m: g["entrez"].tolist() for m, g in mem.groupby("module")}
mod2annot = dict(zip(mem.module, mem.annotation))

print("loading expression (721MB, may take ~1min) ...")
expr = pd.read_csv(ROD / "Expression_data_relative_rodents_Scaled.csv", index_col=0)
# orient so that index=genes, columns=samples
if not set(ann.index[:5]).issubset(set(expr.columns)):
    expr = expr.T
expr.index = expr.index.astype(str)
common = [s for s in ann.index if s in expr.columns]
ann = ann.loc[common]; expr = expr[common]
print(f"expr genes x samples = {expr.shape}; aligned samples = {len(common)}")
groups = ann["Source"].values
n_groups = ann["Source"].nunique()
print(f"datasets (groups) = {n_groups}; species = {ann.Species.value_counts().to_dict()}")

rows = []
for outcome, ycol in TARGETS.items():
    y = pd.to_numeric(ann[ycol], errors="coerce").values
    ok = ~np.isnan(y)
    for mod, genes in mod2genes.items():
        g = [x for x in genes if x in expr.index]
        if len(g) < 5:
            continue
        X = expr.loc[g].T.values[ok]
        yy = y[ok]; gg = groups[ok]
        n_g = len(np.unique(gg)); k = min(5, n_g)
        splits = list(GroupKFold(n_splits=k).split(X, yy, gg))
        base = make_pipeline(
            SimpleImputer(strategy="median"), StandardScaler(),
            ElasticNetCV(l1_ratio=[.5, .9, 1.0], alphas=np.logspace(-3, 1, 20),
                         cv=3, max_iter=4000, n_jobs=-1),
        )
        oof = np.full(len(yy), np.nan)
        for tr, te in splits:
            mdl = clone(base).fit(X[tr], yy[tr])
            oof[te] = mdl.predict(X[te])
        r = pearsonr(oof, yy)[0]; rho = spearmanr(oof, yy)[0]
        final = clone(base).fit(X, yy)
        joblib.dump({"pipeline": final, "genes": g, "module": mod,
                     "annotation": mod2annot[mod], "outcome": outcome},
                    OUT / f"module_{outcome}_{mod}.pkl")
        rows.append({"outcome": outcome, "module": mod, "annotation": mod2annot[mod],
                     "n_genes": len(g), "LOFO_pearson": round(r, 3),
                     "LOFO_spearman": round(rho, 3)})
        print(f"  [{outcome}] {mod:14s} {mod2annot[mod][:32]:32s} n={len(g):3d}  r={r:.3f}")

perf = pd.DataFrame(rows)
perf.to_csv(OUT / "module_clock_performance.csv", index=False)
print("\n===== performance summary (Mortality) =====")
print(perf[perf.outcome == "Mortality"].sort_values("LOFO_pearson", ascending=False).to_string(index=False))
print(f"\nsaved {len(rows)} module clocks to {OUT}")
