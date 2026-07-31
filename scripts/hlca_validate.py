"""Single-cell modality validation: within each healthy lung cell type, do the rodent
module/composite clocks track donor chronological age across 42 donors (ages 20-81)?
Tests whether module-level cross-species transfer RECOVERS at donor x cell-type pseudobulk
resolution (vs the weak bulk-tissue GTEx result)."""
import os
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import joblib

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "python"))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
from tage_predict import _patch_simple_imputer

H = REPO / "data" / "hlca"
MC = REPO / "data" / "module_clocks"
MODELS = REPO / "data" / "clock_models"
RES = REPO / "results"; FIG = RES / "figures"
MIN_DONORS = 12

pb = pd.read_pickle(H / "pseudobulk_counts.pkl")
meta = pd.read_csv(H / "pseudobulk_meta.csv")
meta["key"] = meta["key"].astype(str)

def load_clocks(outcome):
    C = {}
    for f in sorted(MC.glob(f"module_{outcome}_*.pkl")):
        d = joblib.load(f); C[d["annotation"].split("/")[0].strip()] = d
    return C

CH = load_clocks("Chrono"); MO = load_clocks("Mortality")
COMP = {"Chrono": MODELS / "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",
        "Mortality": MODELS / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"}

def predict_module(d, X):
    return d["pipeline"].predict(X.reindex(columns=[str(g) for g in d["genes"]]).values)

def predict_comp(path, X):
    m = joblib.load(path); genes = [str(g) for g in m.feature_names_in_]
    _patch_simple_imputer(m.named_steps["imputation"])
    return m.predict(X.reindex(columns=genes).values)

rows = []
for ct, sub in meta.groupby("celltype"):
    if len(sub) < MIN_DONORS:
        continue
    cols = sub.key.tolist()
    md = sub.set_index("key")
    pp = tp.preprocess(pb[cols], md, species="human", gene_mapping_type="Ensembl")
    X = pp["scaled_diff"]; X.columns = X.columns.map(str)
    age = md.loc[X.index, "age"].values.astype(float)
    for outcome, C in [("Chrono", CH), ("Mortality", MO)]:
        for name, d in C.items():
            rho, p = spearmanr(predict_module(d, X), age)
            rows.append(dict(celltype=ct, outcome=outcome, module=name,
                             n=len(age), spearman=rho, p=p))
        rho, p = spearmanr(predict_comp(COMP[outcome], X), age)
        rows.append(dict(celltype=ct, outcome=outcome, module="COMPOSITE",
                         n=len(age), spearman=rho, p=p))

T = pd.DataFrame(rows)
for oc in ["Chrono", "Mortality"]:
    m = T.outcome == oc
    T.loc[m, "fdr"] = multipletests(T.loc[m, "p"].fillna(1), method="fdr_bh")[1]
T.to_csv(RES / "hlca_validation.csv", index=False)

pd.set_option("display.width", 220)
for oc in ["Mortality", "Chrono"]:
    sub = T[T.outcome == oc]
    comp = sub[sub.module == "COMPOSITE"]
    mods = sub[sub.module != "COMPOSITE"]
    print(f"\n===== {oc} =====")
    print(f"cell types tested: {sub.celltype.nunique()}")
    print(f"COMPOSITE: median rho={comp.spearman.median():.3f}; "
          f"celltypes with rho>0 & FDR<0.10: {((comp.spearman>0)&(comp.fdr<0.10)).sum()}/{len(comp)}")
    sig = mods[(mods.spearman > 0) & (mods.fdr < 0.10)]
    print(f"module x celltype effects positive & FDR<0.10: {len(sig)} / {len(mods)}")
    print("top positive module transfers:")
    print(sig.sort_values("spearman", ascending=False)
          .head(15)[["celltype", "module", "n", "spearman", "fdr"]].round(3).to_string(index=False))

# composite per-celltype summary table
comp_all = T[T.module == "COMPOSITE"].pivot_table(index="celltype", columns="outcome", values="spearman")
print("\n=== composite Spearman(age) per cell type ===")
print(comp_all.round(3).to_string())

# ---------------- figure: composite per-celltype + module recovery heatmap ----------------
mo = T[T.outcome == "Mortality"]
mat = mo[mo.module != "COMPOSITE"].pivot_table(index="celltype", columns="module", values="spearman")
fmat = mo[mo.module != "COMPOSITE"].pivot_table(index="celltype", columns="module", values="fdr").reindex_like(mat)
mat = mat.loc[mat.mean(1).sort_values(ascending=False).index]
mat = mat[mat.mean().sort_values(ascending=False).index]
fmat = fmat.reindex(index=mat.index, columns=mat.columns)
fig, ax = plt.subplots(figsize=(14, 5.5))
v = np.nanpercentile(np.abs(mat.values), 95)
im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=9)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if fmat.values[i, j] < 0.10:
            ax.text(j, i, "*", ha="center", va="center", fontsize=10)
ax.set_title("Single-cell modality: module mortality-tAge vs donor age (Spearman ρ) per cell type\n"
             "healthy lung, 42 donors 20-81y; * FDR<0.10; red=tracks age (older→higher)")
plt.colorbar(im, label="Spearman ρ (module tAge vs age)", shrink=0.7)
plt.tight_layout(); fig.savefig(FIG / "hlca_module_recovery.png", dpi=300); plt.close(fig)
print("\nsaved figure: hlca_module_recovery.png")
