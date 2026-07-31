"""Validate retrained module clocks on Klotho, then apply to IPF (module-resolved map)."""
import os
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import joblib
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 40)

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
MC = REPO / "data" / "module_clocks"
RES = REPO / "results"

# load module clocks
CLOCKS = {"Mortality": {}, "Chrono": {}}
for f in sorted(MC.glob("module_*.pkl")):
    d = joblib.load(f)
    CLOCKS[d["outcome"]][f"{d['module']}|{d['annotation']}"] = d


def predict_modules(scaled_diff, outcome="Mortality"):
    X = scaled_diff.copy(); X.columns = X.columns.map(str)
    out = {}
    for name, d in CLOCKS[outcome].items():
        Xi = X.reindex(columns=[str(g) for g in d["genes"]])
        out[name] = d["pipeline"].predict(Xi.values)
    return pd.DataFrame(out, index=scaled_diff.index)


def short(name):  # pretty label
    mod, ann = name.split("|", 1)
    return f"{ann.split('/')[0].strip()}"


def make_labels(fullnames):
    from collections import Counter
    lab = {fn: fn.split("|", 1)[1].split("/")[0].strip() for fn in fullnames}
    cnt = Counter(lab.values())
    for fn, s in list(lab.items()):
        if cnt[s] > 1:
            lab[fn] = f"{s} ({fn.split('|')[0]})"
    return lab


_ALL = list(CLOCKS["Mortality"].keys())
LAB = make_labels(_ALL)


# ================= 1. KLOTHO VALIDATION =================
print("========== KLOTHO re-validation (retrained module clocks) ==========")
EXT = Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "extdata"
counts = pd.read_csv(EXT / "Exprs_example.csv", index_col=0)
meta = pd.read_csv(EXT / "Metadata_example.csv", index_col=0)
for tissue in meta["Tissue"].unique():
    m_t = meta[meta.Tissue == tissue]; c_t = counts[m_t.index]
    pp = tp.preprocess(c_t, m_t, species="mouse", gene_mapping_type="Ensembl",
                       control_group_column="Genotype", control_group_label="WT")
    S = predict_modules(pp["scaled_diff"], "Mortality"); S["G"] = m_t["Genotype"].values
    d = S.groupby("G").mean(numeric_only=True).T
    d["KO-WT"] = d.get("Klotho KO", 0) - d.get("WT", 0)
    d.index = [short(i) for i in d.index]
    print(f"\n--- {tissue}: KO-WT mortality module score (top) ---")
    print(d["KO-WT"].sort_values(ascending=False).head(8).round(3).to_string())

# ================= 2. IPF MODULE-RESOLVED MAP =================
print("\n========== IPF module-resolved analysis ==========")
pb = pd.read_pickle(RES / "pseudobulk_counts.pkl")
meta_pb = pd.read_csv(RES / "pseudobulk_meta.csv")
meta_pb["group"] = meta_pb["group"].astype(str)
meta_pb = meta_pb[meta_pb.n_cells >= 50]

rows = []
for ct, sub in meta_pb.groupby("Manuscript_Identity"):
    sub = sub[sub.Disease_Identity.isin(["IPF", "Control"])]
    n_ipf = (sub.Disease_Identity == "IPF").sum(); n_ctl = (sub.Disease_Identity == "Control").sum()
    if n_ipf < 4 or n_ctl < 4:
        continue
    cols = sub.group.tolist()
    md = sub.set_index("group")
    pp = tp.preprocess(pb[cols], md, species="human", gene_mapping_type="Ensembl",
                       control_group_column="Disease_Identity", control_group_label="Control")
    S = predict_modules(pp["scaled_diff"], "Mortality")
    S["grp"] = md.loc[S.index, "Disease_Identity"].values
    ipf = S[S.grp == "IPF"]; ctl = S[S.grp == "Control"]
    for mod in [c for c in S.columns if c != "grp"]:
        try:
            p = mannwhitneyu(ipf[mod], ctl[mod]).pvalue
        except Exception:
            p = np.nan
        rows.append({"celltype": ct, "module": LAB[mod], "module_full": mod,
                     "delta": ipf[mod].mean() - ctl[mod].mean(), "p": p,
                     "n_ipf": n_ipf, "n_ctl": n_ctl})

R = pd.DataFrame(rows)
R.to_csv(RES / "ipf_module_resolved.csv", index=False)

# heatmap celltype x module (delta), ordered
mat = R.pivot_table(index="celltype", columns="module", values="delta")
# order cell types by mean immune signal, modules by variance
mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]
mat = mat[mat.var().sort_values(ascending=False).index]
pmat = R.pivot_table(index="celltype", columns="module", values="p").reindex_like(mat)

fig, ax = plt.subplots(figsize=(15, 7))
vmax = np.nanpercentile(np.abs(mat.values), 95)
im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=9)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if pmat.values[i, j] < 0.05:
            ax.text(j, i, "*", ha="center", va="center", color="k", fontsize=11)
ax.set_title("IPF module-resolved mortality-aging (Δ IPF−Control per cell type × module)\n* p<0.05 Mann–Whitney; red=aged, blue=younger")
plt.colorbar(im, label="Δ module mortality score", shrink=0.7)
plt.tight_layout(); fig.savefig(RES / "ipf_module_heatmap.png", dpi=130)
print(f"saved heatmap: {RES/'ipf_module_heatmap.png'}")

print("\n=== significant (p<0.05) cell-type × module effects ===")
sig = R[R.p < 0.05].sort_values("p")
print(sig[["celltype", "module", "delta", "p", "n_ipf", "n_ctl"]].round(3).to_string(index=False))

print("\n=== Fibroblast / Alveolar-Mac module breakdown (artifact check) ===")
for ct in ["Fibroblast", "Macrophage_Alveolar", "NK", "T"]:
    if ct in R.celltype.values:
        sub = R[R.celltype == ct].sort_values("delta", ascending=False)
        print(f"\n{ct}:")
        print(sub[["module", "delta", "p"]].round(3).to_string(index=False))
