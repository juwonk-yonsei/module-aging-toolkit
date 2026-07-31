"""Cross-species validity test: do rodent-trained module/composite transcriptomic
clocks predict HUMAN chronological age in healthy GTEx lung (n=578, ages 20-79)?

Design: apply the identical tAge pipeline used for IPF, but center relative to the
youngest donors (reference), then correlate predicted tAge with chronological age.
This is an independent, disease-free benchmark of cross-species transfer.
"""
import os
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
import joblib

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "python"))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
from tage_predict import _patch_simple_imputer

GTEX = REPO / "data" / "gtex"
MC = REPO / "data" / "module_clocks"
MODELS = REPO / "data" / "clock_models"
RES = REPO / "results"
FIG = RES / "figures"; FIG.mkdir(exist_ok=True)

# ---------------- load data ----------------
counts = pd.read_pickle(GTEX / "lung_counts.pkl")
meta = pd.read_csv(GTEX / "lung_meta.csv", index_col=0)
meta = meta.dropna(subset=["age_mid"])
counts = counts[meta.index]
# reference = two youngest brackets (20-39)
meta["ref"] = np.where(meta["age_mid"] <= 35, "young", "old")
print(f"GTEx lung: {counts.shape[1]} samples; young ref n={(meta.ref=='young').sum()}")

# ---------------- preprocess (young-centered relative profiles) ----------------
pp = tp.preprocess(counts, meta, species="human", gene_mapping_type="Ensembl",
                   control_group_column="ref", control_group_label="young")
X = pp["scaled_diff"]                      # samples x gene_list
X.columns = X.columns.map(str)
age = meta.loc[X.index, "age_mid"].values
sexM = meta.loc[X.index, "sexM"].values

# ---------------- technical / agonal covariates (GTEx post-mortem confounds) ----------------
tech = pd.read_csv(GTEX / "lung_tech.tsv", sep="\t", header=None,
                   names=["SAMPID", "isch", "rin"]).set_index("SAMPID")
ph = pd.read_csv(GTEX / "SubjectPhenotypes.txt", sep="\t").set_index("SUBJID")
cov = pd.DataFrame(index=X.index)
cov["isch"] = tech.reindex(X.index)["isch"].astype(float)
cov["rin"] = tech.reindex(X.index)["rin"].astype(float)
subj = pd.Index(X.index).str.split("-").str[:2].str.join("-")
cov["hardy"] = pd.Series(ph.reindex(subj)["DTHHRDY"].values, index=X.index)
cov["sexM"] = sexM
# design matrix for adjustment: isch, rin (median-imputed), hardy dummies, sex
C = pd.get_dummies(cov, columns=["hardy"], dummy_na=True).astype(float)
C = C.fillna(C.median(numeric_only=True))
C.insert(0, "intercept", 1.0)
Cmat = C.values

def adj_spearman(pred, age):
    """Spearman of age vs prediction residualized on technical/agonal covariates."""
    beta, *_ = np.linalg.lstsq(Cmat, pred, rcond=None)
    resid = pred - Cmat @ beta
    return spearmanr(resid, age)

# ---------------- module clocks ----------------
def load_clocks(outcome):
    C = {}
    for f in sorted(MC.glob(f"module_{outcome}_*.pkl")):
        d = joblib.load(f); C[d["annotation"].split("/")[0].strip()] = d
    return C

def predict_module(d, X):
    Xi = X.reindex(columns=[str(g) for g in d["genes"]])
    return d["pipeline"].predict(Xi.values)

def predict_composite(path, X):
    m = joblib.load(path)
    genes = [str(g) for g in m.feature_names_in_]
    _patch_simple_imputer(m.named_steps["imputation"])
    return m.predict(X.reindex(columns=genes).values)

rows = []
for outcome in ["Chrono", "Mortality"]:
    C = load_clocks(outcome)
    for name, d in C.items():
        pred = predict_module(d, X)
        rho, prho = spearmanr(pred, age)
        arho, aprho = adj_spearman(pred, age)
        rows.append(dict(outcome=outcome, module=name, n_genes=len(d["genes"]),
                         spearman=rho, sp_p=prho, adj_spearman=arho, adj_p=aprho))
    comp_path = (MODELS / (f"EN_{'Chronoage' if outcome=='Chrono' else 'Mortality'}"
                           "_Multispecies_Multitissue_scaleddiff.pkl"))
    cpred = predict_composite(comp_path, X)
    rho, prho = spearmanr(cpred, age); arho, aprho = adj_spearman(cpred, age)
    rows.append(dict(outcome=outcome, module="COMPOSITE", n_genes=np.nan,
                     spearman=rho, sp_p=prho, adj_spearman=arho, adj_p=aprho))
    if outcome == "Chrono":
        comp_chrono = cpred
    if outcome == "Mortality":
        comp_mort = cpred

T = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
T["fdr"] = np.nan; T["adj_fdr"] = np.nan
for oc in ["Chrono", "Mortality"]:
    mask = T.outcome == oc
    T.loc[mask, "fdr"] = multipletests(T.loc[mask, "sp_p"].fillna(1), method="fdr_bh")[1]
    T.loc[mask, "adj_fdr"] = multipletests(T.loc[mask, "adj_p"].fillna(1), method="fdr_bh")[1]
T = T.sort_values(["outcome", "adj_spearman"], ascending=[True, False])
T.to_csv(RES / "gtex_validation.csv", index=False)

# ---------------- report ----------------
pd.set_option("display.width", 220)
print("\n=== GTEx lung: clock vs chronological age (raw & covariate-adjusted Spearman) ===")
print("    adjusted for ischemic time, RIN, Hardy death scale, sex")
print(T[["outcome", "module", "n_genes", "spearman", "adj_spearman", "adj_p", "adj_fdr"]]
      .round(4).to_string(index=False))
for oc in ["Chrono", "Mortality"]:
    sub = T[T.outcome == oc]; comp = sub[sub.module == "COMPOSITE"].iloc[0]
    mods = sub[sub.module != "COMPOSITE"]
    nsig = ((mods.adj_spearman > 0) & (mods.adj_fdr < 0.05)).sum()
    print(f"\n[{oc}] composite adj Spearman={comp.adj_spearman:.3f} (p={comp.adj_p:.2e}); "
          f"modules positive & adj-FDR<0.05: {nsig}/23")

# ---------------- figures (focus on Mortality clock = the one used for IPF) ----------------
mort = T[T.outcome == "Mortality"]
comp_m = mort[mort.module == "COMPOSITE"].iloc[0]

# (A) composite mortality scatter
fig, ax = plt.subplots(figsize=(6.2, 6))
jitter = age + np.random.RandomState(0).uniform(-2, 2, len(age))
ax.scatter(jitter, comp_mort, s=14, alpha=0.5, color="#2c7fb8")
b1, b0 = np.polyfit(age, comp_mort, 1)
xs = np.array([age.min(), age.max()]); ax.plot(xs, b0 + b1 * xs, "r-", lw=2)
ax.set_xlabel("Chronological age (GTEx bracket midpoint, yrs)")
ax.set_ylabel("Predicted composite mortality tAge (a.u.)")
ax.set_title(f"Rodent composite mortality clock tracks HUMAN age in GTEx lung\n"
             f"Spearman ρ={comp_m.spearman:.2f}; covariate-adjusted ρ={comp_m.adj_spearman:.2f} "
             f"(p={comp_m.adj_p:.1e}), n={len(age)}")
plt.tight_layout(); fig.savefig(FIG / "gtex_composite_scatter.png", dpi=300); plt.close(fig)

# (B) per-module adjusted Spearman bar (Mortality)
m2 = mort[mort.module != "COMPOSITE"].sort_values("adj_spearman")
colors = ["#b2182b" if (r > 0 and f < 0.05) else ("#2166ac" if (r < 0 and f < 0.05) else "#cccccc")
          for r, f in zip(m2.adj_spearman, m2.adj_fdr)]
fig, ax = plt.subplots(figsize=(7, 8))
y = np.arange(len(m2))
ax.barh(y, m2.adj_spearman, color=colors, edgecolor="#333", lw=0.4)
ax.set_yticks(y); ax.set_yticklabels(m2.module, fontsize=8)
ax.axvline(0, color="k", lw=0.7)
ax.axvline(comp_m.adj_spearman, color="#2c7fb8", ls="--", lw=1.5,
           label=f"composite adj ρ={comp_m.adj_spearman:.2f}")
ax.set_xlabel("Covariate-adjusted Spearman ρ (module mortality-tAge vs age)")
ax.set_title("Cross-species transfer per module in human lung (GTEx)\n"
             "red=+ & FDR<0.05, blue=− & FDR<0.05, grey=n.s.")
ax.legend(loc="lower right", fontsize=9)
plt.tight_layout(); fig.savefig(FIG / "gtex_module_transfer.png", dpi=300); plt.close(fig)
print("\nsaved figures: gtex_composite_scatter.png, gtex_module_transfer.png")
