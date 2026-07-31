"""Generalized disease-vs-control module-aging contrast for a new tissue/disease.

Reads data/<tag>/pseudobulk_counts.pkl + pseudobulk_meta.csv (cols: key, donor,
celltype, n_cells, disease, grp[Control/<DISEASE>], sex, age), scores module +
composite mortality clocks per cell type (relative to Control), then fits per
cell-type x module models (grp + sex [+ age if available]) with partial-regression
permutation p and BH-FDR, plus a donor-clustered global test. Writes stats CSV,
heatmap, and a summary md.

Usage: python disease_module_pipeline.py <tag> <DISEASE_LABEL>
"""
import os
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import joblib

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
sys.path.insert(0, str(Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "python"))
import tage_prep as tp
from tage_predict import _patch_simple_imputer
MC = REPO / "data" / "module_clocks"
CM = REPO / "data" / "clock_models"
RES = REPO / "results"; FIG = RES / "figures"
rng = np.random.default_rng(0)
MIN_N = 4

CLOCKS = {}
for f in sorted(MC.glob("module_Mortality_*.pkl")):
    d = joblib.load(f); CLOCKS[f"{d['module']}|{d['annotation']}"] = d
comp = joblib.load(CM / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
_patch_simple_imputer(comp.named_steps["imputation"])
comp_genes = [str(g) for g in comp.feature_names_in_]


def make_labels(fullnames):
    from collections import Counter
    lab = {fn: fn.split("|", 1)[1].split("/")[0].strip() for fn in fullnames}
    cnt = Counter(lab.values())
    for fn, s in list(lab.items()):
        if cnt[s] > 1:
            lab[fn] = f"{s} ({fn.split('|')[0]})"
    return lab


LAB = make_labels(list(CLOCKS.keys()))


def score_all(scaled_diff):
    out = {}
    X = scaled_diff.copy(); X.columns = X.columns.map(str)
    for name, d in CLOCKS.items():
        out[LAB[name]] = d["pipeline"].predict(X.reindex(columns=[str(g) for g in d["genes"]]).values)
    out["Composite"] = comp.predict(X.reindex(columns=comp_genes).values) * 122.5
    return pd.DataFrame(out, index=scaled_diff.index)


def perm_p(y, d_case, X0, B=2000):
    ok = np.isfinite(y)
    y, d, X0 = y[ok], d_case[ok], X0[ok]
    H = X0 @ np.linalg.pinv(X0.T @ X0) @ X0.T
    M = np.eye(len(y)) - H
    ry, rd = M @ y, M @ d
    obs = (rd @ ry) / (rd @ rd)
    null = np.empty(B)
    for b in range(B):
        rdp = M @ rng.permutation(d)
        null[b] = (rdp @ ry) / (rdp @ rdp)
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (B + 1)


def main(tag, disease):
    DAT = REPO / "data" / tag
    pb = pd.read_pickle(DAT / "pseudobulk_counts.pkl")
    meta = pd.read_csv(DAT / "pseudobulk_meta.csv"); meta["key"] = meta["key"].astype(str)
    meta = meta[meta.n_cells >= 50]
    age_frac = meta["age"].notna().mean()
    has_age = age_frac > 0.4
    print(f"[{tag}] {disease} vs Control | pseudobulks={len(meta)} donors={meta.donor.nunique()} "
          f"age_available={has_age} (parsed {age_frac*100:.0f}%)")

    # per-sample module scores (relative to Control within each cell type)
    recs = []
    for ct, sub in meta.groupby("celltype"):
        sub = sub[sub.grp.isin([disease, "Control"])]
        n = sub.grp.value_counts()
        if n.get(disease, 0) < MIN_N or n.get("Control", 0) < MIN_N:
            continue
        md = sub.set_index("key")
        pp = tp.preprocess(pb[sub.key.tolist()], md, species="human", gene_mapping_type="Ensembl",
                           control_group_column="grp", control_group_label="Control")
        S = score_all(pp["scaled_diff"])
        S["donor"] = md.loc[S.index, "donor"].values
        S["celltype"] = ct
        S["grp"] = md.loc[S.index, "grp"].values
        S["sexM"] = (md.loc[S.index, "sex"].astype(str).str.lower().str[0] == "m").astype(float).values
        S["age"] = md.loc[S.index, "age"].values
        recs.append(S)
    if not recs:
        print(f"[{tag}] no cell types with enough samples"); return
    D = pd.concat(recs, ignore_index=True)
    if has_age:                       # impute partially-missing donor ages (median) to keep N
        D["age"] = D["age"].fillna(D["age"].median())
    modules = list(LAB.values()) + ["Composite"]
    real_modules = list(LAB.values())
    D.to_csv(RES / f"{tag}_persample_module_scores.csv", index=False)

    covs = "grp + sexM" + (" + age" if has_age else "")
    rows = []
    for ct, g in D.groupby("celltype"):
        cc = g[g.grp.isin(["Control", disease])].copy()
        cc = cc.dropna(subset=["sexM"] + (["age"] if has_age else []))
        if (cc.grp == disease).sum() < MIN_N or (cc.grp == "Control").sum() < MIN_N:
            continue
        cc["grp"] = pd.Categorical(cc["grp"], ["Control", disease])
        cols0 = [np.ones(len(cc)), cc["sexM"].values] + ([cc["age"].values] if has_age else [])
        X0 = np.column_stack(cols0)
        d_case = (cc["grp"] == disease).astype(float).values
        for mod in modules:
            naive = cc.groupby("grp")[mod].mean()
            nd = naive.get(disease, np.nan) - naive.get("Control", np.nan)
            try:
                m = smf.ols(f"Q('{mod}') ~ C(grp) + sexM" + (" + age" if has_age else ""), cc).fit()
                key = [k for k in m.params.index if "grp" in k][0]
                beta, p = m.params[key], m.pvalues[key]
            except Exception:
                beta = p = np.nan
            try:
                _, pp_ = perm_p(cc[mod].values.astype(float), d_case, X0, B=2000)
            except Exception:
                pp_ = np.nan
            rows.append(dict(celltype=ct, module=mod, naive_delta=nd, beta=beta, p=p,
                             perm_p=pp_, n_dis=int((cc.grp == disease).sum()),
                             n_ctl=int((cc.grp == "Control").sum())))
    R = pd.DataFrame(rows)
    mreal = R[~R.module.eq("Composite")]
    R["fdr"] = np.nan
    R.loc[mreal.index, "fdr"] = multipletests(mreal["p"].fillna(1), method="fdr_bh")[1]
    R.to_csv(RES / f"{tag}_module_stats.csv", index=False)
    print(f"[{tag}] tests={len(mreal)} nominal p<0.05={int((mreal.p<0.05).sum())} FDR<0.10={int((R.fdr<0.10).sum())}")
    sig = R[(R.fdr < 0.10) & (~R.module.eq("Composite"))].sort_values("p")
    print("top FDR<0.10:")
    print(sig[["celltype", "module", "beta", "p", "perm_p", "fdr"]].head(12).round(3).to_string(index=False))

    # global donor-clustered
    grows = []
    DD = D[D.grp.isin(["Control", disease])].dropna(subset=["sexM"] + (["age"] if has_age else [])).copy()
    DD["grp"] = pd.Categorical(DD["grp"], ["Control", disease])
    for mod in modules:
        try:
            md = smf.ols(f"Q('{mod}') ~ C(grp) + sexM + C(celltype)" + (" + age" if has_age else ""),
                         DD).fit(cov_type="cluster", cov_kwds={"groups": DD["donor"]})
            key = [k for k in md.params.index if "grp" in k][0]
            grows.append(dict(module=mod, beta=md.params[key], p=md.pvalues[key]))
        except Exception:
            grows.append(dict(module=mod, beta=np.nan, p=np.nan))
    G = pd.DataFrame(grows); G["fdr"] = multipletests(G["p"].fillna(1), method="fdr_bh")[1]
    G.to_csv(RES / f"{tag}_module_global.csv", index=False)
    print("\nglobal (donor-clustered) top:")
    print(G.sort_values("p").head(8).round(4).to_string(index=False))

    # heatmap
    real = R[~R.module.eq("Composite")].copy()
    mat = real.pivot_table(index="celltype", columns="module", values="beta")
    fmat = real.pivot_table(index="celltype", columns="module", values="fdr").reindex_like(mat)
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]
    mat = mat[mat.var().sort_values(ascending=False).index]
    fmat = fmat.reindex(index=mat.index, columns=mat.columns)
    fig, ax = plt.subplots(figsize=(15, max(3.5, 0.5 * mat.shape[0] + 1.5)))
    vmax = np.nanpercentile(np.abs(mat.values), 95) if np.isfinite(mat.values).any() else 1
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=9)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            f = fmat.values[i, j]
            if f < 0.10:
                ax.text(j, i, "*" if f >= 0.05 else "**", ha="center", va="center", fontsize=9)
    adj = "age+sex adj" if has_age else "sex adj (age n/a)"
    ax.set_title(f"{disease} module aging, {adj} (β {disease}−Control)\n* FDR<0.10, ** FDR<0.05; red=aged, blue=younger")
    plt.colorbar(im, label="β (adjusted)", shrink=0.7); plt.tight_layout()
    fig.savefig(FIG / f"{tag}_module_heatmap.png", dpi=300); plt.close(fig)

    # summary
    lines = [
        f"# {disease} module-resolved aging ({tag}), {adj}",
        "",
        f"- pseudobulks: {len(meta)} ({disease} donors + Control); cell types tested: {R.celltype.nunique()}.",
        f"- cell-type × module tests: {len(mreal)}; nominal p<0.05: {int((mreal.p<0.05).sum())}; "
        f"**BH-FDR<0.10: {int((R.fdr<0.10).sum())}**.",
        "",
        "## Top accelerated module × cell-type effects (FDR<0.10)",
        (sig.head(10)[["celltype", "module", "beta", "p", "fdr"]].round(3).to_string(index=False)
         if len(sig) else "(none reached FDR<0.10)"),
        "",
        "## Global module effects (donor-clustered, FDR<0.10)",
        G[G.fdr < 0.10].sort_values("p")[["module", "beta", "p", "fdr"]].round(4).to_string(index=False)
        if (G.fdr < 0.10).any() else "(none)",
    ]
    (RES / f"{tag}_summary.md").write_text("\n".join(lines))
    print(f"\n[{tag}] saved: {tag}_module_stats.csv, {tag}_module_global.csv, "
          f"figures/{tag}_module_heatmap.png, {tag}_summary.md")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
