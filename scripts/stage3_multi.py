"""Stage 3 (extended): does attention localize a known disease-associated
sub-state in MULTIPLE cell types (not just macrophage)?

For each coarse cell type we define a marker axis A (disease-expanded / effector
pole) vs B (homeostatic pole), train attention-MIL (IPF vs Control, out-of-fold),
and correlate per-instance attention weight with marker expression, per disease.
Positive control = in IPF, attention up-weights the A pole and down-weights B,
more strongly than in Control (separation = rho_A - rho_B).

  Macrophage : A mono-derived/profibrotic  vs  B resident-alveolar
  Monocyte   : A classical/inflammatory     vs  B non-classical/patrolling
  Tcell      : A cytotoxic/effector         vs  B naive/memory
  NK         : A cytotoxic                   vs  B regulatory/immature

Outputs: results/mil/stage3_multi.csv, results/mil/stage3_multi_summary.csv,
         results/figures/stage3_multi.png
"""
import os
os.environ.setdefault("CUDA_DEVICE_ORDER", "PCI_BUS_ID")
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "2")
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from scipy.stats import spearmanr

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import mil_model as M
from mil_score import MODULES

IPFD = REPO / "data" / "ipf"
MIL = REPO / "results" / "mil"
FIG = REPO / "results" / "figures"

ADAMS_COARSE = {"Macrophage": "Macrophage", "Macrophage_Alveolar": "Macrophage",
                "cMonocyte": "Monocyte", "ncMonocyte": "Monocyte",
                "NK": "NK", "T": "Tcell", "T_Cytotoxic": "Tcell"}

CONFIG = {
    "Macrophage": {"A_label": "mono-derived/profibrotic",
                   "A": ["SPP1", "FCN1", "VCAN", "CTHRC1", "INHBA", "MERTK", "SPARC"],
                   "B_label": "resident-alveolar", "B": ["FABP4", "MARCO", "MRC1", "PPARG"]},
    "Monocyte": {"A_label": "classical/inflammatory",
                 "A": ["S100A8", "S100A9", "S100A12", "FCN1", "VCAN", "CD14"],
                 "B_label": "non-classical/patrolling", "B": ["FCGR3A", "CX3CR1", "LST1"]},
    "Tcell": {"A_label": "cytotoxic/effector",
              "A": ["GZMB", "GZMK", "GNLY", "NKG7", "PRF1", "CCL5"],
              "B_label": "naive/memory", "B": ["CCR7", "SELL", "TCF7", "IL7R", "LEF1"]},
    "NK": {"A_label": "cytotoxic",
           "A": ["GNLY", "GZMB", "PRF1", "NKG7", "FGFBP2", "FCGR3A"],
           "B_label": "regulatory/immature", "B": ["XCL1", "GZMK", "SELL", "IL7R"]},
}


def sym2ens():
    g = pd.read_csv(IPFD / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    g.columns = [c.strip('"') for c in g.columns]
    for c in g.select_dtypes("object"):
        g[c] = g[c].str.strip('"')
    return dict(zip(g["HGNC_EnsemblAlt_GeneID"], g["Ensembl_GeneID"]))


S2E = sym2ens()


def load_ct(coarse, counts_all, sc_all):
    fine = [f for f, c in ADAMS_COARSE.items() if c == coarse]
    sc = sc_all[sc_all.celltype.isin(fine) & sc_all.Disease_Identity.isin(["IPF", "Control"])].copy()
    if not len(sc):
        return None, None
    cnt = counts_all[sc.index]
    lib = cnt.sum(0).values + 1e-9
    cfg = CONFIG[coarse]
    expr = {}
    for name in cfg["A"] + cfg["B"]:
        ens = S2E.get(name)
        if ens in cnt.index:
            expr[name] = np.log1p(cnt.loc[ens].values / lib * 1e6)
    return sc, pd.DataFrame(expr, index=sc.index)


def oof_attention(sc, seed=0):
    donors = sc.Subject_Identity.values
    uniq = pd.Index(sorted(sc.Subject_Identity.unique()))
    dlabel = sc.groupby("Subject_Identity").Disease_Identity.first().reindex(uniq)
    y = (dlabel == "IPF").astype(float).values
    if y.sum() < 4 or (1 - y).sum() < 4:
        return None
    X = sc[MODULES].values.astype(np.float32)
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    rows_of = {d: np.where(donors == d)[0] for d in uniq}
    w = np.full(len(sc), np.nan)
    k = min(5, int(y.sum()), int((1 - y).sum()))
    for tr, va in StratifiedKFold(n_splits=k, shuffle=True, random_state=seed).split(np.arange(len(uniq)), y):
        bags_tr = [X[rows_of[d]] for d in uniq[tr]]
        model, _ = M.train_eval(bags_tr, y[tr], bags_tr, y[tr], in_dim=len(MODULES),
                                pool="attention", seed=seed, max_epochs=200, patience=25)
        for d in uniq[va]:
            r = rows_of[d]
            w[r] = M.bag_attention(model, X[r]) * len(r)
    return w


def zscore(v):
    return (v - np.nanmean(v)) / (np.nanstd(v) + 1e-9)


def main():
    sc_all = pd.read_csv(MIL / "inst_scores_adams_bs60.csv", index_col=0)
    sc_all.index = sc_all.index.astype(str)
    counts_all = pd.read_pickle(MIL / "inst_counts_adams_bs60.pkl")
    counts_all.columns = counts_all.columns.map(str)

    per_marker, summary = [], []
    panels = {}
    for coarse, cfg in CONFIG.items():
        sc, expr = load_ct(coarse, counts_all, sc_all)
        if sc is None or expr is None or expr.shape[1] < 3:
            print(f"skip {coarse}"); continue
        w = oof_attention(sc)
        if w is None:
            print(f"skip {coarse}: too few donors"); continue
        sc = sc.assign(att=w)
        disp = np.nanstd(w); frac_hi = np.mean(w > 1.5)
        rec = {"celltype": coarse}
        for grp in ["IPF", "Control"]:
            m = sc.Disease_Identity.values == grp
            aw = sc.att.values[m]
            rA, rB = [], []
            for name in expr.columns:
                rho = spearmanr(aw, expr[name].values[m])[0]
                per_marker.append({"celltype": coarse, "group": grp, "marker": name,
                                   "pole": "A" if name in cfg["A"] else "B", "rho": rho})
                (rA if name in cfg["A"] else rB).append(rho)
            rec[f"{grp}_A"] = np.nanmean(rA); rec[f"{grp}_B"] = np.nanmean(rB)
            rec[f"{grp}_sep"] = np.nanmean(rA) - np.nanmean(rB)
        rec["att_std"] = disp; rec["frac_hi"] = frac_hi
        summary.append(rec)
        panels[coarse] = (sc, expr, cfg)
        print(f"{coarse:11s} attn std={disp:.2f} hi={frac_hi:.0%} | "
              f"IPF sep={rec['IPF_sep']:+.3f} (A {rec['IPF_A']:+.2f}/B {rec['IPF_B']:+.2f})  "
              f"Control sep={rec['Control_sep']:+.3f}")

    PM = pd.DataFrame(per_marker); PM.to_csv(MIL / "stage3_multi.csv", index=False)
    SM = pd.DataFrame(summary); SM.to_csv(MIL / "stage3_multi_summary.csv", index=False)

    # ---- figure: one panel per cell type (per-marker rho, IPF vs Control)
    n = len(panels)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    axes = axes.ravel()
    for ax, (coarse, (sc, expr, cfg)) in zip(axes, panels.items()):
        order = [m for m in cfg["A"] if m in expr.columns] + [m for m in cfg["B"] if m in expr.columns]
        x = np.arange(len(order)); wdt = 0.38
        for i, grp in enumerate(["IPF", "Control"]):
            vals = [PM[(PM.celltype == coarse) & (PM.group == grp) & (PM.marker == mm)].rho.values[0]
                    for mm in order]
            ax.bar(x + (i - 0.5) * wdt, vals, wdt, label=grp,
                   color="#c0392b" if grp == "IPF" else "#7f8c8d")
        ax.axhline(0, color="k", lw=0.6)
        nA = sum(m in expr.columns for m in cfg["A"])
        ax.axvline(nA - 0.5, color="k", ls="--", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8)
        ax.set_title(f"{coarse}:  A={cfg['A_label']}  |  B={cfg['B_label']}", fontsize=9)
        ax.set_ylabel("Spearman(attn, marker)")
        ax.legend(fontsize=8)
    for ax in axes[len(panels):]:
        ax.axis("off")
    fig.suptitle("Stage 3 (multi cell type): does attention localize the disease-expanded sub-state?",
                 fontsize=12)
    plt.tight_layout(); fig.savefig(FIG / "stage3_multi.png", dpi=130); plt.close(fig)
    print(f"\nsaved figure: {FIG/'stage3_multi.png'}")
    print(f"saved tables: {MIL/'stage3_multi.csv'}, {MIL/'stage3_multi_summary.csv'}")


if __name__ == "__main__":
    main()
