"""Stage 3: interpretability positive control + attention-meaningfulness check.

Trains the attention-MIL (IPF vs Control) on Adams MACROPHAGE donor bags
(out-of-fold), extracts per-instance attention weights, and asks:

  Q1 (meaningful?): are attention weights informative, i.e. far from uniform?
  Q2 (biology):    within IPF bags, do high-attention instances carry the known
                   mono-derived / profibrotic program (SPP1, FCN1, VCAN, CTHRC1,
                   INHBA, MERTK up) and lose the resident-alveolar program
                   (FABP4, MARCO, MRC1)?  Contrast with Control bags.

A positive Q2 means attention localizes the biologically expected disease
sub-state even though it did not rescue cross-cohort replication (Stage 2) --
i.e. attention is an interpretive tool, and the myeloid non-replication is a
genuine cross-sampling/biology gap, not a broken model.

Outputs: results/mil/stage3_interpret.csv, results/figures/stage3_interpret.png
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

MIL = REPO / "results" / "mil"
FIG = REPO / "results" / "figures"

MONO_DERIVED = {"SPP1": "ENSG00000118785", "FCN1": "ENSG00000085265",
                "VCAN": "ENSG00000038427", "CTHRC1": "ENSG00000164932",
                "INHBA": "ENSG00000122641", "MERTK": "ENSG00000153208"}
RESIDENT = {"FABP4": "ENSG00000170323", "MARCO": "ENSG00000019169",
            "MRC1": "ENSG00000260314"}
FINE = ["Macrophage", "Macrophage_Alveolar"]


def load():
    sc = pd.read_csv(MIL / "inst_scores_adams_bs60.csv", index_col=0)
    sc = sc[sc.celltype.isin(FINE) & sc.Disease_Identity.isin(["IPF", "Control"])].copy()
    sc.index = sc.index.astype(str)
    counts = pd.read_pickle(MIL / "inst_counts_adams_bs60.pkl")
    counts.columns = counts.columns.map(str)
    counts = counts[sc.index]                       # genes x instances (aligned)
    lib = counts.sum(0).values + 1e-9
    marker_expr = {}
    for name, ens in {**MONO_DERIVED, **RESIDENT}.items():
        if ens in counts.index:
            cpm = counts.loc[ens].values / lib * 1e6
            marker_expr[name] = np.log1p(cpm)
    return sc, pd.DataFrame(marker_expr, index=sc.index)


def oof_attention(sc, seed=0):
    donors = sc.Subject_Identity.values
    uniq = pd.Index(sorted(sc.Subject_Identity.unique()))
    dlabel = sc.groupby("Subject_Identity").Disease_Identity.first().reindex(uniq)
    y = (dlabel == "IPF").astype(float).values
    X = sc[MODULES].values.astype(np.float32)
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    rows_of = {d: np.where(donors == d)[0] for d in uniq}
    w = np.full(len(sc), np.nan)
    k = min(5, int(y.sum()), int((1 - y).sum()))
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    for tr, va in skf.split(np.arange(len(uniq)), y):
        bags_tr = [X[rows_of[d]] for d in uniq[tr]]
        model, _ = M.train_eval(bags_tr, y[tr], bags_tr, y[tr], in_dim=len(MODULES),
                                pool="attention", seed=seed, max_epochs=200, patience=25)
        for d in uniq[va]:
            r = rows_of[d]
            aw = M.bag_attention(model, X[r])
            w[r] = aw * len(r)          # normalize: 1.0 == uniform weight
    return w, dlabel.reindex(uniq)


def main():
    sc, mk = load()
    print(f"macrophage instances={len(sc)}  markers={list(mk.columns)}")
    w, _ = oof_attention(sc)
    sc = sc.assign(att=w)

    # Q1: attention non-uniformity (relative weight; 1.0 = uniform)
    disp = np.nanstd(w)
    frac_hi = np.mean(w > 1.5)
    print(f"\n[Q1 meaningfulness] relative attention weight: std={disp:.2f}, "
          f"fraction >1.5x uniform = {frac_hi:.1%}  (0 => uniform/noise)")

    # Q2: correlation of attention with markers, per disease
    rows = []
    for grp in ["IPF", "Control"]:
        m = sc.Disease_Identity.values == grp
        for name in mk.columns:
            rho, p = spearmanr(sc.att.values[m], mk[name].values[m])
            klass = "mono-derived" if name in MONO_DERIVED else "resident"
            rows.append({"group": grp, "marker": name, "class": klass, "rho": rho, "p": p})
    R = pd.DataFrame(rows)
    R.to_csv(MIL / "stage3_interpret.csv", index=False)
    piv = R.pivot_table(index=["class", "marker"], columns="group", values="rho")
    print("\n[Q2 biology] Spearman(attention weight, marker expression):")
    print(piv.round(3).to_string())

    ipf_mono = R[(R.group == "IPF") & (R["class"] == "mono-derived")].rho.mean()
    ipf_res = R[(R.group == "IPF") & (R["class"] == "resident")].rho.mean()
    print(f"\n  IPF: mean rho mono-derived={ipf_mono:+.3f}  resident={ipf_res:+.3f}  "
          f"separation={ipf_mono - ipf_res:+.3f}")

    # ---- figure
    fig, ax = plt.subplots(figsize=(9, 5.5))
    order = list(MONO_DERIVED) + list(RESIDENT)
    x = np.arange(len(order)); w2 = 0.36
    for i, grp in enumerate(["IPF", "Control"]):
        vals = [R[(R.group == grp) & (R.marker == mm)].rho.values[0] for mm in order]
        ax.bar(x + (i - 0.5) * w2, vals, w2, label=grp,
               color="#c0392b" if grp == "IPF" else "#7f8c8d")
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(len(MONO_DERIVED) - 0.5, color="k", ls="--", lw=0.8)
    ax.text(len(MONO_DERIVED) / 2 - 0.5, ax.get_ylim()[1] * 0.9, "mono-derived / profibrotic",
            ha="center", fontsize=9)
    ax.text(len(MONO_DERIVED) + len(RESIDENT) / 2 - 0.5, ax.get_ylim()[1] * 0.9,
            "resident alveolar", ha="center", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel("Spearman(attention weight, marker expr)")
    ax.set_title("Stage 3: attention localizes the mono-derived/profibrotic\nmacrophage sub-state in IPF")
    ax.legend()
    plt.tight_layout(); fig.savefig(FIG / "stage3_interpret.png", dpi=130); plt.close(fig)
    print(f"\nsaved figure: {FIG/'stage3_interpret.png'}")
    print(f"saved table:  {MIL/'stage3_interpret.csv'}")


if __name__ == "__main__":
    main()
