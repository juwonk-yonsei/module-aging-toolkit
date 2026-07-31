"""Stage 2 (decisive): does ATTENTION pooling rescue cross-cohort replication
that MEAN pooling misses?

For each cohort (Adams GSE136831, Habermann GSE135893) and each shared coarse
cell type, we derive a donor-level 23-module IPF-vs-Control effect vector two ways:
  - mean-pool : donor score = plain mean over its instance sub-pseudobulks
                (= the existing pseudobulk approach).
  - attention : donor score = OUT-OF-FOLD attention-weighted mean of instances,
                where attention is trained (IPF vs Control) on OTHER donors.
Then we measure Adams<->Habermann concordance (Pearson r, sign %) of the effect
vectors, per cell type, for each pooling.

Cross-cohort concordance is the fair metric: two INDEPENDENTLY trained attention
models cannot produce correlated spurious effects, so any concordance gain is real.

Prior result (mean/pseudobulk): NK & T replicate; Macrophage/Monocyte/Ciliated
(myeloid/epithelial) do not. Success = attention lifts the failing types while
preserving the replicating ones.

Outputs: results/mil/stage2_replication.csv, results/figures/stage2_replication.png
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
from scipy.stats import pearsonr

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import mil_model as M
from mil_score import MODULES

MIL = REPO / "results" / "mil"
FIG = REPO / "results" / "figures"

# fine celltype -> coarse
ADAMS_COARSE = {"Macrophage": "Macrophage", "Macrophage_Alveolar": "Macrophage",
                "cMonocyte": "Monocyte", "ncMonocyte": "Monocyte", "NK": "NK",
                "T": "Tcell", "T_Cytotoxic": "Tcell", "Ciliated": "Ciliated"}
HAB_COARSE = {"Macrophages": "Macrophage", "Monocytes": "Monocyte", "NK Cells": "NK",
              "T Cells": "Tcell", "Ciliated": "Ciliated"}
# prior (pseudobulk) replication verdict, for annotation
PRIOR = {"NK": "replicated", "Tcell": "replicated",
         "Macrophage": "failed", "Monocyte": "failed", "Ciliated": "failed"}


def load(tag, coarse_map):
    df = pd.read_csv(MIL / f"inst_scores_{tag}.csv", index_col=0)
    df = df[df.Disease_Identity.isin(["IPF", "Control"])].copy()
    df["coarse"] = df.celltype.map(coarse_map)
    return df.dropna(subset=["coarse"])


import argparse


def mean_effect(df_ct):
    donor = df_ct.groupby(["Subject_Identity", "Disease_Identity"])[MODULES].mean().reset_index()
    ipf = donor[donor.Disease_Identity == "IPF"][MODULES].mean()
    ctl = donor[donor.Disease_Identity == "Control"][MODULES].mean()
    return (ipf - ctl).values, donor


def attention_effect(df_ct, seed=0):
    """Out-of-fold attention-weighted donor effect vector (23,)."""
    donors = df_ct.Subject_Identity.values
    uniq = pd.Index(sorted(df_ct.Subject_Identity.unique()))
    dlabel = df_ct.groupby("Subject_Identity").Disease_Identity.first().reindex(uniq)
    y = (dlabel == "IPF").astype(float).values
    n_ipf = int(y.sum()); n_ctl = int((1 - y).sum())
    if n_ipf < 4 or n_ctl < 4:
        return None
    Xraw = df_ct[MODULES].values.astype(np.float32)
    mu, sd = Xraw.mean(0), Xraw.std(0) + 1e-8
    Xstd = (Xraw - mu) / sd
    bag_rows = {d: np.where(donors == d)[0] for d in uniq}

    donor_vec = {}
    k = min(5, n_ipf, n_ctl)
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    for tr, va in skf.split(np.arange(len(uniq)), y):
        tr_d = uniq[tr]; va_d = uniq[va]
        bags_tr = [Xstd[bag_rows[d]] for d in tr_d]
        model, _ = M.train_eval(bags_tr, y[tr], bags_tr, y[tr],
                                in_dim=len(MODULES), pool="attention", seed=seed,
                                max_epochs=200, patience=25)
        for d in va_d:
            rows = bag_rows[d]
            aw = M.bag_attention(model, Xstd[rows])          # weights from std features
            w = aw / (aw.sum() + 1e-12)
            donor_vec[d] = (w[:, None] * Xraw[rows]).sum(0)   # weighted RAW module scores
    V = pd.DataFrame(donor_vec).T                             # donors x 23 modules
    lab = dlabel.reindex(V.index).values
    ipf = V.loc[lab == "IPF"].mean()
    ctl = V.loc[lab == "Control"].mean()
    return (ipf - ctl).values.astype(float)


def concordance(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan, np.nan
    r = pearsonr(a[ok], b[ok])[0]
    sign = np.mean(np.sign(a[ok]) == np.sign(b[ok]))
    return r, sign


def main(adams_tag="adams_bs60", hab_tag="hab_bs60", suffix=""):
    A = load(adams_tag, ADAMS_COARSE)
    H = load(hab_tag, HAB_COARSE)
    shared = [c for c in ["NK", "Tcell", "Macrophage", "Monocyte", "Ciliated"]
              if c in set(A.coarse) and c in set(H.coarse)]
    print("shared coarse cell types:", shared)

    rows = []
    for ct in shared:
        a_mean, _ = mean_effect(A[A.coarse == ct])
        h_mean, _ = mean_effect(H[H.coarse == ct])
        a_att = attention_effect(A[A.coarse == ct])
        h_att = attention_effect(H[H.coarse == ct])
        r_mean, s_mean = concordance(a_mean, h_mean)
        if a_att is not None and h_att is not None:
            r_att, s_att = concordance(a_att, h_att)
        else:
            r_att, s_att = np.nan, np.nan
        rows.append({"celltype": ct, "prior": PRIOR.get(ct, "?"),
                     "r_mean": r_mean, "sign_mean": s_mean,
                     "r_att": r_att, "sign_att": s_att,
                     "delta_r": r_att - r_mean})
        print(f"  {ct:11s} [{PRIOR.get(ct,'?'):9s}]  mean r={r_mean:+.2f} (sign {s_mean:.0%})   "
              f"attention r={r_att:+.2f} (sign {s_att:.0%})   Δr={r_att-r_mean:+.2f}")

    R = pd.DataFrame(rows)
    R.to_csv(MIL / f"stage2_replication{suffix}.csv", index=False)

    # ---- figure: paired r per cell type
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(R)); w = 0.36
    ax.bar(x - w / 2, R.r_mean, w, label="mean-pool (pseudobulk)", color="#2c3e50")
    ax.bar(x + w / 2, R.r_att, w, label="attention-pool", color="#c0392b")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels([f"{c}\n({p})" for c, p in zip(R.celltype, R.prior)])
    ax.set_ylabel("Adams <-> Habermann module-effect concordance (Pearson r)")
    ax.set_title("Stage 2: does attention pooling rescue cross-cohort replication?")
    ax.legend()
    plt.tight_layout(); fig.savefig(FIG / f"stage2_replication{suffix}.png", dpi=130); plt.close(fig)
    print(f"\nsaved figure: {FIG/('stage2_replication'+suffix+'.png')}")
    print(f"saved table:  {MIL/('stage2_replication'+suffix+'.csv')}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--adams_tag", default="adams_bs60")
    ap.add_argument("--hab_tag", default="hab_bs60")
    ap.add_argument("--suffix", default="")
    a = ap.parse_args()
    main(a.adams_tag, a.hab_tag, a.suffix)
