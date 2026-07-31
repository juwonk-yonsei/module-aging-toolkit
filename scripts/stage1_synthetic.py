"""Stage 0 (heterogeneity precondition) + Stage 1 (semi-synthetic ground-truth).

Uses REAL macrophage instance module-scores as background (real scRNA
heterogeneity/noise preserved), then injects a KNOWN signal and asks whether
attention-pooling beats mean-pooling, and when.

Two regimes (this distinction is the scientific point):
  * DIFFUSE   : disease effect added to a RANDOM, unidentifiable subset of
                instances -> attention has no feature to locate them by.
                Expectation: no attention advantage (fair negative control).
  * SUBSTATE  : a minority sub-state S is made identifiable (shifted on marker
                modules in BOTH case & control); the disease effect is applied
                ONLY to S instances of case bags. Mean-pool dilutes the effect by
                the sub-state fraction p_S; attention can focus on S.
                Expectation: attention advantage that GROWS as S gets rarer.

Ground-truth for attention recovery = the signal-bearing instances.

Outputs: results/mil/stage1_results.csv, results/figures/stage1_synthetic.png
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
from sklearn.metrics import roc_auc_score

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import mil_model as M
from mil_score import MODULES

MIL = REPO / "results" / "mil"
FIG = REPO / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

CELLTYPE = "Macrophage"
DELTA = 2.0            # disease effect size within signal-bearing instances (SD)
MARKER_SHIFT = 2.5     # how identifiable the sub-state is (SD on marker modules)
PS_GRID = [0.05, 0.1, 0.2, 0.35, 0.5]   # sub-state rarity
N_SEEDS = 8
N_FOLDS = 5
MIN_INST = 6           # drop donors with too few instances (unstable bags)


def load_macrophage():
    df = pd.read_csv(MIL / "inst_scores_adams_bs60.csv", index_col=0)
    df = df[df.celltype == CELLTYPE].copy()
    X = df[MODULES].values.astype(np.float32)
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    donors = df["Subject_Identity"].values
    return df, X, donors


def stage0(X, donors):
    tmp = pd.DataFrame(X, columns=MODULES); tmp["donor"] = donors
    within = tmp.groupby("donor")[MODULES].var().mean().mean()
    overall = tmp[MODULES].var().mean().mean()
    n_per = tmp.groupby("donor").size()
    print("=== Stage 0: macrophage instance heterogeneity ===")
    print(f"donors={tmp.donor.nunique()}  instances={len(tmp)}  "
          f"median instances/donor={int(n_per.median())} (min {n_per.min()}, max {n_per.max()})")
    print(f"mean within-donor variance / overall variance = "
          f"{within:.3f} / {overall:.3f} = {within/overall:.2%}")
    print("  -> substantial within-donor spread: attention has structure to exploit\n")


def make_bags(X, donors):
    bag_idx, bag_donor = [], []
    for d, g in pd.Series(range(len(donors))).groupby(donors, sort=True):
        if len(g) >= MIN_INST:
            bag_idx.append(g.values); bag_donor.append(d)
    return bag_idx, bag_donor


def build_condition(X, bag_idx, regime, p_s, delta, seed):
    rng = np.random.default_rng(seed)
    n_bags = len(bag_idx)
    y = np.zeros(n_bags, dtype=np.float32)
    y[rng.permutation(n_bags)[: n_bags // 2]] = 1.0

    D = X.shape[1]
    marker = np.array([2, 7, 12])          # fixed marker modules (identify S)
    effect = 17                            # disease-effect module (distinct)
    bags, gt = [], []
    for bi, idx in enumerate(bag_idx):
        b = X[idx].copy()
        n = len(idx)
        if regime == "substate":
            is_s = rng.random(n) < p_s
            b[is_s, marker[0]] += MARKER_SHIFT
            b[is_s, marker[1]] += MARKER_SHIFT
            b[is_s, marker[2]] += MARKER_SHIFT
            if y[bi] == 1.0:
                b[is_s, effect] += delta
            flags = is_s & (y[bi] == 1.0)
        else:  # diffuse: effect on a random unidentifiable subset of case bags
            flags = np.zeros(n, dtype=bool)
            if y[bi] == 1.0:
                k = max(1, int(round(p_s * n)))
                sel = rng.choice(n, size=k, replace=False)
                b[sel, effect] += delta
                flags[sel] = True
        bags.append(b.astype(np.float32)); gt.append(flags)
    return bags, y, gt


def evaluate(bags, y, gt, in_dim, seed, want_recovery):
    out = {}
    for pool in ["attention", "mean"]:
        skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
        va_probs = np.full(len(y), np.nan)
        aurocs = []
        for tr, va in skf.split(np.arange(len(y)), y):
            model, vp = M.train_eval([bags[i] for i in tr], y[tr],
                                     [bags[i] for i in va], y[va],
                                     in_dim=in_dim, pool=pool, seed=seed)
            va_probs[va] = vp
            if pool == "attention" and want_recovery:
                for i in va:
                    if gt[i].sum() > 0 and (~gt[i]).sum() > 0:
                        aw = M.bag_attention(model, bags[i])
                        aurocs.append(roc_auc_score(gt[i], aw))
        out[pool] = {"bag_auc": roc_auc_score(y, va_probs),
                     "att_auroc": float(np.mean(aurocs)) if aurocs else np.nan}
    return out


def main():
    df, X, donors = load_macrophage()
    stage0(X, donors)
    bag_idx, bag_donor = make_bags(X, donors)
    print(f"=== Stage 1: semi-synthetic (macrophage, {len(bag_idx)} donor-bags "
          f">= {MIN_INST} instances, delta={DELTA} SD) ===")

    rows = []
    # main sweep: SUBSTATE regime across sub-state rarity
    print("\n[SUBSTATE regime: disease effect confined to an identifiable minority]")
    for p_s in PS_GRID:
        for seed in range(N_SEEDS):
            res = evaluate(*build_condition(X, bag_idx, "substate", p_s, DELTA, seed),
                           in_dim=X.shape[1], seed=seed, want_recovery=True)
            for pool in ["attention", "mean"]:
                rows.append({"regime": "substate", "p_s": p_s, "seed": seed, "pool": pool,
                             **res[pool]})
        cur = pd.DataFrame([r for r in rows if r["regime"] == "substate" and r["p_s"] == p_s])
        a = cur[cur.pool == "attention"]["bag_auc"]; m = cur[cur.pool == "mean"]["bag_auc"]
        awr = cur[cur.pool == "attention"]["att_auroc"]
        print(f"  p_S={p_s:>4}: AUC attention={a.mean():.3f}+-{a.std():.3f}  "
              f"mean={m.mean():.3f}+-{m.std():.3f}  delta={a.mean()-m.mean():+.3f}  "
              f"| attn-recovery AUROC={awr.mean():.3f}")

    # control: DIFFUSE regime (unidentifiable subset) at a matched fraction
    print("\n[DIFFUSE regime (control): effect on random unidentifiable subset]")
    for p_s in [0.1, 0.2, 0.35]:
        for seed in range(N_SEEDS):
            res = evaluate(*build_condition(X, bag_idx, "diffuse", p_s, DELTA, seed),
                           in_dim=X.shape[1], seed=seed, want_recovery=True)
            for pool in ["attention", "mean"]:
                rows.append({"regime": "diffuse", "p_s": p_s, "seed": seed, "pool": pool,
                             **res[pool]})
        cur = pd.DataFrame([r for r in rows if r["regime"] == "diffuse" and r["p_s"] == p_s])
        a = cur[cur.pool == "attention"]["bag_auc"]; m = cur[cur.pool == "mean"]["bag_auc"]
        print(f"  f={p_s:>4}: AUC attention={a.mean():.3f}+-{a.std():.3f}  "
              f"mean={m.mean():.3f}+-{m.std():.3f}  delta={a.mean()-m.mean():+.3f}")

    R = pd.DataFrame(rows)
    R.to_csv(MIL / "stage1_results.csv", index=False)

    # ---- figure
    sub = R[R.regime == "substate"]
    g = sub.groupby(["p_s", "pool"])["bag_auc"].agg(["mean", "std"]).reset_index()
    gr = sub[sub.pool == "attention"].groupby("p_s")["att_auroc"].agg(["mean", "std"]).reset_index()
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for pool, c in [("attention", "#c0392b"), ("mean", "#2c3e50")]:
        s = g[g.pool == pool]
        ax[0].errorbar(s["p_s"], s["mean"], yerr=s["std"], marker="o", capsize=3,
                       color=c, label=f"{pool}-pool")
    ax[0].axhline(0.5, ls="--", color="grey", lw=0.8)
    ax[0].set_xlabel("sub-state fraction p_S (rarer ->)"); ax[0].invert_xaxis()
    ax[0].set_ylabel("held-out bag-classification AUC")
    ax[0].set_title("Stage 1 SUBSTATE: attention wins when\nsignal is a rare identifiable sub-state")
    ax[0].legend()
    ax[1].errorbar(gr["p_s"], gr["mean"], yerr=gr["std"], marker="s", capsize=3, color="#c0392b")
    ax[1].axhline(0.5, ls="--", color="grey", lw=0.8); ax[1].invert_xaxis()
    ax[1].set_xlabel("sub-state fraction p_S"); ax[1].set_ylabel("attention-weight AUROC (recover sub-state)")
    ax[1].set_title("Attention localizes the signal-bearing sub-state")
    plt.tight_layout(); fig.savefig(FIG / "stage1_synthetic.png", dpi=130); plt.close(fig)
    print(f"\nsaved figure: {FIG/'stage1_synthetic.png'}")
    print(f"saved table:  {MIL/'stage1_results.csv'}")


if __name__ == "__main__":
    main()
