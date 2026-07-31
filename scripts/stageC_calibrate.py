"""Step 1-2: is myeloid non-replication driven by donor-count POWER or by a real
cohort/sampling difference? And are the 'replicating' lymphoid types even
significant at this n?

Uses cached donor-level instance module scores (mean-pool = the pseudobulk
approach). For each coarse cell type:

  (A) Cross-cohort observed r  : Adams<->Habermann IPF-vs-Control module-effect
      concordance, with a donor-level disease-permutation null (both cohorts
      shuffled independently) -> empirical p. Re-confirms which types (e.g. NK)
      are actually significant.

  (B) Within-Adams SPLIT-HALF ceiling : randomly partition Adams donors (disease-
      stratified) into two disjoint halves (~Habermann-scale n each), correlate
      their effect vectors. This is the BEST achievable reproducibility for the
      SAME biology/protocol at this n. If the cross-cohort r sits inside the
      split-half spread -> non-replication is a POWER (n) limit, not a cohort
      difference. If split-half is high but cross-cohort is far below -> a real
      cohort/sampling effect.

Outputs: results/mil/stageC_calibrate.csv, results/figures/stageC_calibrate.png
"""
import os
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import pearsonr
REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
from mil_score import MODULES
from stage2_replication import ADAMS_COARSE, HAB_COARSE, load, PRIOR
MIL = REPO / "results" / "mil"
FIG = REPO / "results" / "figures"
NPERM = 2000
NSPLIT = 500


def donor_table(df_ct):
    """donor x 23 module means, plus disease Series aligned to that index."""
    Dm = df_ct.groupby("Subject_Identity")[MODULES].mean()
    dis = df_ct.groupby("Subject_Identity").Disease_Identity.first().reindex(Dm.index)
    return Dm, dis


def effect(Dm, dis_vec):
    ipf = Dm[dis_vec == "IPF"].mean().values
    ctl = Dm[dis_vec == "Control"].mean().values
    return ipf - ctl


def conc(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return pearsonr(a[ok], b[ok])[0] if ok.sum() >= 3 else np.nan


def main():
    A = load("adams_bs60", ADAMS_COARSE)
    H = load("hab_bs60", HAB_COARSE)
    shared = [c for c in ["NK", "Tcell", "Macrophage", "Monocyte", "Ciliated"]
              if c in set(A.coarse) and c in set(H.coarse)]
    rng = np.random.default_rng(0)
    rows = []
    for ct in shared:
        DmA, disA = donor_table(A[A.coarse == ct])
        DmH, disH = donor_table(H[H.coarse == ct])
        r_obs = conc(effect(DmA, disA.values), effect(DmH, disH.values))

        # (A) cross-cohort disease-permutation null
        vA, vH = disA.values, disH.values
        null = np.empty(NPERM)
        for i in range(NPERM):
            null[i] = conc(effect(DmA, rng.permutation(vA)), effect(DmH, rng.permutation(vH)))
        null = null[np.isfinite(null)]
        p_cross = (1 + np.sum(null >= r_obs)) / (1 + len(null))

        # (B) within-Adams split-half ceiling (disease-stratified)
        idx = np.arange(len(disA)); ipf_i = idx[vA == "IPF"]; ctl_i = idx[vA == "Control"]
        sh = []
        for _ in range(NSPLIT):
            rng.shuffle(ipf_i); rng.shuffle(ctl_i)
            h1 = np.concatenate([ipf_i[:len(ipf_i) // 2], ctl_i[:len(ctl_i) // 2]])
            h2 = np.concatenate([ipf_i[len(ipf_i) // 2:], ctl_i[len(ctl_i) // 2:]])
            e1 = effect(DmA.iloc[h1], pd.Series(vA).iloc[h1].values)
            e2 = effect(DmA.iloc[h2], pd.Series(vA).iloc[h2].values)
            sh.append(conc(e1, e2))
        sh = np.array([x for x in sh if np.isfinite(x)])

        rows.append({"celltype": ct, "prior": PRIOR.get(ct, "?"),
                     "r_cross_obs": r_obs, "p_cross": p_cross,
                     "null_sd": null.std(),
                     "splithalf_med": np.median(sh),
                     "splithalf_lo": np.percentile(sh, 10),
                     "splithalf_hi": np.percentile(sh, 90),
                     "n_IPF_A": int((vA == "IPF").sum()), "n_Ctl_A": int((vA == "Control").sum()),
                     "n_IPF_H": int((vH == "IPF").sum()), "n_Ctl_H": int((vH == "Control").sum())})
        print(f"  {ct:11s} [{PRIOR.get(ct,'?'):9s}]  cross r={r_obs:+.2f} (p={p_cross:.3f}, null sd {null.std():.2f})   "
              f"Adams split-half r med={np.median(sh):+.2f} [{np.percentile(sh,10):+.2f},{np.percentile(sh,90):+.2f}]")

    R = pd.DataFrame(rows)
    R.to_csv(MIL / "stageC_calibrate.csv", index=False)

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    x = np.arange(len(R))
    ax.bar(x, R.splithalf_med, 0.55, color="#7fb3d5",
           yerr=[R.splithalf_med - R.splithalf_lo, R.splithalf_hi - R.splithalf_med],
           capsize=4, label="Adams split-half ceiling (same biology, ~Hab n)")
    ax.scatter(x, R.r_cross_obs, color="#c0392b", zorder=5, s=70,
               label="Adams<->Habermann cross-cohort (observed)")
    for xi, (_, r) in zip(x, R.iterrows()):
        ax.annotate(("p=%.2f" % r.p_cross), (xi, r.r_cross_obs),
                    textcoords="offset points", xytext=(8, 0), fontsize=8, color="#c0392b")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels([f"{c}\n({p})" for c, p in zip(R.celltype, R.prior)])
    ax.set_ylabel("module-effect concordance (Pearson r)")
    ax.set_title("Step 1-2: power (split-half ceiling) vs cross-cohort replication")
    ax.legend(fontsize=8)
    plt.tight_layout(); fig.savefig(FIG / "stageC_calibrate.png", dpi=130); plt.close(fig)
    print(f"\nsaved figure: {FIG/'stageC_calibrate.png'}")
    print(f"saved table:  {MIL/'stageC_calibrate.csv'}")


if __name__ == "__main__":
    main()
