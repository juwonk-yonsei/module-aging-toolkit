"""Supplementary coverage + imputation-robustness analysis (reviewer request).

CORRECTED measurement: a clock feature is "detected" only if it survives the
preprocessing detection filter (count>=10 in >=20% of samples) AND the
human->mouse ortholog map, i.e. it has a real (non-NaN) value entering the
clock. Features that are absent are median-imputed by the clock's SimpleImputer.

Because that imputation uses the clock's FIXED training median (a per-gene
constant, identical for every test donor), each imputed feature adds the same
offset to every donor's module score and therefore CANCELS in the within-study
IPF-vs-control contrast on which all claims rest. We verify this empirically by
recomputing every cell-type x module effect from detected features only and
comparing to the full (imputed-included) effect.

Outputs:
  - results/module_coverage.csv              (per-module detection + robustness)
  - results/module_coverage_bycelltype.csv   (module x cell type detection)
  - manuscript/figures_submission/FigS1.png
  - prints summary numbers for the manuscript text
"""
import os
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import joblib
from scipy.stats import pearsonr

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
sys.path.insert(0, str(Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "python"))
import tage_prep as tp

BASE = REPO
MC = BASE / "data" / "module_clocks"
CM = BASE / "data" / "clock_models"
RES = BASE / "results"
OUT = BASE / "figures"

CB_BLUE = "#0072B2"; CB_VERM = "#D55E00"
MIN_CELLS = 50
HOT = "inflamm|innate immun|respirat|mitochond|oxphos"


def module_effect(pipe, X):
    """Return (effect_full, effect_detected) helpers via the clock's linear form.
    X: samples x clock-genes DataFrame (NaN = undetected -> median-imputed)."""
    imp = pipe.named_steps["simpleimputer"]
    sca = pipe.named_steps["standardscaler"]
    est = pipe.named_steps["elasticnetcv"]
    Xi = imp.transform(X.values)
    Xs = sca.transform(Xi)
    contrib = Xs * est.coef_[None, :]           # samples x genes
    detected = ~X.isna().all(axis=0).values     # per gene: has real values
    return contrib, detected


def main():
    pbdf = pd.read_pickle(RES / "pseudobulk_counts.pkl")
    gmeta = pd.read_csv(RES / "pseudobulk_meta.csv")
    gmeta = gmeta[gmeta["n_cells"] >= MIN_CELLS].copy()

    comp = joblib.load(CM / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    comp_genes = [str(g) for g in comp.feature_names_in_]

    clocks = {}
    for f in sorted(MC.glob("module_Mortality_*.pkl")):
        d = joblib.load(f)
        clocks[d["annotation"]] = ([str(g) for g in d["genes"]], d["pipeline"])

    print("reproducing per-cell-type preprocessing (real detection) ...")
    bycell, comp_det, robo = [], {}, []
    n_ct = 0
    for ct, sub in gmeta.groupby("Manuscript_Identity"):
        dis = sub["Disease_Identity"]
        if (dis == "Control").sum() < 3 or (dis == "IPF").sum() < 3:
            continue
        n_ct += 1
        cols = sub["group"].astype(str).tolist()
        m_ct = sub.set_index(sub["group"].astype(str))
        pp = tp.preprocess(pbdf[cols], m_ct, species="human", gene_mapping_type="Ensembl",
                           control_group_column="Disease_Identity", control_group_label="Control")
        Xall = pp["scaled_diff"]                      # samples x gene_list (NaN where undetected)
        lab = m_ct.loc[Xall.index, "Disease_Identity"].values
        ipf, ctl = (lab == "IPF"), (lab == "Control")

        comp_det[ct] = float(1 - Xall.reindex(columns=comp_genes).isna().all(axis=0).mean())

        for ann, (genes, pipe) in clocks.items():
            X = Xall.reindex(columns=genes)
            det_frac = float(1 - X.isna().all(axis=0).mean())
            contrib, detected = module_effect(pipe, X)
            score_full = contrib.sum(axis=1)
            score_det = contrib[:, detected].sum(axis=1)
            eff_full = score_full[ipf].mean() - score_full[ctl].mean()
            eff_det = score_det[ipf].mean() - score_det[ctl].mean()
            bycell.append({"celltype": ct, "annotation": ann, "n_genes": len(genes),
                           "detected_frac": det_frac})
            robo.append({"celltype": ct, "annotation": ann,
                         "eff_full": eff_full, "eff_detected": eff_det})

    print(f"  cell types analysed: {n_ct}")
    BC = pd.DataFrame(bycell)
    RB = pd.DataFrame(robo)
    BC.to_csv(RES / "module_coverage_bycelltype.csv", index=False)

    agg = BC.groupby("annotation").agg(
        n_genes=("n_genes", "first"),
        detected_mean=("detected_frac", "mean"),
        detected_min=("detected_frac", "min"),
        detected_max=("detected_frac", "max")).reset_index()
    agg["imputed_mean"] = 1 - agg["detected_mean"]
    perf = pd.read_csv(MC / "module_clock_performance.csv")
    perf = perf[perf.outcome == "Mortality"][["annotation", "LOFO_pearson"]]
    agg = agg.merge(perf, on="annotation", how="left").sort_values("detected_mean").reset_index(drop=True)
    agg.to_csv(RES / "module_coverage.csv", index=False)

    comp_mean = np.mean(list(comp_det.values()))
    comp_med = np.median(list(comp_det.values()))
    mod_med = agg.detected_mean.median()
    # robustness: full vs detected-only effect
    r_rob, _ = pearsonr(RB.eff_full, RB.eff_detected)
    sign_match = np.mean(np.sign(RB.eff_full) == np.sign(RB.eff_detected))
    rel_gap = np.abs(RB.eff_full - RB.eff_detected) / (np.abs(RB.eff_full) + 1e-9)

    print("\n=== per-module detection (mean across cell types) ===")
    for _, r in agg.iterrows():
        print(f"  {r.annotation[:40]:40s} n={r.n_genes:4d} det={r.detected_mean*100:5.1f}% "
              f"({r.detected_min*100:4.0f}-{r.detected_max*100:4.0f}%)  LOFO r={r.LOFO_pearson:.2f}")
    print(f"\ncomposite detected: mean {comp_mean*100:.1f}% / median {comp_med*100:.1f}% "
          f"(range {min(comp_det.values())*100:.0f}-{max(comp_det.values())*100:.0f}%)")
    print(f"module detected (median of per-module means): {mod_med*100:.1f}%")
    print(f"ROBUSTNESS full-vs-detected-only effect: Pearson r={r_rob:.4f}, "
          f"sign concordance={sign_match*100:.1f}%, median |rel gap|={np.median(rel_gap):.2e}")

    # ---------------- Supplementary Figure S1 ----------------
    def shortlab(s, n=32):
        s = s.replace("/ ", "/")
        return s if len(s) <= n else s[:n - 1] + "\u2026"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 7),
                                   gridspec_kw={"width_ratios": [1.3, 1]})
    order = agg.sort_values("detected_mean")
    y = np.arange(len(order))
    hot = order.annotation.str.contains(HOT, case=False, na=False).values
    colors = [CB_VERM if h else CB_BLUE for h in hot]
    dm = order.detected_mean.values * 100
    lo = order.detected_min.values * 100
    hi = order.detected_max.values * 100
    ax1.barh(y, dm, color=colors, zorder=2)
    ax1.errorbar(dm, y, xerr=[dm - lo, hi - dm], fmt="none", ecolor="#555",
                 elinewidth=0.8, capsize=2, zorder=3)
    ax1.axvline(mod_med * 100, color="k", ls="--", lw=0.9)
    ax1.set_yticks(y); ax1.set_yticklabels(
        [f"{shortlab(a)}  (n={int(g)})" for a, g in zip(order.annotation, order.n_genes)], fontsize=7)
    ax1.set_xlim(0, 100)
    ax1.set_xlabel("Features detected after ortholog mapping (%),\n"
                   "mean across 13 cell types (bar) and range (whiskers)", fontsize=9)
    ax1.set_title("Per-module coverage in human single-cell pseudobulks", fontsize=10)
    ax1.legend(handles=[Patch(color=CB_VERM, label="key predictive modules\n(inflammation, respiration)"),
                        Patch(color=CB_BLUE, label="other modules"),
                        plt.Line2D([0], [0], color="k", ls="--", label=f"median {mod_med*100:.0f}%")],
               fontsize=7.5, loc="lower right")
    ax1.text(-0.02, 1.03, "a", transform=ax1.transAxes, fontsize=15, fontweight="bold",
             va="bottom", ha="right")

    lim = float(np.nanmax(np.abs(np.r_[RB.eff_full, RB.eff_detected])) * 1.1)
    ax2.plot([-lim, lim], [-lim, lim], color="#999", lw=1, ls="--", zorder=1)
    ax2.scatter(RB.eff_full, RB.eff_detected, s=26, color=CB_BLUE, alpha=0.75,
                edgecolor="k", lw=0.3, zorder=3)
    ax2.set_xlim(-lim, lim); ax2.set_ylim(-lim, lim); ax2.set_aspect("equal")
    ax2.set_xlabel("IPF\u2212control module effect\n(all features; imputed included)", fontsize=9)
    ax2.set_ylabel("IPF\u2212control module effect\n(detected features only)", fontsize=9)
    ax2.set_title("Effects are carried by detected features", fontsize=10)
    ax2.text(0.03, 0.97,
             f"Pearson r = {r_rob:.3f}\nsign concordance = {sign_match*100:.0f}%\n"
             f"median-imputed features add a\nconstant offset that cancels\nin the contrast",
             transform=ax2.transAxes, fontsize=8, va="top",
             bbox=dict(boxstyle="round", fc="white", ec="#bbb"))
    ax2.text(-0.02, 1.03, "b", transform=ax2.transAxes, fontsize=15, fontweight="bold",
             va="bottom", ha="right")

    fig.tight_layout()
    fig.savefig(OUT / "FigS1.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nsaved {OUT/'FigS1.png'}, {RES/'module_coverage.csv'}, "
          f"{RES/'module_coverage_bycelltype.csv'}")

    # numbers block for the manuscript text
    infl = agg[agg.annotation.str.contains("inflamm|innate immun", case=False, na=False)]
    resp = agg[agg.annotation.str.contains("^respiration", case=False, na=False)]
    print("\n--- text numbers ---")
    print(f"composite median detection per cell type: {comp_med*100:.0f}%")
    print(f"module median detection: {mod_med*100:.0f}%")
    if len(infl):
        print(f"inflammation detected {infl.iloc[0].detected_mean*100:.0f}%")
    if len(resp):
        print(f"respiration detected {resp.iloc[0].detected_mean*100:.0f}%")
    print(f"robustness r={r_rob:.3f}, sign {sign_match*100:.0f}%")


if __name__ == "__main__":
    main()
