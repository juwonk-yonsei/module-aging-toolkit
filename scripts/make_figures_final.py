"""Regenerate ALL manuscript figures as submission-ready, multi-panel figures.

Fixes applied (per figure review):
  - remove titles baked into images (GeroScience: no titles inside illustrations)
  - assemble multi-panel figures with (a)(b)(c) panel letters
  - complete Fig4 (a per-celltype, b split-half ceiling, c bulk recovery) and add Fig5 (attention)
  - de-clutter / un-truncate labels (Fig1b, Fig3a, Fig3c)
  - colorblind-safe palette for Fig1c (Okabe-Ito), avoid blue/purple
  - upgrade Fig1a pipeline schematic to depict the novel analysis/validation layer

All panels are re-plotted from saved result tables (no heavy recompute) except Fig1c,
which scores the 23 module clocks on the bundled tAge example data (fast).

Outputs -> manuscript/figures_submission/Fig1.png ... Fig5.png  (300 dpi)
"""
import os
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch, Circle
from scipy.stats import pearsonr

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
BASE = REPO
RES = BASE / "results"
MC = BASE / "data" / "module_clocks"
MIL = RES / "mil"
OUT = BASE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Okabe-Ito colorblind-safe palette
CB_BLUE = "#0072B2"; CB_ORANGE = "#E69F00"; CB_RED = "#D55E00"
CB_GREY = "#999999"; CB_SKY = "#56B4E9"; CB_VERM = "#c0392b"

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "savefig.dpi": 300})


def panel_letter(ax, letter, dx=-0.12, dy=1.02, size=15):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="right")


def short(s, n=26):
    """Truncate a long annotation at a word/'/' boundary (no mid-word cut)."""
    s = str(s)
    if len(s) <= n:
        return s
    cut = s[:n]
    for sep in ["/ ", " ", "/"]:
        j = cut.rfind(sep)
        if j >= n * 0.5:
            return cut[:j].rstrip(" /,")
    return cut.rstrip(" /,")


def draw_pipeline(ax):
    """Redesigned Fig1a schematic (used by both Fig1 panel a and standalone Fig1a).

    Aligned two-tier grid: grey 'core pipeline (reused)' band on top, amber
    'this work' band below. The two are joined by an explicit hand-off node
    (the pipeline's actual output) and a single labeled connector, so the
    relationship reads as: pipeline output -> analyzed & validated by this work.
    Numbered stage badges; colorblind-safe (blue/amber).
    """
    ax.axis("off"); ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(0, 16); ax.set_ylim(0, 10.6)

    C_EDGE = "#37597b"; A_EDGE = "#c8922e"
    # ---- containers (top: reused core pipeline; bottom: this work)
    ax.add_patch(FancyBboxPatch((0.15, 6.2), 15.7, 4.2, boxstyle="round,pad=0.02",
                 fc="#f6f8fa", ec="#d7dce2", lw=1.3))
    ax.add_patch(FancyBboxPatch((0.15, 1.0), 15.7, 3.0, boxstyle="round,pad=0.02",
                 fc="#fdf4e3", ec=A_EDGE, lw=1.4))
    ax.text(0.5, 10.02, "CORE PIPELINE  ·  retrained from public data (reused)",
            ha="left", va="center", fontsize=9.5, fontweight="bold", color="#66707d")
    # 'this work' pill (top-left of amber band)
    ax.add_patch(FancyBboxPatch((0.5, 3.46), 1.95, 0.5, boxstyle="round,pad=0.03",
                 fc="#f0b64a", ec=A_EDGE, lw=1.0))
    ax.text(1.475, 3.71, "THIS WORK", ha="center", va="center", fontsize=8.6,
            fontweight="bold", color="#5a3d08")

    # ---- top row: 5 stage boxes on a common grid
    core = [
        ("Public rodent\nmeta-dataset", "4,539 samples · 96 datasets", "#eaf3fb"),
        ("Retrain 23\nmodule clocks", "ElasticNet · group-OOF CV", "#dcebf7"),
        ("Human scRNA-seq\npseudobulk", "IPF / COPD · GSE136831", "#cfe3f3"),
        ("Cross-species\nmodule scoring", "per cell type", "#c2dcef"),
        ("Age/sex-adj\nstatistics", "perm · FDR · COPD\n+ replication · GSE135893", "#b5d4eb"),
    ]
    w, h, gap, x0, ty = 2.5, 2.15, 0.675, 0.4, 6.85
    centers = [x0 + i * (w + gap) + w / 2 for i in range(5)]
    for i, (title, detail, c) in enumerate(core):
        x = x0 + i * (w + gap); cx = x + w / 2
        ax.add_patch(FancyBboxPatch((x, ty), w, h, boxstyle="round,pad=0.04",
                     fc=c, ec=C_EDGE, lw=1.2))
        ax.text(cx, ty + h - 0.62, title, ha="center", va="center", fontsize=9.2,
                fontweight="bold", color="#16232f")
        ax.text(cx, ty + h - 1.42, detail, ha="center", va="center", fontsize=7.4,
                color="#33465a")
        ax.add_patch(Circle((x + 0.30, ty + h), 0.26, fc="#2c5f88",
                     ec="white", lw=1.2, zorder=5))
        ax.text(x + 0.30, ty + h, str(i + 1), ha="center", va="center",
                fontsize=8.6, color="white", fontweight="bold", zorder=6)
    # arrows between top boxes
    for i in range(4):
        xr = x0 + i * (w + gap) + w
        ax.add_patch(FancyArrowPatch((xr + 0.08, ty + h / 2), (xr + gap - 0.08, ty + h / 2),
                     arrowstyle="-|>", mutation_scale=15, color=C_EDGE, lw=1.7))

    # ---- shared hand-off node: the pipeline's actual output, consumed by 'this work'
    jx0, jx1, jy0, jy1 = 1.4, 14.6, 5.05, 6.0
    for cx in centers:  # every stage contributes to the shared output
        ax.plot([cx, cx], [ty, jy1], color="#9aa3ad", lw=1.2, zorder=1)
    ax.add_patch(FancyBboxPatch((jx0, jy0), jx1 - jx0, jy1 - jy0,
                 boxstyle="round,pad=0.03", fc="#eef1f4", ec="#8a96a3", lw=1.3, zorder=2))
    ax.text(jx0 + 0.35, jy1 - 0.24, "PIPELINE OUTPUT", ha="left", va="center",
            fontsize=7, fontweight="bold", color="#7a828c", zorder=3)
    ax.text((jx0 + jx1) / 2, jy0 + 0.34,
            "Per-cell-type module aging scores  +  age/sex-adjusted statistics",
            ha="center", va="center", fontsize=9, fontweight="bold",
            color="#33414f", zorder=3)
    # single labeled connector: output is analyzed & validated by THIS WORK
    ax.add_patch(FancyArrowPatch((8.0, jy0 - 0.03), (8.0, 4.0), arrowstyle="-|>",
                 mutation_scale=24, color="#a9760f", lw=3.0, zorder=2))
    ax.text(8.5, 4.52, "analyzed &\nvalidated by", ha="left", va="center",
            fontsize=8.2, style="italic", color="#8a5e0c")
    bn_w, bn_h, bn_gap, bn_y = 4.4, 1.6, 1.0, 1.3

    # ---- bottom row: 3 novel components
    novel = [
        ("Composite-masking\ndiagnostic", "opposing modules cancel"),
        ("Attribution ladder", "technical → power →\ncomposition → resolution"),
        ("Attention sub-state\nlocalizer", "interpretable MIL"),
    ]
    for i, (title, detail) in enumerate(novel):
        x = x0 + i * (bn_w + bn_gap); cx = x + bn_w / 2
        ax.add_patch(FancyBboxPatch((x, bn_y), bn_w, bn_h, boxstyle="round,pad=0.04",
                     fc="#f9e6c0", ec=A_EDGE, lw=1.3))
        ax.text(cx, bn_y + bn_h - 0.5, title, ha="center", va="center", fontsize=9.4,
                fontweight="bold", color="#4a3410")
        ax.text(cx, bn_y + bn_h - 1.12, detail, ha="center", va="center", fontsize=7.8,
                color="#6b5417")


# ============================ FIGURE 1 ============================
def figure1():
    # top row hosts the fixed-aspect (16:10.6) pipeline schematic; give it enough height so its
    # equal-aspect box is ~9.7in wide (matching the standalone Fig1a) and text does not overflow
    fig = plt.figure(figsize=(13.0, 14.3))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.45, 2.3], hspace=0.16, wspace=0.32)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    # ---- (a) upgraded pipeline schematic (core pipeline + novel analysis layer)
    draw_pipeline(ax_a)
    panel_letter(ax_a, "a", dx=-0.03, dy=0.98)

    # ---- (b) module-clock performance
    perf = pd.read_csv(MC / "module_clock_performance.csv")
    mo = perf[perf.outcome == "Mortality"].set_index("annotation")["LOFO_pearson"]
    ch = perf[perf.outcome == "Chrono"].set_index("annotation")["LOFO_pearson"]
    order = mo.sort_values().index
    y = np.arange(len(order)); h = 0.4
    ax_b.barh(y + h/2, mo.loc[order], height=h, color=CB_VERM, label="Mortality")
    ax_b.barh(y - h/2, ch.loc[order], height=h, color=CB_BLUE, label="Chronological")
    ax_b.set_yticks(y); ax_b.set_yticklabels([short(a, 30) for a in order], fontsize=7)
    ax_b.set_xlabel("Group-out-of-fold Pearson r (predicted vs observed)", fontsize=9)
    ax_b.axvline(0, color="k", lw=0.6); ax_b.legend(loc="lower right", fontsize=9)
    ax_b.set_xlim(0, 0.7); ax_b.tick_params(axis="x", labelsize=8)
    panel_letter(ax_b, "b", dx=-0.30)

    # ---- (c) Klotho-KO validation (colorblind-safe)
    import joblib, tage_prep as tp
    CLOCKS = {}
    for f in sorted(MC.glob("module_Mortality_*.pkl")):
        d = joblib.load(f); CLOCKS[d["annotation"].split("/")[0].strip()] = d
    EXT = Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "extdata"
    counts = pd.read_csv(EXT / "Exprs_example.csv", index_col=0)
    meta = pd.read_csv(EXT / "Metadata_example.csv", index_col=0)

    def module_scores(sd):
        X = sd.copy(); X.columns = X.columns.map(str)
        return pd.DataFrame({n: d["pipeline"].predict(
            X.reindex(columns=[str(g) for g in d["genes"]]).values) for n, d in CLOCKS.items()},
            index=sd.index)

    kod = {}
    for tissue in meta["Tissue"].unique():
        mt = meta[meta.Tissue == tissue]; ct = counts[mt.index]
        pp = tp.preprocess(ct, mt, species="mouse", gene_mapping_type="Ensembl",
                           control_group_column="Genotype", control_group_label="WT")
        S = module_scores(pp["scaled_diff"]); S["G"] = mt["Genotype"].values
        g = S.groupby("G").mean(numeric_only=True).T
        kod[tissue] = g.get("Klotho KO", 0) - g.get("WT", 0)
    K = pd.DataFrame(kod)
    K = K.loc[K.mean(axis=1).sort_values().index]
    yy = np.arange(len(K)); h = 0.4
    cb_pair = [CB_ORANGE, CB_BLUE]
    for i, tissue in enumerate(K.columns):
        ax_c.barh(yy + (h/2 if i == 0 else -h/2), K[tissue], height=h,
                  label=tissue, color=cb_pair[i % 2])
    ax_c.set_yticks(yy); ax_c.set_yticklabels([short(a, 28) for a in K.index], fontsize=7)
    ax_c.axvline(0, color="k", lw=0.6)
    ax_c.set_xlabel("Δ module mortality score (Klotho-KO − WT)", fontsize=9)
    ax_c.legend(fontsize=9); ax_c.tick_params(axis="x", labelsize=8)
    panel_letter(ax_c, "c", dx=-0.28)

    fig.savefig(OUT / "Fig1.png", dpi=300, bbox_inches="tight"); plt.close(fig)
    print("saved Fig1.png")


# ==================== FIGURE 1 (split panels) ====================
def figure1_split():
    """Emit Fig1a/b/c as separate files (no panel letters) for manual assembly."""
    # ---- (a) pipeline schematic
    fig, ax_a = plt.subplots(figsize=(13, 8.2))
    draw_pipeline(ax_a)
    fig.savefig(OUT / "Fig1a.png", dpi=300, bbox_inches="tight"); plt.close(fig)

    # ---- (b) module-clock performance
    perf = pd.read_csv(MC / "module_clock_performance.csv")
    mo = perf[perf.outcome == "Mortality"].set_index("annotation")["LOFO_pearson"]
    ch = perf[perf.outcome == "Chrono"].set_index("annotation")["LOFO_pearson"]
    order = mo.sort_values().index
    y = np.arange(len(order)); h = 0.4
    fig, ax_b = plt.subplots(figsize=(7.2, 8))
    ax_b.barh(y + h/2, mo.loc[order], height=h, color=CB_VERM, label="Mortality")
    ax_b.barh(y - h/2, ch.loc[order], height=h, color=CB_BLUE, label="Chronological")
    ax_b.set_yticks(y); ax_b.set_yticklabels([short(a, 34) for a in order], fontsize=8)
    ax_b.set_xlabel("Group-out-of-fold Pearson r (predicted vs observed)", fontsize=9)
    ax_b.axvline(0, color="k", lw=0.6); ax_b.legend(loc="lower right", fontsize=9)
    ax_b.set_xlim(0, 0.7)
    fig.savefig(OUT / "Fig1b.png", dpi=300, bbox_inches="tight"); plt.close(fig)

    # ---- (c) Klotho-KO validation (colorblind-safe)
    import joblib, tage_prep as tp
    CLOCKS = {}
    for f in sorted(MC.glob("module_Mortality_*.pkl")):
        d = joblib.load(f); CLOCKS[d["annotation"].split("/")[0].strip()] = d
    EXT = Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "extdata"
    counts = pd.read_csv(EXT / "Exprs_example.csv", index_col=0)
    meta = pd.read_csv(EXT / "Metadata_example.csv", index_col=0)

    def module_scores(sd):
        X = sd.copy(); X.columns = X.columns.map(str)
        return pd.DataFrame({n: d["pipeline"].predict(
            X.reindex(columns=[str(g) for g in d["genes"]]).values) for n, d in CLOCKS.items()},
            index=sd.index)

    kod = {}
    for tissue in meta["Tissue"].unique():
        mt = meta[meta.Tissue == tissue]; ct = counts[mt.index]
        pp = tp.preprocess(ct, mt, species="mouse", gene_mapping_type="Ensembl",
                           control_group_column="Genotype", control_group_label="WT")
        S = module_scores(pp["scaled_diff"]); S["G"] = mt["Genotype"].values
        g = S.groupby("G").mean(numeric_only=True).T
        kod[tissue] = g.get("Klotho KO", 0) - g.get("WT", 0)
    K = pd.DataFrame(kod)
    K = K.loc[K.mean(axis=1).sort_values().index]
    fig, ax_c = plt.subplots(figsize=(6.5, 8))
    yy = np.arange(len(K)); h = 0.4
    cb_pair = [CB_ORANGE, CB_BLUE]
    for i, tissue in enumerate(K.columns):
        ax_c.barh(yy + (h/2 if i == 0 else -h/2), K[tissue], height=h,
                  label=tissue, color=cb_pair[i % 2])
    ax_c.set_yticks(yy); ax_c.set_yticklabels([short(a, 30) for a in K.index], fontsize=8)
    ax_c.axvline(0, color="k", lw=0.6)
    ax_c.set_xlabel("Δ module mortality score (Klotho-KO − WT)", fontsize=9)
    ax_c.legend(fontsize=9)
    fig.savefig(OUT / "Fig1c.png", dpi=300, bbox_inches="tight"); plt.close(fig)
    print("saved Fig1a.png, Fig1b.png, Fig1c.png (separate panels)")


# ============================ FIGURE 2 ============================
def figure2():
    R = pd.read_csv(RES / "ipf_module_stats_adjusted.csv")
    real = R[R.module != "Composite"].copy()
    mat = real.pivot_table(index="celltype", columns="module", values="beta_IPF")
    fmat = real.pivot_table(index="celltype", columns="module", values="fdr_IPF")
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]
    mat = mat[mat.var().sort_values(ascending=False).index]
    fmat = fmat.reindex(index=mat.index, columns=mat.columns)

    fig, ax = plt.subplots(figsize=(15, 6.2))
    vmax = np.nanpercentile(np.abs(mat.values), 95)
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=9)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            f = fmat.values[i, j]
            if f < 0.10:
                ax.text(j, i, "*" if f >= 0.05 else "**", ha="center", va="center", fontsize=9)
    cb = plt.colorbar(im, label="β  (IPF − Control, age/sex adjusted)", shrink=0.75)
    ax.set_xlabel(""); ax.set_ylabel("")
    plt.tight_layout()
    # cell-type × module heatmap is Figure 3 (cited after the composite-masking figure)
    fig.savefig(OUT / "Fig3.png", dpi=300, bbox_inches="tight"); plt.close(fig)
    print("saved Fig3.png (heatmap)")


# ============================ FIGURE 3 ============================
def figure3():
    fig = plt.figure(figsize=(13.5, 8.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.28, 1.0], hspace=0.32, wspace=0.30)
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])

    # ---- (a) composite masking
    G = pd.read_csv(RES / "module_global_robust.csv")
    comp = G[G.module == "Composite"].iloc[0]
    M = G[G.module != "Composite"].sort_values("beta")
    colors = []
    for _, r in M.iterrows():
        if r.fdr < 0.05:
            colors.append("#b2182b" if r.beta > 0 else "#2166ac")
        elif r.fdr < 0.10:
            colors.append("#ef8a62" if r.beta > 0 else "#67a9cf")
        else:
            colors.append("#cccccc")
    y = np.arange(len(M))
    ax_a.barh(y, M.beta, color=colors, edgecolor="#333", lw=0.4)
    ax_a.set_yticks(y); ax_a.set_yticklabels([short(m, 30) for m in M.module], fontsize=8)
    xr = float(np.nanmax(np.abs(M.beta)))
    ax_a.set_xlim(-xr * 1.45, xr * 1.35)
    for i, (_, r) in enumerate(M.iterrows()):
        if r.fdr < 0.10:
            ax_a.text(r.beta + (0.012 * xr if r.beta > 0 else -0.012 * xr), i,
                      f"FDR={r.fdr:.3f}", va="center",
                      ha="left" if r.beta > 0 else "right", fontsize=6.5, color="#222")
    ax_a.axvline(0, color="k", lw=0.7)
    ax_a.set_xlabel("Global IPF−Control effect (β; donor-clustered, cell-type/age/sex adjusted)", fontsize=9)
    ax_a.text(0.97, 0.03,
              f"Composite mortality clock:\nβ = {comp.beta:.1f} months, p = {comp.p:.2f}  (n.s.)",
              transform=ax_a.transAxes, ha="right", va="bottom", fontsize=9.5,
              bbox=dict(boxstyle="round,pad=0.4", fc="#f7f7f7", ec="#999"))
    leg = [Patch(fc="#b2182b", label="up, FDR<0.05"), Patch(fc="#2166ac", label="down, FDR<0.05"),
           Patch(fc="#ef8a62", label="up, FDR<0.10"), Patch(fc="#67a9cf", label="down, FDR<0.10"),
           Patch(fc="#cccccc", label="n.s.")]
    ax_a.legend(handles=leg, loc="upper left", fontsize=8)
    panel_letter(ax_a, "a", dx=-0.16)

    R = pd.read_csv(RES / "ipf_module_stats_adjusted.csv")
    real = R[R.module != "Composite"].copy()

    # ---- (b) age-confound robustness
    r_na = real[["naive_delta", "beta_IPF"]].corr().iloc[0, 1]
    ax_b.scatter(real.naive_delta, real.beta_IPF, s=16, alpha=0.6, color=CB_BLUE)
    lim = float(np.nanmax(np.abs(np.r_[real.naive_delta, real.beta_IPF])) * 1.05)
    ax_b.plot([-lim, lim], [-lim, lim], "k--", lw=0.7)
    ax_b.set_xlim(-lim, lim); ax_b.set_ylim(-lim, lim)
    ax_b.set_xlabel("naive Δ (unadjusted)", fontsize=9)
    ax_b.set_ylabel("β IPF (age/sex adjusted)", fontsize=9)
    ax_b.text(0.04, 0.94, f"r = {r_na:.2f}", transform=ax_b.transAxes, fontsize=10, va="top")
    panel_letter(ax_b, "b", dx=-0.18)

    # ---- (c) IPF vs COPD specificity (de-cluttered labels)
    cc = real.dropna(subset=["beta_COPD", "beta_IPF"]).copy()
    sigipf = cc.fdr_IPF < 0.10
    ax_c.scatter(cc.beta_COPD[~sigipf], cc.beta_IPF[~sigipf], s=16, color="#bbb", label="n.s.")
    ax_c.scatter(cc.beta_COPD[sigipf], cc.beta_IPF[sigipf], s=30, color=CB_VERM, label="IPF FDR<0.10")
    lim = float(np.nanmax(np.abs(np.r_[cc.beta_COPD, cc.beta_IPF])) * 1.08)
    ax_c.plot([-lim, lim], [-lim, lim], "k--", lw=0.7, alpha=0.6)
    ax_c.axhline(0, color="k", lw=0.5); ax_c.axvline(0, color="k", lw=0.5)
    # label only the most IPF-specific significant points (far from diagonal),
    # staggering label positions + leader lines to avoid overlap
    sig = cc[sigipf].copy()
    sig["offdiag"] = (sig.beta_IPF - sig.beta_COPD).abs()
    # spread labels well clear of the dense lower-left cluster; full module names (no truncation)
    stagger = [(12, 18), (14, -6), (-96, 20), (-96, -18), (16, -22)]
    for k, (_, r) in enumerate(sig.sort_values("offdiag", ascending=False).head(5).iterrows()):
        ax_c.annotate(f"{short(r.celltype,12)}: {short(r.module,20)}",
                      (r.beta_COPD, r.beta_IPF), fontsize=6.8,
                      xytext=stagger[k % len(stagger)], textcoords="offset points",
                      ha="left" if stagger[k % len(stagger)][0] >= 0 else "left",
                      alpha=0.95, arrowprops=dict(arrowstyle="-", lw=0.4, color="#777"))
    ax_c.set_xlabel("β COPD − Control", fontsize=9); ax_c.set_ylabel("β IPF − Control", fontsize=9)
    ax_c.legend(fontsize=8, loc="lower right"); ax_c.set_xlim(-lim, lim); ax_c.set_ylim(-lim, lim)
    panel_letter(ax_c, "c", dx=-0.18)

    # composite-masking / age / specificity is Figure 2 (first figure cited in Results)
    fig.savefig(OUT / "Fig2.png", dpi=300, bbox_inches="tight"); plt.close(fig)
    print("saved Fig2.png (composite masking)")


# ============================ FIGURE 4 ============================
def figure4():
    fig = plt.figure(figsize=(14, 8.4))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.15], hspace=0.42)

    # ---- (a) per-cell-type cross-cohort replication
    cmp = pd.read_csv(RES / "replication_coarse_compare.csv")
    order = cmp.groupby("coarse").apply(lambda g: pearsonr(g.beta_adams, g.beta_hab)[0]).sort_values(ascending=False)
    cells = order.index.tolist()
    top = outer[0].subgridspec(1, len(cells), wspace=0.42)
    for k, co in enumerate(cells):
        ax = fig.add_subplot(top[0, k])
        g = cmp[cmp.coarse == co]
        r = pearsonr(g.beta_adams, g.beta_hab)[0]
        sc = np.mean(np.sign(g.beta_adams) == np.sign(g.beta_hab))
        ok = r > 0.3
        ax.scatter(g.beta_adams, g.beta_hab, s=20, color=CB_VERM if ok else "#888", alpha=0.8)
        lim = float(np.nanmax(np.abs(np.r_[g.beta_adams, g.beta_hab])) * 1.15)
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.7, alpha=0.6)
        ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.text(0.5, 1.06, f"{co}\nr={r:+.2f}, sign {sc:.0%}", transform=ax.transAxes,
                ha="center", va="bottom", fontsize=8.5, color=CB_VERM if ok else "#333")
        ax.set_xlabel("Adams β", fontsize=8)
        if k == 0:
            ax.set_ylabel("Habermann β", fontsize=9)
            panel_letter(ax, "a", dx=-0.42, dy=1.12)
        ax.tick_params(labelsize=7)

    bot = outer[1].subgridspec(1, 2, wspace=0.28)
    # ---- (b) split-half ceiling vs cross-cohort, for the NON-replicating cell types only
    # (NK/T are the replicated positive controls shown in panel a; panel b isolates the
    #  myeloid/epithelial cell types whose cross-cohort failure we attribute to power vs sampling)
    ax_b = fig.add_subplot(bot[0, 0])
    C = pd.read_csv(MIL / "stageC_calibrate.csv")
    C = C[C.prior == "failed"].reset_index(drop=True)
    x = np.arange(len(C))
    ax_b.bar(x, C.splithalf_med, 0.55, color=CB_SKY,
             yerr=[C.splithalf_med - C.splithalf_lo, C.splithalf_hi - C.splithalf_med],
             capsize=4, label="within-cohort split-half ceiling")
    ax_b.scatter(x, C.r_cross_obs, color=CB_VERM, zorder=5, s=64,
                 label="Adams↔Habermann cross-cohort")
    for xi, (_, r) in zip(x, C.iterrows()):
        ax_b.annotate(f"p={r.p_cross:.2f}", (xi, r.r_cross_obs), textcoords="offset points",
                      xytext=(7, 0), fontsize=7.5, color=CB_VERM)
    ax_b.axhline(0, color="k", lw=0.6)
    ax_b.set_xticks(x); ax_b.set_xticklabels(
        [f"{c}\n(non-replicating)" for c in C.celltype], fontsize=8)
    ax_b.set_ylabel("module-effect concordance (Pearson r)", fontsize=9)
    ax_b.legend(fontsize=8, loc="upper left")
    panel_letter(ax_b, "b", dx=-0.16)

    # ---- (c) bulk recovery across 3 cohorts
    ax_c = fig.add_subplot(bot[0, 1])
    D = pd.read_csv(MIL / "stageD_bulk.csv")
    x = np.arange(len(D))
    ax_c.bar(x, D.r, 0.6, color=[CB_VERM if p < 0.05 else "#95a5a6" for p in D.p_perm])
    for xi, row in zip(x, D.itertuples()):
        ax_c.annotate(f"p={row.p_perm:.3f}", (xi, row.r), textcoords="offset points",
                      xytext=(0, 4), ha="center", fontsize=8)
    ax_c.axhline(0, color="k", lw=0.6)
    ax_c.set_xticks(x); ax_c.set_xticklabels(
        [p.replace(" vs ", "\nvs ") for p in D.pair], rotation=0, fontsize=7.5)
    ax_c.set_ylabel("whole-tissue module-effect concordance (r)", fontsize=9)
    ax_c.set_ylim(0, max(0.6, float(D.r.max()) * 1.2))
    panel_letter(ax_c, "c", dx=-0.16)

    fig.savefig(OUT / "Fig4.png", dpi=300, bbox_inches="tight"); plt.close(fig)
    print("saved Fig4.png")


# ============================ FIGURE 5 ============================
def figure5():
    PM = pd.read_csv(MIL / "stage3_multi.csv")
    ct_order = [c for c in ["Macrophage", "Monocyte", "Tcell", "NK"] if c in PM.celltype.unique()]
    letters = ["a", "b", "c", "d"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.6)); axes = axes.ravel()
    for ax, co, let in zip(axes, ct_order, letters):
        sub = PM[PM.celltype == co]
        # marker order: pole A then B, preserving appearance order
        seen, order = set(), []
        for pole in ["A", "B"]:
            for m in sub[sub.pole == pole].marker:
                if m not in seen:
                    seen.add(m); order.append((m, pole))
        markers = [m for m, _ in order]
        nA = sum(1 for _, p in order if p == "A")
        x = np.arange(len(markers)); wdt = 0.38
        for i, grp in enumerate(["IPF", "Control"]):
            vals = [sub[(sub.group == grp) & (sub.marker == m)].rho.values
                    for m in markers]
            vals = [v[0] if len(v) else np.nan for v in vals]
            ax.bar(x + (i - 0.5) * wdt, vals, wdt, label=grp,
                   color=CB_VERM if grp == "IPF" else "#7f8c8d")
        ax.axhline(0, color="k", lw=0.6)
        ax.axvline(nA - 0.5, color="k", ls="--", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(markers, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Spearman(attention, marker)", fontsize=9)
        ax.text(0.02, 1.03, f"{co}   (left of dashed = disease-expanded pole)",
                transform=ax.transAxes, fontsize=8.5, va="bottom")
        ax.legend(fontsize=8, loc="best")
        panel_letter(ax, let, dx=-0.11)
    for ax in axes[len(ct_order):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(OUT / "Fig5.png", dpi=300, bbox_inches="tight"); plt.close(fig)
    print("saved Fig5.png")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--split-fig1", action="store_true",
                    help="also emit Fig1a/b/c as separate panel files")
    args, _ = ap.parse_known_args()
    figure1()
    if args.split_fig1:
        figure1_split()
    figure2()
    figure3()
    figure4()
    figure5()
    print("\nAll figures written to", OUT)
