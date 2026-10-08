"""Figures for the revised manuscript (main Fig. 1-6, Supplementary Fig. S1-S6).

Reads the tables written by rev_01 ... rev_15 (results/revision) and the MIL tables
(results/mil); writes PNG (600 dpi) and PDF to figures_revision/. The submission figures in
figures/ are left untouched.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Patch  # noqa: E402

import rev_common as rc  # noqa: E402

FIG = rc.TOOLKIT / "figures_revision"
FIG.mkdir(exist_ok=True)
O = rc.OUT
MIL = rc.RES / "mil"

plt.rcParams.update({
    "font.family": ["Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.titlesize": 7.5,
    "axes.labelsize": 7, "xtick.labelsize": 6.3, "ytick.labelsize": 6.3, "legend.fontsize": 6.3,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False,
    "axes.spines.right": False, "savefig.dpi": 600, "pdf.fonttype": 42, "legend.frameon": False})

MODULE_ORDER = [
    "Innate immunity", "Interferon signaling", "Adaptive immunity",
    "Respiration", "OxPhos", "Heme metabolism",
    "Cell cycle", "Chromatin modification (1)", "Chromatin modification (2)", "mRNA splicing", "Translation",
    "Protein processing in ER", "Heat stress response", "Nrf2 signaling", "Apoptosis",
    "Lipid met", "Cholesterol met (darkred)", "Cholesterol met (salmon)", "Fatty acid met", "Amino acid met",
    "ECM organization", "VEGF signaling", "Muscle contraction"]
LABEL = {
    "Innate immunity": "Innate immunity/inflammation", "Interferon signaling": "Interferon signaling",
    "Adaptive immunity": "Adaptive immunity/T cell", "Respiration": "Mito. respiration/translation",
    "OxPhos": "OxPhos/TCA cycle", "Heme metabolism": "Heme metabolism/ROS",
    "Cell cycle": "Cell cycle/DNA replication", "Chromatin modification (1)": "Chromatin modification 1",
    "Chromatin modification (2)": "Chromatin modification 2", "mRNA splicing": "mRNA splicing",
    "Translation": "Translation/ribosome", "Protein processing in ER": "ER protein processing/UPR",
    "Heat stress response": "Heat stress response", "Nrf2 signaling": "Nrf2 signaling/proteasome",
    "Apoptosis": "Apoptosis/proteasome", "Lipid met": "Lipid metabolism/PPAR",
    "Cholesterol met (darkred)": "Cholesterol metabolism/mTOR",
    "Cholesterol met (salmon)": "Cholesterol metabolism/platelet", "Fatty acid met": "Fatty acid β-oxidation",
    "Amino acid met": "Amino acid/xenobiotic metabolism", "ECM organization": "ECM organization/EMT",
    "VEGF signaling": "VEGF signaling", "Muscle contraction": "Muscle contraction/glycolysis",
    "Composite": "Composite clock", "Non-module genes": "Non-module genes"}
CT_LABEL = {"B": "B cell", "Ciliated": "Ciliated", "Fibroblast": "Fibroblast", "Macrophage": "Macrophage",
            "Macrophage_Alveolar": "Alveolar macrophage", "NK": "NK cell", "T": "T cell",
            "T_Cytotoxic": "Cytotoxic T cell", "cDC2": "cDC2", "cMonocyte": "Classical monocyte",
            "ncMonocyte": "Non-classical monocyte"}
CT_SHORT = {"B": "B", "Ciliated": "Ciliated", "Fibroblast": "Fibro", "Macrophage": "Mac",
            "Macrophage_Alveolar": "AM", "NK": "NK", "T": "T", "T_Cytotoxic": "CTL", "cDC2": "cDC2",
            "cMonocyte": "cMono", "ncMonocyte": "ncMono"}
CT_ORDER = ["Macrophage", "Macrophage_Alveolar", "cMonocyte", "ncMonocyte", "cDC2", "NK", "T", "T_Cytotoxic",
            "B", "Ciliated", "Fibroblast"]
RED, BLUE, GREY, ORANGE, DARK = "#c0392b", "#2c6fb0", "#9aa3ab", "#e08a1e", "#2b2f33"


def read(name, base=O):
    return pd.read_csv(base / name, keep_default_na=False, na_values=[""])


def panel(ax, letter, x=-0.12, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ============================================================================ Fig. 1
def fig1():
    fig, ax = plt.subplots(figsize=(7.1, 4.9))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 72)
    ax.axis("off")
    col = {"pub": ("#eceff1", "#7d8790"), "re": ("#dceaf7", "#2c6fb0"), "new": ("#fde9cf", "#c46f00")}

    def box(x, y, w, h, kind, title, body, num=None):
        fc, ec = col[kind]
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.2",
                                    fc=fc, ec=ec, lw=0.9))
        ax.text(x + w / 2, y + h - 1.9, title, ha="center", va="top", fontsize=6.4, fontweight="bold", color=DARK)
        ax.text(x + w / 2, y + h - 5.6, body, ha="center", va="top", fontsize=5.3, color="#33393f",
                linespacing=1.3)
        if num is not None:
            ax.add_patch(plt.Circle((x + 1.0, y + h - 0.2), 1.45, fc=ec, ec="white", lw=0.8, zorder=5))
            ax.text(x + 1.0, y + h - 0.2, str(num), ha="center", va="center", fontsize=5.8, color="white",
                    fontweight="bold", zorder=6)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color="#4a5560", mutation_scale=8))

    def label(y, text):
        ax.text(99.5, y, text, fontsize=6.4, color="#4a5560", style="italic", va="bottom", ha="right")

    w, h = 22.5, 16
    xs = [1, 26, 51, 76]
    y1, y2, y3 = 47, 25.5, 3.5
    label(y1 + h + 1.2, "rodent domain")
    box(xs[0], y1, w, h, "pub", "Rodent training data",
        "4,539 relative profiles,\n96 datasets (mouse 3,876;\nrat 663); targets Δlog10\nhazard and Δ(age/lifespan)", 1)
    box(xs[1], y1, w, h, "pub", "Predefined modules",
        "23 co-expression modules,\n1,995 disjoint genes,\nused verbatim\n(Tyshkovskiy et al.)", 2)
    box(xs[2], y1, w, h, "re", "Module clocks",
        "elastic net per module\nand target; fitted pipelines\nreleased; reach reported\naccuracy", 3)
    box(xs[3], y1, w, h, "new", "Rodent-domain checks",
        "dataset-grouped and species-\nstratified CV; gene-set and\nrandom-gene benchmarks;\nKlotho-KO, exact tests", 4)
    for a, b in zip(xs[:3], xs[1:]):
        arrow(a + w + 0.4, y1 + h / 2, b - 0.4, y1 + h / 2)
    label(y2 + h + 1.2, "human single-cell application")
    box(xs[0], y2, w, h, "pub", "Ortholog mapping,\npreprocessing",
        "\nhuman → mouse Entrez, RLE,\nper-sample scaling, control-\nmedian subtraction (tAge)", 5)
    box(xs[1], y2, w, h, "new", "Pseudobulk scoring",
        "donor × cell-type pseudobulks\n(≥ 50 cells; 11 cell types);\nGSE136831: 32 IPF,\n28 control, 18 COPD", 6)
    box(xs[2], y2, w, h, "new", "Disease contrasts",
        "IPF and COPD vs control;\nHC3, permutation, mixed model,\nwild bootstrap; direct IPF–COPD;\ncompetitive null gene sets", 7)
    box(xs[3], y2, w, h, "new", "Human validity checks",
        "within-cell-type age (HLCA,\ncontrol donors); GTEx lung;\ndepth and cell-number\nmatching; age overlap", 8)
    for a, b in zip(xs[:3], xs[1:]):
        arrow(a + w + 0.4, y2 + h / 2, b - 0.4, y2 + h / 2)
    xm = xs[2] + w / 2
    ax.plot([xm, xm, xs[0] + w / 2], [y1 - 0.6, y2 + h + 2.4, y2 + h + 2.4], color="#4a5560", lw=0.9)
    arrow(xs[0] + w / 2, y2 + h + 2.45, xs[0] + w / 2, y2 + h + 0.6)
    label(y3 + h + 0.2, "diagnostics introduced here")
    box(xs[0], y3, 35, h - 1, "new", "Composite decomposition",
        "composite clock split into its own module-gene\ncomponents; total signal S, net N,\nmasking index MI = (S − N)/S;\npermutation null and donor bootstrap", 9)
    box(38.5, y3, 35, h - 1, "new", "Replication diagnostics",
        "two-sided permutation concordance; specific effects;\nsplit-half reliability; composition and\nsub-state × disease models; checklist with\nsimulated operating characteristics", 10)
    box(76, y3, 22.5, h - 1, "new", "Supplementary",
        "attention-based MIL\n(negative result: attention\nrecovers the annotation)", None)
    arrow(xs[1] + w / 2, y2 - 0.6, 18.5, y3 + h - 0.4)
    arrow(xs[2] + w / 2, y2 - 0.6, 56, y3 + h - 0.4)
    handles = [Patch(fc=col[k][0], ec=col[k][1], label=lab) for k, lab in
               [("pub", "previously published (adopted)"), ("re", "re-implemented / retrained here"),
                ("new", "introduced in this study")]]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.03), ncol=3, fontsize=6.0,
              handlelength=1.4, columnspacing=1.2)
    save(fig, "Fig1")


# ============================================================================ Fig. 2
def fig2():
    cv = read("rev02_species_cv.csv")
    cv = cv[cv.kind == "module"]
    bench = read("rev02_benchmarks.csv")
    union = bench[bench.kind == "union_modules"].set_index("outcome").r_oof_all
    allg = read("rev02b_allgene_benchmark.csv").set_index("outcome").r_oof_all \
        if (O / "rev02b_allgene_benchmark.csv").exists() else None
    fig = plt.figure(figsize=(7.1, 6.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1], height_ratios=[1, 0.9], hspace=0.45, wspace=0.55)

    ax = fig.add_subplot(gs[:, 0])
    m = cv[cv.outcome == "Mortality"].set_index("module").reindex(MODULE_ORDER)
    c = cv[cv.outcome == "Chrono"].set_index("module").reindex(MODULE_ORDER)
    y = np.arange(len(MODULE_ORDER))[::-1]
    ax.hlines(y, m.r_random_mean, m.r_random_max, color="#c9ced3", lw=3.2, zorder=1)
    ax.scatter(m.r_oof_all, y, s=16, color=RED, zorder=3, label="mortality clock")
    ax.scatter(c.r_oof_all, y, s=16, facecolor="white", edgecolor=BLUE, lw=0.9, zorder=3,
               label="chronological clock")
    ax.axvline(union["Mortality"], color=RED, lw=0.7, ls="--")
    ax.text(union["Mortality"] + 0.012, 9, f"all 1,995 module genes (r = {union['Mortality']:.2f})", color=RED,
            fontsize=5.4, rotation=90, ha="left", va="center")
    if allg is not None:
        ax.axvline(allg["Mortality"], color=DARK, lw=0.7, ls=":")
        ax.text(allg["Mortality"] + 0.012, 9, f"all 18,286 genes (r = {allg['Mortality']:.2f})", color=DARK,
                fontsize=5.4, rotation=90, ha="left", va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[k] for k in MODULE_ORDER])
    ax.set_xlabel("Out-of-fold Pearson r\n(5-fold, whole datasets held out)")
    ax.set_xlim(-0.05, 0.95)
    ax.set_ylim(-0.8, len(y) + 1.6)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=RED, ms=4, label="mortality clock"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec=BLUE, ms=4, label="chronological clock"),
                       Line2D([], [], color="#c9ced3", lw=3.2, label="random genes, mean–max (10 draws)")],
              loc="upper left", fontsize=5.6, bbox_to_anchor=(-0.02, 1.01))
    panel(ax, "a", x=-0.9, y=1.0)

    ax = fig.add_subplot(gs[0, 1:])
    cols = [("r_oof_all", "pooled CV,\nall samples"), ("r_oof_mouse", "pooled CV,\nmouse"),
            ("r_oof_rat", "pooled CV,\nrat"), ("r_mouse_to_rat", "train mouse,\ntest rat"),
            ("r_rat_to_mouse", "train rat,\ntest mouse"), ("r_within_rat_lodo", "rat only,\nLODO")]
    data = [m[k].dropna().values for k, _ in cols]
    bp = ax.boxplot(data, widths=0.5, patch_artist=True, showfliers=False,
                    medianprops=dict(color=DARK, lw=1), boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6),
                    capprops=dict(lw=0.6))
    for p in bp["boxes"]:
        p.set_facecolor("#f4d6d2")
    rng = np.random.default_rng(1)
    for i, d in enumerate(data):
        ax.scatter(i + 1 + rng.uniform(-0.15, 0.15, len(d)), d, s=5, color=RED, alpha=0.7, lw=0, zorder=3)
        ax.text(i + 1, 0.93, f"{np.median(d):.2f}", ha="center", fontsize=5.8)
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_xticks(range(1, len(cols) + 1))
    ax.set_xticklabels([lab for _, lab in cols], fontsize=5.8)
    ax.set_ylabel("Pearson r (mortality clocks)")
    ax.set_ylim(-0.35, 1.0)
    ax.set_title("Species-stratified performance (23 modules; median above)", fontsize=7)
    panel(ax, "b", x=-0.1)

    ax = fig.add_subplot(gs[1, 1])
    k = read("rev04_klotho_modules.csv")
    kk = k[k.tissue == "Kidney"].set_index("module")
    km = k[k.tissue == "Skeletal muscle"].set_index("module")
    mods = MODULE_ORDER
    sig_both = [(kk.loc[x, "fdr_perm"] < .1) and (km.loc[x, "fdr_perm"] < .1) for x in mods]
    ax.axhline(0, color=GREY, lw=0.5)
    ax.axvline(0, color=GREY, lw=0.5)
    ax.errorbar(kk.loc[mods, "delta"], km.loc[mods, "delta"],
                xerr=[kk.loc[mods, "delta"] - kk.loc[mods, "ci_lo"], kk.loc[mods, "ci_hi"] - kk.loc[mods, "delta"]],
                yerr=[km.loc[mods, "delta"] - km.loc[mods, "ci_lo"], km.loc[mods, "ci_hi"] - km.loc[mods, "delta"]],
                fmt="none", ecolor="#d5d9dd", elinewidth=0.5, zorder=1)
    ax.scatter(kk.loc[mods, "delta"], km.loc[mods, "delta"], s=12,
               c=[RED if s else "#7d8790" for s in sig_both], zorder=3, lw=0)
    ax.scatter(kk.loc["Composite", "delta"], km.loc["Composite", "delta"], marker="*", s=70, color=ORANGE,
               edgecolor=DARK, lw=0.4, zorder=4)
    ax.annotate("composite", (kk.loc["Composite", "delta"], km.loc["Composite", "delta"]), xytext=(-8, -11),
                textcoords="offset points", fontsize=5.6)
    ax.annotate("adaptive\nimmunity", (kk.loc["Adaptive immunity", "delta"], km.loc["Adaptive immunity", "delta"]),
                xytext=(-34, 4), textcoords="offset points", fontsize=5.6, color=RED)
    s = read("rev04_klotho_summary.csv").set_index("tissue")
    ax.text(0.02, 0.98, f"modules > 0: kidney {int(s.loc['Kidney','n_positive'])}/23 "
                        f"(p = {s.loc['Kidney','p_count_exact']:.2f})\nmuscle {int(s.loc['Skeletal muscle','n_positive'])}/23 "
                        f"(p = {s.loc['Skeletal muscle','p_count_exact']:.2f}); both "
                        f"{int(s.loc['Kidney','n_positive_both_tissues'])}/23",
            transform=ax.transAxes, fontsize=5.3, va="top")
    ax.set_xlabel("Kidney: KO − WT (Δlog10 hazard)")
    ax.set_ylabel("Muscle: KO − WT (Δlog10 hazard)")
    ax.set_title("Klotho-KO vs WT (n = 6 vs 6 per tissue)", fontsize=7)
    panel(ax, "c", x=-0.3)

    ax = fig.add_subplot(gs[1, 2])
    h = read("rev05_hlca_age.csv")
    a = read("rev05_adams_control_age.csv")
    h, a = h[h.outcome == "Mortality"], a[a.outcome == "Mortality"]
    sets = [(h[h.module != "Composite"].rho, "HLCA\nunadjusted"), (h[h.module != "Composite"].rho_adj, "HLCA\nadjusted"),
            (a[a.module != "Composite"].rho_adj, "GSE136831\ncontrols, adj.")]
    vp = ax.violinplot([s_.dropna() for s_, _ in sets], showmedians=True, widths=0.75)
    for b in vp["bodies"]:
        b.set_facecolor("#d6e4f2")
        b.set_edgecolor(BLUE)
        b.set_alpha(1)
        b.set_linewidth(0.5)
    for key in ["cmedians", "cmins", "cmaxes", "cbars"]:
        vp[key].set_color(BLUE)
        vp[key].set_linewidth(0.6)
    fd = [(h[h.module != "Composite"].fdr < .1).sum(), (h[h.module != "Composite"].fdr_adj < .1).sum(),
          (a[a.module != "Composite"].fdr_adj < .1).sum()]
    for i, f in enumerate(fd):
        ax.text(i + 1, -1.17, f"FDR < 0.1\nn = {f}", ha="center", fontsize=5.2)
    g = read("rev05b_gtex_validation.csv")
    gm = g[(g.outcome == "Mortality") & (g.module == "Composite")]
    if len(gm):
        ax.axhline(gm.adj_spearman.iloc[0], color=ORANGE, lw=0.8, ls="--")
        ax.text(0.55, 1.0, f"– – GTEx lung bulk, composite ρ = {gm.adj_spearman.iloc[0]:.2f}",
                fontsize=5.3, color=ORANGE, ha="left", va="bottom")
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([lab for _, lab in sets], fontsize=5.8)
    ax.set_ylim(-1.25, 1.1)
    ax.set_ylabel("Spearman ρ with donor age\n(module × cell type)")
    ax.set_title("Within-cell-type age association", fontsize=7)
    panel(ax, "d", x=-0.3)
    save(fig, "Fig2")


# ============================================================================ Fig. 3
def fig3():
    G = read("rev07_masking_global.csv")
    summ = G[G.component == "__summary__"].iloc[0]
    G = G[G.component != "__summary__"].copy()
    fig = plt.figure(figsize=(7.1, 5.0))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.3, 0.85, 1.15], hspace=0.75, wspace=1.05)

    ax = fig.add_subplot(gs[:, 0])
    G["fdr"] = G.fdr.astype(float)
    mods = G[G.component != "Non-module genes"].sort_values("delta")
    rows = list(mods.itertuples()) + [G[G.component == "Non-module genes"].iloc[0]]
    labels = [LABEL.get(r.component, r.component) for r in rows] + ["Composite (net)"]
    vals = [r.delta for r in rows] + [summ.delta_composite]
    lo = [r.lo for r in rows] + [summ.comp_lo]
    hi = [r.hi for r in rows] + [summ.comp_hi]
    fdr = [r.fdr for r in rows] + [np.nan]
    y = np.arange(len(vals))[::-1]
    colors = [(RED if v > 0 else BLUE) if (f < 0.10) else "#c3c8cd" for v, f in zip(vals, fdr)]
    colors[-1] = ORANGE
    colors[-2] = "#7d8790"
    ax.barh(y, vals, color=colors, height=0.7)
    ax.errorbar(vals, y, xerr=[np.array(vals) - np.array(lo), np.array(hi) - np.array(vals)], fmt="none",
                ecolor=DARK, elinewidth=0.5, capsize=1.2)
    ax.axvline(0, color=DARK, lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.9)
    ax.set_xlabel("Contribution to composite IPF − control\n(Δlog10 hazard; 95% donor-bootstrap CI)")
    ax.set_xlim(-0.17, 0.17)
    ax.legend(handles=[Patch(color=RED, label="higher in IPF, FDR < 0.10"), Patch(color=BLUE, label="lower in IPF, FDR < 0.10"),
                       Patch(color="#c3c8cd", label="FDR ≥ 0.10")], loc="upper center", fontsize=5.4,
              bbox_to_anchor=(0.45, -0.15), ncol=2, columnspacing=0.8)
    ax.set_title("Composite clock decomposed into\nits own module-gene components", fontsize=7)
    panel(ax, "a", x=-0.78)

    ax = fig.add_subplot(gs[0, 1])
    fmt = (lambda p: "p<0.001" if p < 0.001 else f"p={p:.2f}")
    items = [(summ.S, summ.S_lo, summ.S_hi, summ.S_null_median), (summ.N, np.nan, np.nan, np.nan),
             (summ.MI, summ.MI_lo, summ.MI_hi, summ.MI_null_median)]
    for i, (v, l, hh, nul) in enumerate(items):
        ax.errorbar([i], [v], yerr=[[v - l], [hh - v]] if np.isfinite(l) else None, fmt="o", color=DARK, ms=3.5,
                    elinewidth=0.8, capsize=2)
        if np.isfinite(nul):
            ax.scatter([i + 0.22], [nul], marker="D", s=12, color=GREY, zorder=3)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels([f"S\n{fmt(summ.p_S)}", f"N\n{fmt(summ.p_N)}", f"MI\n{fmt(summ.p_MI)}"], fontsize=5.0)
    ax.set_xlim(-0.5, 2.6)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("S, N (Δlog10 hazard);\nMI (unitless)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=DARK, ms=3.5, label="observed, 95% CI"),
                       Line2D([], [], marker="D", ls="", color=GREY, ms=3, label="null median")],
              loc="upper left", fontsize=5.2, bbox_to_anchor=(-0.04, 1.04))
    ax.set_title("Masking statistics (global)", fontsize=7)
    panel(ax, "b", x=-0.95)

    ax = fig.add_subplot(gs[1, 1])
    C = read("rev07_masking_celltype.csv").set_index("celltype").reindex(CT_ORDER)
    y = np.arange(len(C))[::-1]
    ax.hlines(y, C.S_lo, C.S_hi, color="#c3c8cd", lw=1.2)
    ax.scatter(C.S, y, s=12, color=[RED if f < .1 else DARK for f in C.fdr_S], zorder=3)
    ax.scatter(C.S_null_median, y, marker="D", s=8, color=GREY, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([CT_LABEL[c] for c in C.index], fontsize=5.4)
    ax.set_xlabel("Total component signal S (Δlog10 hazard)")
    ax.set_title("Per cell type\n(red: FDR < 0.10 vs null)", fontsize=7)
    panel(ax, "c", x=-0.95)

    ax = fig.add_subplot(gs[:, 2])
    g = read("rev08_global_stats.csv").set_index("module")
    order = MODULE_ORDER + ["Composite"]
    g = g.reindex(order)
    y = np.arange(len(g))[::-1]
    sigc = (g.fdr_cr1_normal < .1).values
    sigm = (g.fdr_mixed < .1).values
    sigw = (g.fdr_wild_cluster < .1).values
    ax.axvline(0, color=GREY, lw=0.5)
    ax.hlines(y, g.ci_lo_mixed, g.ci_hi_mixed, color="#c3c8cd", lw=1.0)
    ax.scatter(g.beta_mixed, y, s=11, color=[RED if w else (ORANGE if (c or mm) else DARK)
                                             for c, mm, w in zip(sigc, sigm, sigw)], zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[k] for k in order], fontsize=5.6)
    ax.set_xlabel("Global IPF − control, mixed model\n(Δlog10 hazard; 95% CI)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=RED, ms=3.5, label="FDR < 0.10: CR1, mixed, wild bootstrap"),
                       Line2D([], [], marker="o", ls="", color=ORANGE, ms=3.5, label="FDR < 0.10: CR1 or mixed only")],
              loc="lower left", fontsize=5.0, bbox_to_anchor=(-0.05, 1.0))
    ax.set_title("Separately trained module clocks", fontsize=7, pad=22)
    panel(ax, "d", x=-1.0, y=1.06)
    save(fig, "Fig3")


# ============================================================================ Fig. 4
def fig4():
    st = read("rev07_standardized_effects.csv")
    stats = read("rev08_celltype_stats.csv")
    piv = st.pivot(index="celltype", columns="score", values="std_effect").reindex(CT_ORDER)[MODULE_ORDER + ["Composite"]]
    fdr = stats.pivot(index="celltype", columns="module", values="fdr_IPF").reindex(CT_ORDER)[MODULE_ORDER]
    rob = stats.assign(rob=(stats.fdr_IPF < .1) & (stats.fdr_IPF_hc3 < .1) & (stats.fdr_IPF_perm < .1)) \
        .pivot(index="celltype", columns="module", values="rob").reindex(CT_ORDER)[MODULE_ORDER]
    fig = plt.figure(figsize=(7.1, 7.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.25], hspace=0.62)
    ax = fig.add_subplot(gs[0])
    v = 2.0
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    for i in range(fdr.shape[0]):
        for j in range(fdr.shape[1]):
            if fdr.iat[i, j] < 0.10:
                ax.scatter(j, i, s=14, marker="o", facecolor=DARK if rob.iat[i, j] else "white",
                           edgecolor="white" if rob.iat[i, j] else DARK, lw=0.6)
    ax.axvline(len(MODULE_ORDER) - 0.5, color="white", lw=2)
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([LABEL[c] for c in piv.columns], rotation=55, ha="right", fontsize=5.7)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([CT_LABEL[c] for c in piv.index], fontsize=6)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("Standardized IPF − control effect\n(β / residual SD; age, sex adjusted)", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=DARK, mec=DARK, ms=3.5,
                              label="FDR < 0.10 by classical, HC3 and permutation tests"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec=DARK, ms=3.5,
                              label="FDR < 0.10 by classical test only")],
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.8)
    panel(ax, "a", x=-0.13, y=1.08)

    ax = fig.add_subplot(gs[1])
    h = stats[stats.fdr_IPF < 0.10].copy()
    h["ct_rank"] = h.celltype.map({c: i for i, c in enumerate(CT_ORDER)})
    h = h.sort_values(["ct_rank", "beta_IPF"])
    q = 1.96
    x = np.arange(len(h))
    ax.axhline(0, color=GREY, lw=0.5)
    ax.errorbar(x - 0.2, h.beta_IPF, yerr=[h.beta_IPF - h.ci_lo, h.ci_hi - h.beta_IPF], fmt="o", ms=2.8,
                color=RED, elinewidth=0.6, label="IPF − control")
    ax.errorbar(x, h.beta_COPD, fmt="s", ms=2.6, color="#7d8790", label="COPD − control")
    spec = (h.fdr_IPF_vs_COPD < .1).values
    ax.errorbar(x + 0.2, h.IPF_minus_COPD, yerr=q * h.se_IPF_minus_COPD, fmt="none", ecolor=BLUE, elinewidth=0.6)
    ax.scatter(x + 0.2, h.IPF_minus_COPD, marker="D", s=9, facecolor=[BLUE if s_ else "white" for s_ in spec],
               edgecolor=BLUE, lw=0.6, zorder=3, label="IPF − COPD (filled: FDR < 0.10)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{CT_LABEL[c]} · {LABEL[m_]}" for c, m_ in zip(h.celltype, h.module)], rotation=65,
                       ha="right", fontsize=5.2)
    ax.set_ylabel("Effect (Δlog10 hazard, 95% CI)")
    ax.set_xlim(-0.8, len(h) - 0.2)
    ax.legend(loc="upper left", ncol=3, fontsize=5.8, bbox_to_anchor=(0, 1.1))
    ax.set_title("The 30 FDR < 0.10 IPF − control effects with the direct IPF − COPD contrast", fontsize=7, pad=14)
    panel(ax, "b", x=-0.08, y=1.1)
    save(fig, "Fig4")


# ============================================================================ Fig. 5
def fig5():
    from scipy.stats import spearmanr
    fig = plt.figure(figsize=(7.1, 5.4))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.5, height_ratios=[1, 0.95])
    ax = fig.add_subplot(gs[0, 0])
    hc = read("rev06_hits_competitive.csv")
    xr, yp = -np.log10(hc.p_comp_random), -np.log10(hc.p_comp_perm)
    spec = hc.module_specific.astype(str).str.lower().eq("true").values
    ax.scatter(xr, yp, s=12, c=[RED if s_ else "#7d8790" for s_ in spec], lw=0, zorder=3)
    lab_pos = {0: (-4, 0, "right"), 1: (-4, 0, "right"), 2: (4, 2, "left"), 3: (-4, 3, "right"), 4: (4, -4, "left")}
    k = 0
    for xv, yv, r_, s_ in sorted(zip(xr, yp, hc.itertuples(), spec), key=lambda t: -t[1]):
        if s_:
            dx, dy, ha = lab_pos[k]
            sep = " · " if (ha == "left" or k == 0) else " ·\n"
            ax.annotate(f"{CT_SHORT[r_.celltype]}{sep}{LABEL[r_.module].split('/')[0]}", (xv, yv), xytext=(dx, dy),
                        textcoords="offset points", fontsize=4.6, color=RED, va="center", ha=ha)
            k += 1
    ax.axvline(-np.log10(.05), color=GREY, lw=0.6, ls="--")
    ax.axhline(-np.log10(.05), color=GREY, lw=0.6, ls="--")
    ax.set_xlabel("−log10 p vs size-matched random genes")
    ax.set_ylabel("−log10 p vs permuted coefficients")
    nr, npm = int((hc.p_comp_random < .05).sum()), int((hc.p_comp_perm < .05).sum())
    ax.set_title(f"Competitive nulls, 30 effects\n(p < 0.05: random genes {nr}, permuted {npm})", fontsize=6.5)
    ax.set_xlim(0, 4.2)
    ax.set_ylim(0, 3.9)
    panel(ax, "a", x=-0.3)

    ax = fig.add_subplot(gs[0, 1])
    r = read("rev06_rodent_geneset_cv.csv")
    for oc, col_, mk in [("Mortality", RED, "o"), ("Chrono", BLUE, "s")]:
        d = r[r.outcome == oc]
        ax.scatter(d.r_oof_geneset_signed, d.r_oof_clock, s=9, color=col_, marker=mk, lw=0,
                   label="mortality" if oc == "Mortality" else "chronological")
    ax.plot([0, .65], [0, .65], color=GREY, lw=0.6, ls="--")
    ax.set_xlabel("Signed mean z of module genes (r)")
    ax.set_ylabel("Module clock (r)")
    ax.set_title(f"Rodent CV: clock vs gene-set score\n(median Δr = {r.clock_minus_geneset.median():+.2f})", fontsize=6.5)
    ax.legend(fontsize=5.4, loc="lower right")
    panel(ax, "b", x=-0.3)

    ax = fig.add_subplot(gs[0, 2])
    b = read("rev06_ipf_baselines.csv")
    rs, ru = spearmanr(b.t_geneset_signed, b.t_clock)[0], spearmanr(b.t_geneset, b.t_clock)[0]
    ax.scatter(b.t_geneset_signed, b.t_clock, s=4, color=BLUE, lw=0, alpha=0.7, label=f"signed (ρ = {rs:.2f})")
    ax.scatter(b.t_geneset, b.t_clock, s=4, color=ORANGE, lw=0, alpha=0.6, label=f"unsigned (ρ = {ru:.2f})")
    ax.axhline(0, color=GREY, lw=0.4)
    ax.axvline(0, color=GREY, lw=0.4)
    ax.set_xlabel("Gene-set score t (IPF − control)")
    ax.set_ylabel("Module clock t")
    ax.set_title("IPF, 253 cell-type × module tests", fontsize=6.5)
    ax.legend(fontsize=5.2, loc="upper left", handletextpad=0.1)
    panel(ax, "c", x=-0.3)

    ax = fig.add_subplot(gs[1, :])
    s = read("rev09_sensitivity_summary.csv").set_index("analysis")
    order = [("age_overlap", "controls ≥ 54 y only"), ("tech_covariates", "+ cells, depth, detected fraction"),
             ("cpm", "CPM instead of RLE"), ("min_cells_30", "≥ 30 cells per pseudobulk"),
             ("min_cells_100", "≥ 100 cells per pseudobulk"), ("depth_group_matched", "depth matched between groups"),
             ("cells_group_matched", "cell number matched between groups"),
             ("depth_equalized", "depth: all thinned to minimum"), ("depth_loss_control", "depth: same loss, groups unmatched"),
             ("cells_equalized", "cells: one 60-cell sample per donor"), ("cells_loss_control", "cells: same loss, groups unmatched")]
    y = np.arange(len(order))[::-1]
    avail = s.loc[[k for k, _ in order], "main_hits_available"].values
    same = s.loc[[k for k, _ in order], "main_hits_same_sign"].values
    p05 = s.loc[[k for k, _ in order], "main_hits_p05"].values
    ratio = s.loc[[k for k, _ in order], "main_hits_median_beta_ratio"].values
    ax.barh(y, same, color="#d6e4f2", height=0.7, label="same sign as main analysis")
    ax.barh(y, p05, color=BLUE, height=0.7, label="same sign and p < 0.05")
    for yy, r_ in zip(y, ratio):
        ax.text(30.8, yy, f"{r_:.2f}", va="center", fontsize=5.6)
    ax.text(30.8, len(order) - 0.55, "median β ratio\n(analysis / main)", fontsize=5.4, va="bottom")
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab in order], fontsize=6)
    ax.set_xlim(0, 36)
    ax.set_xticks([0, 5, 10, 15, 20, 25, 30])
    ax.set_xlabel("Number of the 30 main FDR < 0.10 effects")
    ax.legend(loc="lower left", fontsize=5.6, bbox_to_anchor=(0.33, 1.0), ncol=2)
    ax.set_title("Sensitivity analyses", fontsize=7, loc="left")
    panel(ax, "d", x=-0.2)
    save(fig, "Fig5")


# ============================================================================ Fig. 6
def coarse_fractions():
    from rev_10_replication import COARSE_A, COARSE_H
    out = []
    for name, f, ctc, cmap in [("GSE136831\n(whole-lung\ndissociation)", "pseudobulk_meta.csv", "Manuscript_Identity", COARSE_A),
                               ("GSE135893\n(three sites)", "hab_pseudobulk_meta.csv", "celltype", COARSE_H)]:
        m = pd.read_csv(rc.RES / f)
        m = m[m.Disease_Identity.isin(["IPF", "Control"])].copy()
        m["coarse"] = m[ctc].map(cmap).fillna("Other")
        tot = m.groupby("Subject_Identity").n_cells.sum()
        mac = m[m.coarse == "Macrophage"].groupby("Subject_Identity").n_cells.sum().reindex(tot.index).fillna(0)
        dis = m.groupby("Subject_Identity").Disease_Identity.first()
        out.append(pd.DataFrame({"cohort": name, "disease": dis, "frac": mac / tot}))
    return pd.concat(out)


def fig6():
    fig = plt.figure(figsize=(7.1, 5.6))
    top = fig.add_gridspec(1, 3, wspace=0.55, width_ratios=[1.6, 0.8, 0.85], top=0.93, bottom=0.6)
    bot = fig.add_gridspec(1, 2, wspace=1.0, width_ratios=[1.3, 1], top=0.43, bottom=0.07)
    ax = fig.add_subplot(top[0, 0])
    C = read("rev15_checklist_realdata.csv")
    order = ["NK", "T cell", "Monocyte", "Macrophage", "Ciliated"]
    C = C.set_index("coarse").reindex(order)
    x = np.arange(len(C))
    ax.bar(x - 0.18, C.r_cross, width=0.34, color=[RED if (p < .05 and r > 0) else "#c3c8cd"
                                                  for p, r in zip(C.p_cross, C.r_cross)], label="cross-cohort r")
    ax.bar(x + 0.18, C.r_splithalf, width=0.34, color="#d6e4f2", edgecolor=BLUE, lw=0.5,
           label="split-half r (GSE136831)")
    ax.scatter(x + 0.18, C.tau_splithalf, marker="_", s=80, color=BLUE, lw=1.2, zorder=3,
               label="split-half null, 95th pct")
    for xi, r_ in zip(x, C.itertuples()):
        ax.text(xi - 0.18, (r_.r_cross + 0.03) if r_.r_cross > 0 else (r_.r_cross - 0.03), f"p={r_.p_cross:.2f}",
                ha="center", va="bottom" if r_.r_cross > 0 else "top", fontsize=4.8)
    short = {"replicated": "replicated", "inconclusive (power)": "low\nreliability",
             "between-cohort difference (unresolved)": "between-\ncohort", "composition": "composition"}
    ax.set_xticks(x)
    ax.set_xticklabels(C.index, fontsize=5.3)
    for xi, a in zip(x, C.attribution):
        ax.text(xi, -0.98, short[a], ha="center", va="top", fontsize=4.9, style="italic",
                color=RED if a == "replicated" else DARK)
    ax.axhline(0, color=DARK, lw=0.5)
    ax.set_ylim(-0.75, 1.0)
    ax.set_ylabel("Pearson r of 23-module\neffect vectors")
    ax.legend(fontsize=5.0, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.0), columnspacing=0.8)
    ax.set_title("Replication and checklist outcome", fontsize=7, pad=24)
    panel(ax, "a", x=-0.3, y=1.12)

    ax = fig.add_subplot(top[0, 1])
    F = coarse_fractions()
    cohorts = F.cohort.unique()
    pos = []
    for i, co in enumerate(cohorts):
        for j, dis in enumerate(["Control", "IPF"]):
            d = F[(F.cohort == co) & (F.disease == dis)].frac
            p_ = i * 2.6 + j
            pos.append(p_)
            ax.boxplot([d], positions=[p_], widths=0.6, showfliers=False, patch_artist=True,
                       boxprops=dict(fc="#eceff1" if dis == "Control" else "#f4d6d2", lw=0.5),
                       medianprops=dict(color=DARK, lw=0.9), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
            ax.scatter(p_ + np.random.default_rng(i * 2 + j).uniform(-0.15, 0.15, len(d)), d, s=4,
                       color=GREY if dis == "Control" else RED, lw=0, zorder=3)
    ax.set_xticks([0.5, 3.1])
    ax.set_xticklabels(cohorts, fontsize=5.4)
    ax.set_ylabel("Macrophage share of all cells")
    ax.set_ylim(0, 1.05)
    ax.legend(handles=[Patch(fc="#eceff1", ec=DARK, lw=0.4, label="control"), Patch(fc="#f4d6d2", ec=DARK, lw=0.4, label="IPF")],
              fontsize=5.2, loc="upper right", bbox_to_anchor=(1.05, 1.03))
    ax.set_title("Composition differs\nbetween cohorts", fontsize=7)
    panel(ax, "b", x=-0.5, y=1.12)

    ax = fig.add_subplot(top[0, 2])
    bc = read("rev10_bulk_concordance.csv")
    lab = ["136831 vs\n135893", "136831 vs\n134692", "135893 vs\n134692"]
    x = np.arange(len(bc))
    ax.bar(x, bc.r, color=["#c3c8cd" if p >= .05 else RED for p in bc.p_two_sided], width=0.55)
    for xi, r_ in zip(x, bc.itertuples()):
        ax.text(xi, r_.r + 0.03, f"p={r_.p_two_sided:.3f}", ha="center", fontsize=5.0)
    ax.set_xticks(x)
    ax.set_xticklabels(lab, fontsize=5.2)
    ax.set_ylim(0, 0.75)
    ax.set_xlabel("GSE pairs (whole-donor or bulk)", fontsize=6)
    ax.set_ylabel("Pearson r")
    ax.set_title("Whole-tissue concordance\n(two-sided permutation p)", fontsize=7)
    panel(ax, "c", x=-0.42, y=1.12)

    ax = fig.add_subplot(bot[0, 1])
    I = read("rev11_interaction.csv")
    I = I[I.celltype == "Macrophage"].set_index("module").reindex(MODULE_ORDER)
    y = np.arange(len(I))[::-1]
    sig = (I.fdr_interaction < .1).values
    ax.hlines(y, I.beta_IPF_resHi, I.beta_IPF_monoHi, color=["#7d8790" if s_ else "#dde1e4" for s_ in sig], lw=0.9)
    ax.scatter(I.beta_IPF_resHi, y, s=9, color=BLUE, zorder=3, label="resident-marker-high")
    ax.scatter(I.beta_IPF_monoHi, y, s=9, color=ORANGE, zorder=3, label="monocyte-marker-high")
    for yy, s_ in zip(y, sig):
        if s_:
            ax.text(0.95, yy, "*", fontsize=7, va="center", ha="center")
    ax.axvline(0, color=GREY, lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[m_] for m_ in I.index], fontsize=5.0)
    ax.set_xlim(-1.0, 1.05)
    ax.set_xlabel("IPF − control within sub-state (Δlog10 hazard)")
    ax.legend(fontsize=5.0, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, columnspacing=0.6,
              handletextpad=0.1)
    ax.set_title("Macrophage sub-state × disease (*: interaction FDR < 0.10)", fontsize=6.5, pad=14)
    panel(ax, "e", x=-1.05, y=1.07)

    ax = fig.add_subplot(bot[0, 0])
    M = read("rev12_checklist_confusion.csv").set_index("scenario")
    cols = ["replicated", "inconclusive (power)", "composition", "between-cohort difference (unresolved)"]
    rows = [("shared", "shared effect"), ("power", "weak shared effect"), ("composition", "sub-state composition"),
            ("composition_shift", "composition shift only"), ("confound", "cohort-specific confound"),
            ("null", "no disease effect")]
    mat = M.loc[[k for k, _ in rows], cols].values.astype(float)
    im = ax.imshow(mat, cmap="Greys", vmin=0, vmax=1, aspect="auto")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6,
                    color="white" if mat[i, j] > 0.55 else DARK)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(["replicated", "inconclusive\n(low reliability)", "composition", "between-cohort\ndifference"],
                       fontsize=5.8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for _, lab in rows], fontsize=6)
    ax.set_xlabel("Checklist attribution")
    ax.set_ylabel("Simulated truth")
    ax.set_title("Checklist operating characteristics\n(300 simulated cohort pairs per scenario)", fontsize=6.8)
    panel(ax, "d", x=-0.38, y=1.07)
    save(fig, "Fig6")


# ============================================================================ Supplementary
def figS1():
    fig = plt.figure(figsize=(7.1, 4.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.04], wspace=0.08, left=0.2, right=0.93, top=0.86, bottom=0.25)
    axs = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
    mf = read("rev03_meff.csv")
    titles = ["Rodent training data\n(in-sample predictions)", "IPF pseudobulks\n(centred within cell type)"]
    for i, (ax, f) in enumerate(zip(axs, ["rev03_corr_training.csv", "rev03_corr_ipf.csv"])):
        c = pd.read_csv(O / f, index_col=0).reindex(index=MODULE_ORDER, columns=MODULE_ORDER)
        im = ax.imshow(c.values, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(23))
        ax.set_yticks(range(23))
        ax.set_xticklabels([LABEL[k] for k in MODULE_ORDER], rotation=75, ha="right", fontsize=4.6)
        ax.set_yticklabels([LABEL[k] for k in MODULE_ORDER] if i == 0 else [], fontsize=4.6)
        r_ = mf.iloc[i]
        ax.set_title(f"{titles[i]}\nmedian |ρ| = {r_.median_abs_r:.2f}; effective modules "
                     f"{r_.meff_liji:.0f}–{r_.meff_nyholt:.0f}", fontsize=6.2)
        ax.spines[:].set_visible(False)
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[2]))
    cb.set_label("Spearman ρ", fontsize=6)
    panel(axs[0], "a", x=-0.08, y=1.2)
    panel(axs[1], "b", x=-0.08, y=1.2)
    save(fig, "FigS1")


def figS2():
    D = read("rev15_detected_by_module.csv")
    P = D.pivot(index="module", columns="celltype", values="frac").reindex(index=MODULE_ORDER, columns=CT_ORDER)
    F = read("rev15_detected_features.csv")
    fig = plt.figure(figsize=(7.1, 3.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1], wspace=0.75, top=0.78, bottom=0.2)
    ax = fig.add_subplot(gs[0])
    im = ax.imshow(P.values, cmap="viridis", vmin=0.4, vmax=1, aspect="auto")
    ax.set_yticks(range(len(P)))
    ax.set_yticklabels([LABEL[k] for k in P.index], fontsize=5.4)
    ax.set_xticks(range(P.shape[1]))
    ax.set_xticklabels([CT_LABEL[c] for c in P.columns], rotation=55, ha="right", fontsize=5.6)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Fraction of module genes detected", fontsize=6)
    panel(ax, "a", x=-0.62, y=1.3)
    ax = fig.add_subplot(gs[1])
    F2 = F[F.celltype != "detected in all cell types"].set_index("celltype").reindex(CT_ORDER)
    y = np.arange(len(F2))[::-1]
    ax.barh(y + 0.18, F2.module_frac, height=0.34, color=BLUE, label="module-clock genes (1,995)")
    ax.barh(y - 0.18, F2.composite_frac, height=0.34, color=ORANGE, label="composite-clock genes (10,487)")
    allr = F[F.celltype == "detected in all cell types"].iloc[0]
    ax.axvline(allr.module_frac, color=BLUE, ls="--", lw=0.7,
               label=f"module genes detected in all 11 ({int(allr.module_genes_detected):,})")
    ax.axvline(allr.composite_frac, color=ORANGE, ls="--", lw=0.7,
               label=f"composite genes detected in all 11 ({int(allr.composite_genes_detected):,})")
    ax.set_yticks(y)
    ax.set_yticklabels([CT_LABEL[c] for c in F2.index], fontsize=5.8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction detected (not imputed)")
    ax.legend(fontsize=5.2, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=1)
    panel(ax, "b", x=-0.6, y=1.3)
    save(fig, "FigS2")


def figS3():
    R = read("rev14_published_reproduction.csv")
    Cf = read("rev14_coef_agreement.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.9), gridspec_kw=dict(wspace=0.45))
    ax = axs[0]
    r = R[(R.table == "5B module") & (R.module != "All module genes")]
    conv = {"as_released": "as released", "gene_standardized": "gene-standardized",
            "within_module_sample_z": "within-module z"}
    for i, (k, lab) in enumerate(conv.items()):
        d = r[r.convention == k]
        ax.scatter(np.full(len(d), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(d)), d.r_insample, s=6,
                   color=BLUE, lw=0)
    rep = r[r.convention == "as_released"].r_reported.dropna()
    ax.scatter(np.full(len(rep), 3) + np.random.default_rng(9).uniform(-0.15, 0.15, len(rep)), rep, s=6, color=GREY, lw=0)
    comp = R[(R.table == "5A composite") & (R.convention == "as_released")]
    ax.scatter([4, 4.2], comp.r_insample, marker="*", s=40, color=ORANGE)
    ax.set_xticks(range(5))
    ax.set_xticklabels(["published\nmodule coef.,\nas released", "gene-\nstandardized", "within-\nmodule z",
                        "reported\n(Supp. Table 4B)", "published\ncomposite\ncoef."], fontsize=5.4)
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_ylabel("In-sample Pearson r on the\nreleased training matrix")
    ax.set_title("Published coefficients do not reproduce module-clock accuracy", fontsize=6.5)
    panel(ax, "a", x=-0.25)
    ax = axs[1]
    for oc, col_, mk in [("Mortality", RED, "o"), ("Chrono", BLUE, "s")]:
        d = Cf[Cf.outcome == oc].merge(R[(R.table == "5B module") & (R.convention == "as_released") & (R.outcome == oc)]
                                       [["module", "r_reported"]], left_on="color", right_on="module", how="left")
        ax.scatter(d.r_reported, d.r_cv_retrained, s=9, color=col_, marker=mk, lw=0,
                   label="mortality" if oc == "Mortality" else "chronological")
    ax.plot([0, .9], [0, .9], color=GREY, ls="--", lw=0.6)
    ax.set_xlabel("Reported r of published module clock")
    ax.set_ylabel("Retrained clock, dataset-grouped CV r")
    ax.legend(fontsize=5.6)
    ax.set_title("Retrained clocks reach the reported accuracy", fontsize=6.5)
    panel(ax, "b", x=-0.25)
    save(fig, "FigS3")


def figS4():
    B = read("rev10_bulk_direction.csv")
    D = read("rev10_bulk_decomposition.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.7), gridspec_kw=dict(width_ratios=[1.6, 1], wspace=0.4))
    ax = axs[0]
    mods = B.module.unique()
    cohorts = B.cohort.unique()
    w = 0.25
    for i, co in enumerate(cohorts):
        d = B[B.cohort == co].set_index("module").reindex(mods)
        ax.bar(np.arange(len(mods)) + (i - 1) * w, d.delta, width=w, label=co,
               color=[RED, ORANGE, BLUE][i], alpha=0.85)
    g = B.drop_duplicates("module").set_index("module").reindex(mods)
    ax.scatter(np.arange(len(mods)), g.beta_global_A, marker="_", s=120, color=DARK, lw=1.4, zorder=4,
               label="GSE136831 cell-type-resolved global effect")
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_xticks(range(len(mods)))
    ax.set_xticklabels([LABEL[m_] for m_ in mods], rotation=25, ha="right", fontsize=5.6)
    ax.set_ylabel("IPF − control (Δlog10 hazard)")
    ax.legend(fontsize=5.2, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Direction of the global effects in whole-donor and bulk data", fontsize=6.5, pad=26)
    panel(ax, "a", x=-0.12, y=1.2)
    ax = axs[1]
    x = np.arange(2)
    for i, co in enumerate(["GSE136831", "GSE135893"]):
        d = D[D.cohort == co].set_index("component")
        ax.bar(x + (i - 0.5) * 0.35, d.loc[["composition_only", "expression_only"], "r_with_observed"], width=0.35,
               color=[RED, ORANGE][i], label=co)
    ax.set_xticks(x)
    ax.set_xticklabels(["composition only\n(control profiles,\nobserved proportions)",
                        "expression only\n(observed profiles,\ncontrol proportions)"], fontsize=5.5)
    ax.set_ylabel("r with observed whole-donor effects")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=5.4, loc="upper left")
    ax.set_title("Whole-donor signal decomposition", fontsize=6.5, pad=26)
    panel(ax, "b", x=-0.3, y=1.2)
    save(fig, "FigS4")


def figS5():
    s1 = pd.read_csv(MIL / "stage1_results.csv")
    auc = read("rev13_mil_auc.csv")
    wl = pd.read_csv(MIL / "stage3_within_label_raw.csv")
    fig, axs = plt.subplots(1, 4, figsize=(7.1, 2.4), gridspec_kw=dict(wspace=0.72))
    ax = axs[0]
    for pool, col_ in [("attention", RED), ("mean", GREY)]:
        d = s1[(s1.regime == "substate") & (s1.pool == pool)].groupby("p_s").bag_auc.agg(["mean", "std"])
        ax.errorbar(d.index, d["mean"], yerr=d["std"], color=col_, marker="o", ms=2.5, lw=0.8, capsize=1.5,
                    label=f"{pool} pooling")
    ax.set_xlabel("Sub-state fraction")
    ax.set_ylabel("Bag AUC (semi-synthetic)")
    ax.legend(fontsize=5.2)
    panel(ax, "a", x=-0.45)
    ax = axs[1]
    for reg, col_ in [("substate", RED), ("diffuse", BLUE)]:
        d = s1[(s1.regime == reg) & (s1.pool == "attention")].groupby("p_s").att_auroc.agg(["mean", "std"])
        ax.errorbar(d.index, d["mean"], yerr=d["std"], color=col_, marker="o", ms=2.5, lw=0.8, capsize=1.5, label=reg)
    ax.set_xlabel("Signal-bearing fraction")
    ax.set_ylabel("Attention AUROC for\nsignal-bearing instances")
    ax.legend(fontsize=5.2)
    panel(ax, "b", x=-0.45)
    ax = axs[2]
    order = ["attention", "mean", "logistic"]
    ax.boxplot([auc[auc.model == k].auc for k in order], widths=0.5, showfliers=False,
               medianprops=dict(color=DARK), boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    for i, k in enumerate(order):
        d = auc[auc.model == k].auc
        ax.scatter(np.full(len(d), i + 1) + np.random.default_rng(i).uniform(-0.12, 0.12, len(d)), d, s=5, color=RED, lw=0)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["attention\nMIL", "mean-pool\nMIL", "logistic,\ndonor\nmeans"], fontsize=5.2)
    ax.set_ylabel("Out-of-fold donor AUC\n(IPF vs control; 10 seeds)")
    ax.set_ylim(0.85, 1.0)
    panel(ax, "c", x=-0.45)
    ax = axs[3]
    w = wl[wl.group == "IPF"]
    vals = [w.sep_pooled, w.sep_within_Macrophage, w.sep_within_Macrophage_Alveolar]
    ax.boxplot(vals, widths=0.5, showfliers=False, medianprops=dict(color=DARK), boxprops=dict(lw=0.5),
               whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["pooled", "within\nmacro-\nphage", "within\nalveolar\nmacro."], fontsize=5.2)
    ax.set_ylabel("Attention separation\n(monocyte − resident markers)")
    ax.set_title(f"IPF donors; cell-type\nannotation AUROC {w.auroc_annotation.mean():.2f}", fontsize=6)
    panel(ax, "d", x=-0.45)
    save(fig, "FigS5")


def figS6():
    HL = {"ATI": "AT1", "ATII": "AT2", "B": "B cell", "Ciliated": "Ciliated", "DC": "Dendritic cell",
          "Endothelial": "Endothelial", "Macrophage": "Macrophage", "Monocyte": "Monocyte", "NK": "NK cell",
          "Stromal": "Stromal", "T": "T cell"}
    cols = ["Composite"] + MODULE_ORDER
    fig = plt.figure(figsize=(7.1, 7.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.85, 1.25], width_ratios=[1, 0.025], hspace=0.95, wspace=0.03,
                          left=0.17, right=0.92, top=0.95, bottom=0.06)
    specs = [("rev05_hlca_age.csv", "HLCA healthy lung, within cell type (11 cell types, 14–38 donors each, 20–81 y)",
              HL, "a"),
             ("rev05_adams_control_age.csv", "GSE136831 control donors, within cell type (28 donors, 20–80 y)",
              CT_LABEL, "b")]
    for k, (f, title, names, letter) in enumerate(specs):
        T = read(f)
        T = T[T.outcome == "Mortality"]
        M = T.pivot_table(index="celltype", columns="module", values="rho").reindex(columns=cols)
        Fq = T.pivot_table(index="celltype", columns="module", values="fdr").reindex(index=M.index, columns=cols)
        if letter == "b":
            M = M.reindex([c for c in CT_ORDER if c in M.index])
            Fq = Fq.reindex(M.index)
        ax = fig.add_subplot(gs[k, 0])
        im = ax.imshow(M.values, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if Fq.values[i, j] < 0.10:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=7)
        ax.axvline(0.5, color=DARK, lw=0.8)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([LABEL[c] for c in cols], rotation=60, ha="right", fontsize=5.0)
        ax.set_yticks(range(M.shape[0]))
        ax.set_yticklabels([names.get(c, c) for c in M.index], fontsize=5.6)
        ax.spines[:].set_visible(False)
        ax.set_title(title, fontsize=6.6, loc="left")
        panel(ax, letter, x=-0.16, y=1.02)
        cb = fig.colorbar(im, cax=fig.add_subplot(gs[k, 1]))
        cb.set_label("Spearman ρ with age", fontsize=5.6)
    G = read("rev05b_gtex_validation.csv")
    G = G[G.outcome == "Mortality"].set_index("module").reindex(cols)
    ax = fig.add_subplot(gs[2, 0])
    y = np.arange(len(G))[::-1]
    ax.barh(y, G.adj_spearman, color=[ORANGE if m_ == "Composite" else ("#c3c8cd" if q >= 0.05 else BLUE)
                                     for m_, q in zip(G.index, G.adj_fdr.fillna(0))], height=0.7)
    ax.axvline(0, color=DARK, lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[c] for c in G.index], fontsize=5.2)
    ax.set_xlabel("Spearman ρ with age, adjusted for ischaemic time, RIN, Hardy scale and sex")
    c = G.loc["Composite"]
    ax.set_title(f"GTEx v8 lung bulk tissue (n = {int(c.n)}): composite ρ = {c.adj_spearman:.2f} "
                 f"(p = {c.adj_p:.1e}); no module clock at FDR < 0.10", fontsize=6.6, loc="left")
    panel(ax, "c", x=-0.16, y=1.02)
    save(fig, "FigS6")


if __name__ == "__main__":
    for f in [fig1, fig2, fig3, fig4, fig5, fig6, figS1, figS2, figS3, figS4, figS5, figS6]:
        f()
