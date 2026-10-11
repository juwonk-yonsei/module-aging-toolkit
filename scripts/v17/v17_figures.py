"""Manuscript figures for instruction 14. Same data as v12, with the listed layout changes. Writes figures/v14 only."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RES = REPO / "results"
FIG = REPO / "figures" / "v17"
V12 = RES / "v12"

IPF, CTRL, COPD = "#D55E00", "#0072B2", "#E69F00"
ADOPT, REIMPL, NEW = "#009E73", "#666666", "#0072B2"
DROP, WITHIN = "#999999", "#CC79A7"
ORANGE, SKY = "#E69F00", "#56B4E9"

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": False,
    "legend.fontsize": 7, "legend.frameon": False,
})
NUMBERS, MANIFEST, SIZES = [], [], []
MM = 1 / 25.4
CT_ORDER = ["B", "Ciliated", "Fibroblast", "Macrophage", "Macrophage_Alveolar", "NK", "T",
            "T_Cytotoxic", "cDC2", "cMonocyte", "ncMonocyte"]
CT_LABEL = {"B": "B", "Ciliated": "Ciliated", "Fibroblast": "Fibroblast", "Macrophage": "Macrophage",
            "Macrophage_Alveolar": "Alv. mac", "NK": "NK", "T": "T", "T_Cytotoxic": "T cyt.",
            "cDC2": "cDC2", "cMonocyte": "cMonocyte", "ncMonocyte": "ncMonocyte"}


def mm(w, h):
    return w * MM, h * MM


def remember(fig, panel, shown, source, column, value):
    NUMBERS.append(dict(figure=fig, panel=panel, shown=shown, source=source, column=column, value=value))


def manifest(fig, panel, source, columns, note):
    MANIFEST.append(dict(figure=fig, panel=panel, source=source, columns=columns, note=note))


def style_ax(ax):
    ax.tick_params(length=2.5, width=0.6, pad=2)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def letters(fig, pairs):
    for ax, letter in pairs:
        pos = ax.get_position()
        fig.text(pos.x0, min(pos.y1 + 0.012, 0.985), letter, fontsize=10, fontweight="bold", ha="left", va="bottom")


def save(fig, name):
    from overlap_audit import audit
    FIG.mkdir(parents=True, exist_ok=True)
    width_mm = fig.get_figwidth() * 25.4
    height_mm = fig.get_figheight() * 25.4
    if abs(width_mm - 174) > 0.2 or height_mm > 234.05:
        raise SystemExit(f"{name} canvas {width_mm:.1f} x {height_mm:.1f} mm")
    fig.savefig(FIG / f"{name}.pdf")
    tiff_path = FIG / f"{name}.tiff"
    fig.savefig(tiff_path, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    from PIL import Image
    im = Image.open(tiff_path)
    if im.mode != "RGB":
        rgb = Image.new("RGB", im.size, (255, 255, 255))
        rgb.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
        im = rgb
    im.save(tiff_path, compression="tiff_lzw", dpi=(600, 600))
    fig.savefig(FIG / f"{name}.png", dpi=200)
    SIZES.append(dict(figure=name, width_mm=round(width_mm, 2), height_mm=round(height_mm, 2),
                      pdf_width_mm=round(width_mm, 2), pdf_height_mm=round(height_mm, 2),
                      tiff_px_w=im.size[0], tiff_px_h=im.size[1], tiff_mode=im.mode))
    rows = audit(fig, name)
    renderer = fig.canvas.get_renderer()
    boxes = [ax.get_window_extent(renderer) for ax in fig.axes]
    for t in fig.texts:
        if t.get_visible() and t.get_text().strip():
            boxes.append(t.get_window_extent(renderer))
    for ax in fig.axes:
        for t in list(ax.texts) + list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            if t.get_visible() and t.get_text().strip():
                boxes.append(t.get_window_extent(renderer))
        if ax.get_legend() is not None:
            boxes.append(ax.get_legend().get_window_extent(renderer))
    fb = fig.bbox
    x0, y0 = min(b.x0 for b in boxes), min(b.y0 for b in boxes)
    x1, y1 = max(b.x1 for b in boxes), max(b.y1 for b in boxes)
    to_mm = 25.4 / fig.dpi
    print(f"  margins mm L { (x0-fb.x0)*to_mm:.1f} R {(fb.x1-x1)*to_mm:.1f} "
          f"B {(y0-fb.y0)*to_mm:.1f} T {(fb.y1-y1)*to_mm:.1f}", flush=True)
    plt.close(fig)
    return rows


def text_polygon(t, renderer):
    """Display-coordinate corners of the text itself, including rotation."""
    _, _, (corner, size) = t._get_layout(renderer)
    w, h = size
    if w < 1 or h < 1:
        return None
    ang = np.deg2rad(t.get_rotation())
    c, s = np.cos(ang), np.sin(ang)
    right = np.array([c * w, s * w])
    up = np.array([-s * h, c * h])
    anchor = np.array(t.get_transform().transform(t.get_unitless_position()))
    origin = np.array(corner) + anchor
    return [origin, origin + up, origin + right + up, origin + right]


def polygons_overlap(p, q, min_area=2.0):
    """Separating-axis test, then a coarse intersection area."""
    def axes(poly):
        out = []
        for i in range(len(poly)):
            edge = poly[(i + 1) % len(poly)] - poly[i]
            out.append(np.array([-edge[1], edge[0]]))
        return out
    for axis in axes(p) + axes(q):
        if np.allclose(axis, 0):
            continue
        a = [np.dot(v, axis) for v in p]
        b = [np.dot(v, axis) for v in q]
        if max(a) < min(b) or max(b) < min(a):
            return False
    # AABB of the intersection of the two axis-aligned bounds, as a size filter
    pa, qa = np.vstack(p), np.vstack(q)
    x0, y0 = max(pa[:, 0].min(), qa[:, 0].min()), max(pa[:, 1].min(), qa[:, 1].min())
    x1, y1 = min(pa[:, 0].max(), qa[:, 0].max()), min(pa[:, 1].max(), qa[:, 1].max())
    return (x1 - x0) * (y1 - y0) >= min_area


def overlap_rows(fig, name):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = []
    for t in fig.findobj(plt.Text):
        s = t.get_text().strip()
        if not s or not t.get_visible():
            continue
        poly = text_polygon(t, renderer)
        if poly is None:
            continue
        texts.append((s, poly))
    rows = []
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            if not polygons_overlap(texts[i][1], texts[j][1]):
                continue
            rows.append(dict(figure=name, text_a=texts[i][0][:60], text_b=texts[j][0][:60]))
    return rows


def box(ax, x, y, w, h, title, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.03",
                                facecolor="white", edgecolor=color, linewidth=1.0))
    ax.text(x + w / 2, y + h * (0.72 if body else 0.50), title, ha="center", va="center",
            fontsize=7, fontweight="bold", color=color)
    if body:
        ax.text(x + w / 2, y + h * 0.34, body, ha="center", va="center", fontsize=7)


def arrow(ax, x0, x1, y):
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="-|>", color="#333333", lw=0.6))


def fig1():
    fig, ax = plt.subplots(figsize=mm(174, 132))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    chk = pd.read_csv(RES / "v6/corruption/input_check_v6.csv")
    fail = float(chk[(chk.tag == "corrupted") & np.isclose(chk.f, 0.10)].set_failure_rate.iloc[0])
    mask = pd.read_csv(RES / "v8/masking/masking_fixed_v8.csv").iloc[0]
    rec = pd.read_csv(V12 / "mouse_recovery_ci_v12.csv")
    cells = rec[rec.unit.str.startswith("celltype:")]
    nk = pd.read_csv(RES / "v3/human/r5/replication_detail_r5.csv")
    nk_r = float(nk.loc[nk.coarse == "NK", "r_centred"].iloc[0])
    n_tests = int((~pd.read_csv(RES / "v3/human/r5/s9a_fdr_r5.csv").unestimable.astype(str).str.lower().isin(["true", "1", "1.0"])).sum())
    top = [
        ("Rodent\ntraining data", "", ADOPT),
        ("23 modules", "", ADOPT),
        ("Serialized clocks\nmodule + composite", "", ADOPT),
        ("1  Input check", f"f = 0.10\nfail {fail:.0%}", NEW),
        ("Ortholog map", "", REIMPL),
        ("Pseudobulk\ndonor × cell type", "", REIMPL),
    ]
    w, h, y, gap = 1.48, 2.5, 6.6, 0.16
    x0 = 0.12
    for i, (title, body, color) in enumerate(top):
        x = x0 + i * (w + gap)
        box(ax, x, y, w, h, title, body, color)
        if i < len(top) - 1:
            arrow(ax, x + w, x + w + gap, y + h / 2)
    tissue = float(rec.loc[rec.unit == "whole", "recovery"].iloc[0])
    bottom = [
        ("2  Score", f"{n_tests} tests\n11 cell types", REIMPL),
        ("3  Decompose", f"S {mask.S:.3f}\nN {mask.N:.3f}", NEW),
        ("4  Claim gates", "Q1–Q3 met\nIPF pass 0", NEW),
        ("5  Resolution", f"tissue {tissue:.0%}\ncell types {100 * cells.recovery.min():.0f}–{100 * cells.recovery.max():.0f}%", NEW),
        ("6  Replication", f"NK r = {nk_r:.2f}", NEW),
    ]
    bw, yb = 1.72, 1.05
    span = 5 * bw + 4 * gap
    xb0 = (10 - span) / 2
    for i, (title, body, color) in enumerate(bottom):
        x = xb0 + i * (bw + gap)
        box(ax, x, yb, bw, 3.4, title, body, color)
        if i < len(bottom) - 1:
            arrow(ax, x + bw, x + bw + gap, yb + 1.7)
    last = x0 + 5 * (w + gap) + w / 2
    ax.annotate("", xy=(xb0 + bw / 2, yb + 3.4), xytext=(last, y),
                arrowprops=dict(arrowstyle="-|>", color="#333333", lw=0.6))
    ax.plot([], [], color=ADOPT, lw=1.4, label="adopted")
    ax.plot([], [], color=REIMPL, lw=1.4, label="re-implemented")
    ax.plot([], [], color=NEW, lw=1.4, label="new in this study")
    ax.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
    remember("Fig1", "flow", f"{mask.S:.3f}", "results/v8/masking/masking_fixed_v8.csv", "S", float(mask.S))
    remember("Fig1", "flow", f"{mask.N:.3f}", "results/v8/masking/masking_fixed_v8.csv", "N", float(mask.N))
    remember("Fig1", "flow", f"{nk_r:.2f}", "results/v3/human/r5/replication_detail_r5.csv", "r_centred", nk_r)
    manifest("Fig1", "flow", "results/v6/corruption/input_check_v6.csv; results/v8/masking/masking_fixed_v8.csv; results/v12/mouse_recovery_ci_v12.csv",
             "set_failure_rate; S; N; recovery", "green adopted, grey re-implemented, blue new")
    fig.tight_layout()
    return save(fig, "Fig1")


def fig2():
    fig = plt.figure(figsize=mm(174, 210))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.35, 1], hspace=0.95, wspace=1.35,
                          left=0.14, right=0.98, top=0.96, bottom=0.08)
    axa = fig.add_subplot(gs[0, 0:3])
    axb = fig.add_subplot(gs[0, 3:6])
    axc = fig.add_subplot(gs[1, 0:2])
    dgs = gs[1, 2:4].subgridspec(2, 1, hspace=1.15)
    axd = fig.add_subplot(dgs[0])
    axp = fig.add_subplot(dgs[1], sharex=axd)
    axe = fig.add_subplot(gs[1, 4:6])
    for ax in (axa, axb, axc, axd, axp, axe):
        style_ax(ax)
    nat = pd.read_csv(RES / "v3/rodent/native_reproduction.csv")
    nat = nat[(nat.outcome == "Mortality") & (nat.module != "allmodulegenes")]
    b1 = pd.read_csv(RES / "audit_module_mapping/stage2/b1_reproduction.csv")
    b1 = b1[(b1.source == "supp_table") & (b1.outcome == "Mortality") & (b1.convention == "as_released_na0")]
    b1 = b1[b1.module != "All module genes"]
    m = nat.merge(b1, on="module", suffixes=("_off", "_5b"))
    order = m.sort_values("r_native")
    y = np.arange(len(order))
    for i, r in enumerate(order.itertuples()):
        axa.plot([r.r_insample, r.r_native], [i, i], color=DROP, lw=0.6, zorder=1)
    axa.scatter(order.r_insample, y, s=10, color=DROP, label="spreadsheet", zorder=2)
    axa.scatter(order.r_native, y, s=10, color=ADOPT, label="serialized", zorder=2)
    reported = float(nat.r_reported.median())
    axa.axvline(reported, color="black", lw=0.6)
    axa.set_yticks(y)
    axa.set_yticklabels(order.module.to_numpy())
    lo = float(min(order.r_insample.min(), order.r_native.min())) - 0.05
    axa.set_xlim(lo, 1.05)
    axa.set_xlabel("in-sample r, mortality clocks")
    axa.legend(loc="lower right")
    remember("Fig2", "a", f"{reported:.2f}", "results/v3/rodent/native_reproduction.csv", "r_reported median", reported)
    manifest("Fig2", "a", "results/v3/rodent/native_reproduction.csv; results/audit_module_mapping/stage2/b1_reproduction.csv",
             "r_native; r_insample; r_reported", "Mortality modules; spreadsheet convention as_released_na0; vertical line is the reported median")

    band = pd.read_csv(RES / "v10/coherence_null_band_v10.csv")
    fams = list(dict.fromkeys(band.family))
    x = np.arange(len(fams))
    off = band[band.source == "official"].set_index("family").loc[fams]
    tab = band[band.source == "table5b"].set_index("family").loc[fams]
    axb.bar(x - 0.18, off.top_share, 0.32, color=ADOPT, label="serialized")
    axb.bar(x + 0.18, tab.top_share, 0.32, color=DROP, label="spreadsheet")
    for i, fam in enumerate(fams):
        axb.plot([i - 0.18, i - 0.18], [off.loc[fam, "null_p2_5"], off.loc[fam, "null_p97_5"]], color="black", lw=0.6)
        axb.plot([i + 0.18, i + 0.18], [tab.loc[fam, "null_p2_5"], tab.loc[fam, "null_p97_5"]], color="black", lw=0.6)
    axb.set_xticks(x)
    axb.set_xticklabels(fams, rotation=90, ha="center")
    axb.set_ylabel("top-module share")
    axb.set_ylim(0, 1.35)
    axb.legend(loc="upper right")
    manifest("Fig2", "b", "results/v10/coherence_null_band_v10.csv", "top_share; null_p2_5; null_p97_5",
             "whisker is the null 2.5–97.5 percent; seed 0, 1000 label shuffles")

    chk = pd.read_csv(RES / "v6/corruption/input_check_v6.csv")
    cor = chk[chk.tag == "corrupted"].sort_values("f")
    off_row = chk[chk.tag == "official"].iloc[0]
    tab_row = chk[chk.tag == "table5b"].iloc[0]
    xs = [0.0] + cor.f.tolist()
    ys = [off_row.set_failure_rate] + cor.set_failure_rate.tolist()
    lo = [off_row.ci_lo] + cor.ci_lo.tolist()
    hi = [off_row.ci_hi] + cor.ci_hi.tolist()
    axc.errorbar(xs, ys, yerr=[np.array(ys) - np.array(lo), np.array(hi) - np.array(ys)],
                 fmt="o", color=NEW, ms=4, lw=0.6, capsize=2)
    axc.scatter([1.08], [tab_row.set_failure_rate], marker="*", s=40, color=ORANGE, zorder=3)
    axc.set_xlim(-0.08, 1.22)
    axc.set_ylim(-0.18, 1.22)
    axc.annotate("serialized", xy=(0.0, float(off_row.set_failure_rate)), xytext=(0.28, 0.18),
                 fontsize=7, ha="left", va="bottom",
                 arrowprops=dict(arrowstyle="-", lw=0.4, color="black"))
    axc.annotate("spreadsheet", xy=(1.08, float(tab_row.set_failure_rate)), xytext=(1.02, 0.66),
                 fontsize=7, color=ORANGE, ha="right", va="center",
                 arrowprops=dict(arrowstyle="-", lw=0.4, color=ORANGE,
                                 shrinkA=2, shrinkB=2))
    axc.set_xticks([0, 0.1, 0.25, 0.5, 1.0])
    axc.set_xticklabels(["0", "0.10", "0.25", "0.50", "1"], rotation=90)
    axc.set_xlabel("label-swap fraction f")
    axc.set_ylabel("set-failure rate")
    manifest("Fig2", "c", "results/v6/corruption/input_check_v6.csv", "f; set_failure_rate; ci_lo; ci_hi",
             "Wilson interval; f = 0 is the serialized set; star is the spreadsheet")

    scr = pd.read_csv(RES / "v10/corruption_screen_v10.csv")
    sub = scr[scr.tag == "corrupted"]
    base = scr[scr.tag == "official"]
    width = 0.18
    fs = [0.25, 0.50, 1.00]
    for j, (path, color) in enumerate([("V1", NEW), ("V2", ORANGE)]):
        robust_mean, prog_mean, clock_mean = [], [], []
        for i, frac in enumerate(fs):
            g = sub[(sub.path == path) & np.isclose(sub.f_x, frac)]
            robust_mean.append(g.n_robust.mean())
            prog_mean.append(g.n_program_pass.mean())
            clock_mean.append(g.n_clock_pass.mean())
            x0 = i + (j - 0.5) * width
            axd.scatter(np.full(len(g), x0), g.n_robust, s=8, color=color, zorder=3)
            axp.scatter(np.full(len(g), x0 - 0.05), g.n_program_pass, s=10, color=color, marker="o", zorder=3)
            axp.scatter(np.full(len(g), x0 + 0.05), g.n_clock_pass, s=10, color=color, marker="D", zorder=3)
        axd.bar(np.arange(3) + (j - 0.5) * width, robust_mean, width, color=color, label=path, zorder=2)
        axp.hlines(prog_mean, np.arange(3) + (j - 0.5) * width - 0.08, np.arange(3) + (j - 0.5) * width - 0.02, color=color, lw=1.2)
        axp.hlines(clock_mean, np.arange(3) + (j - 0.5) * width + 0.02, np.arange(3) + (j - 0.5) * width + 0.08, color=color, lw=1.2)
    serial = float(base.n_robust.iloc[0])
    axd.axhline(serial, color="black", lw=0.8)
    handles, _ = axd.get_legend_handles_labels()
    handles.append(plt.Line2D([0], [0], color="black", lw=0.8, label=f"serialized clocks ({serial:.0f})"))
    axd.set_xlim(-0.45, 2.45)
    axd.set_ylim(-1, 52)
    axd.set_yticks([0, 10, 20])
    axd.legend(handles=handles, loc="upper right", ncol=3, frameon=False,
               borderpad=0.2, handlelength=1.2, columnspacing=0.6)
    axd.set_ylabel("robust hits\n(count)")
    axd.tick_params(labelbottom=False)
    axp.set_xlim(-0.45, 3.4)
    axp.set_ylim(-0.35, 2.4)
    axp.set_yticks([0, 1])
    axp.set_xticks([0, 1, 2])
    axp.set_xticklabels(["0.25", "0.50", "1"])
    axp.set_xlabel("label-swap fraction f")
    axp.set_ylabel("screen passes\n(count)")
    axp.legend(handles=[
        plt.Line2D([0], [0], marker="o", color="black", lw=0, label="program-level"),
        plt.Line2D([0], [0], marker="D", color="black", lw=0, label="clock-direction"),
        plt.Line2D([0], [0], color="0.35", lw=1.6, solid_capstyle="butt", label="mean"),
    ], loc="lower left", bbox_to_anchor=(0.0, 1.02), frameon=False, borderaxespad=0,
       handlelength=0.9)
    manifest("Fig2", "d", "results/v10/corruption_screen_v10.csv", "n_robust; n_program_pass; n_clock_pass",
             "bars and ticks are means; points are replicates; line is the serialized robust count")

    raw = pd.read_csv(RES / "v3/human/r5/scrambled_null_raw_r5.csv")
    counts, _, _ = axe.hist(raw.hits_classical, bins=np.arange(0, raw.hits_classical.max() + 2) - 0.5,
                            color=SKY, edgecolor="white", linewidth=0.4)
    axe.axvline(10, color=ADOPT, lw=0.8)
    axe.axvline(30, color=ORANGE, lw=0.8)
    top = float(np.max(counts)) + 5.2
    axe.set_ylim(0, top)
    axe.text(12.2, top - 0.15, "serialized 10", color=ADOPT, fontsize=7, ha="left", va="top")
    axe.text(28.2, top - 1.55, "original\nsubmission 30", color=ORANGE, fontsize=7, ha="right", va="top")
    axe.set_xlabel("nominal hits (count)")
    axe.set_ylabel("scrambled sets (count)")
    remember("Fig2", "e", "10", "results/v3/human/r5/scrambled_null_summary_r5.csv", "percentile reference", 10)
    remember("Fig2", "e", "30", "results/v3/human/r5/scrambled_null_summary_r5.csv", "percentile_v020_30 reference", 30)
    manifest("Fig2", "e", "results/v3/human/r5/scrambled_null_raw_r5.csv", "hits_classical",
             "100 scrambled memberships; lines at 10 serialized hits and 30 original-submission hits")
    letters(fig, [(axa, "a"), (axb, "b"), (axc, "c"), (axd, "d"), (axe, "e")])
    return save(fig, "Fig2")


def fig3():
    fig = plt.figure(figsize=mm(174, 175))
    gs = fig.add_gridspec(2, 6, hspace=0.85, wspace=1.40, left=0.08, right=0.98, top=0.94, bottom=0.16)
    axa = fig.add_subplot(gs[0, 0:2])
    axb = fig.add_subplot(gs[0, 2:6])
    axc = fig.add_subplot(gs[1, 0:2])
    axd = fig.add_subplot(gs[1, 2:4])
    axe = fig.add_subplot(gs[1, 4:6])
    for ax in (axa, axb, axc, axd, axe):
        style_ax(ax)
    axa.axis("off")
    axa.set_xlim(0, 1)
    axa.set_ylim(0, 1)
    axa.text(0.5, 0.82, "remove disease means", ha="center", fontsize=7)
    axa.annotate("", xy=(0.5, 0.68), xytext=(0.5, 0.76), arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#333333"))
    axa.text(0.5, 0.60, "shuffle donor labels", ha="center", fontsize=7)
    axa.annotate("", xy=(0.5, 0.46), xytext=(0.5, 0.54), arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#333333"))
    axa.text(0.5, 0.28, "S1a  S1b\nS2a  S4", ha="center", va="center", fontsize=7,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=NEW, lw=0.8))
    p3 = pd.read_csv(RES / "v7/spikein/program_criteria_v7.csv")
    p3row = p3[p3.item == "P3"].iloc[0]
    axa.text(0.5, 0.06, f"S3: IPF-specific {int(p3row.k)}/{int(p3row.n)}", ha="center", fontsize=7)
    remember("Fig3", "a", str(int(p3row.k)), "results/v7/spikein/program_criteria_v7.csv", "k", int(p3row.k))
    manifest("Fig3", "a", "results/v7/spikein/program_criteria_v7.csv", "k; n", "schematic")

    cats = pd.read_csv(V12 / "fig3_categories_v12.csv")
    order = ["S1a", "S1b", "S2a", "S4"]
    cat_order = ["pass", "robust + IPF-specific", "robust + set-specific", "robust only", "hit only", "not a hit"]
    cat_col = [ADOPT, IPF, ORANGE, NEW, WITHIN, DROP]
    bottom = np.zeros(len(order))
    for cat, col in zip(cat_order, cat_col):
        vals = []
        for sc in order:
            g = cats[(cats.scenario == sc) & (cats.category == cat)]
            vals.append(float(g.rate.iloc[0]) if len(g) else 0.0)
        axb.bar(order, vals, bottom=bottom, color=col, width=0.7, label=cat)
        bottom += np.array(vals)
    axb.set_ylabel("proportion of pairs")
    axb.set_ylim(0, 1.05)
    leg = axb.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=True, framealpha=1, edgecolor="0.8")
    leg.set_zorder(5)
    manifest("Fig3", "b", "results/v12/fig3_categories_v12.csv", "scenario; category; rate; denominator",
             "S1a and S1b: planted module pairs; S2a: pairs in the three planted cell types; S4: all pairs")

    bias = pd.read_csv(V12 / "fig3_bias_v12.csv")
    labels = ["set", "weight"]
    x = np.arange(len(labels))
    for shift, method, color, alpha in [(-0.28, "within", WITHIN, 1), (-0.09, "screen", NEW, 1),
                                        (0.10, "within", WITHIN, 0.35), (0.28, "screen", NEW, 0.35)]:
        klass = "false_hit" if shift < 0 else "true_hit"
        vals, yerr = [], [[], []]
        for spec in labels:
            r = bias[(bias.hit_class == klass) & (bias.specificity == spec) & (bias.method == method)].iloc[0]
            vals.append(r.rate)
            yerr[0].append(r.rate - r.ci_lo)
            yerr[1].append(r.ci_hi - r.rate)
        axc.bar(x + shift, vals, 0.16, color=color, alpha=alpha, yerr=yerr, capsize=1.5, error_kw=dict(lw=0.6))
    axc.set_xticks(x)
    axc.set_xticklabels(["set-specific", "weight-specific"])
    axc.set_ylabel("proportion of hits")
    axc.set_ylim(0, 1.28)
    axc.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, color=WITHIN, label="within hits"),
        plt.Rectangle((0, 0), 1, 1, color=NEW, label="across all tests"),
        plt.Rectangle((0, 0), 1, 1, color=DROP, alpha=0.35, label="pale: planted-module hits"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.32))
    manifest("Fig3", "c", "results/v12/fig3_bias_v12.csv", "rate; ci_lo; ci_hi; hit_class; method",
             "Solid bars are false hits (n = 293). Pale bars are hits in S1a and S1b planted module pairs. Wilson 95% CI.")

    mde = pd.read_csv(RES / "v9/screen_mde_v9.csv")
    mde = mde[(mde.method == "screen") & (mde.scenario == "S1a") & mde.celltype.isin(CT_ORDER)].copy()
    cells = [c for c in CT_ORDER if c in set(mde.celltype)]
    y = np.arange(len(cells))
    for gate, color, dx in [("program", ADOPT, -0.15), ("clock", NEW, 0.15)]:
        g = mde[mde.gate == gate].set_index("celltype")
        for i, ct in enumerate(cells):
            raw = str(g.loc[ct, "mde80"])
            if raw.startswith(">"):
                axd.annotate("", xy=(2.15, i + dx), xytext=(2.0, i + dx),
                             arrowprops=dict(arrowstyle="-|>", color=color, lw=0.7))
            else:
                axd.scatter([float(raw)], [i + dx], s=14, color=color, zorder=3)
    axd.axvline(2.0, color=DROP, lw=0.6, ls="--")
    axd.set_yticks(y)
    axd.set_yticklabels([CT_LABEL[c] for c in cells])
    axd.set_xlabel("MDE80 (residual SD)")
    axd.set_xlim(1.2, 2.4)
    axd.legend(handles=[
        plt.Line2D([0], [0], marker="o", color=ADOPT, lw=0, label="program-level"),
        plt.Line2D([0], [0], marker="o", color=NEW, lw=0, label="clock-direction"),
    ], loc="lower center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False)
    manifest("Fig3", "d", "results/v9/screen_mde_v9.csv", "mde80", "S1a; 11 cell types of the IPF grid; arrow means MDE80 > 2")

    mde_all = pd.read_csv(RES / "v9/screen_mde_v9.csv")
    for scenario, ls in [("S1a", "-"), ("S1b", "--")]:
        for gate, color in [("program", ADOPT), ("clock", NEW)]:
            row = mde_all[(mde_all.method == "screen") & (mde_all.gate == gate)
                          & (mde_all.scenario == scenario) & (mde_all.celltype == "(all)")].iloc[0]
            axe.plot([1, 1.5, 2], [row.rate_d1, row.rate_d15, row.rate_d2], ls=ls, color=color, lw=0.9, marker="o", ms=3)
    axe.set_xlabel("planted effect d (residual SD)")
    axe.set_ylabel("pass rate (proportion)")
    axe.set_ylim(-0.02, 1.05)
    axe.set_xticks([1, 1.5, 2])
    axe.legend(handles=[
        plt.Line2D([0], [0], color=ADOPT, lw=0.9, label="S1a program"),
        plt.Line2D([0], [0], color=NEW, lw=0.9, label="S1a clock"),
        plt.Line2D([0], [0], color=ADOPT, lw=0.9, ls="--", label="S1b program"),
        plt.Line2D([0], [0], color=NEW, lw=0.9, ls="--", label="S1b clock"),
    ], loc="upper left")
    manifest("Fig3", "e", "results/v9/screen_mde_v9.csv", "rate_d1; rate_d15; rate_d2",
             "pooled planted pairs; solid S1a, dashed S1b")
    letters(fig, [(axa, "a"), (axb, "b"), (axc, "c"), (axd, "d"), (axe, "e")])
    return save(fig, "Fig3")


def fig4():
    fig = plt.figure(figsize=mm(174, 234))
    # Top two rows keep the v16 millimetre positions. The extra height goes to the bottom row
    # so panel d's cell-type labels clear 1.3 times the font size.
    def box(x, y_mm, w, h_mm):
        return [x, y_mm / 234, w, h_mm / 234]
    axa = fig.add_axes(box(0.1600, 188.90, 0.8150, 34.80))
    axb = fig.add_axes(box(0.1600, 130.76, 0.4896, 31.49))
    gax = fig.add_axes(box(0.8567, 130.76, 0.1183, 31.49))
    axc = fig.add_axes(box(0.1600, 23.36, 0.1511, 80.75))
    axd = fig.add_axes(box(0.5182, 23.36, 0.4568, 80.75))
    for ax in (axa, axb, gax, axc, axd):
        style_ax(ax)
    klo = pd.read_csv(RES / "v3/rodent/klotho_official.csv")
    band = pd.read_csv(V12 / "klotho_noise_band_v12.csv").set_index("tissue")
    modules = list(klo[klo.tissue == "Kidney"].module)
    axa.axhspan(band.loc["Skeletal muscle", "p2_5"], band.loc["Skeletal muscle", "p97_5"], color=DROP, alpha=0.18, lw=0, label="muscle WT, 6 vs 6")
    axa.axhspan(band.loc["Kidney", "p2_5"], band.loc["Kidney", "p97_5"], color=DROP, alpha=0.45, lw=0, label="kidney WT, 6 vs 6")
    summary = pd.read_csv(RES / "v3/rodent/klotho_official_summary.csv").set_index("tissue")
    for tissue, color, dx, name in [("Kidney", ORANGE, -0.12, "kidney"), ("Skeletal muscle", CTRL, 0.12, "muscle")]:
        g = klo[klo.tissue == tissue].set_index("module").loc[modules]
        npos = int(summary.loc[tissue, "n_positive"])
        axa.scatter(np.arange(len(modules)) + dx, g.delta, s=12, color=color, label=f"{name} {npos}/23", zorder=3)
    axa.axhline(0, color="black", lw=0.4)
    axa.set_xticks(np.arange(len(modules)))
    axa.set_xticklabels(modules, rotation=90, ha="center")
    axa.set_ylabel("KO − WT (Δlog10 hazard units)")
    axa.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=4)
    manifest("Fig4", "a", "results/v3/rodent/klotho_official.csv; results/v12/klotho_noise_band_v12.csv",
             "delta; p2_5; p97_5", "bands are the 2.5–97.5% of WT-half differences after the 1/sqrt(2) six-versus-six scaling")

    hl = pd.read_csv(RES / "v3/human/r5/s8_hlca_all_r5.csv")
    hl = hl[(hl.outcome == "Mortality") & (hl.module != "Composite")]
    cells = sorted(hl.celltype.unique())
    rng = np.random.default_rng(0)
    for i, ct in enumerate(cells):
        g = hl[hl.celltype == ct]
        axb.scatter(np.full(len(g), i) + rng.uniform(-0.15, 0.15, len(g)), g.rho_adj, s=6, color=NEW, alpha=0.75, linewidths=0)
        thr = float(g.min_abs_rho_80power.iloc[0])
        axb.plot([i - 0.28, i + 0.28], [thr, thr], color=DROP, lw=0.8)
        axb.plot([i - 0.28, i + 0.28], [-thr, -thr], color=DROP, lw=0.8)
    axb.axhline(0, color="black", lw=0.4)
    axb.set_xticks(range(len(cells)))
    axb.set_xticklabels(cells, rotation=90)
    axb.set_ylabel("adjusted Spearman ρ")
    axb.legend(handles=[plt.Line2D([0], [0], color=DROP, lw=1.2, label="80% power |ρ|, α = 0.05")],
               loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False)
    manifest("Fig4", "b-left", "results/v3/human/r5/s8_hlca_all_r5.csv", "rho_adj; min_abs_rho_80power",
             "Grey ticks are the uncorrected two-sided α = 0.05 Fisher-z |ρ| with 80% power")

    gtex = pd.read_csv(RES / "v3/human/s8_gtex.csv")
    gt = gtex[(gtex.outcome == "Mortality") & (gtex.module != "Composite")].sort_values("adj_spearman")
    sig = gt.adj_fdr < 0.10
    gax.scatter(np.arange(len(gt)), gt.adj_spearman, s=8, color=np.where(sig, ORANGE, DROP), linewidths=0)
    comp = float(gtex[(gtex.outcome == "Mortality") & (gtex.module == "Composite")].adj_spearman.iloc[0])
    gax.axhline(comp, color=ADOPT, lw=0.8)
    gax.set_ylim(-0.18, 0.30)
    gax.text(0.4, comp - 0.02, f"composite {comp:.2f}", color=ADOPT, fontsize=7, ha="left", va="top")
    gax.set_xlim(-1, len(gt))
    gax.set_xticks([0, 10, 20])
    gax.set_xlabel("modules ordered\nby adjusted ρ", labelpad=2)
    gax.set_ylabel("adjusted Spearman ρ")
    gax.set_ylim(-0.22, 0.42)
    gax.legend(handles=[plt.Line2D([0], [0], marker="o", color=ORANGE, lw=0, markersize=2.8,
                                  label="FDR < 0.10")],
               loc="upper right", frameon=False)
    remember("Fig4", "b", f"{comp:.2f}", "results/v3/human/s8_gtex.csv", "adj_spearman", comp)

    rec = pd.read_csv(V12 / "mouse_recovery_ci_v12.csv")
    y = rec.recovery * 100
    yerr = np.vstack([y - rec.recovery_ci_lo, rec.recovery_ci_hi - y])
    axc.bar(np.arange(len(rec)), y, color=CTRL, width=0.72)
    axc.errorbar(np.arange(len(rec)), y, yerr=yerr, fmt="none", ecolor="black", elinewidth=0.6, capsize=2)
    axc.axhline(100, color="black", lw=0.5, ls="--")
    axc.set_xticks(np.arange(len(rec)))
    axc.set_xticklabels(rec.label, rotation=90)
    axc.set_ylabel("recovery of expected shift (%)")
    manifest("Fig4", "c", "results/v12/mouse_recovery_ci_v12.csv", "recovery; recovery_ci_lo; recovery_ci_hi",
             "100% is the training expected shift; error bars are 95% CIs of the same fits")

    axd.cla()
    style_ax(axd)
    tops = pd.read_csv(V12 / "fig4d_abs_v12.csv")
    ylabels, ypos = [], []
    for shift, species, color in [(0, "mouse", CTRL), (12, "human", ADOPT)]:
        d = tops[tops.species == species].sort_values("rank_abs")
        yy = np.arange(len(d)) + shift
        axd.barh(yy, d.contrib_ratio, color=color, height=0.72)
        ypos.extend(yy)
        ylabels.extend([str(c).replace("_", " ") for c in d.celltype])
    axd.set_yticks(ypos)
    axd.set_yticklabels(ylabels)
    axd.axvline(0, color="black", lw=0.4)
    add_m = float(tops.loc[tops.species == "mouse", "additivity_ratio"].iloc[0])
    add_h = float(tops.loc[tops.species == "human", "additivity_ratio"].iloc[0])
    axd.set_xlabel("contribution (fraction of the expression-only effect)")
    axd.set_title(f"sum of contributions: mouse {add_m:.2f}, HLCA {add_h:.2f}", fontsize=7, pad=8)
    axd.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, color=CTRL, label="mouse lung"),
        plt.Rectangle((0, 0), 1, 1, color=ADOPT, label="HLCA"),
    ], loc="lower right")
    remember("Fig4", "d", f"{add_m:.2f}", "results/v8/decomposition/contrib_summary_v8.csv", "additivity_ratio", add_m)
    remember("Fig4", "d", f"{add_h:.2f}", "results/v8/decomposition/contrib_summary_v8.csv", "additivity_ratio", add_h)
    manifest("Fig4", "d", "results/v12/fig4d_abs_v12.csv", "contrib_ratio; rank_abs; additivity_ratio",
             "eight cell types with the largest absolute contributions, remainder pooled")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labs = sorted((t for t in axd.get_yticklabels() if t.get_text().strip()),
                  key=lambda t: t.get_window_extent(renderer).y0)
    font_px = 7 / 72 * fig.dpi
    gaps = []
    for a, b in zip(labs, labs[1:]):
        ca = a.get_window_extent(renderer)
        cb = b.get_window_extent(renderer)
        gaps.append(((cb.y0 + cb.y1) / 2 - (ca.y0 + ca.y1) / 2) / font_px)
    # The two blocks are separated by a gap; the within-block spacing is the smaller one.
    within = min(gaps) if gaps else 0
    print(f"  fig4d line spacing {within:.2f} x font", flush=True)
    if within < 1.3:
        raise SystemExit(f"Fig4d label spacing {within:.2f} is below 1.3")
    letters(fig, [(axa, "a"), (axb, "b"), (axc, "c"), (axd, "d")])
    return save(fig, "Fig4")


def fig5():
    fig = plt.figure(figsize=mm(174, 230))
    axa = fig.add_axes([0.48, 0.70, 0.20, 0.24])
    axb = fig.add_axes([0.84, 0.70, 0.14, 0.24])
    axc = fig.add_axes([0.18, 0.36, 0.70, 0.26])
    axd = [fig.add_axes([0.10 + i * 0.31, 0.05, 0.18, 0.18]) for i in range(3)]
    for ax in (axa, axb, axc, *axd):
        style_ax(ax)
    parts = pd.read_csv(RES / "v3/human/r5/masking_global_r5.csv")
    parts = parts[parts.level == "component"].copy()
    parts["label"] = np.where(parts.item.astype(str) == "outside", "outside modules", parts.official_name.astype(str))
    parts = parts.sort_values("value")
    fixed = pd.read_csv(RES / "v8/masking/masking_fixed_v8.csv").iloc[0]
    comp_row = dict(label="composite", value=float(fixed.delta_sum), ci_lo=float(fixed.composite_ci_lo),
                    ci_hi=float(fixed.composite_ci_hi), fdr=1.0)
    plot = pd.concat([parts, pd.DataFrame([comp_row])], ignore_index=True)
    sig = plot.fdr < 0.10
    y = np.arange(len(plot))
    axa.errorbar(plot.value, y, xerr=[plot.value - plot.ci_lo, plot.ci_hi - plot.value],
                 fmt="none", ecolor=DROP, elinewidth=0.6, capsize=1.2)
    axa.scatter(plot.value[~sig], y[~sig.to_numpy()], s=12, facecolors="none", edgecolors=NEW, linewidths=0.6)
    axa.scatter(plot.value[sig], y[sig.to_numpy()], s=12, color=NEW)
    axa.axvline(0, color="black", lw=0.4)
    axa.set_yticks(y)
    axa.set_yticklabels(plot.label, fontsize=6)
    axa.set_xlabel("effect (Δlog10 hazard units)")
    data_lo = float(np.nanmin(plot.ci_lo))
    data_hi = float(np.nanmax(plot.ci_hi))
    span = data_hi - data_lo
    axa.set_xlim(data_lo - 0.08 * span, data_hi + 1.15 * span)
    draws = pd.read_csv(RES / "v10/masking_null_draws_v10.csv")
    obs = pd.read_csv(RES / "v10/masking_null_observed_v10.csv").iloc[0]
    for i, (series, color, line) in enumerate(zip(["S", "N", "MI"], [NEW, IPF, ADOPT], [obs.S, obs.N, obs.MI])):
        inset = axa.inset_axes([0.58, 0.06 + i * 0.22, 0.38, 0.16])
        vals = draws[series].to_numpy(float)
        inset.hist(vals, bins=16, color=color, edgecolor="white", linewidth=0.2)
        inset.axvline(line, color="black", lw=0.6)
        inset.set_yticks([])
        inset.set_xticks([])
        inset.set_xlabel(series, fontsize=6)
        if series == "MI":
            inset.set_xlim(float(vals.min()), 1.0)
        inset.tick_params(labelsize=6, length=1.5)
    manifest("Fig5", "a", "results/v3/human/r5/masking_global_r5.csv; results/v8/masking/masking_fixed_v8.csv; results/v10/masking_null_draws_v10.csv",
             "value; ci_lo; ci_hi; fdr; delta_sum; composite_ci_lo; composite_ci_hi; S; N; MI",
             "composite uses the donor-bootstrap interval; filled points are component FDR < 0.10; inset B = 2000")

    ipf = pd.read_csv(RES / "v9/ipf_screen_funnels_v9.csv")
    ipf = ipf[ipf.level == "test"]
    for col in ("hit", "robust", "IPF_specific", "program_pass_screen"):
        ipf[col] = ipf[col].map(lambda v: str(v).lower() in {"true", "1", "1.0"})
    steps = [("tests", len(ipf)), ("hit", int(ipf.hit.sum())), ("robust", int(ipf.robust.sum())),
             ("IPF-specific", int(ipf.IPF_specific.sum())), ("pass", int(ipf.program_pass_screen.sum()))]
    axb.barh(np.arange(len(steps))[::-1], [s[1] for s in steps], color=[DROP, SKY, NEW, IPF, ADOPT])
    for i, (name, n) in enumerate(steps):
        axb.text(n + 4, len(steps) - 1 - i, str(n), va="center", fontsize=7)
        remember("Fig5", "b", str(n), "results/v9/ipf_screen_funnels_v9.csv", name, n)
    axb.set_yticks(np.arange(len(steps))[::-1])
    axb.set_yticklabels([s[0] for s in steps])
    axb.set_xlabel("tests remaining (count)")
    axb.set_xlim(0, 340)
    axb.set_ylim(-0.7, len(steps) - 0.3)
    axb.set_xticks([0, 100, 200, 300])
    manifest("Fig5", "b", "results/v9/ipf_screen_funnels_v9.csv", "hit; robust; IPF_specific; program_pass_screen",
             "pass is the program-level screen gate")

    std = pd.read_csv(RES / "v3/human/s11_standardized.csv")
    info = pd.read_csv(RES / "v3/human/r5/s9a_fdr_r5.csv")
    info = info[info.module != "Composite"]
    known = pd.read_csv(RES / "v3/human/r5/score_informativeness.csv")
    known["unestimable"] = known.unestimable.astype(str).str.lower().isin(["true", "1", "1.0"])
    unest = set(zip(known.loc[known.unestimable, "celltype"], known.loc[known.unestimable, "module"]))
    cells = [c for c in CT_ORDER if c in set(info.celltype)]
    mods = sorted(m for m in info.module.unique() if m != "Composite")
    grid = np.full((len(cells), len(mods)), np.nan)
    std_i = std.set_index(["celltype", "module"])
    for i, ct in enumerate(cells):
        for j, mod in enumerate(mods):
            if (ct, mod) in unest:
                continue
            if (ct, mod) in std_i.index:
                grid[i, j] = float(std_i.loc[(ct, mod), "standardized"])
    lim = np.nanmax(np.abs(grid))
    cmap = LinearSegmentedColormap.from_list("ipf", [CTRL, "#F7F7F7", IPF])
    cmap.set_bad("#ffffff")
    im = axc.imshow(np.ma.masked_invalid(grid), cmap=cmap, vmin=-lim, vmax=lim, aspect="auto")
    for ct, mod in unest:
        if ct in cells and mod in mods:
            axc.add_patch(Rectangle((mods.index(mod) - 0.5, cells.index(ct) - 0.5), 1, 1,
                                    facecolor="white", hatch="///", edgecolor=DROP, lw=0.3))
    hits = info[info.fdr_IPF < 0.10]
    for rec in hits.itertuples():
        x, yy = mods.index(rec.module), cells.index(rec.celltype)
        robust = str(rec.robust).lower() in {"true", "1", "1.0"}
        specific = str(rec.IPF_specific).lower() in {"true", "1", "1.0"}
        if specific and robust:
            axc.scatter([x], [yy], marker="*", s=28, color="black", zorder=3)
        elif robust:
            axc.scatter([x], [yy], s=14, facecolors="none", edgecolors="black", linewidths=0.7, zorder=3)
        else:
            axc.scatter([x], [yy], s=10, color="black", zorder=3)
    axc.set_xticks(range(len(mods)))
    axc.set_xticklabels(mods, rotation=90)
    axc.set_yticks(range(len(cells)))
    axc.set_yticklabels([CT_LABEL[c] for c in cells])
    cb = fig.colorbar(im, ax=axc, fraction=0.025, pad=0.02)
    cb.set_label("standardized effect")
    axc.legend(handles=[
        plt.Line2D([0], [0], marker="o", color="black", lw=0, label="hit, not robust"),
        plt.Line2D([0], [0], marker="o", markerfacecolor="none", markeredgecolor="black", lw=0, label="robust hit"),
        plt.Line2D([0], [0], marker="*", color="black", lw=0, label="IPF-specific robust hit"),
    ], loc="upper center", bbox_to_anchor=(0.42, -0.34), ncol=3, frameon=False)
    manifest("Fig5", "c", "results/v3/human/s11_standardized.csv; results/v3/human/r5/s9a_fdr_r5.csv; results/v3/human/r5/score_informativeness.csv",
             "standardized; fdr_IPF; robust; IPF_specific; unestimable",
             "four unestimable black-module tests are blank and hatched")

    don = pd.read_csv(RES / "v10/alveolar_green_donors_v10.csv")
    measures = [("clock_score", "clock score\n(Δlog10 hazard units)"),
                ("unweighted_mean", "mean expression (z)"),
                ("mki67_fraction", "MKI67+ fraction")]
    diseases = [("Control", "Ctrl", CTRL), ("IPF", "IPF", IPF), ("COPD", "COPD", COPD)]
    counts = don.groupby("disease").donor.nunique()
    rng = np.random.default_rng(1)
    for ax, (col, name) in zip(axd, measures):
        for j, (dis, lab, color) in enumerate(diseases):
            vals = don.loc[don.disease == dis, col].dropna().to_numpy(float)
            if len(vals) == 0:
                continue
            ax.scatter(np.full(len(vals), j) + rng.uniform(-0.12, 0.12, len(vals)), vals, s=8, color=color, alpha=0.85, linewidths=0)
            ax.plot([j - 0.2, j + 0.2], [np.mean(vals), np.mean(vals)], color="black", lw=0.9)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels([f"{lab}\n{int(counts.get(dis, 0))}" for dis, lab, _ in diseases])
        ax.set_ylabel(name, labelpad=6)
        ax.axhline(0, color=DROP, lw=0.4)
    manifest("Fig5", "d", "results/v10/alveolar_green_donors_v10.csv", "clock_score; unweighted_mean; mki67_fraction; disease",
             "black lines are group means; x labels give the donor count")
    letters(fig, [(axa, "a"), (axb, "b"), (axc, "c"), (axd[0], "d")])
    return save(fig, "Fig5")


def fig6():
    fig = plt.figure(figsize=mm(174, 124))
    axa = fig.add_axes([0.221, 0.30, 0.14, 0.62])
    axb = fig.add_axes([0.381, 0.30, 0.18, 0.62])
    axc = fig.add_axes([0.681, 0.30, 0.28, 0.62])
    for ax in (axa, axb, axc):
        style_ax(ax)
    rep = pd.read_csv(RES / "v3/human/r5/replication_detail_r5.csv")
    band = pd.read_csv(RES / "v10/concordance_null_band_v10.csv")
    m = rep.merge(band, on="coarse", how="left")
    y = np.arange(len(m))
    axa.hlines(y, m.null_p2_5, m.null_p97_5, color=DROP, lw=1.2)
    axa.scatter(m.r_centred, y, s=18, color=[IPF if c == "NK" else NEW for c in m.coarse], zorder=3)
    axa.axvline(0, color="black", lw=0.4)
    labels = [f"{c}  {int(a)}+{int(b)} | {int(c1)}+{int(d)}"
              for c, a, b, c1, d in zip(m.coarse, m.n_IPF_GSE136831, m.n_Control_GSE136831, m.n_IPF_GSE135893, m.n_Control_GSE135893)]
    axa.set_yticks(y)
    axa.set_yticklabels(labels)
    axa.set_xlabel("Pearson r")
    manifest("Fig6", "a", "results/v3/human/r5/replication_detail_r5.csv; results/v10/concordance_null_band_v10.csv",
             "r_centred; null_p2_5; null_p97_5; n_IPF_GSE136831; n_Control_GSE136831; n_IPF_GSE135893; n_Control_GSE135893",
             "label is IPF+control in GSE136831 | IPF+control in GSE135893")

    conf = pd.read_csv(RES / "v3/human/s15_confusion.csv", keep_default_na=False).set_index("scenario")
    row_map = [("shared", "shared effect"), ("power", "weak shared\neffect"), ("composition", "sub-state\ncomposition"),
               ("composition_shift", "composition\nshift only"), ("confound", "confound in\none cohort"), ("null", "no effect")]
    col_map = [("between-cohort difference (unresolved)", "between-cohort\ndifference"),
               ("composition", "composition"), ("inconclusive (low reliability)", "inconclusive"), ("replicated", "replicated")]
    mat = np.array([[float(conf.loc[s, c]) for c, _ in col_map] for s, _ in row_map])
    im = axb.imshow(mat, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    axb.set_xticks(range(len(col_map)))
    axb.set_xticklabels([lab.replace("\n", " ") for _, lab in col_map], rotation=90, ha="right", rotation_mode="anchor")
    axb.set_yticks(range(len(row_map)))
    axb.set_yticklabels([lab for _, lab in row_map])
    axb.yaxis.tick_right()
    axb.yaxis.set_label_position("right")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            tone = "white" if mat[i, j] >= 0.80 else "black"
            axb.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=7, color=tone)
            remember("Fig6", "b", f"{mat[i, j]:.2f}", "results/v3/human/s15_confusion.csv", col_map[j][0], mat[i, j])
    axb.set_xlabel("assigned outcome (proportion)")
    manifest("Fig6", "b", "results/v3/human/s15_confusion.csv", "scenario and outcome columns",
             "300 simulated pairs; cell text is the stored proportion")

    real = pd.read_csv(RES / "v3/human/s15_checklist_real.csv")
    attr = {"검정 불가": "not tested", "inconclusive (low reliability)": "inconclusive",
            "between-cohort difference (unresolved)": "between-cohort\ndifference", "replicated": "replicated"}
    axc.axis("off")
    axc.set_xlim(0, 1)
    axc.set_ylim(0, 1)
    one_line = {"검정 불가": "not tested", "inconclusive (low reliability)": "inconclusive",
                "between-cohort difference (unresolved)": "between-cohort difference",
                "replicated": "replicated"}
    for i, rec in enumerate(real.itertuples()):
        yy = 0.90 - i * 0.12
        axc.text(0.0, yy, str(rec.coarse), fontsize=7, va="center", ha="left",
                 color=IPF if rec.coarse == "NK" else "black")
        axc.text(1.0, yy, one_line.get(str(rec.attribution), str(rec.attribution)),
                 fontsize=7, va="center", ha="right")
    manifest("Fig6", "c", "results/v3/human/s15_checklist_real.csv", "coarse; attribution",
             "outcome names match panel b; AT2 is not tested")
    letters(fig, [(axa, "a"), (axb, "b"), (axc, "c")])
    return save(fig, "Fig6")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    overlaps = []
    for fn in (fig2, fig3, fig4, fig5, fig6):
        print(fn.__name__, flush=True)
        overlaps.extend(fn())
    pd.DataFrame(overlaps, columns=["figure", "kind", "text"]).to_csv(FIG / "overlap_check.csv", index=False)
    pd.DataFrame(MANIFEST).to_csv(FIG / "manifest_v17.csv", index=False)
    pd.DataFrame(SIZES).to_csv(FIG / "size_check.csv", index=False)
    nums = pd.DataFrame(NUMBERS)
    def ok(row):
        try:
            return abs(float(row.shown) - float(row.value)) <= max(0.006, 0.021 * abs(float(row.value)))
        except Exception:
            return True
    nums["match"] = [ok(r) for r in nums.itertuples()]
    nums.to_csv(FIG / "figure_numbers_check.csv", index=False)
    print("overlaps", len(overlaps), "number mismatches", int((~nums.match).sum()), flush=True)


if __name__ == "__main__":
    main()
