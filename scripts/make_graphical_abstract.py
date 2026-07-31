"""Graphical (visual) abstract for the GeroScience submission.

Renders a clean, colorblind-safe schematic (~2:1, 300 dpi) summarizing the central message:
composite clocks mask opposing programs -> module- & cell-type resolution reveals them ->
attribution ladder + attention localizer make the readout reproducible and interpretable.

Output: manuscript/figures_submission/graphical_abstract.png
"""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Wedge, Circle

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "figures"; OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#0072B2"; ORANGE = "#E69F00"; VERM = "#D55E00"; GREY = "#8a8f98"
INK = "#222222"

fig = plt.figure(figsize=(8.0, 4.0)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, 16); ax.set_ylim(0, 8)

# ---- Title band
ax.text(8, 7.55, "Module- and cell-type-resolved transcriptomic aging clocks",
        ha="center", va="center", fontsize=13.5, fontweight="bold", color=INK)
ax.text(8, 6.98, "expose aging programs that a single composite score hides",
        ha="center", va="center", fontsize=10.5, color="#555")

def box(x, y, w, h, fc, ec="#444", lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec=ec, lw=lw))

def arrow(x0, y0, x1, y1, color=INK, lw=2.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                 mutation_scale=20, color=color, lw=lw))

# ================= Stage 1: composite clock = "no change" =================
box(0.4, 2.7, 3.5, 3.3, "#eef3fb")
ax.text(2.15, 5.6, "Composite\naging clock", ha="center", va="center",
        fontsize=10.5, fontweight="bold", color=INK)
# gauge
cx, cy, r = 2.15, 4.1, 1.05
ax.add_patch(Wedge((cx, cy), r, 0, 180, width=0.26, fc="#cfd8e3", ec="none"))
ax.add_patch(Wedge((cx, cy), r, 60, 120, width=0.26, fc=GREY, ec="none"))
ang = np.deg2rad(90)  # needle at center = no change
ax.plot([cx, cx + 0.82 * r * np.cos(ang)], [cy, cy + 0.82 * r * np.sin(ang)],
        color=INK, lw=2.4)
ax.add_patch(Circle((cx, cy), 0.07, fc=INK, ec="none"))
ax.text(cx, cy - 0.55, "no change", ha="center", va="center", fontsize=9, style="italic",
        color="#444")
ax.text(2.15, 3.02, "IPF:  β = −2.9,  p = 0.76  (n.s.)", ha="center", va="center",
        fontsize=8.2, color="#555")

arrow(4.05, 4.35, 5.15, 4.35)
ax.text(4.6, 4.75, "resolve", ha="center", fontsize=8.5, color="#666")

# ================= Stage 2: decompose into modules x cell types =================
box(5.25, 2.7, 5.0, 3.3, "#fbf3e6")
ax.text(7.75, 5.62, "23 pathway modules  ×  cell types", ha="center", va="center",
        fontsize=10, fontweight="bold", color=INK)
# opposing module bars around zero
labels = ["Adaptive\nimmunity", "Lipid\nmet", "Chromatin", "Interferon", "VEGF"]
vals = [0.9, 0.6, 0.5, -0.75, -0.6]
bx = np.linspace(5.9, 9.6, len(vals)); base = 4.15; sc = 0.95
for x, v in zip(bx, vals):
    c = VERM if v > 0 else BLUE
    ax.add_patch(plt.Rectangle((x - 0.28, base), 0.56, v * sc, fc=c, ec="none"))
ax.plot([5.6, 9.9], [base, base], color=INK, lw=0.9)
ax.text(9.95, base + 0.9 * sc, "↑ accelerated", fontsize=7.6, color=VERM, va="center")
ax.text(9.95, base - 0.75 * sc, "↓ decelerated", fontsize=7.6, color=BLUE, va="center")
ax.text(7.75, 3.0, "opposing programs cancel in the composite", ha="center",
        fontsize=8.4, style="italic", color="#555")

arrow(10.4, 4.35, 11.5, 4.35)
ax.text(10.95, 4.75, "validate", ha="center", fontsize=8.5, color="#666")

# ================= Stage 3: reproducibility + interpretability =================
box(11.6, 4.35, 4.0, 1.65, "#eaf3ec")
ax.text(13.6, 5.66, "Attribution ladder", ha="center", fontsize=9.3, fontweight="bold",
        color=INK)
ax.text(13.6, 4.98, "technical → power →\ncomposition → resolution", ha="center",
        va="center", fontsize=8.0, color="#444")

box(11.6, 2.7, 4.0, 1.5, "#eef3fb")
ax.text(13.6, 3.86, "Attention localizer", ha="center", fontsize=9.3, fontweight="bold",
        color=INK)
ax.text(13.6, 3.2, "pinpoints the SPP1⁺ profibrotic\nsub-state carrying the signal",
        ha="center", va="center", fontsize=8.0, color="#444")

# ================= Bottom banner =================
box(0.4, 0.5, 15.2, 1.5, "#f4f5f7", ec="#ccc", lw=1.0)
ax.text(8, 1.55, "Proof-of-principle:  idiopathic pulmonary fibrosis (single-cell RNA-seq)  ·  "
        "COPD-specificity contrast  ·  external + bulk replication",
        ha="center", va="center", fontsize=8.8, color=INK)
ax.text(8, 0.92, "Reusable, clock-agnostic framework — retrained module clocks & pipeline "
        "released as a resource for single-cell geroscience",
        ha="center", va="center", fontsize=8.4, style="italic", color="#555")

fig.savefig(OUT / "graphical_abstract.png", dpi=300, facecolor="white")
plt.close(fig)
print("saved", OUT / "graphical_abstract.png")
