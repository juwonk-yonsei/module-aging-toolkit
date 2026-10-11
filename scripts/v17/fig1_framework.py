"""Figure 1: the framework (mock-up for the revision).

Layout, 174 mm wide:
  top row      inputs: rodent clocks (adopted) and human lung data (re-implemented), plus the colour key
  two rows     the six steps as cards: what each step asks, its key result as a small glyph and in words,
               and the known truth it was validated against

All numbers come from FACTS. In the repository version, FACTS is filled from the result tables
(source noted beside each entry), so the figure script reads only the result tables.
"""
import re
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
OUT = str(REPO / "figures" / "v17" / "Fig1")

# ----------------------------------------------------------------------------- facts
FACTS = {
    # rodent inputs (Supplementary Table S-M1 / S1)
    'n_profiles': '4,539', 'n_datasets': 96,
    # step 1, input check (S17_input, S17_corrupt; Fig. 2c)
    'f_curve': [(0.0, 0.0), (0.10, 0.90), (0.25, 1.0), (0.50, 1.0), (1.0, 1.0)],
    'f10_failed': '18 of 20',
    # step 2, scoring (S09_ipf, S09_unestimable; S01)
    'n_celltypes': 11, 'n_modules': 23, 'n_planned': 253, 'n_estimable': 249,
    'unestimable': [('T', 'black'), ('T cyt.', 'black'), ('cMonocyte', 'black'), ('ncMonocyte', 'black')],
    'celltypes': ['B', 'Ciliated', 'Fibroblast', 'Macrophage', 'Alv. mac', 'NK', 'T', 'T cyt.', 'cDC2',
                  'cMonocyte', 'ncMonocyte'],
    'modules_first': 'black',
    'gse136831': '32 IPF, 28 control, 18 COPD donors',
    # step 3, decomposition (S11_global)
    'S': 0.221, 'N': 0.023, 'n_comp_fdr': 13, 'n_comp': 24,
    # step 4, claim gates (S09_funnel; S16)
    'funnel': [(249, 'tests'), (10, 'hits'), (9, 'robust'), (2, 'IPF-specific'), (0, 'passed')],
    'null_pass': '0 of 19,334',
    # step 5, resolution (S18_recovery_ci)
    'whole_lung': 58, 'celltype_lo': 2, 'celltype_hi': 15,
    # step 6, replication (S13_replication; S15)
    'cross_r': [('NK', 0.67), ('AT2', 0.37), ('T', 0.25), ('Ciliated', 0.15), ('Macrophage', 0.12),
                ('Monocyte', 0.00), ('B', -0.05)],
    'nk_p': '0.008', 'n_sim_pairs': 300,
}

def compare_facts():
    """Compare the drawn FACTS with stored tables. The figure keeps FACTS."""
    import pandas as pd
    ic = pd.read_csv(REPO / "results/v6/corruption/input_check_v6.csv")
    off = ic[ic.tag == "official"].iloc[0]
    cor = ic[ic.tag == "corrupted"].set_index("f")
    curve = [(0.0, float(off.set_failure_rate)),
             (0.10, float(cor.loc[0.10, "set_failure_rate"])),
             (0.25, float(cor.loc[0.25, "set_failure_rate"])),
             (0.50, float(cor.loc[0.50, "set_failure_rate"])),
             (1.0, float(cor.loc[1.00, "set_failure_rate"]))]
    known = pd.read_csv(REPO / "results/v3/human/r5/score_informativeness.csv")
    known["unestimable"] = known.unestimable.astype(str).str.lower().isin(["true", "1", "1.0"])
    black = known[(known.unestimable) & (known.module == "black")]
    fx = pd.read_csv(REPO / "results/v8/masking/masking_fixed_v8.csv").iloc[0]
    comp = pd.read_csv(REPO / "results/v3/human/r5/masking_global_r5.csv")
    comp = comp[comp.level == "component"]
    funnel = pd.read_csv(REPO / "results/v9/ipf_screen_funnels_v9.csv")
    funnel = funnel[funnel.level == "test"]
    for col in ("hit", "robust", "IPF_specific", "program_pass_screen"):
        funnel[col] = funnel[col].astype(str).str.lower().isin(["true", "1", "1.0"])
    crit = pd.read_csv(REPO / "results/v9/screen_criteria_v9.csv")
    q1 = crit[(crit.gate == "program") & (crit.item == "Q1")].iloc[0]
    rec = pd.read_csv(REPO / "results/v12/mouse_recovery_ci_v12.csv")
    ct = rec[rec.unit.str.startswith("celltype:")]
    rep = pd.read_csv(REPO / "results/v3/human/r5/replication_detail_r5.csv").set_index("coarse")
    order = [("NK", "NK"), ("AT2", "AT2"), ("T", "T cell"), ("Ciliated", "Ciliated"),
             ("Macrophage", "Macrophage"), ("Monocyte", "Monocyte"), ("B", "B cell")]
    pairs = pd.read_csv(REPO / "results/v3/human/s15_simulation_raw.csv")
    coh = pd.read_csv(REPO / "results/v3/human/s1_cohorts.csv")
    adams = coh[coh.cohort.str.startswith("GSE136831")].set_index("group")
    rows = [
        ("n_profiles", FACTS["n_profiles"], "4,539", "training annotation rows 4539"),
        ("n_datasets", FACTS["n_datasets"], 96, "training annotation Source nunique 96"),
        ("f_curve", str(FACTS["f_curve"]), str([(a, round(b, 2)) for a, b in curve]),
         "results/v6/corruption/input_check_v6.csv"),
        ("f10_failed", FACTS["f10_failed"], f"{int(cor.loc[0.10, 'set_failure_k'])} of {int(cor.loc[0.10, 'n_reps'])}",
         "results/v6/corruption/input_check_v6.csv"),
        ("n_estimable", FACTS["n_estimable"], int((~known.unestimable).sum()) if "unestimable" in known else len(funnel),
         "results/v9/ipf_screen_funnels_v9.csv"),
        ("unestimable_n", len(FACTS["unestimable"]), int(black.shape[0]),
         "results/v3/human/r5/score_informativeness.csv"),
        ("gse136831", FACTS["gse136831"],
         f"{int(adams.loc['IPF','donors'])} IPF, {int(adams.loc['Control','donors'])} control, {int(adams.loc['COPD','donors'])} COPD donors",
         "results/v3/human/s1_cohorts.csv"),
        ("S", FACTS["S"], round(float(fx.S), 3), "results/v8/masking/masking_fixed_v8.csv"),
        ("N", FACTS["N"], round(float(fx.N), 3), "results/v8/masking/masking_fixed_v8.csv"),
        ("n_comp_fdr", FACTS["n_comp_fdr"], int((comp.fdr < 0.10).sum()),
         "results/v3/human/r5/masking_global_r5.csv"),
        ("funnel", str([a for a, _ in FACTS["funnel"]]),
         str([len(funnel), int(funnel.hit.sum()), int(funnel.robust.sum()),
              int(funnel.IPF_specific.sum()), int(funnel.program_pass_screen.sum())]),
         "results/v9/ipf_screen_funnels_v9.csv"),
        ("null_pass", FACTS["null_pass"], f"{int(q1.k)} of {int(q1.n):,}",
         "results/v9/screen_criteria_v9.csv"),
        ("whole_lung", FACTS["whole_lung"], int(round(100 * float(rec.loc[rec.unit == "whole", "recovery"].iloc[0]))),
         "results/v12/mouse_recovery_ci_v12.csv"),
        ("celltype_lo", FACTS["celltype_lo"], int(round(100 * float(ct.recovery.min()))),
         "results/v12/mouse_recovery_ci_v12.csv"),
        ("celltype_hi", FACTS["celltype_hi"], int(round(100 * float(ct.recovery.max()))),
         "results/v12/mouse_recovery_ci_v12.csv"),
        ("cross_r", str(FACTS["cross_r"]),
         str([(lab, round(float(rep.loc[key, "r_centred"]), 2)) for lab, key in order]),
         "results/v3/human/r5/replication_detail_r5.csv"),
        ("nk_p", FACTS["nk_p"], f"{float(rep.loc['NK', 'p_two_sided']):.3f}",
         "results/v3/human/r5/replication_detail_r5.csv"),
        ("n_sim_pairs", FACTS["n_sim_pairs"], int(pairs.iloc[:, 0].value_counts().iloc[0]),
         "results/v3/human/s15_simulation_raw.csv"),
    ]
    out = pd.DataFrame(rows, columns=["key", "facts", "stored", "file"])
    out["match"] = out.facts.astype(str) == out.stored.astype(str)
    out.to_csv(REPO / "figures" / "v17" / "fig1_facts.csv", index=False)
    print("Fig1 FACTS mismatches", int((~out.match).sum()))


compare_facts()

# ----------------------------------------------------------------------------- style
W, H = 174.0, 124.0                       # mm
MM = 1 / 25.4
PT2MM = 25.4 / 72
GREEN, GREY, BLUE, ACCENT = '#009E73', '#505050', '#0072B2', '#D55E00'
TINT = {GREEN: '#EAF6F1', GREY: '#F1F1F1', BLUE: '#EAF2F8'}
INK, INK2, MUTED, RULE = '#1A1A1A', '#4D4D4D', '#8C8C8C', '#C9C9C9'
FONT = ['Liberation Sans', 'DejaVu Sans']
plt.rcParams.update({'font.family': FONT, 'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.major.size': 1.6, 'ytick.major.size': 1.6, 'xtick.major.pad': 1.2,
                     'ytick.major.pad': 1.2, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'axes.labelsize': 6, 'axes.labelpad': 1.2, 'axes.edgecolor': INK2,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.labelcolor': INK2})

fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis('off')

def text_w(s, size, weight='normal', style='normal'):
    """Rendered width of s in mm."""
    fp = FontProperties(family=FONT, size=size, weight=weight, style=style)
    return TextPath((0, 0), s, size=size, prop=fp).get_extents().width * PT2MM

def wrap(s, width_mm, size, weight='normal', style='normal'):
    lines, cur = [], ''
    for word in s.split(' '):
        trial = (cur + ' ' + word).strip()
        if text_w(trial, size, weight, style) <= width_mm or not cur:
            cur = trial
        else:
            lines.append(cur); cur = word
    lines.append(cur)
    return '\n'.join(lines)

def put(x, y, s, size=7, color=INK, weight='normal', style='normal', ha='left', va='top', ls=1.18, **kw):
    return ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, fontstyle=style, ha=ha, va=va,
                   linespacing=ls, **kw)

def rbox(x, y, w, h, fc, ec, lw=0.6, r=1.2, z=1):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw,
                       zorder=z, mutation_aspect=1)
    ax.add_patch(p); return p

def arrow(x0, y0, x1, y1, color=INK2, lw=0.8, style='-|>', conn='arc3'):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=6, lw=lw, color=color,
                        connectionstyle=conn, shrinkA=0, shrinkB=0, zorder=5)
    ax.add_patch(a)

def check(x, y, size=2.2, color=BLUE, lw=0.9):
    """A check mark centred at (x, y), drawn as a path so that the figure uses one font."""
    h = size / 2
    ax.plot([x - h, x - h * 0.25, x + h], [y, y + h * 0.75, y - h * 0.8], color=color, lw=lw,
            solid_capstyle='round', solid_joinstyle='round', zorder=6)

def inset(x, y, w, h):
    """An axes placed in mm (x, y = top-left)."""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])

def tidy(a):
    for s in ('top', 'right'):
        a.spines[s].set_visible(False)
    a.patch.set_alpha(0)

# ----------------------------------------------------------------------------- inputs row
MARGIN, GAP = 1.0, 4.0
COLW = (W - 2 * MARGIN - 2 * GAP) / 3
COLX = [MARGIN + i * (COLW + GAP) for i in range(3)]
IN_Y, IN_H = 1.2, 17.0

def input_box(x, w, title, sub, col):
    rbox(x, IN_Y, w, IN_H, TINT[col], col, lw=0.7)
    put(x + w / 2, IN_Y + 2.6, title, size=6.5, weight='bold', ha='center', va='top')
    if sub:
        put(x + w / 2, IN_Y + IN_H - 2.0, sub, size=6.5, color=INK2, ha='center', va='bottom')

# rodent stream (adopted), above step 1
bw = (COLW - 2 * 2.6) / 3
xs = [COLX[0] + i * (bw + 2.6) for i in range(3)]
input_box(xs[0], bw, 'Rodent\ntraining data', f"{FACTS['n_profiles']} profiles\n{FACTS['n_datasets']} datasets", GREEN)
input_box(xs[1], bw, 'Modules', f"{FACTS['n_modules']}, defined by\nco-expression", GREEN)
input_box(xs[2], bw, 'Clocks', 'module and\ncomposite,\nas released', GREEN)
for i in range(2):
    arrow(xs[i] + bw + 0.3, IN_Y + IN_H / 2, xs[i + 1] - 0.3, IN_Y + IN_H / 2)

# human stream (re-implemented), above step 2
hw = (COLW - 2.6) / 2
hx = [COLX[1], COLX[1] + hw + 2.6]
input_box(hx[0], hw, 'Human lung\nscRNA-seq', 'IPF, COPD and control;\ntwo cohorts', GREY)
input_box(hx[1], hw, 'Ortholog mapping,\npseudobulks', 'donor × cell type', GREY)
arrow(hx[0] + hw + 0.3, IN_Y + IN_H / 2, hx[1] - 0.3, IN_Y + IN_H / 2)

# key, above step 3
kx, ky = COLX[2] + 3.0, IN_Y + 1.6
for i, (col, lab) in enumerate([(GREEN, 'Adopted from previous work'), (GREY, 'Re-implemented'),
                                (BLUE, 'New in this study')]):
    yy = ky + i * 3.6
    rbox(kx, yy, 4.2, 2.6, TINT[col], col, lw=0.7, r=0.6)
    put(kx + 6.0, yy + 1.3, lab, size=7, va='center')
yy = ky + 3 * 3.6
check(kx + 2.1, yy + 1.3)
put(kx + 6.0, yy + 1.3, 'Tested against known truth', size=7, va='center')

# ----------------------------------------------------------------------------- step cards
CARD_H, HEAD = 44.0, 6.2
ROW_Y = [IN_Y + IN_H + 7.0, IN_Y + IN_H + 7.0 + CARD_H + 7.5]

def card(col, row, num, title, question, color, result, foot, foot_check):
    x, y = COLX[col], ROW_Y[row]
    rbox(x, y, COLW, CARD_H, TINT[color], color, lw=0.7, z=1)
    rbox(x, y, COLW, HEAD, color, color, lw=0.7, z=2)
    ax.add_patch(Rectangle((x, y + HEAD - 1.4), COLW, 1.4, fc=color, ec=color, lw=0.7, zorder=2))
    ax.add_patch(Circle((x + 3.6, y + HEAD / 2), 2.15, fc='white', ec='none', zorder=3))
    put(x + 3.6, y + HEAD / 2 + 0.1, str(num), size=7.5, color=color, weight='bold', ha='center', va='center',
        zorder=4)
    put(x + 7.0, y + HEAD / 2 + 0.1, title, size=8, color='white', weight='bold', va='center', zorder=4)
    put(x + 2.4, y + HEAD + 1.6, wrap(question, COLW - 4.8, 7, style='italic'), size=7, color=INK2,
        style='italic')
    result = re.sub(r' ([=<]) ', '\u00a0\\1\u00a0', result)     # keep 'N = 0.02', 'FDR < 0.10' on one line
    put(x + 29.0, y + HEAD + 9.4, wrap(result, COLW - 29.0 - 2.2, 7), size=7)
    fy = y + CARD_H - 5.4
    ax.plot([x + 2.4, x + COLW - 2.4], [fy, fy], color=RULE, lw=0.5, zorder=3)
    if foot_check:
        check(x + 3.5, fy + 2.7)
        put(x + 5.6, fy + 2.7, foot, size=6.5, color=INK2, va='center')
    else:
        put(x + 2.4, fy + 2.7, foot, size=6.5, color=INK2, va='center')
    return x, y

GX, GY, GW, GH = 2.6, HEAD + 9.6, 23.0, 21.0      # glyph box inside a card (offsets, mm)

# 1  input check
x, y = card(0, 0, 1, 'Input check', 'Are the module gene labels correct?', BLUE,
            'Serialized clocks pass; the published spreadsheet labels fail. With 10% of genes '
            f"relabelled, {FACTS['f10_failed']} sets fail.",
            'memberships with relabelled genes', True)
a = inset(x + GX + 4.2, y + GY, GW - 4.2, GH - 6.0)
fx, fv = zip(*FACTS['f_curve'])
a.plot(fx, fv, color=BLUE, lw=1.1, zorder=2)
a.scatter(fx, fv, s=9, color=BLUE, edgecolor='white', linewidth=0.6, zorder=3)
a.scatter([1.16], [1.0], s=26, marker='*', color=ACCENT, edgecolor='white', linewidth=0.4, zorder=4,
          clip_on=False)
a.annotate('spreadsheet', (1.16, 1.0), xytext=(1.16, 0.55), fontsize=6, color=INK2, ha='right', va='center',
           arrowprops=dict(arrowstyle='-', lw=0.4, color=MUTED, shrinkA=0, shrinkB=1.5))
a.set_xlim(-0.06, 1.24); a.set_ylim(-0.08, 1.12)
a.set_xticks([0, 0.5, 1]); a.set_xticklabels(['0', '0.5', '1'])
a.set_yticks([0, 1]); a.set_yticklabels(['0', '1'])
a.set_xlabel('fraction relabelled'); a.set_ylabel('sets failing', labelpad=0.5)
tidy(a)

# 2  score
x, y = card(1, 0, 2, 'Score', 'Every clock scored in donor × cell-type pseudobulks', GREY,
            f"{FACTS['n_planned']} planned tests; {FACTS['n_planned'] - FACTS['n_estimable']} (crosses) are not "
            f"estimable, leaving {FACTS['n_estimable']}.",
            f"GSE136831: {FACTS['gse136831']}", False)
a = inset(x + GX, y + GY + 0.4, GW + 1.0, GH - 8.0)
nr, nc = FACTS['n_celltypes'], FACTS['n_modules']
ct = FACTS['celltypes']
for i in range(nr):
    for j in range(nc):
        hatch = j == 0 and any(c == ct[i] for c, _ in FACTS['unestimable'])
        if not hatch:
            a.add_patch(Rectangle((j + 0.1, i + 0.1), 0.8, 0.8, fc='#C4C4C4', ec='none'))
        else:                                       # non-estimable tests: a cross in an empty cell
            for xx in ([j + 0.2, j + 0.8], [j + 0.8, j + 0.2]):
                a.plot(xx, [i + 0.2, i + 0.8], color=ACCENT, lw=0.7, solid_capstyle='round')
a.set_xlim(0, nc); a.set_ylim(nr, 0); a.set_aspect('equal', anchor='N')
a.set_xticks([]); a.set_yticks([])
for sp in a.spines.values():
    sp.set_visible(False)
a.patch.set_alpha(0)
put(x + GX + (GW + 1.0) / 2, y + GY + GH - 6.4, f'{nr} cell types × {nc} modules', size=6.5, color=INK2,
    ha='center')

# 3  decompose
x, y = card(2, 0, 3, 'Decompose', 'Does an unchanged composite hide opposing module signals?', BLUE,
            f"In IPF the composite barely moves (N = {FACTS['N']:.2f}) while its components do "
            f"(S = {FACTS['S']:.2f}; {FACTS['n_comp_fdr']} of {FACTS['n_comp']} at FDR < 0.10).",
            'permutation null and donor bootstrap', False)
a = inset(x + GX + 8.4, y + GY + 1.0, GW - 8.4, GH - 6.4)
vals = [FACTS['S'], FACTS['N']]
a.barh([0, 1], vals, height=0.52, color=BLUE, zorder=2)
for k, v in enumerate(vals):
    a.text(v + 0.008, k, f'{v:.3f}', fontsize=6, color=INK, va='center', ha='left')
a.set_yticks([0, 1]); a.set_yticklabels(['S total', 'N net'])
a.set_ylim(1.6, -0.6); a.set_xlim(0, 0.40)
a.set_xticks([0, 0.2, 0.4]); a.set_xticklabels(['0', '0.2', '0.4'])
a.set_xlabel('Δlog10 hazard units')
a.tick_params(axis='y', length=0)
a.spines['left'].set_visible(False)
tidy(a)

# 4  claim gates
x, y = card(0, 1, 4, 'Claim gates', 'Is a cell-type × module effect robust, specific to IPF and to the module?',
            BLUE,
            'Robust inference, a COPD contrast and competitive nulls applied across all tests. '
            'No IPF effect passed all gates.',
            f"spike-ins: {FACTS['null_pass']} null pairs passed", True)
fx0, fy0 = x + GX + 0.4, y + GY - 0.4
for k, (n, lab) in enumerate(FACTS['funnel']):
    yy = fy0 + k * 4.2
    last = k == len(FACTS['funnel']) - 1
    put(fx0 + 6.4, yy, f'{n:,}', size=7.5 if last else 7, color=BLUE if last else INK, weight='bold',
        ha='right', va='top')
    put(fx0 + 7.8, yy, lab, size=7, color=BLUE if last else INK2, weight='bold' if last else 'normal',
        va='top')

# 5  resolution
x, y = card(1, 1, 5, 'Resolution', 'Does the clock carry age signal at the resolution it is read?', BLUE,
            f"Aging mouse lung: the composite recovers {FACTS['whole_lung']}% of the expected shift in whole "
            f"lung, {FACTS['celltype_lo']}–{FACTS['celltype_hi']}% within cell types.",
            '3- versus 24-month mice', True)
a = inset(x + GX + 9.0, y + GY + 1.0, GW - 9.0, GH - 6.4)
a.barh([0], [FACTS['whole_lung']], height=0.52, color=BLUE, zorder=2)
a.plot([FACTS['celltype_lo'], FACTS['celltype_hi']], [1, 1], color=BLUE, lw=1.2, solid_capstyle='butt', zorder=2)
for v in (FACTS['celltype_lo'], FACTS['celltype_hi']):
    a.plot([v, v], [0.78, 1.22], color=BLUE, lw=1.0, zorder=2)
a.axvline(100, color=MUTED, lw=0.6, ls=(0, (2, 1.5)), zorder=1)
a.text(FACTS['whole_lung'] + 3, 0, f"{FACTS['whole_lung']}%", fontsize=6, va='center', color=INK)
a.text(FACTS['celltype_hi'] + 3, 1, f"{FACTS['celltype_lo']}–{FACTS['celltype_hi']}%", fontsize=6,
       va='center', color=INK)
a.set_yticks([0, 1]); a.set_yticklabels(['whole lung', 'cell types'])
a.set_ylim(1.6, -0.6); a.set_xlim(0, 108)
a.set_xticks([0, 50, 100]); a.set_xticklabels(['0', '50', '100%'])
a.set_xlabel('expected shift recovered')
a.tick_params(axis='y', length=0)
a.spines['left'].set_visible(False)
tidy(a)

# 6  replication
x, y = card(2, 1, 6, 'Replication', 'Does the effect pattern hold in a second cohort?', BLUE,
            f"Only the NK effect pattern agrees across cohorts (r = {FACTS['cross_r'][0][1]:.2f}, "
            f"p = {FACTS['nk_p']}); no individual effect replicates.",
            f"{FACTS['n_sim_pairs']} simulated cohort pairs per scenario", True)
a = inset(x + GX + 1.0, y + GY + 0.2, GW - 1.0, GH - 5.0)
a.axvline(0, color=RULE, lw=0.5, zorder=1)
for k, (n, r) in enumerate(FACTS['cross_r']):
    nk = n == 'NK'
    a.plot([0, r], [k, k], color=ACCENT if nk else RULE, lw=0.7, zorder=1)
    a.scatter([r], [k], s=13, color=ACCENT if nk else MUTED, edgecolor='white', linewidth=0.5, zorder=3)
    if nk:
        a.text(r + 0.05, k, 'NK', fontsize=6.5, color=INK, va='center', fontweight='bold')
a.text(0.42, 3.4, 'six other\ncell types', fontsize=6, color=INK2, va='center', linespacing=1.1)
a.set_yticks([])
a.set_ylim(len(FACTS['cross_r']) - 0.4, -0.6); a.set_xlim(-0.15, 0.95)
a.set_xticks([0, 0.5]); a.set_xticklabels(['0', '0.5'])
a.set_xlabel('cross-cohort r')
a.spines['left'].set_visible(False)
tidy(a)

# ----------------------------------------------------------------------------- flow arrows
# inputs into steps 1 and 2
arrow(xs[2] + bw / 2, IN_Y + IN_H + 0.3, xs[2] + bw / 2, ROW_Y[0] - 0.3)
arrow(hx[1] + hw / 2, IN_Y + IN_H + 0.3, hx[1] + hw / 2, ROW_Y[0] - 0.3)
# along the rows
for row in (0, 1):
    yy = ROW_Y[row] + HEAD / 2
    for c in (0, 1):
        arrow(COLX[c] + COLW + 0.4, yy, COLX[c + 1] - 0.4, yy)
# from step 3 to step 4: down, left along the gap, down
y3 = ROW_Y[0] + CARD_H
ymid = y3 + 3.6
xs3 = COLX[2] + COLW / 2
xs4 = COLX[0] + COLW / 2
ax.plot([xs3, xs3, xs4], [y3 + 0.3, ymid, ymid], color=INK2, lw=0.8, solid_capstyle='butt',
        solid_joinstyle='miter', zorder=5)
arrow(xs4, ymid, xs4, ROW_Y[1] - 0.3)

from PIL import Image
fig.savefig(f'{OUT}.png', dpi=200)
fig.savefig(f'{OUT}.pdf')
tif = f'{OUT}.tif'
fig.savefig(tif, dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
im = Image.open(tif)
if im.mode != 'RGB':
    rgb = Image.new('RGB', im.size, (255, 255, 255))
    rgb.paste(im, mask=im.split()[-1] if im.mode == 'RGBA' else None)
    im = rgb
im.save(tif, compression='tiff_lzw', dpi=(600, 600))
print('saved', OUT, f'{W:.0f} x {H:.0f} mm', im.mode, im.size)

import pandas as pd
from overlap_audit import audit
rows = audit(fig, "Fig1")
overlap_path = REPO / "figures" / "v17" / "overlap_check.csv"
prev = pd.read_csv(overlap_path) if overlap_path.exists() else pd.DataFrame(columns=["figure", "kind", "text"])
prev = prev[prev.figure != "Fig1"] if len(prev) else prev
pd.concat([prev, pd.DataFrame(rows)], ignore_index=True).to_csv(overlap_path, index=False)
size_path = REPO / "figures" / "v17" / "size_check.csv"
size_row = dict(figure="Fig1", width_mm=W, height_mm=H, pdf_width_mm=W, pdf_height_mm=H,
                tiff_px_w=im.size[0], tiff_px_h=im.size[1], tiff_mode=im.mode)
sizes = pd.read_csv(size_path) if size_path.exists() else pd.DataFrame()
sizes = sizes[sizes.figure != "Fig1"] if len(sizes) else sizes
pd.concat([sizes, pd.DataFrame([size_row])], ignore_index=True).to_csv(size_path, index=False)
print("Fig1 overlaps", len(rows))

# ----------------------------------------------------------------------------- self-check
if '--check' in sys.argv:
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    items = []
    for a in fig.axes:
        arts = list(a.texts)
        if a is not ax:                              # the canvas axes is off; its ticks are not drawn
            arts += [a.xaxis.label, a.yaxis.label] + a.get_xticklabels() + a.get_yticklabels()
        for t in arts:
            if t.get_visible() and t.get_text().strip():
                items.append((t.get_text().replace('\n', ' ')[:28], t.get_window_extent(r)))
    fb = fig.bbox
    bad = 0
    for i in range(len(items)):
        n1, b1 = items[i]
        if b1.x0 < fb.x0 - 0.5 or b1.x1 > fb.x1 + 0.5 or b1.y0 < fb.y0 - 0.5 or b1.y1 > fb.y1 + 0.5:
            print('outside canvas:', n1); bad += 1
        for j in range(i + 1, len(items)):
            n2, b2 = items[j]
            if b1.overlaps(b2):
                print('text overlap:', repr(n1), '<->', repr(n2)); bad += 1
    print('texts checked:', len(items), 'problems:', bad)
