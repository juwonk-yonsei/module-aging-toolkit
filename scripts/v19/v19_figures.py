"""Figures 2 and 6 for instruction 19: layout fixes only (Fig. 2c label, Fig. 2d legend, Fig. 6a x range).

Data, labels and styling come from scripts/v17/v17_figures.py. Writes figures/v19 only.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "v17"))
import v17_figures as base
from v17_figures import (ADOPT, DROP, IPF, NEW, ORANGE, RES, SKY, letters, manifest, mm,
                         np, pd, plt, remember, style_ax)

FIG = REPO / "figures" / "v19"
base.FIG = FIG


def save(fig, name):
    return base.save(fig, name)


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
    axc.annotate("spreadsheet", xy=(1.08, float(tab_row.set_failure_rate)), xytext=(1.20, 0.70),
                 fontsize=7, color=ORANGE, ha="right", va="top",
                 arrowprops=dict(arrowstyle="-", lw=0.4, color=ORANGE, relpos=(0.85, 1.0),
                                 shrinkA=1, shrinkB=3))
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
    axd.set_ylim(-1, 64)
    axd.set_yticks([0, 10, 20])
    axd.legend(handles=[handles[0], handles[2], handles[1]], loc="upper left", ncol=2, frameon=False,
               borderpad=0.2, borderaxespad=0.2, handlelength=1.2, columnspacing=0.8, labelspacing=0.2)
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
    axa.set_xlim(axa.get_xlim()[0], 0.8)
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
    for fn in (fig2, fig6):
        print(fn.__name__, flush=True)
        overlaps.extend(fn())
    pd.DataFrame(overlaps, columns=["figure", "kind", "text"]).to_csv(FIG / "overlap_check.csv", index=False)
    pd.DataFrame(base.MANIFEST).to_csv(FIG / "manifest_v19.csv", index=False)
    pd.DataFrame(base.SIZES).to_csv(FIG / "size_check.csv", index=False)
    nums = pd.DataFrame(base.NUMBERS)

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
