"""Build the final supplementary workbook (v17) on top of the v16 tables.

Changes: Klotho wild-type noise from the serialized clocks, MDE80 at full
precision, the red-module retrain check, S19_mde restricted to the second
experiment, permutation counts, and the S09_planned source.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "v18"))
import build_v16_tables as v16

OUT = REPO / "supplementary_revision" / "Supplementary_Tables.xlsx"
DS = (1.0, 1.5, 2.0)


def mde80(rates):
    """First d at which the pass rate reaches 0.8, linear between grid points (scripts/v9/v9_aggregate.py)."""
    if any(not np.isfinite(r) for r in rates):
        return np.nan, "not estimable"
    if rates[0] >= 0.8:
        return 1.0, ""
    for i in range(2):
        if rates[i] < 0.8 <= rates[i + 1]:
            span = rates[i + 1] - rates[i]
            return (DS[i] + (0.8 - rates[i]) / span * (DS[i + 1] - DS[i]) if span else DS[i + 1]), ""
    return np.nan, "> 2.0"


def klotho_noise(noise):
    w = pd.read_csv(REPO / "results" / "v3" / "rodent" / "klotho_official_wtwt.csv")
    mods = [c for c in w.columns if c not in ("tissue", "split") and not c.endswith("_div_sqrt2")]
    med = {t: float(np.median(g[mods].abs().to_numpy().ravel())) for t, g in w.groupby("tissue")}
    noise = noise.drop(columns=["median_abs_shift"], errors="ignore")
    noise["median_abs_shift_3v3"] = noise.tissue.map(med)
    noise["median_abs_shift_6v6"] = noise.median_abs_shift_3v3 / np.sqrt(2)
    return noise


def b_candidates(p, b_max=2000):
    out = []
    for b in range(1, b_max + 1):
        k = p * (b + 1) - 1
        if abs(k - round(k)) < 1e-9 and 0 <= round(k) <= b:
            out.append(b)
    return out


def planned():
    r = pd.read_csv(REPO / "results" / "v3" / "human" / "s12_ratios.csv")
    m = r[r.analysis == "main"][["celltype", "module", "coefficient", "p", "fdr"]].copy()
    info = pd.read_csv(REPO / "results" / "v3" / "human" / "r5" / "score_informativeness.csv")
    m = m.merge(info[["celltype", "module", "unestimable"]], on=["celltype", "module"], how="left")
    return pd.DataFrame({
        "celltype": m.celltype, "module": m.module, "coefficient": m.coefficient,
        "p_IPF": m.p, "fdr_IPF_253": m.fdr, "hit_253": m.fdr < 0.10,
        "unestimable": m.unestimable.fillna(False).astype(bool),
    })


def note(readme, sheet, description, source, n_rows, columns):
    row = {"sheet": sheet, "description": description, "source": source, "n_rows": n_rows, "columns": columns}
    if (readme.sheet == sheet).any():
        i = readme.index[readme.sheet == sheet][0]
        for k, v in row.items():
            readme.at[i, k] = v
        return readme
    return pd.concat([readme, pd.DataFrame([row])], ignore_index=True)


def main():
    book, readme = v16.build()
    readme = readme.reset_index(drop=True)

    book["S07_klotho_noise"] = klotho_noise(book["S07_klotho_noise"])
    readme = note(
        readme, "S07_klotho_noise",
        "Wild-type half-split noise of the serialized clocks after scaling by 1/sqrt(2). p2_5 and p97_5 are the 2.5 and 97.5 percentiles of the stored scaled draws; the uncorrected columns are the same percentiles before scaling. median_abs_shift_3v3 is the median absolute module difference between two halves of three wild-type animals, over all splits and the 23 modules; median_abs_shift_6v6 divides it by sqrt(2) to the scale of a six versus six contrast.",
        "results/v12/klotho_noise_band_v12.csv; results/v3/rodent/klotho_official_wtwt.csv",
        len(book["S07_klotho_noise"]), "; ".join(book["S07_klotho_noise"].columns))

    mde = book["S16_mde2"].reset_index(drop=True)
    vals = [mde80(r) for r in mde[["rate_d1.0", "rate_d1.5", "rate_d2.0"]].to_numpy(float)]
    mde["mde80"] = [v for v, _ in vals]
    mde["mde80_text"] = [t for _, t in vals]
    book["S16_mde2"] = mde
    readme = note(
        readme, "S16_mde2",
        "Second spike-in experiment, screen definition. mde80 is the first d at which the pass rate reaches 80%, by linear interpolation between d = 1.0, 1.5 and 2.0, recomputed at full precision from the stored rates (the stored mde80 column is rounded to three decimals). mde80 is blank when 80% is not reached by d = 2 (mde80_text > 2.0) or when a rate is missing (mde80_text not estimable). Rates are fractions.",
        "results/v9/screen_mde_v9.csv method=screen (rates)",
        len(mde), "; ".join(mde.columns))

    s19 = book["S19_mde"]
    s19 = s19[s19.scenario == "S1a"].drop(columns=["mde80", "experiment"], errors="ignore")
    lookup = mde[mde.scenario == "S1a"].set_index(["gate", "celltype"]).mde80
    s19["mde80"] = [lookup.get((g, c), np.nan) for g, c in zip(s19.system, s19.celltype)]
    s19["comparison"] = np.where(s19.abs_standardized >= s19.mde80, "above MDE", "below MDE")
    s19["mde_source"] = "spike-in 2"
    book["S19_mde"] = s19.reset_index(drop=True)
    readme = note(
        readme, "S19_mde",
        "Standardized effect of the two candidates against the S1a MDE80 of the second experiment (S16_mde2). comparison is above MDE when abs_standardized >= mde80.",
        "results/v8/candidates/candidate_power_v8.csv abs_standardized; S16_mde2 mde80",
        len(book["S19_mde"]), "; ".join(book["S19_mde"].columns))

    rc = pd.read_csv(REPO / "results" / "v19" / "retrain_check_v19.csv")
    book["S17_retrain_check"] = rc[["module", "n_genes", "spearman_rho_coefficients", "max_abs_diff_training_mean"]]
    readme = note(
        readme, "S17_retrain_check",
        "Genes of the serialized red mortality module refitted on the full rodent training matrix with the published pipeline and the serialized alpha and l1_ratio. spearman_rho_coefficients compares the refitted and serialized coefficients; max_abs_diff_training_mean is the largest absolute difference of the training means. No gene names or coefficients are listed.",
        "results/v19/retrain_check_v19.csv; scripts/v19/retrain_check.py",
        1, "module; n_genes; spearman_rho_coefficients; max_abs_diff_training_mean")

    rep = book["S13_replication"]
    cands = [b_candidates(float(p)) for p in rep.p_two_sided]
    rep["B_valid"] = [max(c) if c else np.nan for c in cands]
    rep["B_valid_unique"] = [len(c) == 1 for c in cands]
    book["S13_replication"] = rep
    hit = readme[readme.sheet == "S13_replication"]
    readme = note(
        readme, "S13_replication",
        "Module-effect correlation between cohorts. r is the stored correlation. Permutations in which a cell type lacked three donors per group in either cohort were dropped. B_valid is inferred from the stored p: the largest B up to 2,000 for which p = (1+k)/(B+1) with integer k. B_valid_unique is false where several B fit the stored p.",
        str(hit.source.iloc[0]) if len(hit) else "results/v3/human/r5/replication_detail_r5.csv",
        len(rep), "; ".join(rep.columns))

    sub = book["S14_substate"]
    rerun = pd.read_csv(REPO / "results" / "v19" / "substate_rerun_v19.csv")
    if not rerun.p_matches.all():
        raise SystemExit("sub-state rerun does not reproduce the stored p")
    sub["B_valid"] = sub.merge(rerun, on=["split", "celltype", "substate"], how="left").B_valid_y.to_numpy()
    book["S14_substate"] = sub
    readme = note(
        readme, "S14_substate",
        "Correlation of module-effect vectors between matched sub-states. r is Pearson. p_two_sided is the stored two-sided permutation p. Permutations in which a sub-state lacked three donors per group in either cohort were dropped; B is the number kept. B_valid is the number kept in a rerun with the stored seed, which reproduces every stored p.",
        "results/v3/human/r5/substate_concordance_r5.csv; results/v19/substate_rerun_v19.csv",
        len(sub), "; ".join(sub.columns))

    book["S09_planned"] = planned()
    readme = note(
        readme, "S09_planned",
        "All 253 planned tests. fdr_IPF_253 is the stored Benjamini-Hochberg FDR over the 253 tests. hit_253 is that FDR below 0.10. unestimable is the stored score-informativeness flag.",
        "results/v3/human/s12_ratios.csv analysis=main; results/v3/human/r5/score_informativeness.csv",
        len(book["S09_planned"]), "; ".join(book["S09_planned"].columns))

    names = [k for k in book if k != "S17_retrain_check"]
    last17 = max(i for i, k in enumerate(names) if k.startswith("S17_"))
    names.insert(last17 + 1, "S17_retrain_check")
    book = {k: book[k] for k in names}
    rank = {k: i for i, k in enumerate(names)}
    readme = readme.assign(_o=readme.sheet.map(rank).fillna(len(rank))).sort_values("_o", kind="stable").drop(columns="_o")

    v16.write(book, readme, OUT)
    print(mde[mde.mde80.notna() & (mde.scenario == "S1a")][["gate", "celltype", "mde80"]].to_string(index=False))
    print(book["S19_mde"].to_string(index=False))
    print(book["S07_klotho_noise"].to_string(index=False))
    print(rep[["coarse", "p_two_sided", "B_valid", "B_valid_unique"]].to_string(index=False))


if __name__ == "__main__":
    main()
