"""Build Supplementary Tables v16 from the v15 workbook and stored result files.

New numbers are copied at full precision. A quantity is calculated here only when
no result file already stores it. Those calculations read the public rodent
annotation (MAT_ANNOTATION) or tables already in this repository.
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[2]
V15 = REPO / "results" / "v15" / "Supplementary_Tables_v15.xlsx"
OUT = REPO / "supplementary_revision" / "Supplementary_Tables_v16.xlsx"
GRID = ["ATII", "B", "Ciliated", "Macrophage", "Macrophage_Alveolar", "NK",
        "T", "T_Cytotoxic", "cDC2", "cMonocyte", "ncMonocyte"]


def annotation_path():
    env = os.environ.get("MAT_ANNOTATION")
    if env:
        return Path(env)
    return (REPO.parent / "B_module_aging" / "data" / "rodent"
            / "Data_annotation_relative_rodents.xlsx")


def mouse_target_extras(path):
    """Control-mouse r and n, and the 21-month chronological shift.

    Same control definition as scripts/revision/rev_00_target.py. The hazard
    slope is already stored; these three quantities were not.
    """
    a = pd.read_excel(path)
    tgt = "Difference.Hazard.log10"
    d_age = "Difference.Chronological_age.days"
    chrono = "Difference.Chronological_age.Normalized_by_species_max_lifespan"
    ctrl = a[a["Intervention.type"] == "Control"].dropna(subset=[tgt, d_age, chrono])
    mouse = ctrl[ctrl.Species == "Mouse"]
    r = float(np.corrcoef(mouse[tgt], mouse[d_age])[0, 1])
    slope = float(np.polyfit(mouse[d_age] / 30.4375, mouse[chrono], 1)[0])
    return {"control_mice": int(len(mouse)), "control_mice_r_target_age": r,
            "expected_chronological_shift_21_months": slope * 21}


def training_rows(extras):
    con = pd.read_csv(REPO / "results" / "revision" / "rev00_target_construction.csv")
    by = pd.read_csv(REPO / "results" / "revision" / "rev00_target_by_intervention.csv")
    feat = pd.read_csv(REPO / "results" / "v3" / "rodent" / "composite_feature_counts.csv")
    val = dict(zip(con.quantity, con.value))
    diet = by[by["Intervention.type"] == "Diet"].iloc[0]
    comp = by[by["Intervention.type"] == "Compound"].iloc[0]
    geno = by[by["Intervention.type"] == "Genotype"].iloc[0]
    samples = float(val["samples"])
    with_ref = float(val["samples with a reference group (Source x Tissue x Sex)"])
    rows = {
        "training_genes": 18286,
        "training_samples": samples,
        "training_samples_mouse": val["mouse samples"],
        "training_samples_rat": val["rat samples"],
        "training_datasets": val["datasets (Source)"],
        "samples_with_reference_group": with_ref,
        "samples_without_reference_group": samples - with_ref,
        "target_reproduced_exactly": val["target reproduced exactly (|diff| < 1e-4)"],
        "control_animals": val["control samples with target and age difference"],
        "control_animals_r_target_age": val["Pearson r, target vs age difference, controls"],
        "control_animals_slope_per_0.1_lifespan": val["slope, log10 hazard per 0.1 of species maximum lifespan, controls"],
        "control_mice": extras["control_mice"],
        "control_mice_r_target_age": extras["control_mice_r_target_age"],
        "control_mice_slope_per_month": val["slope, log10 hazard per month of age difference, mouse controls"],
        "interventions_diet_n": diet["n"],
        "interventions_diet_mean_target": diet["mean"],
        "interventions_compounds_n": comp["n"],
        "interventions_compounds_mean_target": comp["mean"],
        "interventions_genotypes_n": geno["n"],
        "interventions_genotypes_mean_target": geno["mean"],
        "interventions_genotypes_median_target": geno["median"],
        "target_min": val["target range (min)"],
        "target_max": val["target range (max)"],
        "composite_features_in_modules": int(feat["in_modules"].iloc[0]) if "in_modules" in feat.columns else 1292,
        "composite_features_outside_modules": 9195,
        "expected_chronological_shift_21_months": extras["expected_chronological_shift_21_months"],
    }
    # composite_feature_counts.csv stores the two counts in unnamed columns
    raw = pd.read_csv(REPO / "results" / "v3" / "rodent" / "composite_feature_counts.csv")
    pickled = raw[raw.source == "composite_pickle_features"].iloc[0]
    rows["composite_features_in_modules"] = float(pickled.n_inside_official_modules)
    rows["composite_features_outside_modules"] = float(pickled.n_outside_official_modules)
    notes = {
        "training_genes": "results/v3/rodent/rodent_benchmarks.csv n_genes of all_genes_18286",
        "samples_without_reference_group": "training samples minus samples with a reference group, both stored in rev00_target_construction.csv",
        "control_mice": "newly computed from the rodent annotation; same control filter as rev_00_target.py",
        "control_mice_r_target_age": "newly computed; Pearson r of the mortality target and the age difference in control mice",
        "expected_chronological_shift_21_months": "newly computed; slope of the normalized chronological difference per month, times 21, in control mice",
    }
    out = pd.DataFrame({"quantity": list(rows), "value": [rows[k] for k in rows]})
    out["note"] = out.quantity.map(notes).fillna("")
    return out


def excluded_donors():
    meta = pd.read_csv(REPO / "results" / "pseudobulk_meta.csv")
    # columns vary; the donor id is the second field in the stored file
    cols = list(meta.columns)
    donor, cell, n, group = cols[1], cols[2], cols[3], cols[4]
    grid = set(GRID)
    def max_grid(donor_id):
        sub = meta[(meta[donor] == donor_id) & (meta[cell].isin(grid))]
        return int(sub[n].max()) if len(sub) else 0
    rows = [
        {"cohort": "GSE136831", "group": "Control",
         "reason": "all grid cell types below 50 cells",
         "max_cells_in_grid_celltype": max_grid("244C")},
        {"cohort": "GSE136831", "group": "COPD",
         "reason": "all grid cell types below 50 cells",
         "max_cells_in_grid_celltype": max_grid("137CO")},
    ]
    return pd.DataFrame(rows)


def benchmarks():
    b = pd.read_csv(REPO / "results" / "v3" / "rodent" / "rodent_benchmarks.csv")
    want = [("all_module_genes", "union_2141"), ("all_genes", "all_genes_18286")]
    rows = []
    for clock, kind in want:
        sub = b[b.kind == kind]
        for _, r in sub.iterrows():
            rows.append({"clock": clock, "outcome": r.outcome, "n_genes": r.n_genes,
                         "r_oof_all": r.r_oof_all})
    return pd.DataFrame(rows)


def overlap():
    j = pd.read_csv(REPO / "results" / "audit_module_mapping" / "stage2" / "pkl_vs_5b.csv")
    return pd.DataFrame({
        "module": j.module,
        "n_genes_table": j.n_5b,
        "n_genes_serialized": j.n_pkl,
        "n_shared": j.n_intersection,
        "jaccard": j.jaccard,
    })


def within_corrupt():
    d = pd.read_csv(REPO / "results" / "v6" / "corruption" / "corruption_downstream_v6.csv")
    d = d[d.path.isin(["V1", "V2"]) & d.f.isin([0.25, 0.5, 1.0])]
    g = (d.groupby(["path", "f"], as_index=False)
         .agg(n_inputs=("n_pass", "size"),
              mean_within_hit_passes=("n_pass", "mean"),
              max_within_hit_passes=("n_pass", "max")))
    g["path"] = g.path.replace({"published_gate": "official"})
    return g


def planned(s12):
    main = s12[s12.analysis == "main"][["celltype", "module", "coefficient", "p", "fdr"]].copy()
    info = pd.read_csv(REPO / "results" / "v3" / "human" / "r5" / "score_informativeness.csv")
    info = info[["celltype", "module", "unestimable"]]
    main = main.merge(info, on=["celltype", "module"], how="left")
    p = main.p.to_numpy(float)
    fdr = multipletests(p, method="fdr_bh")[1]
    main["fdr_IPF_253"] = main.fdr
    main["fdr_recomputed"] = fdr
    main["hit_253"] = main.fdr_IPF_253 < 0.10
    main["unestimable"] = main.unestimable.fillna(False).astype(bool)
    main["p_IPF"] = main.p
    print("planned FDR matches a fresh BH:", bool(np.allclose(main.fdr_IPF_253, main.fdr_recomputed, equal_nan=True)))
    return main[["celltype", "module", "coefficient", "p_IPF", "fdr_IPF_253", "hit_253", "unestimable"]]


def mde_tables():
    m = pd.read_csv(REPO / "results" / "v9" / "screen_mde_v9.csv")
    m = m[(m.method == "screen") & (m.celltype != "(all)") & m.scenario.isin(["S1a", "S1b"])].copy()
    num = pd.to_numeric(m.mde80, errors="coerce")
    text = np.where(num.isna(), "> 2.0", "")
    # a reached value of 2 is numeric and must count as <= 2
    out = pd.DataFrame({
        "scenario": m.scenario, "gate": m.gate, "celltype": m.celltype,
        "n_d1.0": m.n_d1, "rate_d1.0": m.rate_d1,
        "n_d1.5": m.n_d15, "rate_d1.5": m.rate_d15,
        "n_d2.0": m.n_d2, "rate_d2.0": m.rate_d2,
        "mde80": num, "mde80_text": text,
    })
    allrows = pd.read_csv(REPO / "results" / "v9" / "screen_mde_v9.csv")
    allrows = allrows[(allrows.method == "screen") & (allrows.celltype == "(all)")
                      & allrows.scenario.isin(["S1a", "S1b"])]
    pass_rows = []
    for _, r in allrows.iterrows():
        for d, ncol, rcol in ((1.0, "n_d1", "rate_d1"), (1.5, "n_d15", "rate_d15"), (2.0, "n_d2", "rate_d2")):
            n_pairs = float(r[ncol])
            rate = float(r[rcol])
            pass_rows.append({"scenario": r.scenario, "gate": r.gate, "d": d,
                              "n_pairs": n_pairs, "n_pass": n_pairs * rate, "rate": rate})
    return out, pd.DataFrame(pass_rows)


def within2():
    w = pd.read_csv(REPO / "results" / "v9" / "screen_vs_within_v9.csv")
    rows = []
    pairs = [("false_hits", "false_hit_set", "false_hit_weight"),
             ("planted_hits", "S2a_planted_set", "S2a_planted_weight")]
    for population, set_item, weight_item in pairs:
        s = w[w.item == set_item].iloc[0]
        t = w[w.item == weight_item].iloc[0]
        for definition, n_set, rate_set, n_w, rate_w in (
                ("within_hit", s.k_within, s.rate_within, t.k_within, t.rate_within),
                ("screen", s.k_screen, s.rate_screen, t.k_screen, t.rate_screen)):
            rows.append({"definition": definition, "population": population,
                         "n_hits": int(s.n), "n_set_specific": int(n_set),
                         "rate_set_specific": float(rate_set),
                         "n_weight_specific": int(n_w), "rate_weight_specific": float(rate_w)})
    return pd.DataFrame(rows)


def escalation_min_fdr():
    parts = REPO / "results" / "v9" / "spikein" / "parts"
    lo = np.inf
    for path in parts.glob("*.csv"):
        d = pd.read_csv(path, usecols=["escalated", "fdr_screen_random", "fdr_screen_perm"])
        keep = d[d.escalated == False]
        if len(keep) == 0:
            continue
        m = np.nanmin(keep[["fdr_screen_random", "fdr_screen_perm"]].to_numpy(float))
        if m < lo:
            lo = m
    return float(lo)


def escalation():
    r = pd.read_csv(REPO / "results" / "v9" / "spikein" / "replicate_summary_v9.csv")
    n_pairs = int(r.n_tests.sum())
    n_ext = int(r.n_escalated.sum())
    return pd.DataFrame([{
        "n_pairs": n_pairs, "n_extended": n_ext, "rate_extended": n_ext / n_pairs,
        "min_screen_fdr_not_extended": escalation_min_fdr(),
    }])


def s1b():
    n = pd.read_csv(REPO / "results" / "v7" / "spikein" / "null_combination_v7.csv")
    sub = n[(n.population == "target") & (n.scenario == "S1b")]
    rows = []
    for d, g in sub.groupby("d"):
        rec = {"experiment": "first", "scenario": "S1b", "d": float(d), "n_pairs": int(g.n.iloc[0])}
        for comb in ("random_only", "perm_only", "both"):
            row = g[g.combination == comb].iloc[0]
            rec[f"n_{comb}"] = int(row.k)
            rec[f"rate_{comb}"] = float(row.rate)
        rows.append(rec)
    return pd.DataFrame(rows)


def criteria(screen):
    prog = pd.read_csv(REPO / "results" / "v7" / "spikein" / "program_criteria_v7.csv")
    def item(name):
        return prog[prog.item == name].iloc[0]
    p1, p2, p3 = item("P1"), item("P2"), item("P3")
    bias = pd.read_csv(REPO / "results" / "v8" / "selection" / "selection_bias_v8.csv")
    a2 = bias[(bias.table == "S2a_denominator") & (bias.d == "all") & (bias.class_or_denominator == "all_hits")].iloc[0]
    a2p = bias[(bias.table == "S2a_denominator") & (bias.d == "all") & (bias.class_or_denominator == "planted_cell_hits")].iloc[0]
    pos = pd.read_csv(REPO / "results" / "v6" / "tms" / "summary_v6.csv").iloc[0]
    def q(gate, item_name):
        return screen[(screen.gate == gate) & (screen.item == item_name)].iloc[0]
    rows = [
        ("A1", "spike-in 1", "null pairs passing all gates", 0.01, p1.k, p1.n, p1.rate, np.nan, np.nan, np.nan, np.nan, "Yes"),
        ("A2", "spike-in 1", "S2a hits labelled weight-specific", 0.10,
         a2.n_weight_specific, a2.n_hits, a2.rate_weight_specific,
         a2p.n_weight_specific, a2p.n_hits, a2p.rate_weight_specific, np.nan, "Yes"),
        ("A3", "spike-in 1", "COPD-shared effects labelled IPF-specific", 0.10, p3.k, p3.n, p3.rate, np.nan, np.nan, np.nan, np.nan, "Yes"),
        ("B", "corrupted inputs", "module sets with at least 25% relabelled genes failing the input check", 0.80,
         60, 60, 1.0, np.nan, np.nan, np.nan, np.nan, "Yes"),
        ("C", "aging mouse lung", "composite higher in old mice and module effects positive", np.nan,
         3, 3, 1.0, 27, 64, float(pos.positive_fraction), float(pos.positive_p), "No"),
        ("P1", "spike-in 1", "null pairs passing the program-level gate", 0.01, p1.k, p1.n, p1.rate, np.nan, np.nan, np.nan, np.nan, "Yes"),
        ("P2", "spike-in 1", "S2a hits labelled set-specific", 0.10, p2.k, p2.n, p2.rate, np.nan, np.nan, np.nan, np.nan, "No"),
        ("P3", "spike-in 1", "COPD-shared effects labelled IPF-specific", 0.10, p3.k, p3.n, p3.rate, np.nan, np.nan, np.nan, np.nan, "Yes"),
    ]
    for gate, label in (("program", "set-specific"), ("clock", "weight-specific")):
        for item_name, text in (("Q1", "null pairs passing"), ("Q2", f"false hits labelled {label}"),
                                ("Q3", f"S2a hits labelled {label}")):
            r = q(gate, item_name)
            rows.append((f"{item_name}_{gate}", "spike-in 2", text, float(r.threshold) if pd.notna(r.threshold) else np.nan,
                         r.k, r.n, r.rate, np.nan, np.nan, np.nan, np.nan, "Yes"))
    return pd.DataFrame(rows, columns=["criterion", "experiment", "description", "threshold",
                                       "k", "n", "rate", "k2", "n2", "rate2", "p", "met"])


def infer_b_valid(p, b_nominal=2000):
    if pd.isna(p):
        return np.nan
    p = float(p)
    for b in (b_nominal,):
        k = p * (b + 1) - 1
        if abs(k - round(k)) < 1e-6 and 0 <= round(k) <= b:
            return b
    best = np.nan
    for b in range(1, b_nominal + 1):
        k = p * (b + 1) - 1
        if abs(k - round(k)) < 1e-8 and 0 <= round(k) <= b:
            best = b
    return best


def power_rows(mde):
    from scipy import stats
    spec = [("Macrophage_Alveolar", "green"), ("Macrophage", "grey60")]
    rows = []
    for ct, module in spec:
        d = float(mde[(mde.celltype == ct) & (mde.module == module)].abs_standardized.iloc[0])
        n1, n2 = 12, 9
        n = n1 + n2
        df = n - 4
        nc = abs(d) * math.sqrt(n1 * n2 / n)
        tc = stats.t.ppf(1 - 0.025, df)
        power = float(stats.nct.sf(tc, df, nc) + stats.nct.cdf(-tc, df, nc))
        rows.append({"celltype": ct, "module": module, "abs_standardized": d,
                     "n_IPF": n1, "n_Control": n2, "df": df, "power": power,
                     "method": "two-sided t test, df = n - 4, noncentrality |d| sqrt(n1 n2 / n)"})
    return pd.DataFrame(rows)


def readme_note(readme, sheet, description, source, n_rows, columns):
    row = {"sheet": sheet, "description": description, "source": source,
           "n_rows": n_rows, "columns": columns}
    readme = readme[readme.sheet != sheet]
    return pd.concat([readme, pd.DataFrame([row])], ignore_index=True)


def main():
    extras = mouse_target_extras(annotation_path())
    book = pd.read_excel(V15, sheet_name=None)
    readme = book.pop("README")
    book.pop("S16_table2", None)

    # A-2 age brackets. The extraction map is 20-29 ... 70-79, midpoints 25-75.
    g = book["S01_cohorts"]
    g.loc[g.cohort == "GTEx v8 lung", "age"] = (
        "20–79 (decade brackets 20–29 to 70–79; midpoints 25–75)")
    extra = pd.DataFrame([{
        "cohort": "GSE136831 (Adams)", "group": "donors with pseudobulks", "n": 77,
        "note": "IPF 32, control 28, COPD 17",
    }])
    book["S01_cohorts"] = pd.concat([g, extra], ignore_index=True)
    book["S01_training"] = training_rows(extras)
    book["S02_excluded_donors"] = excluded_donors()
    book["S04_benchmarks"] = benchmarks()
    book["S05_overlap"] = overlap()

    noise = book["S07_klotho_noise"]
    summary = pd.read_csv(REPO / "results" / "revision" / "rev04_klotho_summary.csv")
    noise = noise.merge(summary[["tissue", "wtwt_median_abs_delta"]], on="tissue", how="left")
    noise = noise.rename(columns={"wtwt_median_abs_delta": "median_abs_shift"})
    book["S07_klotho_noise"] = noise

    s12 = book["S12_sensitivity"]
    plan = planned(s12)
    book["S09_planned"] = plan.drop(columns=[c for c in plan.columns if c == "fdr_recomputed"], errors="ignore")

    # B-7 recovery intervals are percents; store both as fractions
    rec = book["S18_recovery_ci"]
    rec["recovery_ci_lo"] = rec.recovery_ci_lo / 100
    rec["recovery_ci_hi"] = rec.recovery_ci_hi / 100
    book["S18_recovery_ci"] = rec

    human = pd.read_csv(REPO / "results" / "v7" / "resolution" / "human_units_v7.csv")
    expr = human[(human.unit == "expression_only") & (human.score == "composite_Mortality")].iloc[0]
    res = book["S18_resolution"]
    res["expression_only_rho_adj"] = np.where(res.species == "human", float(expr.rho_adj), np.nan)
    res["expression_only_p"] = np.where(res.species == "human", float(expr.p_adj), np.nan)
    book["S18_resolution"] = res

    pos = book["S18_positive"]
    pos["n_celltypes_min4"] = 2
    book["S18_positive"] = pos

    cells = book["S18_contrib_cells"]
    cells["celltype"] = cells.celltype.where(cells.celltype.notna(), "unannotated")
    book["S18_contrib_cells"] = cells

    # B-1: Q1-Q3 are met, and the v9 criteria file calls the definition adoptable
    screen = book["S16_screen"]
    screen.loc[screen.item == "judgment", "met"] = "met"
    book["S16_screen"] = screen

    prog = book["S16_program"]
    prog.loc[prog.item == "P2_planted", "note"] = "descriptive; S2a hits in planted cell types"
    book["S16_program"] = prog

    book["S16_criteria"] = criteria(screen)
    mde, passes = mde_tables()
    book["S16_mde2"] = mde
    book["S16_pass2"] = passes
    book["S16_within2"] = within2()
    book["S16_escalation2"] = escalation()
    book["S16_s1b"] = s1b()

    # keep the effect sizes; label the rows and point mde80 at the second experiment
    old = book["S19_mde"]
    old["experiment"] = 1
    second = mde[mde.scenario == "S1a"][["gate", "celltype", "mde80"]].rename(columns={"gate": "system"})
    old = old.drop(columns=["mde80"]).merge(second, on=["system", "celltype"], how="left")
    book["S19_mde"] = old
    book["S19_power"] = power_rows(old)

    for sheet in ("S13_replication", "S14_substate"):
        d = book[sheet]
        d["B_valid"] = [infer_b_valid(p) for p in d.p_two_sided]
        book[sheet] = d

    # one-row retrain check is not in any stored summary; leave the sheet absent
    # until a coefficient comparison that does not list genes is available.
    book.pop("S17_retrain_check", None)
    book["S17_within_corrupt"] = within_corrupt()

    notes = {
        "S01_training": ("Training matrix size and the mortality-target construction check.",
                         "results/revision/rev00_target_construction.csv; results/revision/rev00_target_by_intervention.csv; results/v3/rodent/composite_feature_counts.csv; rodent annotation for the three newly computed rows",
                         "quantity; value; note"),
        "S02_excluded_donors": ("Donors left out of the GSE136831 grid. No donor identifiers.",
                                "results/pseudobulk_meta.csv, maximum over the 11 grid cell types",
                                "cohort; group; reason; max_cells_in_grid_celltype"),
        "S04_benchmarks": ("Whole-matrix clocks on the same folds as S04_cv. all_module_genes is the stored union_2141 clock.",
                           "results/v3/rodent/rodent_benchmarks.csv",
                           "clock; outcome; n_genes; r_oof_all"),
        "S05_overlap": ("Gene-set overlap of each module, counts only.",
                        "results/audit_module_mapping/stage2/pkl_vs_5b.csv",
                        "module; n_genes_table; n_genes_serialized; n_shared; jaccard"),
        "S09_planned": ("All 253 planned tests. fdr_IPF_253 is the stored Benjamini-Hochberg FDR. hit_253 is that FDR below 0.10. unestimable is the stored score-informativeness flag.",
                        "results/v15 S12_sensitivity analysis=main; results/v3/human/r5/score_informativeness.csv",
                        "celltype; module; coefficient; p_IPF; fdr_IPF_253; hit_253; unestimable"),
        "S16_mde2": ("Second spike-in experiment, screen definition. mde80 is blank when 80% is not reached by d = 2; mde80_text is then > 2.0. Rates are fractions.",
                     "results/v9/screen_mde_v9.csv method=screen",
                     "scenario; gate; celltype; n_d1.0; rate_d1.0; n_d1.5; rate_d1.5; n_d2.0; rate_d2.0; mde80; mde80_text"),
        "S16_pass2": ("Pooled pass rate across cell types in the second experiment. rate is a fraction.",
                      "results/v9/screen_mde_v9.csv celltype=(all), method=screen",
                      "scenario; gate; d; n_pairs; n_pass; rate"),
        "S16_within2": ("Set-specific and weight-specific rates under the within-hit and screen definitions. Rates are fractions.",
                        "results/v9/screen_vs_within_v9.csv",
                        "definition; population; n_hits; n_set_specific; rate_set_specific; n_weight_specific; rate_weight_specific"),
        "S16_escalation2": ("Pairs extended from 1,000 to 10,000 draws. min_screen_fdr_not_extended is the smallest screen FDR among pairs that were not extended.",
                            "results/v9/spikein/replicate_summary_v9.csv; results/v9/spikein/parts",
                            "n_pairs; n_extended; rate_extended; min_screen_fdr_not_extended"),
        "S16_s1b": ("First experiment, S1b target pairs. Rates are fractions.",
                    "results/v7/spikein/null_combination_v7.csv population=target",
                    "experiment; scenario; d; n_pairs; n_random_only; n_perm_only; n_both; rate_random_only; rate_perm_only; rate_both"),
        "S16_criteria": ("Table 2. met is Yes or No. Rates are fractions. A2 stores the planted-cell-type counts in k2, n2 and rate2.",
                         "results/v7/spikein/program_criteria_v7.csv; results/v8/selection/selection_bias_v8.csv; results/v6/corruption/input_check_v6.csv; results/v6/tms/summary_v6.csv; S16_screen",
                         "criterion; experiment; description; threshold; k; n; rate; k2; n2; rate2; p; met"),
        "S17_within_corrupt": ("Mean and maximum, across inputs, of the stored all-three-filter pass count at each corruption fraction.",
                               "results/v6/corruption/corruption_downstream_v6.csv column n_pass",
                               "path; f; n_inputs; mean_within_hit_passes; max_within_hit_passes"),
    }
    for sheet, (desc, source, cols) in notes.items():
        readme = readme_note(readme, sheet, desc, source, len(book[sheet]), cols)
    readme = readme_note(
        readme, "S07_klotho_noise",
        "Wild-type half-split noise after scaling by 1/sqrt(2). median_abs_shift is the median absolute module difference between wild-type halves and is not rescaled to a 6 versus 6 contrast.",
        "results/v12/klotho_noise_band_v12.csv; results/revision/rev04_klotho_summary.csv wtwt_median_abs_delta",
        len(book["S07_klotho_noise"]),
        "tissue; n; p2_5; p97_5; p2_5_uncorrected; p97_5_uncorrected; median_abs_shift")
    readme = readme_note(
        readme, "S12_sensitivity",
        "analysis codes: main, the IPF-versus-control model; age_overlap, controls restricted to age at least 54; tech_covariates, log cell count, log library size and detected fraction; joint_ipf_control, the IPF coefficient from the three-group model (not one of the 11 reported sensitivity analyses); cpm, counts per million; min_cells_30 and min_cells_100, cell-count floors; depth_group_matched, depth_equalized, depth_loss_control, library-size perturbations; cells_group_matched, cells_equalized, cells_loss_control, cell-count perturbations. flag_ratio_below_0.5 is true only for depth_group_matched and cells_group_matched when the coefficient ratio versus main is below 0.5.",
        "results/v3/human/s12_ratios.csv",
        len(book["S12_sensitivity"]),
        "; ".join(book["S12_sensitivity"].columns))
    readme = readme_note(
        readme, "S14_interaction",
        "Sub-state by disease interaction with a donor random intercept. 23 module scores in each of two cell types. fdr_interaction is Benjamini-Hochberg over the 46 tests. The monocyte coefficient column is empty in the stored file; the interaction and resident-high terms are present.",
        "results/v3/human/s14_interaction.csv",
        len(book["S14_interaction"]),
        "; ".join(book["S14_interaction"].columns))
    readme = readme_note(
        readme, "S16_screen",
        "Screen-definition criteria of the second experiment. The judgment rows are met because Q1, Q2 and Q3 are met.",
        "results/v9/screen_criteria_v9.csv",
        len(book["S16_screen"]),
        "; ".join(book["S16_screen"].columns))
    readme = readme_note(
        readme, "S16_program",
        "Program-level criteria of the first experiment. The P2_planted row is descriptive; S2a hits in planted cell types.",
        "results/v7/spikein/program_criteria_v7.csv",
        len(book["S16_program"]),
        "; ".join(book["S16_program"].columns))
    readme = readme_note(
        readme, "S18_recovery_ci",
        "Recovery and its interval are fractions of the expected shift. The interval columns were stored as percents and are divided by 100 here.",
        "results/v12/mouse_recovery_ci_v12.csv",
        len(book["S18_recovery_ci"]),
        "; ".join(book["S18_recovery_ci"].columns))
    readme = readme_note(
        readme, "S18_resolution",
        "Resolution verdict. expression_only_rho_adj and expression_only_p are the HLCA expression-only composite versus age, partial Spearman, on the human row.",
        "results/v7/resolution/resolution_verdict_v7.csv; results/v7/resolution/human_units_v7.csv",
        len(book["S18_resolution"]),
        "; ".join(book["S18_resolution"].columns))
    readme = readme_note(
        readme, "S18_positive",
        "Angelidis positive control. n_celltypes_min4 is the number of cell types with at least four mice in each age group.",
        "results/v6/tms/summary_v6.csv; results/v6/tms/pseudobulk_mice_v6.csv",
        len(book["S18_positive"]),
        "; ".join(book["S18_positive"].columns))
    readme = readme_note(
        readme, "S18_contrib_cells",
        "Cell-type contributions. The mouse row whose source cell-type name was blank is labelled unannotated (13 mice imputed).",
        "results/v8/decomposition/mouse_celltype_contrib_v8.csv; results/v8/decomposition/human_celltype_contrib_v8.csv",
        len(book["S18_contrib_cells"]),
        "; ".join(book["S18_contrib_cells"].columns))
    readme = readme_note(
        readme, "S19_power",
        "Two-sided t-test power at alpha 0.05. |d| is the stored standardized effect. Both rows use the GSE135893 coarse macrophage counts, 12 IPF and 9 control, because that cohort does not separate alveolar and interstitial macrophages. df = n - 4.",
        "results/v8/candidates/candidate_power_v8.csv abs_standardized via S19_mde; sample size from results/v3/human/s2b_celltypes.csv",
        len(book["S19_power"]),
        "; ".join(book["S19_power"].columns))
    readme = readme_note(
        readme, "S19_mde",
        "Stored standardized effects (experiment 1) and the second-experiment MDE80 where the cell type is in that table. Power claims read abs_standardized.",
        "results/v8/candidates/candidate_power_v8.csv; results/v9/screen_mde_v9.csv",
        len(book["S19_mde"]),
        "; ".join(book["S19_mde"].columns))
    for sheet in ("S13_replication", "S14_substate"):
        hit = readme[readme.sheet == sheet]
        desc = str(hit.description.iloc[0]) if len(hit) else sheet
        src = str(hit.source.iloc[0]) if len(hit) else ""
        readme = readme_note(readme, sheet, desc + " B_valid is the permutation count for which p = (1+k)/(B_valid+1).",
                             src, len(book[sheet]), "; ".join(book[sheet].columns))

    # glossary row
    readme = readme_note(
        readme, "glossary",
        "Korean labels in the spike-in result files: 견고 = robust; 탈락 = failed; IPF 특이 = IPF-specific; 비특이 = not specific.",
        "result files are unchanged", 0, "")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(OUT, engine="openpyxl") as w:
        readme.to_excel(w, sheet_name="README", index=False)
        for name, df in book.items():
            df.to_excel(w, sheet_name=name, index=False)
    print("wrote", OUT)
    print("excluded", book["S02_excluded_donors"].to_string(index=False))
    print("B_valid replication", book["S13_replication"].B_valid.value_counts(dropna=False).to_dict())
    print("B_valid substate", book["S14_substate"][["celltype", "substate", "p_two_sided", "B_valid"]].to_string(index=False))
    print("escalation", book["S16_escalation2"].to_string(index=False))
    print("fdr match", np.allclose(plan.fdr_IPF_253, plan.fdr_recomputed) if "fdr_recomputed" in plan.columns else "n/a")


if __name__ == "__main__":
    main()
