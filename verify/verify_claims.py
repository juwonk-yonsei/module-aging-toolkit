#!/usr/bin/env python3
"""Check the numbers printed in the manuscript against the shipped result tables.

Each check reads a value from ``results/`` and compares it with the value printed in the
manuscript at the printed precision: a value printed as 0.57 passes if the table holds a value
that rounds to 0.57. No raw data, GPU or network access is needed, and the script runs in seconds.

Usage
-----
    python verify/verify_claims.py                # revised manuscript (default)
    python verify/verify_claims.py --submission   # original submission (v0.1.0)
    MAT_RESULTS=/path/to/results python verify/verify_claims.py

The exit code is non-zero if any check fails.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS = Path(os.environ.get("MAT_RESULTS", Path(__file__).resolve().parents[1] / "results"))
REVISION = RESULTS / "revision"
MIL = RESULTS / "mil"


def rev(name: str, **kw) -> pd.DataFrame:
    return pd.read_csv(REVISION / name, **kw)


class Checker:
    def __init__(self) -> None:
        self.rows: list[tuple[str, str, str, bool | None]] = []

    def section(self, title: str) -> None:
        self.rows.append((title, "", "", None))

    def printed(self, desc: str, observed, value: float, decimals: int) -> None:
        """Pass if `observed` rounds to `value` printed with `decimals` decimals."""
        obs = float(observed)
        ok = bool(np.isfinite(obs)) and abs(obs - value) <= 0.5 * 10.0 ** -decimals + 1e-9
        self.rows.append((desc, f"{obs:.{decimals + 2}f}", f"{value:.{decimals}f}", ok))

    def sci(self, desc: str, observed, mantissa: int, exponent: int) -> None:
        """Pass if `observed` is printed as mantissa x 10^exponent (one significant digit)."""
        obs = float(observed)
        ok = round(obs / 10.0 ** exponent) == mantissa
        self.rows.append((desc, f"{obs:.2e}", f"{mantissa}e{exponent}", ok))

    def equal(self, desc: str, observed, expected) -> None:
        if isinstance(expected, (int, np.integer)) and not isinstance(expected, bool):
            observed = int(observed)
        elif isinstance(expected, bool):
            observed = bool(observed)
        self.rows.append((desc, str(observed), str(expected), observed == expected))

    def near(self, desc: str, observed, expected: float, tol: float) -> None:
        obs = float(observed)
        self.rows.append((desc, f"{obs:.3f}", f"{expected:.3f}", abs(obs - expected) <= tol))

    def report(self) -> bool:
        checks = [r for r in self.rows if r[3] is not None]
        w = max(len(r[0]) for r in checks)
        print(f"\n{'CLAIM'.ljust(w)}  {'OBSERVED':>12}  {'PRINTED':>12}  RESULT")
        for desc, obs, exp, ok in self.rows:
            if ok is None:
                print(f"\n{desc}\n" + "-" * (w + 42))
            else:
                print(f"{desc.ljust(w)}  {obs:>12}  {exp:>12}  {'PASS' if ok else 'FAIL'}")
        n_fail = sum(1 for r in checks if not r[3])
        print("-" * (w + 42))
        print(f"{len(checks) - n_fail}/{len(checks)} checks passed" + ("" if n_fail == 0 else f"  ({n_fail} FAILED)"))
        return n_fail == 0


# ----------------------------------------------------------------------------- revision
def check_target(c: Checker) -> None:
    c.section("Mortality target (Supplementary Methods S1, Table S-M1)")
    t = rev("rev00_target_construction.csv").set_index("quantity")["value"]
    c.equal("training samples", t["samples"], 4539)
    c.equal("mouse samples", t["mouse samples"], 3876)
    c.equal("rat samples", t["rat samples"], 663)
    c.equal("datasets", t["datasets (Source)"], 96)
    c.equal("samples with a reference group", t["samples with a reference group (Source x Tissue x Sex)"], 3774)
    c.equal("target reproduced exactly", t["target reproduced exactly (|diff| < 1e-4)"], 3466)
    c.printed("fraction reproduced among samples with a reference group",
              t["fraction reproduced exactly among samples with a reference group"], 0.92, 2)
    c.printed("Pearson r, target vs age difference (controls)", t["Pearson r, target vs age difference, controls"], 0.97, 2)
    c.printed("log10 hazard per month (mouse controls)",
              t["slope, log10 hazard per month of age difference, mouse controls"], 0.10, 2)
    c.printed("log10 hazard per 0.1 maximum lifespan", t["slope, log10 hazard per 0.1 of species maximum lifespan, controls"],
              0.50, 2)
    c.printed("target range (max |value|)", t["target range (max)"], 3.63, 2)


def check_cohorts(c: Checker) -> None:
    c.section("Cohorts and cell types (Section 2.4, Supplementary Tables S1-S2)")
    d = rev("rev01_cohort_donors.csv").set_index(["cohort", "group"])
    a = "GSE136831 (Adams)"
    c.equal("GSE136831 IPF donors", d.loc[(a, "IPF"), "donors"], 32)
    c.equal("GSE136831 control donors", d.loc[(a, "Control"), "donors"], 28)
    c.equal("GSE136831 COPD donors", d.loc[(a, "COPD"), "donors"], 18)
    c.equal("IPF age", d.loc[(a, "IPF"), "age"], "65.4 ± 5.4 (54–78)")
    c.equal("control age", d.loc[(a, "Control"), "age"], "45.6 ± 17.7 (20–80)")
    c.equal("COPD age", d.loc[(a, "COPD"), "age"].split(" (")[0], "62.4 ± 4.9")
    c.equal("IPF men", d.loc[(a, "IPF"), "male"], 26)
    c.equal("control men", d.loc[(a, "Control"), "male"], 16)
    c.equal("GSE135893 IPF donors", d.loc[("GSE135893 (Habermann)", "IPF"), "donors"], 12)
    c.equal("GSE135893 control donors", d.loc[("GSE135893 (Habermann)", "Control"), "donors"], 10)
    ct = rev("rev01_celltype_table.csv")
    c.equal("cell types with >= 50-cell pseudobulks", ct.in_coverage_set.sum(), 13)
    c.equal("cell types in the statistical grid", ct.in_statistics_set.sum(), 11)
    c.equal("scored-only cell types have 3 IPF donors", bool((ct.loc[~ct.in_statistics_set, "donors_IPF"] == 3).all()), True)
    det = rev("rev15_detected_features.csv").set_index("celltype").loc["detected in all cell types"]
    c.equal("module genes", det.module_genes_total, 1995)
    c.equal("module genes detected in all 11 cell types", det.module_genes_detected, 1025)
    c.equal("composite genes", det.composite_genes_total, 10487)
    c.equal("composite genes detected in all 11 cell types", det.composite_genes_detected, 3754)
    rle = rev("rev15_rle_support.csv")
    c.equal("RLE support genes, median per cell type", rle.genes_nonzero_all.median(), 5324)
    c.equal("RLE support genes, minimum", rle.genes_nonzero_all.min(), 4669)
    c.equal("RLE support genes, maximum", rle.genes_nonzero_all.max(), 9774)
    dep = rev("rev01_depth_by_disease.csv").set_index(["celltype", "variable"]).loc[("Macrophage", "lib_size")]
    c.printed("macrophage library size, IPF / control", dep.ratio, 2.4, 1)
    c.printed("macrophage library size FDR", dep.fdr, 0.017, 3)


def check_rodent(c: Checker) -> None:
    c.section("Released versus retrained clocks (Section 3.1, Supplementary Methods S2, Table S5)")
    s = rev("rev14_summary.csv").iloc[0]
    rep = rev("rev14_published_reproduction.csv")
    mods = rep[(rep.table == "5B module") & (rep.module != "All module genes")]
    med = mods[mods.outcome == "Mortality"].groupby("convention").r_insample.median()
    c.printed("released module clocks, median in-sample r (as released)", med["as_released"], 0.04, 2)
    c.printed("released module clocks, median r (gene-standardized)", med["gene_standardized"], 0.05, 2)
    c.printed("released module clocks, median r (within-module z)", med["within_module_sample_z"], 0.06, 2)
    c.printed("reported median module r", s.module_r_reported_median_mortality, 0.42, 2)
    c.printed("released module clocks, maximum r (any target or convention)", mods.r_insample.max(), 0.34, 2)
    c.printed("released composite, in-sample r", s.composite_r_insample_mortality, 0.998, 3)
    c.printed("released all-module clock, in-sample r", s.allmodule_r_insample_mortality, 0.22, 2)
    c.printed("all-module clock, reported r", s.allmodule_r_reported_mortality, 0.88, 2)
    c.printed("released vs retrained coefficients, median Spearman", s.coef_spearman_median, 0.04, 2)
    c.printed("released vs retrained IPF t statistics, Spearman", s.ipf_t_spearman, 0.006, 3)
    c.equal("IPF effects at FDR < 0.10 with released coefficients", s.ipf_hits_published_fdr10, 19)
    c.equal("overlap with the 30 retrained effects", s.ipf_hits_overlap, 5)
    c.equal("retrained effects with the same sign", s.ipf_hits_same_sign_published, 14)

    c.section("Rodent cross-validation (Section 3.1, Fig. 2a-b, Supplementary Table S4)")
    cv = rev("rev02_species_cv.csv")
    m, ch = cv[cv.outcome == "Mortality"].set_index("module"), cv[cv.outcome == "Chrono"]
    c.printed("mortality module clocks, minimum r", m.r_oof_all.min(), 0.15, 2)
    c.printed("mortality module clocks, maximum r", m.r_oof_all.max(), 0.54, 2)
    c.printed("mortality module clocks, median r", m.r_oof_all.median(), 0.40, 2)
    c.printed("chronological module clocks, minimum r", ch.r_oof_all.min(), 0.16, 2)
    c.printed("chronological module clocks, maximum r", ch.r_oof_all.max(), 0.60, 2)
    c.printed("chronological module clocks, median r", ch.r_oof_all.median(), 0.44, 2)
    c.equal("mortality modules exceeding all random draws", m.exceeds_all_random.sum(), 8)
    c.equal("chronological modules exceeding all random draws", ch.exceeds_all_random.sum(), 10)
    c.printed("median r, mouse test samples", m.r_oof_mouse.median(), 0.40, 2)
    c.printed("median r, rat test samples", m.r_oof_rat.median(), 0.49, 2)
    c.printed("median r, mouse-trained -> rat", m.r_mouse_to_rat.median(), 0.51, 2)
    c.printed("median r, rat-trained -> mouse", m.r_rat_to_mouse.median(), 0.33, 2)
    c.printed("median r, within-rat leave-one-dataset-out", m.r_within_rat_lodo.median(), 0.39, 2)
    c.printed("protein processing in ER, rat test r", m.loc["Protein processing in ER", "r_oof_rat"], -0.26, 2)
    c.printed("heme metabolism, rat test r", m.loc["Heme metabolism", "r_oof_rat"], -0.18, 2)
    c.printed("lipid metabolism r (Discussion: r ~ 0.3)", m.loc["Lipid met", "r_oof_all"], 0.3, 1)
    c.printed("heat stress response r (Discussion: r ~ 0.3)", m.loc["Heat stress response", "r_oof_all"], 0.3, 1)
    b = rev("rev02_benchmarks.csv")
    u = b[b.kind == "union_modules"].set_index("outcome").r_oof_all
    c.printed("union-of-modules clock r, mortality", u["Mortality"], 0.66, 2)
    c.printed("union-of-modules clock r, chronological", u["Chrono"], 0.75, 2)
    ag = rev("rev02b_allgene_benchmark.csv").set_index("outcome")
    c.equal("all-gene clock, genes", ag.loc["Mortality", "n_genes"], 18286)
    c.printed("all-gene clock r, mortality", ag.loc["Mortality", "r_oof_all"], 0.75, 2)
    c.printed("all-gene clock r, chronological", ag.loc["Chrono", "r_oof_all"], 0.82, 2)
    gs = rev("rev06_rodent_geneset_cv.csv")
    gm = gs[gs.outcome == "Mortality"]
    c.printed("clock minus signed gene-set r, median", gs.clock_minus_geneset.median(), 0.11, 2)
    c.equal("mortality clocks better than their gene-set score", (gm.clock_minus_geneset > 0).sum(), 20)
    c.printed("clock minus gene-set r, minimum", gm.clock_minus_geneset.min(), -0.06, 2)
    c.printed("clock minus gene-set r, maximum", gm.clock_minus_geneset.max(), 0.25, 2)


def check_klotho(c: Checker) -> None:
    c.section("Klotho knockout (Section 3.2, Fig. 2c, Supplementary Table S7)")
    k = rev("rev04_klotho_modules.csv").set_index(["tissue", "module"])
    for tissue, (d, lo, hi, p, pd_) in {"Kidney": (0.65, 0.27, 1.02, 0.009, 3),
                                         "Skeletal muscle": (0.61, 0.12, 1.10, 0.019, 3)}.items():
        r = k.loc[(tissue, "Composite")]
        c.printed(f"{tissue}: composite shift", r.delta, d, 2)
        c.printed(f"{tissue}: composite 95% CI, lower", r.ci_lo, lo, 2)
        c.printed(f"{tissue}: composite 95% CI, upper", r.ci_hi, hi, 2)
        c.printed(f"{tissue}: composite exact p", r.p_perm_exact, p, pd_)
    s = rev("rev04_klotho_summary.csv").set_index("tissue")
    c.equal("kidney modules shifted in the expected direction", s.loc["Kidney", "n_positive"], 9)
    c.equal("muscle modules shifted in the expected direction", s.loc["Skeletal muscle", "n_positive"], 13)
    c.printed("kidney count permutation p", s.loc["Kidney", "p_count_exact"], 0.78, 2)
    c.printed("muscle count permutation p", s.loc["Skeletal muscle", "p_count_exact"], 0.45, 2)
    c.equal("modules positive in both tissues", s.loc["Kidney", "n_positive_both_tissues"], 7)
    c.printed("muscle median |shift|", s.loc["Skeletal muscle", "obs_median_abs_delta"], 0.16, 2)
    c.printed("muscle wild-type split median |shift|", s.loc["Skeletal muscle", "wtwt_median_abs_delta"], 0.16, 2)
    c.printed("kidney median |shift|", s.loc["Kidney", "obs_median_abs_delta"], 0.17, 2)
    c.printed("kidney wild-type split median |shift|", s.loc["Kidney", "wtwt_median_abs_delta"], 0.09, 2)
    mods = k.drop(index="Composite", level="module")
    hits = mods[mods.fdr_perm < 0.05]
    c.equal("modules at FDR < 0.05 (kidney 5, muscle 1)", len(hits), 6)
    for (tissue, module), v in {("Kidney", "Adaptive immunity"): 0.50, ("Skeletal muscle", "Adaptive immunity"): 0.67,
                                ("Kidney", "Apoptosis"): -0.43, ("Kidney", "Nrf2 signaling"): -0.33,
                                ("Kidney", "Fatty acid met"): -0.36, ("Kidney", "Translation"): -0.38}.items():
        c.printed(f"{tissue}: {module}", hits.loc[(tissue, module), "delta"], v, 2)
    c.equal("modules exceeding a competitive null (p < 0.05)",
            int(((mods.p_vs_perm_coef < 0.05) | (mods.p_vs_random_genes < 0.05)).sum()), 0)


def check_human_age(c: Checker) -> None:
    c.section("Human age (Section 3.3, Fig. 2d, Supplementary Table S8)")
    h = rev("rev05_human_age_summary.csv").set_index(["dataset", "outcome"])
    hl = h.loc[("HLCA", "Mortality")]
    c.equal("HLCA tests", hl.n_tests, 253)
    c.printed("HLCA median rho", hl.module_median_rho, -0.17, 2)
    c.equal("HLCA FDR < 0.10 before adjustment", hl.module_pos_fdr10 + hl.module_neg_fdr10, 14)
    c.equal("HLCA FDR < 0.10 before adjustment, negative", hl.module_neg_fdr10, 13)
    c.equal("HLCA FDR < 0.10 after adjustment", hl.module_pos_fdr10_adj + hl.module_neg_fdr10_adj, 0)
    ad = h.loc[("Adams controls", "Mortality")]
    c.printed("GSE136831 controls median rho", ad.module_median_rho, 0.04, 2)
    c.equal("GSE136831 controls FDR < 0.10", ad.module_pos_fdr10 + ad.module_neg_fdr10, 0)
    c.printed("HLCA vs GSE136831 concordance rho", hl.hlca_vs_adams_rho_concordance, 0.09, 2)
    c.printed("HLCA vs GSE136831 concordance p", hl.hlca_vs_adams_concordance_p, 0.22, 2)
    g = rev("rev05b_gtex_validation.csv").set_index(["outcome", "module"])
    c.equal("GTEx donors", g.loc[("Mortality", "Composite"), "n"], 578)
    c.printed("GTEx composite mortality rho (adjusted)", g.loc[("Mortality", "Composite"), "adj_spearman"], 0.21, 2)
    c.sci("GTEx composite mortality p", g.loc[("Mortality", "Composite"), "adj_p"], 6, -7)
    c.printed("GTEx composite chronological rho (adjusted)", g.loc[("Chrono", "Composite"), "adj_spearman"], 0.19, 2)
    c.sci("GTEx composite chronological p", g.loc[("Chrono", "Composite"), "adj_p"], 3, -6)
    gm = g.loc["Mortality"].drop(index="Composite")
    c.equal("GTEx mortality modules scored", len(gm), 23)
    c.equal("GTEx modules FDR < 0.05 before adjustment", (gm.fdr < 0.05).sum(), 14)
    c.equal("... of which positive", ((gm.fdr < 0.05) & (gm.spearman > 0)).sum(), 6)
    c.equal("GTEx modules FDR < 0.05 after adjustment", (gm.adj_fdr < 0.05).sum(), 0)
    top = gm.adj_spearman.idxmax()
    c.equal("GTEx largest adjusted module", top, "Adaptive immunity")
    c.printed("GTEx largest adjusted rho", gm.loc[top, "adj_spearman"], 0.11, 2)
    c.printed("GTEx largest adjusted rho, FDR", gm.loc[top, "adj_fdr"], 0.12, 2)

    c.section("Module correlation (Section 3.3, Supplementary Fig. S1, Table S6)")
    me = rev("rev03_meff.csv")
    tr, ipf = me.iloc[0], me.iloc[1]
    c.printed("training median |rho|", tr.median_abs_r, 0.62, 2)
    c.printed("training effective number (Li-Ji)", tr.meff_liji, 10, 0)
    c.printed("training effective number (Nyholt)", tr.meff_nyholt, 15, 0)
    c.printed("IPF median |rho|", ipf.median_abs_r, 0.13, 2)
    c.printed("IPF effective number (Li-Ji)", ipf.meff_liji, 19, 0)
    c.printed("IPF effective number (Nyholt)", ipf.meff_nyholt, 22, 0)


def check_masking(c: Checker) -> None:
    c.section("Composite decomposition (Section 3.4, Fig. 3, Supplementary Table S11)")
    g = rev("rev08_global_stats.csv").set_index("module")
    comp = g.loc["Composite"]
    c.printed("composite IPF effect (log10 hazard units)", comp.beta, -0.023, 3)
    c.printed("composite 95% CI lower (cluster-robust)", comp.ci_lo, -0.178, 3)
    c.printed("composite 95% CI upper (cluster-robust)", comp.ci_hi, 0.131, 3)
    c.printed("composite p", comp.p_cr1_t, 0.76, 2)
    m = rev("rev07_masking_global.csv").set_index("component")
    s = m.loc["__summary__"]
    c.printed("composite 95% CI lower (donor bootstrap)", s.comp_lo, -0.158, 3)
    c.printed("composite 95% CI upper (donor bootstrap)", s.comp_hi, 0.144, 3)
    c.printed("total component signal S", s.S, 0.226, 3)
    c.printed("S 95% CI lower", s.S_lo, 0.196, 3)
    c.printed("S 95% CI upper", s.S_hi, 0.393, 3)
    c.printed("S permutation median", s.S_null_median, 0.090, 3)
    c.printed("S / permutation median", s.S / s.S_null_median, 2.5, 1)
    c.printed("S permutation p", s.p_S, 0.0005, 4)
    c.printed("net signal N", s.N, 0.023, 3)
    c.printed("N p", s.p_N, 0.64, 2)
    c.printed("masking index MI", s.MI, 0.90, 2)
    c.printed("MI 95% CI lower", s.MI_lo, 0.47, 2)
    c.printed("MI 95% CI upper", s.MI_hi, 0.99, 2)
    c.printed("MI permutation median", s.MI_null_median, 0.59, 2)
    c.printed("MI permutation p", s.p_MI, 0.13, 2)
    comps = m.drop(index="__summary__")
    sig = comps[comps.fdr < 0.10]
    c.equal("components at FDR < 0.10", len(sig), 10)
    for name, v in {"Adaptive immunity": 0.026, "Cell cycle": 0.024, "mRNA splicing": 0.013,
                    "Heat stress response": 0.006, "Respiration": -0.035, "Innate immunity": -0.027,
                    "Chromatin modification (1)": -0.019, "ECM organization": -0.014, "Amino acid met": -0.008,
                    "Muscle contraction": -0.006}.items():
        c.printed(f"component: {name}", sig.loc[name, "delta"], v, 3)
    nm = comps.loc["Non-module genes"]
    c.printed("non-module component", nm.delta, 0.026, 3)
    c.printed("non-module component CI lower", nm.lo, -0.105, 3)
    c.printed("non-module component CI upper", nm.hi, 0.162, 3)
    c.printed("non-module component FDR", nm.fdr, 0.73, 2)
    mc = rev("rev07_masking_celltype.csv")
    c.equal("cell types with S above its null at FDR < 0.10",
            sorted(mc.loc[mc.fdr_S < 0.10, "celltype"]), ["Macrophage", "Macrophage_Alveolar", "T", "cDC2", "ncMonocyte"])

    c.section("Global module effects (Section 3.4, Fig. 3d, Supplementary Table S9b)")
    for name, (b, lo, hi) in {"Interferon signaling": (-0.33, -0.46, -0.21), "VEGF signaling": (-0.27, -0.40, -0.14)}.items():
        r = g.loc[name]
        c.printed(f"{name}: effect", r.beta, b, 2)
        c.printed(f"{name}: 95% CI lower", r.ci_lo, lo, 2)
        c.printed(f"{name}: 95% CI upper", r.ci_hi, hi, 2)
        c.equal(f"{name}: FDR <= 0.007 under all three methods",
                bool(max(r.fdr_cr1_t, r.fdr_wild_cluster, r.fdr_mixed) <= 0.0075), True)
    for name, v in {"Lipid met": 0.31, "Heat stress response": 0.15, "Muscle contraction": -0.19}.items():
        r = g.loc[name]
        c.printed(f"{name}: effect", r.beta, v, 2)
        c.equal(f"{name}: passes CR1 and mixed model, not wild bootstrap",
                bool(r.fdr_cr1_t < 0.10 and r.fdr_mixed < 0.10 and r.fdr_wild_cluster >= 0.10), True)
        c.printed(f"{name}: wild bootstrap FDR", r.fdr_wild_cluster, 0.13, 2)
    t = g.loc["Translation"]
    c.printed("Translation: effect", t.beta, 0.27, 2)
    c.printed("Translation: CR1 FDR", t.fdr_cr1_t, 0.07, 2)
    c.printed("Translation: mixed-model FDR", t.fdr_mixed, 0.23, 2)
    sa = g.loc["Cholesterol met (salmon)"]
    c.printed("Cholesterol (salmon): effect", sa.beta, 0.23, 2)
    c.printed("Cholesterol (salmon): mixed-model FDR", sa.fdr_mixed, 0.06, 2)
    c.equal("Cholesterol (salmon): passes the mixed model only",
            bool(sa.fdr_mixed < 0.10 and sa.fdr_cr1_t >= 0.10 and sa.fdr_wild_cluster >= 0.10), True)


def check_celltype(c: Checker) -> None:
    c.section("Cell-type effects (Section 3.5, Fig. 4, Supplementary Table S9a)")
    s = rev("rev08_summary.csv").iloc[0]
    c.equal("tests", s.n_tests, 253)
    c.equal("FDR < 0.10, classical errors", s.hits_classical, 30)
    c.equal("FDR < 0.10, HC3", s.hits_hc3, 23)
    c.equal("FDR < 0.10, permutation", s.hits_perm, 29)
    c.equal("FDR < 0.10, all three", s.hits_all_three, 23)
    c.equal("IPF-specific (direct IPF-COPD contrast)", s.hits_ipf_specific_direct, 13)
    c.equal("IPF-specific under the previous rule", s.hits_old_rule_copd_ns, 15)
    c.equal("global modules at FDR < 0.10, CR1", s.global_fdr10_cr1, 6)
    c.equal("global modules at FDR < 0.10, wild bootstrap", s.global_fdr10_wild, 2)
    c.equal("global modules at FDR < 0.10, mixed model", s.global_fdr10_mixed, 6)
    st = rev("rev08_celltype_stats.csv").set_index(["celltype", "module"])
    hits = st[st.fdr_IPF < 0.10]
    n_ct = hits.index.get_level_values(0).value_counts()
    c.equal("effects in macrophages", n_ct["Macrophage"], 8)
    c.equal("effects in alveolar macrophages", n_ct["Macrophage_Alveolar"], 7)
    spec = hits[hits.IPF_specific_direct.astype(bool)]
    c.equal("IPF-specific effects in macrophages or alveolar macrophages",
            int(spec.index.get_level_values(0).str.startswith("Macrophage").sum()), 12)
    c.equal("remaining IPF-specific effect", [i for i in spec.index if not i[0].startswith("Macrophage")],
            [("NK", "ECM organization")])
    for (ct, mod), (b, f) in {("NK", "Cell cycle"): (0.86, 0.071), ("NK", "OxPhos"): (0.61, 0.055),
                              ("NK", "ECM organization"): (0.33, 0.084), ("NK", "Heme metabolism"): (0.31, 0.098),
                              ("T", "Interferon signaling"): (-0.95, 0.030)}.items():
        c.printed(f"{ct} {mod}: effect", hits.loc[(ct, mod), "beta_IPF"], b, 2)
        c.printed(f"{ct} {mod}: FDR", hits.loc[(ct, mod), "fdr_IPF"], f, 3)
    nk = hits.loc["NK"]
    c.equal("NK effects robust to HC3", list(nk.index[nk.fdr_IPF_hc3 < 0.10]), ["ECM organization"])
    c.printed("macrophage heat stress, COPD-control", st.loc[("Macrophage", "Heat stress response"), "beta_COPD"], 0.40, 2)
    c.printed("classical monocyte heat stress, COPD-control", st.loc[("cMonocyte", "Heat stress response"), "beta_COPD"],
              0.38, 2)

    c.section("Sensitivity analyses (Section 3.5, Fig. 5d, Supplementary Table S12)")
    v = rev("rev09_sensitivity_summary.csv").set_index("analysis")
    c.equal("age overlap: same sign", v.loc["age_overlap", "main_hits_same_sign"], 30)
    c.equal("age overlap: p < 0.05", v.loc["age_overlap", "main_hits_p05"], 30)
    c.printed("age overlap: coefficient ratio", v.loc["age_overlap", "main_hits_median_beta_ratio"], 1.02, 2)
    c.printed("age overlap: Spearman with main", v.loc["age_overlap", "beta_spearman_vs_main"], 0.96, 2)
    c.equal("technical covariates: same sign", v.loc["tech_covariates", "main_hits_same_sign"], 30)
    c.equal("technical covariates: p < 0.05", v.loc["tech_covariates", "main_hits_p05"], 28)
    c.printed("technical covariates: ratio", v.loc["tech_covariates", "main_hits_median_beta_ratio"], 0.91, 2)
    c.equal("CPM: same sign", v.loc["cpm", "main_hits_same_sign"], 30)
    c.printed("CPM: ratio", v.loc["cpm", "main_hits_median_beta_ratio"], 0.99, 2)
    c.equal("100 cells: testable effects", v.loc["min_cells_100", "main_hits_available"], 28)
    c.equal("100 cells: p < 0.05", v.loc["min_cells_100", "main_hits_p05"], 24)
    c.equal("30 cells: cell types (three added)", v.loc["min_cells_30", "n_celltypes"], 14)
    c.equal("30 cells: same sign", v.loc["min_cells_30", "main_hits_same_sign"], 29)
    c.equal("30 cells: p < 0.05", v.loc["min_cells_30", "main_hits_p05"], 26)
    c.printed("30 cells: ratio", v.loc["min_cells_30", "main_hits_median_beta_ratio"], 0.98, 2)
    c.printed("depth-matched: ratio", v.loc["depth_group_matched", "main_hits_median_beta_ratio"], 1.02, 2)
    c.printed("cell-number-matched: ratio", v.loc["cells_group_matched", "main_hits_median_beta_ratio"], 1.03, 2)
    c.equal("depth-matched: p < 0.05", v.loc["depth_group_matched", "main_hits_p05"], 23)
    c.equal("cell-number-matched: p < 0.05", v.loc["cells_group_matched", "main_hits_p05"], 22)
    eq = v.loc[["depth_equalized", "cells_equalized"]]
    c.equal("equalized: same sign 24-25", sorted(eq.main_hits_same_sign.astype(int)), [24, 25])
    c.printed("equalized: smallest ratio", eq.main_hits_median_beta_ratio.min(), 0.30, 2)
    c.printed("equalized: largest ratio", eq.main_hits_median_beta_ratio.max(), 0.39, 2)
    lc = v.loc[["depth_loss_control", "cells_loss_control"], "main_hits_median_beta_ratio"]
    c.printed("loss controls: smallest ratio", lc.min(), 0.49, 2)
    c.printed("loss controls: largest ratio", lc.max(), 0.51, 2)

    c.section("Baselines and competitive nulls (Section 3.6, Fig. 5a-c, Supplementary Table S10)")
    bl = rev("rev06_ipf_baselines_summary.csv").iloc[0]
    c.equal("unsigned gene-set effects at FDR < 0.10", bl.hits_geneset_fdr10, 34)
    c.equal("signed gene-set effects at FDR < 0.10", bl.hits_geneset_signed_fdr10, 27)
    c.printed("t correlation, clock vs unsigned gene set", bl.corr_t_clock_geneset, 0.06, 2)
    c.printed("t correlation, clock vs signed gene set", bl.corr_t_clock_geneset_signed, 0.62, 2)
    hc = rev("rev06_hits_competitive.csv").set_index(["celltype", "module"])
    c.equal("effects exceeding random genes (p < 0.05; 10,000 draws)", (hc.p_comp_random < 0.05).sum(), 16)
    c.equal("effects exceeding permuted coefficients (p < 0.05)", (hc.p_comp_perm < 0.05).sum(), 10)
    both = sorted(hc.index[hc.module_specific.astype(bool)])
    c.equal("effects passing both nulls at FDR < 0.10",
            both, sorted([("B", "Muscle contraction"), ("ncMonocyte", "Heat stress response"),
                          ("ncMonocyte", "Heme metabolism"), ("T_Cytotoxic", "VEGF signaling"),
                          ("NK", "ECM organization")]))
    c.printed("macrophage VEGF vs random genes", hc.loc[("Macrophage", "VEGF signaling"), "p_comp_random"], 0.012, 3)
    c.printed("macrophage interferon vs random genes", hc.loc[("Macrophage", "Interferon signaling"), "p_comp_random"],
              0.045, 3)
    mp = hc.loc[[("Macrophage", "VEGF signaling"), ("Macrophage", "Interferon signaling")], "p_comp_perm"]
    c.equal("macrophage VEGF and interferon vs permuted coefficients, p 0.05-0.10",
            bool(((mp >= 0.05) & (mp <= 0.10)).all()), True)


def check_replication(c: Checker) -> None:
    c.section("Replication (Section 3.7, Fig. 6a-c, Supplementary Table S13)")
    cc = rev("rev10_coarse_concordance.csv").set_index("coarse")
    for ct, (r, p) in {"NK": (0.57, 0.022), "T cell": (0.55, 0.022), "Monocyte": (0.18, 0.54),
                       "Macrophage": (-0.02, 0.94), "Ciliated": (-0.49, 0.059)}.items():
        c.printed(f"{ct}: cross-cohort r", cc.loc[ct, "r"], r, 2)
        c.printed(f"{ct}: two-sided p", cc.loc[ct, "p_two_sided"], p, len(str(p).split(".")[1]))
    c.printed("NK/T FDR across five cell types", cc.loc["NK", "fdr_two_sided"], 0.056, 3)
    c.printed("NK r after partialling the uniform shift", cc.loc["NK", "r_partial_uniform_shift"], 0.58, 2)
    c.printed("T r after partialling the uniform shift", cc.loc["T cell", "r_partial_uniform_shift"], 0.53, 2)
    c.printed("NK r, Vanderbilt donors only", cc.loc["NK", "r_vanderbilt_only"], 0.57, 2)
    c.printed("T r, Vanderbilt donors only", cc.loc["T cell", "r_vanderbilt_only"], 0.60, 2)
    c.printed("NK sign agreement (centred)", cc.loc["NK", "sign_conc_centred"], 0.65, 2)
    c.printed("T sign agreement (centred)", cc.loc["T cell", "sign_conc_centred"], 0.52, 2)
    se = rev("rev10_specific_effects_summary.csv").iloc[0]
    c.equal("testable FDR < 0.10 effects", se.testable, 27)
    c.equal("same sign in GSE135893", se.same_sign, 16)
    c.equal("permutation median", se.null_same_median, 14)
    c.printed("same-sign count p", se.p_count, 0.26, 2)
    c.equal("effects at p < 0.05 in GSE135893", se.n_p05_H, 6)

    ck = rev("rev15_checklist_realdata.csv").set_index("coarse")
    c.printed("ciliated split-half r", ck.loc["Ciliated", "r_splithalf"], -0.18, 2)
    c.printed("ciliated tau", ck.loc["Ciliated", "tau_splithalf"], 0.35, 2)
    c.printed("macrophage split-half r", ck.loc["Macrophage", "r_splithalf"], 0.81, 2)
    c.printed("monocyte split-half r", ck.loc["Monocyte", "r_splithalf"], 0.49, 2)
    c.printed("macrophage tau", ck.loc["Macrophage", "tau_splithalf"], 0.32, 2)
    c.printed("monocyte tau", ck.loc["Monocyte", "tau_splithalf"], 0.33, 2)
    c.printed("NK split-half r", ck.loc["NK", "r_splithalf"], 0.33, 2)
    c.printed("NK tau", ck.loc["NK", "tau_splithalf"], 0.40, 2)
    c.printed("macrophage sub-state mean r", ck.loc["Macrophage", "r_substate"], 0.15, 2)
    c.printed("monocyte sub-state mean r", ck.loc["Monocyte", "r_substate"], 0.08, 2)
    c.printed("macrophage sub-state Bonferroni p", ck.loc["Macrophage", "p_substate"], 1.00, 2)
    c.printed("monocyte sub-state Bonferroni p", ck.loc["Monocyte", "p_substate"], 0.25, 2)
    c.printed("macrophage random-split mean r", ck.loc["Macrophage", "r_random_split"], 0.14, 2)
    c.printed("monocyte random-split mean r", ck.loc["Monocyte", "r_random_split"], 0.20, 2)
    att = ck.attribution.to_dict()
    c.equal("attribution: ciliated", att["Ciliated"], "inconclusive (power)")
    c.equal("attribution: macrophage", att["Macrophage"], "between-cohort difference (unresolved)")
    c.equal("attribution: monocyte", att["Monocyte"], "between-cohort difference (unresolved)")
    c.equal("attribution: NK", att["NK"], "replicated")
    c.equal("attribution: T cell", att["T cell"], "replicated")

    co = rev("rev10_composition.csv").set_index(["cohort", "cell_type"])
    c.printed("macrophage share, GSE136831 IPF", co.loc[("GSE136831", "Macrophage"), "median_IPF"], 0.61, 2)
    c.printed("macrophage share, GSE136831 control", co.loc[("GSE136831", "Macrophage"), "median_Control"], 0.71, 2)
    c.printed("macrophage share, GSE135893 IPF", co.loc[("GSE135893", "Macrophage"), "median_IPF"], 0.31, 2)
    c.printed("macrophage share, GSE135893 control", co.loc[("GSE135893", "Macrophage"), "median_Control"], 0.28, 2)
    ca = rev("rev11_composition_adjusted.csv")
    fr = ca[ca.kind == "monoHi fraction"].set_index("celltype").loc["Macrophage"]
    c.printed("monocyte-marker-high fraction, IPF", fr.mono_frac_IPF, 0.66, 2)
    c.printed("monocyte-marker-high fraction, control", fr.mono_frac_Control, 0.42, 2)
    c.printed("monocyte-marker-high fraction, p", fr.p_mannwhitney, 0.01, 2)
    adj = ca[ca.kind == "adjusted effect"]
    c.equal("macrophage effects adjusted for the fraction", len(adj), 15)
    c.printed("adjusted coefficient ratio, minimum", adj.ratio.min(), 0.74, 2)
    c.printed("adjusted coefficient ratio, maximum", adj.ratio.max(), 1.36, 2)
    c.equal("all adjusted effects p < 0.05", bool((adj.p_adj_monofrac < 0.05).all()), True)
    ix = rev("rev11_interaction.csv")
    n_int = ix.groupby("celltype").apply(lambda x: (int((x.fdr_interaction < 0.10).sum()), len(x)), include_groups=False)
    c.equal("macrophage sub-state x disease interactions", n_int["Macrophage"], (8, 24))
    c.equal("monocyte sub-state x disease interactions", n_int["Monocyte"], (11, 24))

    bk = rev("rev10_bulk_concordance.csv").set_index("pair")
    for pair, (r, p) in {"GSE136831 whole-donor vs GSE135893 whole-donor": (0.47, 0.071),
                         "GSE136831 whole-donor vs GSE134692 bulk": (0.47, 0.059),
                         "GSE135893 whole-donor vs GSE134692 bulk": (0.51, 0.036)}.items():
        c.printed(f"{pair}: r", bk.loc[pair, "r"], r, 2)
        c.printed(f"{pair}: p", bk.loc[pair, "p_two_sided"], p, 3)
    dc = rev("rev10_bulk_decomposition.csv").set_index(["cohort", "component"]).r_with_observed
    c.printed("composition-only reconstruction, GSE136831", dc[("GSE136831", "composition_only")], 0.83, 2)
    c.printed("composition-only reconstruction, GSE135893", dc[("GSE135893", "composition_only")], 0.71, 2)
    c.printed("expression-only reconstruction, GSE136831", dc[("GSE136831", "expression_only")], 0.89, 2)
    c.printed("expression-only reconstruction, GSE135893", dc[("GSE135893", "expression_only")], 0.74, 2)


def check_checklist_simulation(c: Checker) -> None:
    c.section("Checklist simulation (Section 3.8, Fig. 6d, Supplementary Methods S3, Table S15a)")
    cf = rev("rev12_checklist_confusion.csv", keep_default_na=False).set_index("scenario")
    rep_, comp_, betw, inc = "replicated", "composition", "between-cohort difference (unresolved)", "inconclusive (power)"
    c.printed("shared effect -> replicated", cf.loc["shared", rep_], 1.00, 2)
    c.printed("sub-state composition -> composition", cf.loc["composition", comp_], 0.97, 2)
    c.printed("confound in one cohort -> between-cohort difference", cf.loc["confound", betw], 0.99, 2)
    c.printed("null effect -> inconclusive", cf.loc["null", inc], 0.94, 2)
    c.printed("weak shared effect -> replicated", cf.loc["power", rep_], 0.50, 2)
    c.printed("weak shared effect -> between-cohort difference", cf.loc["power", betw], 0.41, 2)
    c.printed("composition shift -> replicated", cf.loc["composition_shift", rep_], 0.44, 2)
    c.printed("composition shift -> between-cohort difference", cf.loc["composition_shift", betw], 0.44, 2)
    c.printed("simulation tau", cf.tau_splithalf.iloc[0], 0.251, 3)


def check_mil(c: Checker) -> None:
    c.section("Attention MIL, negative result (Section 3.8, Supplementary Results, Table S16)")
    a = rev("rev13_mil_auc_summary.csv").set_index("model")
    for model, (m, sd) in {"attention": (0.95, 0.02), "mean": (0.98, 0.01), "logistic": (0.97, 0.01)}.items():
        c.printed(f"{model}: seed-mean AUC", a.loc[model, "auc_seed_mean"], m, 2)
        c.printed(f"{model}: seed SD", a.loc[model, "auc_seed_sd"], sd, 2)
    c.equal("seeds", a.loc["attention", "n_seeds"], 10)
    c.printed("attention: AUC range, minimum", a.loc["attention", "auc_seed_min"], 0.91, 2)
    c.printed("attention: AUC range, maximum", a.loc["attention", "auc_seed_max"], 0.98, 2)
    for model, (p, lo) in {"attention": (0.96, 0.91), "mean": (0.98, 0.94), "logistic": (0.97, 0.93)}.items():
        c.printed(f"{model}: seed-averaged AUC", a.loc[model, "auc_of_mean_prediction"], p, 2)
        c.printed(f"{model}: 95% CI lower", a.loc[model, "ci_lo"], lo, 2)
        c.printed(f"{model}: 95% CI upper", a.loc[model, "ci_hi"], 1.00, 2)
    w = pd.read_csv(MIL / "stage3_within_label_raw.csv")
    wi = w[w.group == "IPF"]
    c.printed("separation, macrophages pooled (mean)", wi.sep_pooled.mean(), 0.52, 2)
    c.printed("separation, macrophages pooled (SD)", wi.sep_pooled.std(), 0.05, 2)
    c.printed("separation within Macrophage (mean)", wi.sep_within_Macrophage.mean(), -0.03, 2)
    c.printed("separation within Macrophage (SD)", wi.sep_within_Macrophage.std(), 0.03, 2)
    c.printed("separation within Macrophage_Alveolar (mean)", wi.sep_within_Macrophage_Alveolar.mean(), 0.10, 2)
    c.printed("separation within Macrophage_Alveolar (SD)", wi.sep_within_Macrophage_Alveolar.std(), 0.08, 2)
    c.printed("attention AUROC for the annotation label (mean)", wi.auroc_annotation.mean(), 0.88, 2)
    c.printed("attention AUROC for the annotation label (SD)", wi.auroc_annotation.std(), 0.03, 2)
    s1 = pd.read_csv(MIL / "stage1_results.csv")
    bag = s1.groupby(["regime", "p_s", "pool"]).bag_auc.mean()
    att = s1[s1.pool == "attention"].groupby(["regime", "p_s"]).att_auroc.mean()
    c.printed("calibration, sub-state p_S 0.2: attention bag AUC", bag[("substate", 0.2, "attention")], 0.82, 2)
    c.printed("calibration, sub-state p_S 0.2: mean-pool bag AUC", bag[("substate", 0.2, "mean")], 0.65, 2)
    c.printed("calibration, attention AUROC at p_S 0.05", att[("substate", 0.05)], 0.58, 2)
    c.printed("calibration, attention AUROC at p_S 0.5", att[("substate", 0.5)], 0.87, 2)
    c.printed("calibration, diffuse 0.2: attention bag AUC", bag[("diffuse", 0.2, "attention")], 0.73, 2)
    c.printed("calibration, diffuse 0.2: mean-pool bag AUC", bag[("diffuse", 0.2, "mean")], 0.65, 2)
    c.equal("calibration seeds", s1.seed.nunique(), 8)


def check_package(c: Checker) -> None:
    c.section("Reusable checklist (src/module_aging/checklist.py) on Supplementary Table S15b")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from module_aging import checklist

    ck = rev("rev15_checklist_realdata.csv").set_index("coarse")
    for co, r in ck.iterrows():
        sub = {} if pd.isna(r.r_substate) else dict(r_substate=r.r_substate, p_substate=r.p_substate,
                                                     r_random_split=r.r_random_split)
        c.equal(f"{co}: checklist() decision", checklist(r.r_cross, r.p_cross, r.r_splithalf, r.tau_splithalf,
                                                        **sub).decision, r.attribution)


def revision(c: Checker) -> None:
    for f in (check_target, check_cohorts, check_rodent, check_klotho, check_human_age, check_masking,
              check_celltype, check_replication, check_checklist_simulation, check_mil, check_package):
        f(c)


# ----------------------------------------------------------------------------- submission
CORRECTED = """
Values of the original submission that the revision corrects or replaces:
  * composite global effect: printed as "-2.9 months"; it is the composite effect in log10 hazard
    units (-0.023) multiplied by 122.5 and mislabelled as months (revised: -0.023, 95% CI -0.178 to 0.131)
  * cross-cohort r: NK 0.47, T 0.51, monocyte 0.08, macrophage -0.11, ciliated -0.20
    (revised analysis with fine-cell-type fixed effects and two-sided tests: 0.57, 0.55, 0.18, -0.02, -0.49)
  * split-half reliability: macrophage 0.77, monocyte 0.75 (revised: 0.81 and 0.49 with a permutation threshold)
  * bulk "recovery" of replication: withdrawn (composition and expression cannot be separated at bulk level)
  * attention MIL: moved to the Supplementary Information as a negative result
"""


def submission(c: Checker) -> None:
    def load(name: str) -> pd.DataFrame:
        return pd.read_csv(RESULTS / name)

    c.section("Original submission (v0.1.0); see the notes below the table")
    cov = load("module_coverage.csv")
    c.equal("retrained module clocks", len(cov), 23)

    g = load("module_global_robust.csv").set_index("module")
    comp = g.loc["Composite"]
    c.near("composite global effect as printed (-2.9; see note)", comp["beta"], -2.9, tol=0.1)
    c.equal("composite global effect non-significant (p > 0.05)", bool(comp["p"] > 0.05), True)
    c.near("Lipid met module +0.31", g.loc["Lipid met", "beta"], 0.31, tol=0.02)
    c.near("Translation module +0.27", g.loc["Translation", "beta"], 0.27, tol=0.02)
    c.near("Interferon signaling module -0.33", g.loc["Interferon signaling", "beta"], -0.33, tol=0.02)
    c.near("VEGF signaling module -0.27", g.loc["VEGF signaling", "beta"], -0.27, tol=0.02)

    st = load("ipf_module_stats_adjusted.csv")
    c.equal("cell-type x module effects tested", int(st.fdr_IPF.notna().sum()), 253)
    c.equal("effects passing FDR < 0.10", int((st.fdr_IPF < 0.10).sum()), 30)
    c.equal("IPF-specific under the rule COPD p > 0.05", int(((st.fdr_IPF < 0.10) & (st.p_COPD > 0.05)).sum()), 17)

    rep = load("replication_coarse_compare.csv")
    r_by_ct = {ct: float(np.corrcoef(gg.beta_adams, gg.beta_hab)[0, 1]) for ct, gg in rep.groupby("coarse")}
    c.near("cross-cohort r: NK 0.47 (see note)", r_by_ct["NK"], 0.47, tol=0.02)
    c.near("cross-cohort r: T cell 0.51", r_by_ct["Tcell"], 0.51, tol=0.02)
    c.near("cross-cohort r: monocyte 0.08", r_by_ct["Monocyte"], 0.08, tol=0.02)
    c.near("cross-cohort r: macrophage -0.11", r_by_ct["Macrophage"], -0.11, tol=0.02)
    c.near("cross-cohort r: ciliated -0.20", r_by_ct["Ciliated"], -0.20, tol=0.02)

    sc = load("mil/stageC_calibrate.csv").set_index("celltype")
    c.near("split-half macrophage 0.77 (see note)", sc.loc["Macrophage", "splithalf_med"], 0.77, tol=0.02)
    c.near("split-half monocyte 0.75", sc.loc["Monocyte", "splithalf_med"], 0.75, tol=0.02)

    bulk = load("mil/stageD_bulk.csv")
    c.equal("bulk r all in [0.46, 0.51] (withdrawn)", bool(((bulk.r >= 0.46) & (bulk.r <= 0.51)).all()), True)
    c.equal("bulk one-sided permutation p < 0.05 (withdrawn)", bool((bulk.p_perm < 0.05).all()), True)

    syn = load("mil/stage1_results.csv")
    sub = syn[syn.regime == "substate"]
    p02 = sub[np.isclose(sub.p_s, 0.20)]
    c.equal("attention > mean pooling at p_S = 0.2 (bag AUC)",
            bool(p02[p02.pool == "attention"].bag_auc.median() > p02[p02.pool == "mean"].bag_auc.median()), True)
    auroc = sub[sub.pool == "attention"].groupby("p_s").att_auroc.median()
    c.near("attention localization AUROC, maximum ~0.89", float(auroc.max()), 0.89, tol=0.06)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--submission", action="store_true",
                    help="check the numbers of the original submission (v0.1.0) instead of the revision")
    args = ap.parse_args()
    c = Checker()
    (submission if args.submission else revision)(c)
    ok = c.report()
    if args.submission:
        print(CORRECTED)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
