"""Assemble the revision supplementary tables (one XLSX workbook, one sheet per table) from the
CSV outputs of rev_01 ... rev_15, plus the module gene lists with retrained weights.

Output: <toolkit>/supplementary_revision/Supplementary_Tables.xlsx
"""
import importlib.metadata as im
import platform

import joblib
import numpy as np
import pandas as pd

import rev_common as rc

O = rc.OUT
DST = rc.TOOLKIT / "supplementary_revision"
DST.mkdir(exist_ok=True)
META = rc.TAGE_DIR / "inst" / "extdata" / "metadata"


def module_weights():
    memb = pd.read_csv(rc.DATA / "supp" / "module_membership_rodent.csv", dtype={"entrez": str})
    sym = dict(zip(memb.entrez, memb.symbol))
    orth = pd.read_csv(META / "Table_of_orthologs.csv").dropna(subset=["Entrez.Mouse", "Entrez.Human"])
    hsym = pd.read_csv(META / "Gene_table_human.csv").dropna(subset=["Entrez"])
    hsym = dict(zip(hsym.Entrez.astype("int64").astype(str), hsym["Gene.Symbol"]))
    orth["m"] = orth["Entrez.Mouse"].astype("int64").astype(str)
    orth["h"] = orth["Entrez.Human"].astype("int64").astype(str)
    h_of = orth.groupby("m").h.apply(lambda s: ";".join(sorted(set(s))))
    pub = pd.read_csv(rc.DATA / "supp" / "module_coef_rodent.csv", dtype={"entrez": str})
    pub = pub[pub.module != "All module genes"]
    pub["outcome"] = pub.outcome.map({"Chronological": "Chrono", "Mortality": "Mortality"})
    pub = pub.set_index(["outcome", "module", "entrez"]).coef
    cv = pd.read_csv(O / "rev02_species_cv.csv")
    cv = cv[cv.kind == "module"].set_index(["outcome", "module"]).r_oof_all
    genes, models = [], []
    for outcome in ["Mortality", "Chrono"]:
        for name, d in rc.load_module_clocks(outcome).items():
            med, mean, scale, coef, b0 = rc.linear_parts(d["pipeline"])
            est = list(d["pipeline"].named_steps.values())[-1]
            g = [str(x) for x in d["genes"]]
            hu = [h_of.get(x, "") for x in g]
            genes.append(pd.DataFrame({
                "outcome": outcome, "module_colour": d["module"], "module": name, "annotation": d["annotation"],
                "mouse_entrez": g, "mouse_symbol": [sym.get(x, "") for x in g],
                "human_entrez": hu, "human_symbol": [";".join(hsym.get(h, "") for h in s.split(";")) if s else "" for s in hu],
                "coef_standardized": coef, "coef_raw_scale": coef / scale, "train_median": med, "train_mean": mean,
                "train_sd": scale,
                "coef_published_SuppTable5B": [pub.get((outcome, d["module"], x), np.nan) for x in g]}))
            models.append(dict(outcome=outcome, module_colour=d["module"], module=name, annotation=d["annotation"],
                               n_genes=len(g), n_nonzero=int((coef != 0).sum()), alpha=float(est.alpha_),
                               l1_ratio=float(est.l1_ratio_), intercept=b0,
                               r_oof_dataset_grouped=cv.get((outcome, name), np.nan)))
    G = pd.concat(genes, ignore_index=True)
    M = pd.DataFrame(models)
    return G, M


def software():
    pk = ["numpy", "pandas", "scipy", "statsmodels", "scikit-learn", "joblib", "torch", "matplotlib",
          "scanpy", "anndata", "openpyxl"]
    rows = [dict(item="python", value=platform.python_version()), dict(item="platform", value=platform.platform())]
    for p in pk:
        try:
            rows.append(dict(item=p, value=im.version(p)))
        except im.PackageNotFoundError:
            rows.append(dict(item=p, value="not installed in this environment"))
    rows += [dict(item="global random seed (rev_common.SEED)", value=str(rc.SEED)),
             dict(item="module-clock training", value="SimpleImputer(median) -> StandardScaler -> ElasticNetCV("
                  "l1_ratio {0.5, 0.9, 1.0}, alphas logspace(-3, 1, 20), cv=3, max_iter=4000)"),
             dict(item="rodent cross-validation", value="GroupKFold(5) with groups = dataset (Source)"),
             dict(item="permutation nulls", value="Freedman-Lane B = 5,000 for the 253 cell-type tests; B = 2,000 for "
                  "masking metrics and cross-cohort label permutation; competitive nulls 1,000 draws in the 253-test "
                  "screen and in Klotho, 10,000 draws for the 30 FDR < 0.10 effects"),
             dict(item="bootstrap", value="wild cluster bootstrap (Rademacher) B = 9,999; stratified donor bootstrap "
                  "B = 2,000 (masking, MIL AUC)"),
             dict(item="checklist", value="split-half: 200 disease-stratified halves, null 200 label permutations x 20 "
                  "splits; simulation: 300 cohort pairs per scenario, 200 permutations per test"),
             dict(item="MIL", value="10 seeds x StratifiedKFold(5) on donors; seeds 0-9")]
    return pd.DataFrame(rows)


def extra_cohorts():
    return pd.DataFrame([
        dict(cohort="GSE134692 (bulk RNA-seq)", group="IPF / control", donors="46 / 26", age="see GEO",
             male="see GEO", female="see GEO", tissue="lung homogenate"),
        dict(cohort="HLCA core (healthy)", group="healthy", donors="42", age="20–81 y", male="", female="",
             tissue="lung, multiple sampling protocols"),
        dict(cohort="GTEx v8 lung", group="healthy (post-mortem)", donors="578", age="20–79 y (decade brackets)",
             male="", female="", tissue="bulk lung")])


TABLES = [
    ("S1", "Cohorts", "Donors per group, age (mean ± SD, range), sex and tissue source for every cohort.",
     lambda: pd.concat([pd.read_csv(O / "rev01_cohort_donors.csv").astype(str), extra_cohorts()], ignore_index=True)),
    ("S2a", "CellTypes_GSE136831", "Donors and cells per cell type and group in GSE136831; in_statistics_set marks "
     "the 11 cell types with >= 4 donors per group (the 253-test grid = 11 x 23).", "rev01_celltype_table.csv"),
    ("S2b", "CellTypes_GSE135893", "Donors per cell type in GSE135893 after label harmonization.",
     "rev01_hab_celltype_table.csv"),
    ("S2c", "Depth_by_disease", "Cells, library size and detected-gene fraction per pseudobulk, IPF vs control "
     "(Mann-Whitney, BH-FDR).", "rev01_depth_by_disease.csv"),
    ("S3a", "Module_gene_weights", "All genes of the 23 retrained module clocks for both outcomes: mouse Entrez/symbol, "
     "human ortholog, standardized and raw-scale coefficients, training median/mean/SD used for imputation and scaling, "
     "and the published coefficient (Tyshkovskiy et al., Supp. Table 5B). Derived from data distributed under the "
     "MGB Open Access License 1.0: non-commercial academic use only.", None),
    ("S3b", "Module_models", "Model-level parameters of the retrained module clocks (selected alpha, l1 ratio, "
     "intercept, non-zero coefficients) and dataset-grouped out-of-fold r. Same licence terms as S3a.", None),
    ("S4a", "Rodent_CV_species", "Dataset-grouped cross-validation: pooled, species-stratified, mouse->rat, rat->mouse, "
     "rat leave-one-dataset-out, size-matched random-gene benchmarks; union-of-module-genes and all-gene benchmarks.",
     "rev02_species_cv.csv"),
    ("S4b", "Rodent_CV_random_draws", "Each of the 10 size-matched random-gene draws per module and outcome.",
     "rev02_benchmarks.csv"),
    ("S4c", "Rodent_CV_geneset", "Clock vs equal-weight signed gene-set score of the same genes on the same folds.",
     "rev06_rodent_geneset_cv.csv"),
    ("S5a", "Published_reproduction", "Published coefficients applied to the released training matrix under three "
     "scaling conventions vs the reported accuracy.", "rev14_published_reproduction.csv"),
    ("S5b", "Published_vs_retrained_coef", "Agreement between published and retrained coefficients per module.",
     "rev14_coef_agreement.csv"),
    ("S5c", "Published_vs_retrained_IPF", "IPF effects obtained with published vs retrained module clocks.",
     "rev14_ipf_agreement.csv"),
    ("S6a", "Module_redundancy", "Median absolute Spearman correlation between module scores and effective number "
     "of modules (Li-Ji, Nyholt) in the rodent training data and in IPF pseudobulks.", "rev03_meff.csv"),
    ("S6b", "Module_gene_overlap", "Pairwise gene overlap between modules.", "rev03_gene_overlap.csv"),
    ("S6c", "Module_corr_training", "Spearman correlation matrix of module predictions, rodent training data.",
     "rev03_corr_training.csv"),
    ("S6d", "Module_corr_IPF", "Spearman correlation matrix of module scores, IPF pseudobulks (centred within cell "
     "type).", "rev03_corr_ipf.csv"),
    ("S7a", "Klotho_modules", "Klotho-KO vs WT per module and tissue: difference, 95% CI, exact permutation p, FDR, "
     "competitive nulls.", "rev04_klotho_modules.csv"),
    ("S7b", "Klotho_summary", "Klotho-KO: number of modules in the a-priori (positive) direction and exact count test.",
     "rev04_klotho_summary.csv"),
    ("S7c", "Klotho_WT_vs_WT", "WT-vs-WT split null for the Klotho comparison.", "rev04_klotho_wtwt.csv"),
    ("S8a", "Age_HLCA", "Spearman rho with donor age per cell type x module in HLCA healthy lung (raw and adjusted for "
     "technical covariates).", "rev05_hlca_age.csv"),
    ("S8b", "Age_GSE136831_controls", "Spearman rho with donor age per cell type x module in GSE136831 controls.",
     "rev05_adams_control_age.csv"),
    ("S8c", "Age_GTEx", "GTEx v8 lung bulk: composite and 23 module clocks vs age, raw and adjusted for ischaemic time, "
     "RIN, Hardy scale and sex.", "rev05b_gtex_validation.csv"),
    ("S8d", "Age_summary", "Summary of the human age analyses.", "rev05_human_age_summary.csv"),
    ("S9a", "IPF_celltype_stats", "All cell type x module tests: IPF-control (classical, HC3, Freedman-Lane permutation), "
     "COPD-control and direct IPF-COPD contrasts with BH-FDR.", "rev08_celltype_stats.csv"),
    ("S9b", "IPF_global_stats", "Global module effects: cluster-robust (CR1), wild cluster bootstrap and mixed model.",
     "rev08_global_stats.csv"),
    ("S9c", "IPF_stats_summary", "Hit counts under each inference method.", "rev08_summary.csv"),
    ("S10a", "Baselines_all_tests", "Unsigned/signed gene-set scores and competitive nulls for all 253 tests.",
     "rev06_ipf_baselines.csv"),
    ("S10b", "Baselines_30_hits", "Competitive nulls (10,000 size-matched random gene sets; 10,000 coefficient "
     "permutations) for the 30 FDR < 0.10 effects.", "rev06_hits_competitive.csv"),
    ("S10c", "Baselines_summary", "Summary of the baseline comparisons.", "rev06_ipf_baselines_summary.csv"),
    ("S11a", "Masking_global", "Composite decomposition into module and non-module components with bootstrap CIs; "
     "masking metrics S, N, MI with label-permutation nulls (row __summary__).", "rev07_masking_global.csv"),
    ("S11b", "Masking_celltype", "Masking metrics per cell type.", "rev07_masking_celltype.csv"),
    ("S11c", "Masking_components", "Per-cell-type composite components.", "rev07_masking_components.csv"),
    ("S11d", "Standardized_effects", "Standardized effects (beta / residual SD) for all tests.",
     "rev07_standardized_effects.csv"),
    ("S12a", "Sensitivity_hits", "The 30 FDR < 0.10 effects under every sensitivity analysis.",
     "rev09_hits_sensitivity.csv"),
    ("S12b", "Sensitivity_summary", "Summary of each sensitivity analysis.", "rev09_sensitivity_summary.csv"),
    ("S13a", "Replication_vectors", "Cross-cohort correlation of 23-module effect vectors (raw, centred, uniform-shift "
     "partialled; two-sided label-permutation p; Vanderbilt-only).", "rev10_coarse_concordance.csv"),
    ("S13b", "Replication_specific", "Replication of the specific FDR < 0.10 effects in GSE135893.",
     "rev10_specific_effects.csv"),
    ("S13c", "Replication_specific_sum", "Summary of specific-effect replication.", "rev10_specific_effects_summary.csv"),
    ("S13d", "Bulk_concordance", "Whole-donor and bulk concordance between cohorts (two-sided permutation p).",
     "rev10_bulk_concordance.csv"),
    ("S13e", "Bulk_direction", "Direction of the global module effects in whole-donor and bulk data.",
     "rev10_bulk_direction.csv"),
    ("S13f", "Bulk_decomposition", "Composition-only and expression-only reconstructions of whole-donor effects.",
     "rev10_bulk_decomposition.csv"),
    ("S13g", "Cell_composition", "Cell-type proportions per donor and cohort.", "rev10_composition.csv"),
    ("S13h", "Marker_overlap", "Overlap of macrophage identity/sub-state markers with module gene sets.",
     "rev10_marker_overlap.csv"),
    ("S14a", "Substate_interaction", "Macrophage/monocyte sub-state x disease interaction per module.",
     "rev11_interaction.csv"),
    ("S14b", "Substate_adjusted", "Main effects after adjustment for the monocyte-marker-high fraction.",
     "rev11_composition_adjusted.csv"),
    ("S14c", "Substate_concordance", "Cross-cohort concordance of matched sub-states vs random splits.",
     "rev11_substate_concordance.csv"),
    ("S15a", "Checklist_simulation", "Operating characteristics of the replication checklist (300 simulated cohort "
     "pairs per scenario).", "rev12_checklist_confusion.csv"),
    ("S15b", "Checklist_real_data", "Checklist applied to the five coarse cell types (split-half r, null 95th "
     "percentile tau, sub-state and random-split comparisons, attribution).", "rev15_checklist_realdata.csv"),
    ("S15c", "Detected_features", "Detected (non-imputed) composite and module features per cell type and in all "
     "cell types.", "rev15_detected_features.csv"),
    ("S15d", "Detected_by_module", "Detected fraction per module and cell type.", "rev15_detected_by_module.csv"),
    ("S15e", "RLE_support", "Genes retained by the detection filter after ortholog mapping, and the subset with "
     "non-zero counts in every pseudobulk from which the RLE size factors are computed, per cell type.",
     "rev15_rle_support.csv"),
    ("S16a", "MIL_AUC", "Out-of-fold donor AUC per seed for attention MIL, mean-pool MIL and logistic regression.",
     "rev13_mil_auc.csv"),
    ("S16b", "MIL_AUC_summary", "Seed mean/SD/range and stratified donor-bootstrap 95% CI.",
     "rev13_mil_auc_summary.csv"),
    ("S16c", "MIL_within_label", "Attention separation along the monocyte-derived vs resident marker axis, pooled and "
     "within each macrophage annotation, and annotation-recovery AUROC, per seed and donor group.",
     lambda: pd.read_csv(rc.RES / "mil" / "stage3_within_label_raw.csv")),
    ("S16d", "MIL_calibration", "Semi-synthetic calibration: bag AUC and attention AUROC for signal-bearing "
     "instances per regime, signal fraction, pooling and seed.", lambda: pd.read_csv(rc.RES / "mil" / "stage1_results.csv")),
    ("S17", "Software", "Software versions, seeds and fixed settings.", None),
]


def main():
    G, M = module_weights()
    special = {"S3a": G, "S3b": M, "S17": software()}
    ab = O / "rev02b_allgene_benchmark.csv"
    path = DST / "Supplementary_Tables.xlsx"
    readme = []
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        pd.DataFrame().to_excel(xw, sheet_name="README")
        for tid, sheet, desc, src in TABLES:
            if tid in special:
                df = special[tid]
            elif callable(src):
                df = src()
            else:
                df = pd.read_csv(O / src, keep_default_na=False, na_values=[""])
                if tid == "S4a" and ab.exists():
                    df = pd.concat([df, pd.read_csv(ab)], ignore_index=True)
                if src.startswith("rev03_corr"):
                    df = pd.read_csv(O / src, index_col=0).reset_index().rename(columns={"index": "module"})
            name = f"{tid}_{sheet}"[:31]
            df.to_excel(xw, sheet_name=name, index=False)
            readme.append(dict(table=f"Table {tid}", sheet=name, rows=len(df), description=desc))
        R = pd.DataFrame(readme)
        head = pd.DataFrame([dict(table="", sheet="", rows="", description=(
            "Supplementary Tables for 'Opening the composite clock: a single-cell framework with built-in controls for "
            "mortality-associated module signatures, demonstrated in pulmonary fibrosis' (revised manuscript). "
            "Effects are on the clock's native scale (difference in log10 relative mortality hazard for mortality clocks) "
            "unless stated; FDR = Benjamini-Hochberg."))])
        pd.concat([head, R], ignore_index=True).to_excel(xw, sheet_name="README", index=False)
    print(f"wrote {path}  ({len(TABLES)} tables; module genes {len(G)} rows)")
    print(M.groupby("outcome")[["n_genes", "n_nonzero"]].sum())


if __name__ == "__main__":
    main()
