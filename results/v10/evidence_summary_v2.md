# Two public objects for the same module clocks

This note records what the public files support. It does not choose a verdict.

## 1. Supplementary Table 5B

Source workbook sheet `(B) Module rodent clocks` in MOESM7, md5 `990a4d6cc6a9aaa892606f47fa2e1731`. Expression matrix `Expression_data_relative_rodents_Scaled.csv`, Zenodo 18763485, md5 `c283869cf69e35afa200c6a6fc281519` (18,286 genes by 4,539 samples).

From `results/audit_module_mapping/stage2/b1_reproduction.csv` (`source=supp_table`, `convention=as_released_na0`, mortality): the median in-sample Pearson r of the 23 modules is 0.050. The all-module clock has r = 0.222. The same file’s reported column, and `results/v3/rodent/native_reproduction.csv`, place the published module median near 0.42 (median of `r_reported` for mortality modules = 0.419) and the all-module figures at 0.881 (mortality) and 0.904 (chronological).

Sign concordance of non-zero Table 5B coefficients with the gene–target correlation is at chance for the mortality all-module clock (51.1%, permutation p = 0.208) and for 22 of 23 modules. Mean pairwise correlation inside a module is r̄ = 0.021 (median across modules), and no module exceeds z = 3 against size-matched random sets. These two results are the stage-2 compilation in `보고서/02_evidence_summary_en.md`. Gene-family concentration is in `results/v6/corruption/coherence/table5b_mortality_coherence.csv` and `results/v10/coherence_null_band_v10.csv`: all six families have permutation p > 0.45 (seed 0, 1,000 shuffles).

## 2. Official serialized clocks

Zenodo record 22166800, archive `Rodent module clocks.zip`, md5 `be9234547aacfbb2e83f0f1c5e44b3bd`. Loader: `src/module_aging/official_clocks.py` (`load_clocks`, `apply_clock`). Gene coefficients are not copied into this workbook.

From `results/v3/rodent/native_reproduction.csv`: mortality module median in-sample r = 0.443; all-module r = 0.892 (mortality) and 0.913 (chronological), beside reported 0.881 and 0.904. Family concentration in `results/v6/corruption/coherence/official_coherence.csv`: ribosome, OXPHOS, interferon, cell cycle, and collagen sit at the floor p = 0.001 (1/1,001); heat shock top share 0.444, p = 0.014. Within-module coexpression for these official sets is the excess reported with the stage-2 official scores (median pairwise r = 0.385; all 23 modules z > 3), compiled in `보고서/02_evidence_summary_en.md`.

## 3. Gene-count difference

`results/audit_module_mapping/stage2/pkl_vs_5b.csv`. Median Jaccard of the official gene set and the Table 5B non-zero set for the same colour is 0.011. Counts match for 16 of 23 colours. The differences (official minus Table 5B) are brown +64 (255 vs 191), darkgreen +8 (36 vs 28), grey60 +6 (52 vs 46), midnightblue +4 (52 vs 48), darkorange +2 (17 vs 15), and cyan +1 (20 vs 19).

## 4. How to reproduce

| Step | Script | Command | Output |
| --- | --- | --- | --- |
| Spreadsheet scores | `scripts/audit/stage2_alignment.py` | `python scripts/audit/stage2_alignment.py` | `results/audit_module_mapping/stage2/b1_reproduction.csv` |
| Official scores | `scripts/v3/extract_official_clocks.py` | `python scripts/v3/extract_official_clocks.py` | `results/v3/rodent/native_reproduction.csv` |
| Family concentration | `지시서/module_coherence_audit.py` | `python 지시서/module_coherence_audit.py` | `results/v6/corruption/coherence/*_coherence.csv` |

The coherence null bands used here were recomputed with seed 0 and 1,000 permutations and written to `results/v10/coherence_null_band_v10.csv`. They match the stored top shares.

## 5. Questions

1. Are the gene labels in Supplementary Table 5B aligned with the coefficient columns that trained the clocks in Zenodo 22166800?
2. Which object should be cited for Table 4B: the spreadsheet, or the pickled pipelines on Zenodo 22166800?
3. Is a correction of Table 5B planned, and will the membership table behind Table 6B be released?
