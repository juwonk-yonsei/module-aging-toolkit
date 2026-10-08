# Reproducing the manuscript

Two levels of reproduction:

- **Fast (no large data):** check every reported number against the committed result
  tables, and regenerate the figures and supplementary tables from them.
- **Full:** re-run the analyses from the public raw data.

## A. Fast verification (seconds)

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r verify/requirements.txt
python verify/verify_claims.py               # revised manuscript: 354 checks
python verify/verify_claims.py --submission  # original submission (v0.1.0): 21 checks
```

Each check reads a value from `results/revision/` (or `results/` and `results/mil/`) and
compares it with the value printed in the manuscript at the printed precision, grouped
by manuscript section. The exit code is non-zero on any mismatch. The `--submission`
mode ends with a list of the values that the revision corrects or withdraws.

## B. Revised manuscript (`scripts/revision/`)

Prerequisites: `bash setup.sh` (tAge, Zenodo rodent data, retrained clocks), the human
cohorts in `data/` (see `data/README.md`) and the v0.1.0 intermediates in `results/`
(pseudobulks and MIL instances, produced by the pipeline in section C). The scripts read
their inputs from `MAT_ROOT` and write to `MAT_OUT`:

```bash
export MAT_ROOT="$PWD"                         # repository root (data/, results/)
export MAT_OUT="$PWD/results/revision"
export TAGE_DIR="$PWD/third_party/tAge"        # tAge checkout made by setup.sh
cd scripts/revision
for s in rev_00_target rev_01_cohort rev_02_species_cv rev_02b_allgene_benchmark \
         rev_03_module_correlation rev_04_klotho rev_05_human_age rev_05b_gtex rev_06_baselines \
         rev_07_masking rev_08_stats rev_09_confounding rev_10_replication rev_11_substate \
         rev_12_checklist_simulation rev_13_mil rev_14_original_clocks rev_15_checklist_realdata; do
  python $s.py
done
python make_revision_figures.py       # figures_revision/Fig1-6, FigS1-S6 (PNG 600 dpi + PDF)
python make_supplementary_tables.py   # supplementary_revision/Supplementary_Tables.xlsx
```

Run the scripts in this order: later scripts read the outputs of earlier ones (rev_06 and
rev_14 read rev_02; rev_10 and rev_11 read rev_08; rev_12 reads rev_03; rev_15 reads
rev_10 and rev_11). `rev_02` trains 2 × 23 module clocks on 4,539 rodent samples and is
the slowest step; `rev_13` uses PyTorch and runs on CPU or GPU.

| Script | Main outputs (`results/revision/`) | Manuscript |
|---|---|---|
| `rev_00_target.py` | `rev00_target_construction.csv` | Supplementary Methods S1 |
| `rev_01_cohort.py` | `rev01_cohort_donors.csv`, `rev01_celltype_table.csv`, `rev01_depth_by_disease.csv` | Section 2.4; Tables S1–S2 |
| `rev_02_species_cv.py`, `rev_02b_allgene_benchmark.py` | `rev02_species_cv.csv`, `rev02_benchmarks.csv`, `rev02b_allgene_benchmark.csv` | Section 3.1; Fig. 2a–b; Table S4 |
| `rev_03_module_correlation.py` | `rev03_corr_*.csv`, `rev03_meff.csv`, `rev03_gene_overlap.csv` | Section 3.3; Fig. S1; Table S6 |
| `rev_04_klotho.py` | `rev04_klotho_modules.csv`, `rev04_klotho_summary.csv` | Section 3.2; Fig. 2c; Table S7 |
| `rev_05_human_age.py`, `rev_05b_gtex.py` | `rev05_hlca_age.csv`, `rev05_adams_control_age.csv`, `rev05b_gtex_validation.csv` | Section 3.3; Figs. 2d, S6; Table S8 |
| `rev_06_baselines.py` | `rev06_ipf_baselines*.csv`, `rev06_hits_competitive.csv`, `rev06_rodent_geneset_cv.csv` | Sections 2.7, 3.6; Fig. 5a–c; Table S10 |
| `rev_07_masking.py` | `rev07_masking_global.csv`, `rev07_masking_celltype.csv`, `rev07_standardized_effects.csv` | Sections 2.8, 3.4; Fig. 3a–c; Table S11 |
| `rev_08_stats.py` | `rev08_celltype_stats.csv`, `rev08_global_stats.csv`, `rev08_summary.csv` | Sections 2.6, 3.4–3.5; Figs. 3d, 4; Table S9 |
| `rev_09_confounding.py` | `rev09_sensitivity_summary.csv`, `rev09_hits_sensitivity.csv`, `rev09_sensitivity_long.csv` | Sections 2.9, 3.5; Fig. 5d; Table S12 |
| `rev_10_replication.py` | `rev10_coarse_concordance.csv`, `rev10_specific_effects*.csv`, `rev10_bulk_*.csv`, `rev10_composition.csv` | Section 3.7; Fig. 6a–c; Fig. S4; Table S13 |
| `rev_11_substate.py` | `rev11_interaction.csv`, `rev11_composition_adjusted.csv`, `rev11_substate_concordance.csv` | Section 3.7; Fig. 6e; Table S14 |
| `rev_12_checklist_simulation.py` | `rev12_checklist_confusion.csv`, `rev12_checklist_raw.csv` | Section 3.8; Fig. 6d; Methods S3; Table S15a |
| `rev_13_mil.py` (with `scripts/stage3_within_label.py`) | `rev13_mil_auc*.csv`; `results/mil/stage3_within_label*.csv` | Methods S4; Fig. S5; Table S16 |
| `rev_14_original_clocks.py` | `rev14_published_reproduction.csv`, `rev14_coef_agreement.csv`, `rev14_ipf_agreement.csv`, `rev14_summary.csv` | Section 3.1; Methods S2; Fig. S3; Table S5 |
| `rev_15_checklist_realdata.py` | `rev15_checklist_realdata.csv`, `rev15_detected_*.csv`, `rev15_rle_support.csv` | Sections 2.5, 3.7; Table 2; Fig. S2; Table S15b–e |

Seeds are fixed (`rev_common.SEED`); results can differ in the last digits between
platforms, and `rev_13` on a GPU is not bit-reproducible. Files named `_cache_*.pkl` in
`MAT_OUT` are local caches and can be deleted.

`make_revision_figures.py` needs only the committed tables (no tAge checkout, no raw
data), so the figures can be regenerated from a fresh clone with `requirements.txt`
installed. Rendered text can shift by a pixel with different fonts.
`make_supplementary_tables.py` also reads the retrained clocks in `data/module_clocks/`
for Tables S3a–b and therefore needs `setup.sh` first.

## C. Original submission (v0.1.0, `scripts/`)

Prerequisites as above. Set `MAT_ROOT` if the repository is not the working directory.

| Step | Script | Produces | Submission figure |
|---|---|---|---|
| 1 | `train_module_clocks.py` | `data/module_clocks/*.pkl`, `module_clock_performance.csv` | Fig 1b |
| 2 | `validate_modules_klotho.py`, `validate_example.py` | Klotho-KO module shifts | Fig 1c |
| 3 | `gtex_validate.py`, `hlca_validate.py` | cross-tissue human validation | — |
| 4 | `ipf_pipeline.py`, `stats_strengthen.py` | `ipf_module_stats_adjusted.csv`, `module_global_robust.csv` | Fig 2, Fig 3 |
| 5 | `copd_analysis.py`, `disease_module_pipeline.py` | IPF-vs-COPD specificity | Fig 2c |
| 6 | `replication_hab.py`, `diagnose_replication.py` | `replication_coarse_compare.csv` | Fig 4a |
| 7 | `stageA_harmonize.py` → `stageB_perm.py` / `stageB_substate.py` → `stageC_calibrate.py` → `stageD_bulk.py` | attribution ladder (`mil/stageC_calibrate.csv`, `mil/stageD_bulk.csv`) | Fig 4b, 4c |
| 8 | `mil_extract.py`, `mil_extract_hab.py` → `stage1_synthetic.py` | `mil/stage1_results.csv` (semi-synthetic calibration) | Methods S2 |
| 9 | `stage2_replication.py`, `stage3_interpret.py`, `stage3_multi.py` | `mil/stage3_multi.csv` (attention localization) | Fig 5 |
| — | `make_figures_final.py` | submission figures, needs steps 1–9 and the clocks | Fig 1–6 |
| — | `make_supp_coverage.py` | `figures/FigS1.png`, Supplementary Table S1 | Fig S1 |
| — | `make_graphical_abstract.py` | `figures/graphical_abstract.png` | GA |

The committed `figures/*.png` are the exact images of the original submission.
`gtex_validate.py` indexes module clocks by their label, so the two
cholesterol-metabolism modules overwrite each other and 22 module clocks are scored;
`scripts/revision/rev_05b_gtex.py` scores all 23.

## Notes

- MGB-licensed components are non-commercial and academic only
  (`THIRD_PARTY_NOTICES.md`).
- The retrained module clocks are regenerated locally by `train_module_clocks.py`; their
  coefficients are listed in Supplementary Table S3a.
