# Reproducing the manuscript

Two levels of reproduction:

- **Fast (no big data):** verify every headline number and regenerate the
  panels that read only the committed source CSVs.
- **Full:** re-run the whole pipeline from the public raw data.

## A. Fast verification (seconds)

```bash
python verify/verify_claims.py
```

Recomputes each manuscript statistic from `results/*.csv` and prints a
PASS/FAIL table (20 checks: composite β, opposing modules, 253/30 FDR counts,
per-cell-type cross-cohort r, split-half ceilings, bulk recovery, semi-synthetic
attention calibration). Exit code is non-zero on any mismatch.

## B. Full pipeline

Prerequisites: `bash setup.sh` (tAge + Zenodo + retrained clocks) and the human
cohorts in `data/` (see `data/README.md`). Set `MAT_ROOT` if the repo is not the
working root. Run in order:

| Step | Script | Produces | Figure |
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
| — | `make_figures_final.py` | `figures/Fig1..5.png`, needs steps 1–9 + clocks | Fig 1–5 |
| — | `make_supp_coverage.py` | `figures/FigS1.png`, Supplementary Table S1 | Fig S1 |
| — | `make_graphical_abstract.py` | `figures/graphical_abstract.png` | GA |

Example:

```bash
export MAT_ROOT="$PWD"
python scripts/ipf_pipeline.py
python scripts/stats_strengthen.py
python scripts/stageA_harmonize.py && python scripts/stageC_calibrate.py && python scripts/stageD_bulk.py
python scripts/make_figures_final.py
```

## Notes

- Figure regeneration reads the module clocks (`data/module_clocks/`) for Fig 1;
  run `setup.sh` first. Fig 2–5 read only committed `results/*.csv`.
- The committed `figures/*.png` are the exact submission images.
- MGB-licensed components are non-commercial/academic only
  (`THIRD_PARTY_NOTICES.md`).
