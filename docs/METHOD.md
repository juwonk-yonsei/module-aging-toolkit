# Using the method on your own data

This toolkit lets you (1) score single-cell pseudobulks with **module- and
cell-type-resolved** aging clocks and (2) run the **attribution ladder**
(manuscript Table 1) to diagnose *why* a module signal does or does not
replicate across cohorts.

## 1. The attribution ladder (`module_aging`)

The reusable, clock-agnostic core lives in `src/module_aging/`. It works on
plain NumPy/pandas inputs — you do **not** need the transcriptomic clocks to use
it, only per-module effect vectors (and, for the power stage, per-donor module
scores).

```python
from module_aging import (
    effect_concordance, concordance_permutation_p,
    split_half_ceiling, classify,
)

# per-module effect vectors (one value per module) from two cohorts
conc = effect_concordance(eff_cohortA, eff_cohortB)     # Pearson/Spearman/sign
p    = concordance_permutation_p(eff_cohortA, eff_cohortB)

# within-cohort split-half "ceiling" from per-donor module scores + labels
sh = split_half_ceiling(donor_scores_A, labels_A)       # {median, lo, hi}

verdict = classify(conc.pearson_r, p, splithalf_median=sh["median"],
                   r_bulk=..., p_bulk=...)               # Table-1 decision
print(verdict.stage, "->", verdict.conclusion)
```

The ladder stages (applied in order, stop at the first that resolves it):

| Stage | Question | Primitive |
|---|---|---|
| 0 observe | Does the effect vector replicate? | `effect_concordance` + `concordance_permutation_p` |
| 1 technical | Gene-coverage / imputation artifact? | rescore on common gene support, then `effect_concordance` |
| 2 power | Is *n* too small? | `split_half_ceiling` |
| 3 composition | Sub-state mixing? | matched marker-defined sub-states, then `effect_concordance` |
| 4 resolution | Recovers at bulk? | whole-donor pseudobulk / bulk cohort, then `effect_concordance` |

See `examples/reuse_quickstart.py` for a runnable end-to-end demo on synthetic
data.

## 2. Scoring your own pseudobulks with the module clocks

This path needs the tAge dependency and the module clocks (run `setup.sh`
first). The clocks are applied to **donor × cell-type pseudobulk counts**:

```
counts : genes × samples   (raw pseudobulk counts)
meta   : samples × covariates   (must include a control-group column)
```

1. Preprocess with the ported tAge pipeline
   (`third_party/tage_prep.py`: RLE → log → per-sample scale → YuGene →
   ortholog map to mouse Entrez → control-median subtraction).
2. Apply each `data/module_clocks/module_*.pkl` to its module genes to get a
   per-sample module aging score.

`scripts/apply_module_clocks.py` and `scripts/ipf_pipeline.py` show the full
sequence; adapt the input paths to your dataset (or set `MAT_ROOT`).

> **License:** the module clocks / tAge components are MGB Open Access License
> 1.0 (non-commercial, academic). See `THIRD_PARTY_NOTICES.md`.
