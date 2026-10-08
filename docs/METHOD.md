# Using the method on your own data

This toolkit lets you (1) score single-cell pseudobulks with **module- and
cell-type-resolved** clocks and (2) run the **diagnostic checklist**
(revised manuscript, Table 2) to decide how far a cross-cohort non-replication
can be interpreted.

## 1. The diagnostic checklist (`module_aging.checklist`)

The reusable, clock-agnostic core lives in `src/module_aging/`. It needs only NumPy
and works on sample × module score matrices with a case indicator per sample; you do
**not** need the transcriptomic clocks to use it.

```python
from module_aging import cross_cohort_test, split_half_reliability, checklist

# step 1: Pearson r between the two cohorts' effect vectors, two-sided p from
# permuting disease labels within each cohort
r, p = cross_cohort_test(scores_A, case_A, scores_B, case_B,
                         strata_a=fine_celltype_A, strata_b=fine_celltype_B,
                         donors_a=donor_A, donors_b=donor_B)

# step 2: split-half reliability of the primary cohort and its permutation threshold
r_sh, tau = split_half_reliability(scores_A, case_A, strata=fine_celltype_A, donors=donor_A)

# steps 3-4: optional sub-state comparison (mean r over matched sub-states,
# Bonferroni-adjusted p, and the mean r of random splits of the same sizes)
res = checklist(r, p, r_sh, tau, r_substate=..., p_substate=..., r_random_split=...)
print(res.step, res.decision)
```

`strata` adds fixed effects for fine cell types within a coarse cell type, and `donors`
makes permutations and split halves operate on donors rather than on pseudobulks; both
are optional. The steps are applied in order and stop at the first decision:

| Step | Question | Decision |
|---|---|---|
| 1 | Do the effect vectors agree between cohorts? | r > 0 and p < 0.05: replicated |
| 2 | Is the primary cohort's effect vector reliable at this sample size? | split-half r < τ: inconclusive (power) |
| 3 | Does matching cell sub-states restore agreement? | p < 0.05 and r above the random-split r: composition |
| 4 | — | between-cohort difference (unresolved) |

In simulation (Section 3.8, Supplementary Table S15a) the checklist called a shared
effect replicated in 100% of cohort pairs, but it called a weak shared effect a
between-cohort difference in 41% and a pure composition shift replicated in 44%. Treat
its decisions as diagnostic statements about what the data can support, not as proof of
a cause. The analyses in the manuscript use the scripts in `scripts/revision/`
(`rev_10`, `rev_11`, `rev_12`, `rev_15`); `examples/reuse_quickstart.py` is a runnable
demonstration on synthetic data.

The attribution ladder of the original submission (`module_aging.classify` and related
functions) is kept unchanged so that v0.1.0 results can be reproduced.

## 2. Scoring your own pseudobulks with the module clocks

This path needs the tAge dependency and the module clocks (run `setup.sh`
first). The clocks are applied to **donor × cell-type pseudobulk counts**:

```
counts : genes × samples   (raw pseudobulk counts)
meta   : samples × covariates   (must include a control-group column)
```

1. Preprocess with the ported tAge pipeline
   (`third_party/tage_prep.py`: detection filter → ortholog map to mouse Entrez →
   RLE → log10 → per-sample scaling → control-median subtraction).
2. Apply each `data/module_clocks/module_*.pkl` to its module genes to get a
   per-sample module score in Δlog10 hazard units (mortality clocks).

`scripts/apply_module_clocks.py`, `scripts/ipf_pipeline.py` and
`scripts/revision/rev_common.py` show the full sequence; adapt the input paths to your
dataset (or set `MAT_ROOT`). Contrasts between disease groups within a cell type are
unaffected by the imputation of undetected genes; raw scores should not be compared
between cell types.

> **License:** the module clocks / tAge components are MGB Open Access License
> 1.0 (non-commercial, academic). See `THIRD_PARTY_NOTICES.md`.
