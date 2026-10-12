# Verification of the reported numbers

This folder recomputes the numbers reported as results in the revised manuscript and its Supplementary
Information from the released result tables, and compares each with the value printed in the text.

## Files

| File | Content |
|---|---|
| `claims_r2.csv` | One row per number printed in the manuscript and the Supplementary Information, including numbers written as words ("six of the 11 cell types"). References, author information and keywords are excluded. |
| `claims_expr.py` | For each number, the expression that recomputes it from the result tables, and a note where the expression needs one. |
| `verify_claims_r2.py` | Evaluates the expressions and compares each value with the printed one. |
| `verify_report.csv` | The output: one row per number with the recomputed value and the result. |

Columns of `claims_r2.csv`: `id`; `document` (manuscript or SI); `section` (section, figure legend or
table); `loc` (paragraph or table-row label); `form` (digits or word); `printed` (as printed); `value`;
`precision` (decimal places printed); `kind`; `rule` (how the kind was assigned); `context` (the number in
its sentence, marked with ⟦ ⟧).

Each number has one of four kinds:

- **computed**: a result of an analysis. Every computed number has an expression and is recomputed.
- **design**: a setting or threshold fixed in advance (for example B = 2,000 permutations or FDR < 0.10).
  Where the setting is stored in a result table, it is checked too; otherwise it is listed.
- **structural**: section, figure, table, step and reference numbers, version labels and dates. Listed.
- **cited**: a value reported by another study. Checked where the result tables store it; otherwise listed.

## How values are compared

Expressions read only the result tables (`Supplementary_Tables.xlsx`, one sheet per table) through a few
helpers defined in `verify_claims_r2.py`: `v` (the value in the one row selected by the filters), `n`
(number of rows), `med`, `mean`, `mn`, `mx`, `corr` (Pearson or Spearman), `tci` (confidence limits from
a coefficient, its t statistic and residual degrees of freedom), `power_t` (power of a two-sided t test) and
`num` (a number inside a text cell). A recomputed value matches if it rounds to the printed value at the
printed number of decimal places (half up or half to even); values printed in scientific notation are
compared at the printed number of significant digits. An expression that returns True or False checks a
bound ("all p < 0.05"); an expression that returns several values requires each to print as the value.

## Run

```
python verify/verify_claims_r2.py verify/claims_r2.csv supplementary_revision/Supplementary_Tables.xlsx verify/verify_report.csv
```

The script prints the number of rows of each kind and result and a summary line for the computed numbers,
and lists any row that fails. The exit code is 0 only if no row fails and every computed number has an
expression. It needs pandas, numpy, openpyxl and scipy, which are in `requirements.txt`.

## Scope

The script checks that the text reports what the result tables contain. It does not check the result
tables themselves; they are written by the analysis scripts of this repository from the public data, and
re-running those scripts requires the data that `setup.sh` downloads. The claims table was built by
extracting every number from the final text of the revision; a number added to the text later would need a
new row.
