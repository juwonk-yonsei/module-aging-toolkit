# module-aging-toolkit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23229671.svg)](https://doi.org/10.5281/zenodo.23229671)

**A single-cell framework with built-in controls for module- and cell-type-resolved
transcriptomic clocks: competitive nulls, composite decomposition, disease contrasts,
sensitivity analyses and a diagnostic checklist for cross-cohort non-replication.**

Companion code for the manuscript *"Opening the composite clock: a single-cell framework
with built-in controls for mortality-associated module signatures, demonstrated in
pulmonary fibrosis"* (Kang, Jeon & Choi; Methods Paper, under revision at *GeroScience*).

> Zenodo concept DOI [10.5281/zenodo.23229671](https://doi.org/10.5281/zenodo.23229671) ·
> v0.2.0 [10.5281/zenodo.23229672](https://doi.org/10.5281/zenodo.23229672) ·
> v0.3.0 has no version DOI yet · Paper DOI: _added on acceptance_

| Version | Corresponds to |
|---|---|
| **v0.3.0** (this version) | revised manuscript using the serialized module clocks: `figures/v17/`, `results/v15/Supplementary_Tables_v15.xlsx`, `verify/claims_r2.csv` |
| v0.2.0 | withdrawn: that release retrained the module clocks |
| v0.1.0 | original submission: `scripts/`, `results/*.csv`, `figures/` |

A composite clock can report no change while its components move in opposite
directions, and a cell-type-resolved clock signal can look informative while failing
null models or replication. The toolkit scores 23 rodent module clocks in donor ×
cell-type pseudobulks and attaches an explicit control to each claim. In the
demonstration in idiopathic pulmonary fibrosis, the scores are interpreted as
mortality-associated module signatures, not as aging rates, because they do not track
donor age within human cell types.

---

## Two ways to use this repo

### Reuse the checklist on your own data
```bash
git clone https://github.com/juwonk-yonsei/module-aging-toolkit.git && cd module-aging-toolkit
pip install -e .            # exposes the `module_aging` package
python examples/reuse_quickstart.py
```
The checklist (manuscript Table 2) works on plain sample × module score matrices and
needs no clocks. See **[docs/METHOD.md](docs/METHOD.md)**.

### Verify the reported numbers
```bash
pip install -r verify/requirements.txt     # numpy and pandas, pinned
python verify/verify_claims.py --r2        # v0.3.0 manuscript numbers that map to a supplementary-table cell
python verify/verify_claims.py             # v0.2.0 revision tables
python verify/verify_claims.py --submission   # original submission (v0.1.0)
```
Each check reads a value from the committed result tables and compares it with the
number printed in the manuscript at the printed precision. No raw data are needed and
the script runs in seconds; the exit code is non-zero on any mismatch.
`--r2` checks the claims in `verify/claims_r2.csv` against `results/v15/Supplementary_Tables_v15.xlsx`.

Full pipeline: **[docs/REPRODUCE.md](docs/REPRODUCE.md)**.

---

## Install

```bash
# option A: pip
pip install -r requirements.txt && pip install -e .
# option B: conda
conda env create -f environment.yml && conda activate module-aging
# Checksums of the author-licensed inputs already on disk (no download, no retraining):
MAT_EXPRESSION=/path/to/Expression_data_relative_rodents_Scaled.csv \
MAT_CLOCKS=/path/to/Rodent_module_clocks.zip \
MAT_TABLE5B=/path/to/MOESM7.xlsx \
bash setup.sh --verify-only
```
`requirements.txt` pins the versions that produced the results (Python 3.11). For a
CPU-only PyTorch, install `torch==2.5.1` from `https://download.pytorch.org/whl/cpu`
first. The v0.2.0 path that retrained module clocks is withdrawn. Human cohorts
(GEO, HLCA, GTEx) stay with their providers; see **[data/README.md](data/README.md)**.

## Repository layout

```
src/module_aging/        diagnostic checklist (checklist.py) and the v0.1.0
                         attribution ladder                               (MIT)
scripts/                 original pipeline and figure code (v0.1.0)       (MIT)
scripts/revision/        revision analyses rev_00 ... rev_15, figure and
                         supplementary-table scripts                      (MIT)
verify/                  manuscript-number checks + pinned requirements   (MIT)
examples/                runnable quickstart on synthetic data            (MIT)
results/                 v0.1.0 result tables (results/*.csv, results/mil/)
results/revision/        revision result tables (rev00_* ... rev15_*)
figures/                 submission figures (Fig1-6, FigS1, graphical abstract)
figures_revision/        revised figures (Fig1-6, FigS1-S6; PNG 600 dpi + PDF)
supplementary_revision/  Supplementary Tables S1-S17 (one workbook)
third_party/             tage_prep.py port + LICENSE.MGB                  (MGB licence)
data/                    (git-ignored) fetched by setup.sh; see data/README.md
docs/                    METHOD.md (reuse) and REPRODUCE.md (verify, rerun)
```

## Licensing (important)

- **Our code** (`src/`, `scripts/`, `verify/`, `examples/`) is **MIT**.
- The **transcriptomic clocks, tAge API and rodent training data** are provided by
  Mass General Brigham under the **MGB Open Access License 1.0** (**non-commercial,
  academic use only**, with a share-back obligation). `third_party/tage_prep.py` (our
  port of the tAge preprocessing) is a derivative of those Materials and is distributed
  under the same licence.
- The serialized module clocks, their coefficients and the training matrix are not in
  this repository. `scripts/setup_v030.py --verify-only` checks the checksums of copies
  obtained from the authors' Zenodo records. The v0.2.0 retrained clocks are withdrawn.
- Because the pipeline depends on these components, **the pipeline as a whole may only
  be used for non-commercial academic purposes.** Details:
  **[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)**.

## Citing

See **[CITATION.cff](CITATION.cff)**. Please also cite Tyshkovskiy et al. 2026
(Nature; tAge clocks, DOI 10.1038/s41586-026-10542-3).
