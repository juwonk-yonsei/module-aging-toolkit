# module-aging-toolkit

**A single-cell framework with built-in controls for module- and cell-type-resolved
transcriptomic clocks: competitive nulls, composite decomposition, disease contrasts,
sensitivity analyses and a diagnostic checklist for cross-cohort non-replication.**

Companion code for the manuscript *"Opening the composite clock: a single-cell framework
with built-in controls for mortality-associated module signatures, demonstrated in
pulmonary fibrosis"* (Kang, Jeon & Choi; Methods Paper, under revision at *GeroScience*).

> Zenodo archive of this version: _DOI added on release_ · Paper DOI: _added on acceptance_

| Version | Corresponds to |
|---|---|
| **v0.2.0** (this version) | revised manuscript: `scripts/revision/`, `results/revision/`, `figures_revision/`, `supplementary_revision/` |
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
python verify/verify_claims.py             # revised manuscript
python verify/verify_claims.py --submission   # original submission (v0.1.0)
```
Each check reads a value from the committed result tables and compares it with the
number printed in the manuscript at the printed precision. No raw data are needed and
the script runs in seconds; the exit code is non-zero on any mismatch.

```
354/354 checks passed   # target construction, rodent CV, Klotho, human age, composite
                        # decomposition, inference, sensitivity, nulls, replication,
                        # checklist simulation, MIL, and the packaged checklist
```

The `--submission` mode reproduces the 21 checks of v0.1.0 and lists the values that the
revision corrects, including the composite effect printed as "−2.9 months", which was
the effect in log10 hazard units (−0.023) multiplied by 122.5.

Full pipeline: **[docs/REPRODUCE.md](docs/REPRODUCE.md)**.

---

## Install

```bash
# option A: pip
pip install -r requirements.txt && pip install -e .
# option B: conda
conda env create -f environment.yml && conda activate module-aging
# then fetch third-party dependencies and regenerate the module clocks:
bash setup.sh
```
`requirements.txt` pins the versions that produced the results (Python 3.11). For a
CPU-only PyTorch, install `torch==2.5.1` from `https://download.pytorch.org/whl/cpu`
first. `setup.sh` clones the tAge dependency, downloads the rodent training data (Zenodo
18763485) and regenerates the retrained clocks locally. Human cohorts (GEO, HLCA, GTEx)
are downloaded separately; see **[data/README.md](data/README.md)**.

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
- The retrained module clocks are derived from those data. Their fitted model files are
  **not** distributed; `scripts/train_module_clocks.py` regenerates them locally. Their
  coefficients, training medians and scaling parameters are listed in Supplementary
  Tables S3a and S3b (sheets `S3a_Module_gene_weights` and `S3b_Module_models` of
  `supplementary_revision/Supplementary_Tables.xlsx`) under the MGB licence terms.
- Because the pipeline depends on these components, **the pipeline as a whole may only
  be used for non-commercial academic purposes.** Details:
  **[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)**.

## Citing

See **[CITATION.cff](CITATION.cff)**. Please also cite Tyshkovskiy et al. 2026
(Nature; tAge clocks, DOI 10.1038/s41586-026-10542-3).
