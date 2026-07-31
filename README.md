# module-aging-toolkit

**Module- and cell-type-resolved application and validation of transcriptomic
aging clocks in single-cell data — with a reusable "attribution ladder" for
diagnosing cross-cohort non-replication.**

Companion code for the manuscript *"A module- and cell-type-resolved framework
for applying and validating transcriptomic aging clocks in single-cell data: a
proof-of-concept in idiopathic pulmonary fibrosis"* (Kang & Choi, submitted to
*GeroScience*).

> Paper / DOI: _add on acceptance_ · Zenodo archive: _add on release_

A single composite aging score can hide **opposing** biological programs. This
toolkit decomposes aging into pathway-specific **module clocks** scored per
**cell type**, and provides a clock-agnostic **attribution ladder** (Table 1)
that localizes *why* a signal replicates or fails across cohorts (technical /
power / composition / resolution).

---

## Two ways to use this repo

### 🔁 REUSE — apply the method to your data
```bash
git clone https://github.com/juwonk-yonsei/module-aging-toolkit.git && cd module-aging-toolkit
pip install -e .            # exposes the `module_aging` package
python examples/reuse_quickstart.py
```
The attribution ladder runs on plain per-module effect vectors (no clocks
required). See **[docs/METHOD.md](docs/METHOD.md)**.

### ✅ REPRODUCE / VERIFY — check our results
```bash
python verify/verify_claims.py
```
Recomputes every headline number from the committed `results/*.csv` and prints a
PASS/FAIL table — a machine-checkable audit that the reported statistics are
genuine code outputs. Full pipeline: **[docs/REPRODUCE.md](docs/REPRODUCE.md)**.

```
20/20 checks passed   # composite β=-2.9 (n.s.), 253/30 FDR, cross-cohort r,
                      # split-half ceilings, bulk recovery, MIL calibration ...
```

---

## Install

```bash
# option A: pip
pip install -r requirements.txt && pip install -e .
# option B: conda
conda env create -f environment.yml && conda activate module-aging
# then fetch third-party deps + regenerate module clocks:
bash setup.sh
```
`setup.sh` clones the tAge dependency, downloads the module-clock training data
(Zenodo 18763485), and regenerates the retrained clocks locally. Human cohorts
(GEO/GTEx) are downloaded separately — see **[data/README.md](data/README.md)**.

## Repository layout

```
src/module_aging/     attribution ladder + concordance primitives   (MIT)
scripts/              full analysis + figure pipeline (01→figures)   (MIT)
verify/verify_claims.py   manuscript-number audit                    (MIT)
examples/             runnable quickstart on synthetic data          (MIT)
results/              source CSVs to regenerate figures & run verify
figures/              submission figures (Fig1–5, FigS1, graphical abstract)
third_party/          tage_prep.py port + LICENSE.MGB (MGB license)
data/                 (git-ignored) fetched by setup.sh; see data/README.md
docs/                 METHOD.md (reuse) · REPRODUCE.md (verify)
```

## Licensing (important)

- **Our code** (`src/`, `scripts/`, `verify/`, `examples/`) is **MIT**.
- The **transcriptomic clocks, tAge API, and rodent training data** are provided
  by Mass General Brigham under the **MGB Open Access License 1.0**
  (**non-commercial, academic use only**, with a share-back obligation).
  `third_party/tage_prep.py` (our port of tAge preprocessing) is a derivative of
  those Materials and is distributed under the same MGB license.
- Because the pipeline depends on those components, **the pipeline as a whole may
  only be used for non-commercial academic purposes.** Full details:
  **[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)**.

## Citing

See **[CITATION.cff](CITATION.cff)**. Please also cite Tyshkovskiy et al. 2026
(Nature; tAge clocks, DOI 10.1038/s41586-026-10542-3).
