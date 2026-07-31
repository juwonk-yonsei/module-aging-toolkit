# Third-party components and licenses

This toolkit's **own** code (`src/module_aging/`, `scripts/` analysis & figure
code, `verify/`, `examples/`) is released under the MIT License (see `LICENSE`).

It depends on external components that carry their **own, more restrictive
licenses**. These are *not* redistributed inside this repository (except the
clearly marked port below); `setup.sh` fetches them into your environment. By
using this pipeline you agree to those upstream licenses.

---

## 1. tAge — transcriptomic clock API and module clocks

- **Source:** https://github.com/Gladyshev-Lab/tAge
- **Data / models:** Zenodo record **18763485**
  (https://zenodo.org/records/18763485) — processed rodent expression
  meta-dataset (18,286 genes × 4,539 samples, 96 datasets) and transcriptomic
  clock models (`.pkl`).
- **Publication:** Tyshkovskiy A, et al. *Universal transcriptomic hallmarks of
  mammalian ageing and mortality.* Nature (2026).
  https://doi.org/10.1038/s41586-026-10542-3
- **License:** **MGB OPEN ACCESS LICENSE 1.0** (Mass General Brigham). A copy is
  included at `third_party/LICENSE.MGB`. Key terms:
  - **Non-commercial, non-revenue-generating, academic use only** (§3).
  - You may modify and redistribute, but **only** non-commercially, and you must
    **deliver your modifications back** to the Licensor / their GitHub, and you
    **may not remove authorship/notices** (§4, §5(e)).
  - **No transfer** of rights/obligations to third parties (§5(f), §10).

Because our retrained module clocks and our preprocessing port are derived from
MGB-licensed data/code, they **cannot** be relicensed under MIT.

### What we do about it
- **tAge python API (`tage_predict`) and `inst/extdata`, `inst/python`:** not
  redistributed. `setup.sh` clones the upstream repo.
- **Module clock `.pkl` models (retrained by us on the Zenodo data via
  `scripts/train_module_clocks.py`):** not redistributed. `setup.sh` downloads
  the Zenodo data and re-runs training locally, so each user generates the
  models from the MGB-licensed source under the MGB license.
- **`third_party/tage_prep.py`:** our Python port of tAge's R preprocessing
  (`preprocessing.R`). It is a **derivative work of MGB Materials**, so it is
  distributed here under the **MGB Open Access License 1.0** (not MIT), with the
  original attribution retained in its header. Per §4 we will also contribute
  this port back to the upstream tAge project.

## 2. Public datasets (obtain from the original providers)

| Dataset | Accession | Provider terms |
|---|---|---|
| Human scRNA-seq (primary, Adams) | GEO **GSE136831** | NCBI GEO |
| Human scRNA-seq (replication, Habermann) | GEO **GSE135893** | NCBI GEO |
| Human bulk RNA-seq (IPF) | GEO **GSE134692** | NCBI GEO |
| Cross-tissue human validation | **GTEx v8** lung | GTEx / dbGaP terms |

See `data/README.md` for retrieval instructions. These are not redistributed.

## 3. Python libraries

NumPy, pandas, SciPy, scikit-learn, statsmodels, matplotlib, PyTorch, joblib,
Pillow, Scanpy, cellxgene-census — each under its own OSI-approved license.
