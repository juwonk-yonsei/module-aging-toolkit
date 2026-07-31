# Data

Raw data are **not** stored in this repository (they are large and publicly
hosted). `setup.sh` retrieves the module-clock training data from Zenodo and
regenerates the clocks; the human cohorts below must be downloaded from their
original providers under those providers' terms.

Expected layout after download (all git-ignored):

```
data/
├── rodent/     Expression_data_relative_rodents_Scaled.csv   (Zenodo 18763485)
│               Data_annotation_relative_rodents.xlsx
├── supp/       module_membership_rodent.csv
├── clock_models/   EN_*.pkl        base tAge clocks   (Zenodo 18763485)
├── module_clocks/  module_*.pkl    retrained locally by setup.sh
├── ipf/        GSE136831 (Adams)  primary human scRNA-seq pseudobulks
├── ipf2/       GSE135893 (Habermann) replication scRNA-seq pseudobulks
├── ipf_bulk/   GSE134692          independent bulk RNA-seq
├── gtex/       GTEx v8 lung        cross-tissue human validation
└── hlca/       Human Lung Cell Atlas (via cellxgene-census)
```

## Accessions

| Purpose | Accession | Source |
|---|---|---|
| Module-clock training data + base clocks | Zenodo **18763485** | https://zenodo.org/records/18763485 |
| Primary human scRNA-seq (Adams) | GEO **GSE136831** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831 |
| Replication scRNA-seq (Habermann) | GEO **GSE135893** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893 |
| Independent bulk RNA-seq | GEO **GSE134692** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134692 |
| Cross-tissue human validation | **GTEx v8** lung | https://gtexportal.org (bulk); dbGaP for raw |

## Notes

- The **module-clock training data and base clocks (Zenodo 18763485)** are
  distributed under the **MGB Open Access License 1.0** (non-commercial /
  academic). See `../THIRD_PARTY_NOTICES.md`.
- Human single-cell counts are aggregated to **donor × cell-type pseudobulks**
  (≥ 50 cells retained). HLCA/GTEx pull helpers:
  `scripts/hlca_pull_pseudobulk.py`, `scripts/gtex_extract_lung.py`.
- Only the **source CSVs** needed to regenerate figures and to run
  `verify/verify_claims.py` are committed (under `../results/`).
