"""Python port of tAge preprocessing (from Gladyshev-Lab/tAge R package).

LICENSE: This file is a derivative work of tAge (Gladyshev Lab, Mass General
Brigham) and is therefore distributed under the MGB OPEN ACCESS LICENSE 1.0
(non-commercial/academic; see ../third_party/LICENSE.MGB and
../THIRD_PARTY_NOTICES.md), NOT the MIT license that covers the rest of this
toolkit. Original authorship/attribution to the tAge project is retained.

Reproduces tAge_preprocessing() so transcriptomic clocks (sklearn .pkl) can be
applied fully in Python (scanpy ecosystem), without R.

Pipeline (matches R preprocessing.R):
  filter_genes -> map_genes(->mouse Entrez) -> RLE -> log10(x+1)
  -> scale (per-SAMPLE z-score, i.e. per column across genes)
  -> YuGene (on scaled) -> align_to_gene_list -> control_subtraction
Outputs `scaled_diff` and `yugene_diff` matrices (samples x gene_list).

Notes on fidelity:
  * R `scale(expr)` standardizes COLUMNS (samples), since expr is genes x samples.
  * YuGene is applied to the scaled matrix, per column (sample), exactly as in R.
  * RLE = edgeR calcNormFactors(method="RLE") = median-of-ratios vs per-gene
    geometric mean, factors rescaled to geomean 1; effective lib = lib*factor;
    norm = counts/eff*1e7.
"""
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd

# Resolve tAge extdata. Override with TAGE_EXTDATA, else TAGE_DIR/inst/extdata,
# else the toolkit default third_party/tAge/inst/extdata (populated by setup.sh).
_REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
_TAGE = Path(os.environ.get("TAGE_DIR", str(_REPO / "third_party" / "tAge")))
_EXT = Path(os.environ.get("TAGE_EXTDATA", str(_TAGE / "inst" / "extdata")))
_META = _EXT / "metadata"
_GENE_LIST = _EXT / "Gene_list_all_4.6.txt"


def load_gene_list() -> list[str]:
    return [ln.strip() for ln in _GENE_LIST.read_text().splitlines() if ln.strip()]


# ---------- gene mapping -> mouse Entrez ----------
def _load_gene_mapping(species: str, gene_mapping_type: str) -> dict:
    gt = pd.read_csv(_META / f"Gene_table_{species}.csv")
    gt = gt.dropna(subset=["Entrez"]).drop_duplicates(subset=[gene_mapping_type])
    return dict(zip(gt[gene_mapping_type].astype(str), gt["Entrez"].astype("int64")))


def map_genes(counts: pd.DataFrame, species: str, gene_mapping_type: str) -> pd.DataFrame:
    """counts: genes x samples (index = Ensembl or Gene.Symbol). Returns mouse-Entrez x samples."""
    gmap = _load_gene_mapping(species, gene_mapping_type)
    mapped = pd.Series(pd.Index(counts.index.astype(str)).map(gmap), index=counts.index)
    valid = mapped.notna()
    df = counts.loc[valid.values].copy()
    df.index = mapped[valid].astype("int64").values
    df = df.groupby(level=0).sum()  # aggregate duplicate Entrez by sum

    if species not in ("mouse", "monkey"):
        orth = pd.read_csv(_META / "Table_of_orthologs.csv")
        key = f"Entrez.{species.capitalize()}"
        orth = orth.dropna(subset=["Entrez.Mouse"]).drop_duplicates(subset=[key])
        omap = dict(zip(orth[key].astype("int64"), orth["Entrez.Mouse"].astype("int64")))
        newidx = pd.Series(pd.Index(df.index).map(omap), index=df.index)
        v = newidx.notna()
        df = df.loc[v.values]
        df.index = newidx[v].astype("int64").values
        df = df[~df.index.duplicated(keep="first")]

    df.index = df.index.astype("int64").astype(str)
    return df


# ---------- normalization steps ----------
def filter_genes(counts: pd.DataFrame, count_threshold=10, percent_threshold=20) -> pd.DataFrame:
    keep = (counts >= count_threshold).sum(axis=1) >= counts.shape[1] * percent_threshold / 100.0
    return counts.loc[keep.values]


def rle_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    X = counts.values.astype(float)  # genes x samples
    with np.errstate(divide="ignore"):
        logX = np.log(X)
    ref_log = logX.mean(axis=1)            # -inf if gene has any zero
    finite_gene = np.isfinite(ref_log)
    ratios = logX[finite_gene] - ref_log[finite_gene][:, None]
    ratios = np.where(np.isfinite(ratios), ratios, np.nan)
    factors = np.exp(np.nanmedian(ratios, axis=0))
    factors = factors / np.exp(np.mean(np.log(factors)))  # geomean 1 (edgeR)
    lib = X.sum(axis=0)
    eff = lib * factors
    norm = X / eff[None, :] * 1e7
    return pd.DataFrame(norm, index=counts.index, columns=counts.columns)


def log_transform(df: pd.DataFrame) -> pd.DataFrame:
    return np.log10(df + 1.0)


def scale_per_sample(df: pd.DataFrame) -> pd.DataFrame:
    """R scale() on genes x samples => standardize each COLUMN (sample) across genes, ddof=1."""
    X = df.values.astype(float)
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    return pd.DataFrame((X - mu) / sd, index=df.index, columns=df.columns)


def yugene(df: pd.DataFrame) -> pd.DataFrame:
    """Port of tAge YuGene(), applied per column (sample)."""
    X = df.values.astype(float)
    X = X - np.nanmin(X, axis=0, keepdims=True)  # shift each column to >=0
    out = np.empty_like(X)
    for j in range(X.shape[1]):
        col = X[:, j]
        order = np.argsort(-col, kind="mergesort")  # decreasing, stable
        svals = col[order]
        total = svals.sum()
        if total == 0:
            out[:, j] = 1.0
            continue
        cumprop = np.cumsum(svals) / total
        # carry forward for duplicates
        for i in range(1, len(cumprop)):
            if svals[i] == svals[i - 1]:
                cumprop[i] = cumprop[i - 1]
        final = 1.0 - cumprop
        res = np.empty_like(col)
        res[order] = final
        out[:, j] = res
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def align_to_gene_list(df: pd.DataFrame, gene_list: list[str]) -> pd.DataFrame:
    return df.reindex(index=gene_list)  # missing -> NaN


def control_subtraction(df: pd.DataFrame, control_mask: np.ndarray | None) -> pd.DataFrame:
    """Subtract per-gene (row) median of control samples. df: genes x samples."""
    if control_mask is None or control_mask.sum() == 0:
        xc = df.median(axis=1, skipna=True)
    else:
        xc = df.loc[:, control_mask].median(axis=1, skipna=True)
    return df.sub(xc, axis=0)


# ---------- full pipeline ----------
def preprocess(counts: pd.DataFrame, meta: pd.DataFrame, species: str,
               gene_mapping_type: str = "Ensembl",
               control_group_column: str | None = None,
               control_group_label: str | None = None,
               count_threshold=10, percent_threshold=20):
    """counts: genes x samples (raw). meta: samples x vars (index == counts.columns).
    Returns dict with 'scaled_diff' and 'yugene_diff' as samples x gene_list DataFrames."""
    gene_list = load_gene_list()
    counts = counts.loc[:, meta.index]  # align sample order

    f = filter_genes(counts, count_threshold, percent_threshold)
    m = map_genes(f, species, gene_mapping_type)
    r = rle_normalize(m)
    lg = log_transform(r)
    sc = scale_per_sample(lg)
    yg = yugene(sc)

    sc_a = align_to_gene_list(sc, gene_list)
    yg_a = align_to_gene_list(yg, gene_list)

    if control_group_column is not None and control_group_label is not None:
        cmask = (meta[control_group_column].astype(str) == str(control_group_label)).values
    else:
        cmask = None

    sc_diff = control_subtraction(sc_a, cmask)
    yg_diff = control_subtraction(yg_a, cmask)

    # return samples x genes
    return {"scaled_diff": sc_diff.T, "yugene_diff": yg_diff.T}
