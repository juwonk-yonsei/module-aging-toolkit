"""Shared helpers for the revision analyses (revised manuscript; toolkit v0.2.0).

Inputs are read from MAT_ROOT (data/ and results/ produced by the main pipeline);
outputs are written to MAT_OUT (default MAT_ROOT/results/revision).
"""
from __future__ import annotations

import os
import sys
import warnings
from collections import Counter
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

TOOLKIT = Path(__file__).resolve().parents[2]
REPO = Path(os.environ.get("MAT_ROOT", str(TOOLKIT)))
OUT = Path(os.environ.get("MAT_OUT", str(REPO / "results" / "revision")))
FIG = OUT / "figs"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
TAGE_DIR = Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge")))

for p in (TOOLKIT / "scripts", TOOLKIT / "third_party", TAGE_DIR / "inst" / "python"):
    sys.path.insert(0, str(p))
import tage_prep as tp  # noqa: E402

DATA = REPO / "data"
RES = REPO / "results"
MC = DATA / "module_clocks"
CM = DATA / "clock_models"
HUMAN_ADJ = 122.5          # tAge species factor for human (maximum lifespan, years)
SEED = 20260808

OKABE = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]


# ----------------------------------------------------------------------------- clocks
def make_labels(fullnames):
    lab = {fn: fn.split("|", 1)[1].split("/")[0].strip() for fn in fullnames}
    cnt = Counter(lab.values())
    for fn, s in list(lab.items()):
        if cnt[s] > 1:
            lab[fn] = f"{s} ({fn.split('|')[0]})"
    return lab


def load_module_clocks(outcome="Mortality"):
    clocks = {}
    for f in sorted(MC.glob(f"module_{outcome}_*.pkl")):
        d = joblib.load(f)
        clocks[f"{d['module']}|{d['annotation']}"] = d
    lab = make_labels(list(clocks))
    return {lab[k]: v for k, v in clocks.items()}


def load_composite(outcome="Mortality"):
    from tage_predict import _patch_simple_imputer  # tAge checkout (setup.sh); not needed for figures
    name = {"Mortality": "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl",
            "Chrono": "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl"}[outcome]
    m = joblib.load(CM / name)
    _patch_simple_imputer(m.named_steps["imputation"])
    return m, [str(g) for g in m.feature_names_in_]


def linear_parts(pipeline):
    """Decompose a fitted SimpleImputer->StandardScaler->ElasticNet pipeline into
    (median, mean, scale, coef, intercept) so it can be applied as a linear score."""
    steps = pipeline.named_steps
    imp = [s for n, s in steps.items() if "imput" in n.lower()][0]
    sca = [s for n, s in steps.items() if "scal" in n.lower()][0]
    sel = [s for s in steps.values() if hasattr(s, "get_support")]
    est = list(steps.values())[-1]
    med = np.asarray(imp.statistics_, float)
    mean = np.asarray(sca.mean_, float) if sca.mean_ is not None else np.zeros_like(med)
    scale = np.asarray(sca.scale_, float) if sca.scale_ is not None else np.ones_like(med)
    coef = np.zeros_like(med)
    keep = sel[0].get_support() if sel else np.ones(len(med), bool)
    coef[keep] = np.asarray(est.coef_, float)
    return med, mean, scale, coef, float(est.intercept_)


def linear_score(X, median, mean, scale, coef, intercept):
    """X: samples x genes (NaN = undetected). Same arithmetic as pipeline.predict."""
    Xi = np.where(np.isnan(X), median[None, :], X)
    return ((Xi - mean[None, :]) / scale[None, :]) @ coef + intercept


def train_gene_stats():
    """Per-gene median / mean / sd in the rodent training matrix, computed the same way
    as SimpleImputer(median) -> StandardScaler. Used to standardize random genes."""
    cache = OUT / "_cache_train_gene_stats.pkl"
    if cache.exists():
        return pd.read_pickle(cache)
    from rev_02_species_cv import load_training
    _, expr = load_training()                      # genes x samples
    E, genes = expr.values.T, expr.index
    med = np.nanmedian(E, axis=0)
    Ei = np.where(np.isnan(E), med[None, :], E)
    st = pd.DataFrame({"median": med, "mean": Ei.mean(0), "scale": Ei.std(0)}, index=genes.astype(str))
    st = st[np.isfinite(st["median"]) & (st["scale"] > 0)]
    st.to_pickle(cache)
    return st


class NullScorer:
    """Null versions of a module clock applied to one scaled_diff matrix.

    perm_coef(n):     the module's coefficients shuffled across its own genes.
    random_genes(n):  size-matched random genes carrying the module's coefficients,
                      standardized with their own training statistics.
    geneset(signed):  equal-weight mean of the standardized module genes; signed=True
                      uses the sign of the trained coefficient.
    """

    def __init__(self, scaled_diff, universe_stats):
        X = scaled_diff.copy()
        X.columns = X.columns.map(str)
        self.X = X
        st = universe_stats.loc[[g for g in universe_stats.index if g in set(X.columns)]]
        Xu = X[st.index].values
        Xu = np.where(np.isnan(Xu), st["median"].values[None, :], Xu)
        self.Zu = (Xu - st["mean"].values[None, :]) / st["scale"].values[None, :]
        self.universe = np.array(st.index)

    def _module_z(self, clock):
        med, mu, sc, coef, b = linear_parts(clock["pipeline"])
        Xm = self.X.reindex(columns=[str(g) for g in clock["genes"]]).values
        Xm = np.where(np.isnan(Xm), med[None, :], Xm)
        return (Xm - mu[None, :]) / sc[None, :], coef, b

    def perm_coef(self, clock, n, rng):
        Z, coef, b = self._module_z(clock)
        C = np.stack([rng.permutation(coef) for _ in range(n)], axis=1)
        return Z @ C + b

    def random_genes(self, clock, n, rng):
        _, coef, b = self._module_z(clock)
        out = np.empty((self.Zu.shape[0], n))
        for i in range(n):
            idx = rng.choice(self.Zu.shape[1], size=len(coef), replace=False)
            out[:, i] = self.Zu[:, idx] @ coef + b
        return out

    def geneset(self, clock, signed=False):
        Z, coef, _ = self._module_z(clock)
        w = np.sign(coef) if signed else np.ones_like(coef)
        if signed:
            keep = w != 0
            return (Z[:, keep] * w[keep]).mean(1)
        return Z.mean(1)


def score_modules(scaled_diff, clocks=None, composite=None, composite_adj=False):
    """Return samples x (23 modules + Composite). Composite on the raw model scale
    (Δlog10 hazard units, same as the modules) unless composite_adj=True."""
    if clocks is None:
        clocks = load_module_clocks()
    X = scaled_diff.copy()
    X.columns = X.columns.map(str)
    out = {}
    for lab, d in clocks.items():
        out[lab] = d["pipeline"].predict(X.reindex(columns=[str(g) for g in d["genes"]]).values)
    if composite is not False:
        comp, cg = composite if composite is not None else load_composite()
        v = comp.predict(X.reindex(columns=cg).values)
        out["Composite"] = v * (HUMAN_ADJ if composite_adj else 1.0)
    return pd.DataFrame(out, index=scaled_diff.index)


# ----------------------------------------------------------------------------- data
def load_adams(min_cells=50):
    pb = pd.read_pickle(RES / "pseudobulk_counts.pkl")
    meta = pd.read_csv(RES / "pseudobulk_meta.csv")
    meta["group"] = meta["group"].astype(str)
    cov = pd.read_csv(DATA / "ipf" / "subject_covariates.csv").set_index("Subject_Identity")
    meta = meta.merge(cov[["age", "Sex", "disease"]], left_on="Subject_Identity",
                      right_index=True, how="left")
    meta["sexM"] = (meta["Sex"] == "M").astype(float)
    meta["lib_size"] = pb[meta["group"]].sum(0).values
    return pb, meta[meta.n_cells >= min_cells].copy()


def testable_celltypes(meta, min_donors=4):
    keep = []
    for ct, sub in meta.groupby("Manuscript_Identity"):
        n = sub.Disease_Identity.value_counts()
        if n.get("IPF", 0) >= min_donors and n.get("Control", 0) >= min_donors:
            keep.append(ct)
    return keep


def preprocess_ct(pb, meta, ct, diseases=("IPF", "Control", "COPD"), counts=None,
                  normalization="RLE"):
    sub = meta[(meta.Manuscript_Identity == ct) & meta.Disease_Identity.isin(diseases)]
    md = sub.set_index("group")
    cnt = pb[sub.group.tolist()] if counts is None else counts[sub.group.tolist()]
    if normalization == "RLE":
        pp = tp.preprocess(cnt, md, species="human", gene_mapping_type="Ensembl",
                           control_group_column="Disease_Identity", control_group_label="Control")
        X = pp["scaled_diff"]
    else:
        X = preprocess_cpm(cnt, md)
    return X, md


def preprocess_cpm(counts, md):
    """tAge preprocessing with library-size (CPM-like) scaling instead of RLE."""
    gene_list = tp.load_gene_list()
    f = tp.filter_genes(counts.loc[:, md.index])
    m = tp.map_genes(f, "human", "Ensembl")
    norm = m / m.sum(0) * 1e7
    sc = tp.scale_per_sample(tp.log_transform(norm))
    sc = tp.align_to_gene_list(sc, gene_list)
    cmask = (md["Disease_Identity"] == "Control").values
    return tp.control_subtraction(sc, cmask).T


def adams_scores(min_cells=50, normalization="RLE", pb=None, meta=None, extra_scorers=None):
    """Per-sample module scores for all testable Adams cell types (long table)."""
    if pb is None:
        pb, meta = load_adams(min_cells)
    clocks = load_module_clocks()
    comp = load_composite()
    recs = []
    for ct in testable_celltypes(meta):
        X, md = preprocess_ct(pb, meta, ct, normalization=normalization)
        S = score_modules(X, clocks, comp)
        if extra_scorers:
            for name, fn in extra_scorers.items():
                S = S.join(fn(X).add_prefix(f"{name}::"))
        S = S.join(md[["Subject_Identity", "Disease_Identity", "age", "sexM", "n_cells",
                       "lib_size"]])
        S["celltype"] = ct
        recs.append(S)
    return pd.concat(recs)


def module_names(df):
    return [c for c in df.columns if c not in
            ("Subject_Identity", "Disease_Identity", "age", "sexM", "n_cells", "lib_size",
             "celltype", "Composite", "group") and "::" not in c]


# ----------------------------------------------------------------------------- stats
def ols(y, X):
    """OLS with classical SE. Returns beta, se, t, p, df for every column of X."""
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[ok], X[ok]
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ X.T @ y
    resid = y - X @ b
    df = n - k
    s2 = resid @ resid / df
    se = np.sqrt(np.diag(XtX_inv) * s2)
    t = b / se
    p = 2 * stats.t.sf(np.abs(t), df)
    return b, se, t, p, df


def disease_effect(df, y, case="IPF", ref="Control", covars=("age", "sexM")):
    sub = df[df.Disease_Identity.isin([case, ref])].dropna(subset=list(covars) + [y])
    X = np.column_stack([np.ones(len(sub)), (sub.Disease_Identity == case).astype(float)]
                        + [sub[c].values for c in covars])
    b, se, t, p, dof = ols(sub[y].values.astype(float), X)
    q = stats.t.ppf(0.975, dof)
    return dict(beta=b[1], se=se[1], p=p[1], lo=b[1] - q * se[1], hi=b[1] + q * se[1],
                df=dof, n_case=int((sub.Disease_Identity == case).sum()),
                n_ref=int((sub.Disease_Identity == ref).sum()))


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def save(df, name, index=False):
    path = OUT / name
    df.to_csv(path, index=index)
    print(f"saved {path}")
    return path


def setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                         "figure.dpi": 150, "savefig.dpi": 600, "font.family": "DejaVu Sans",
                         "axes.spines.top": False, "axes.spines.right": False})
    return plt
