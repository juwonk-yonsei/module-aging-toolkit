"""Cross-species robustness and performance benchmarks of the retrained module clocks
(Section 3.1, Fig. 2a-b, Supplementary Table S4).

For every module x outcome (Mortality, Chrono):
  * dataset-grouped 5-fold CV exactly as in training -> out-of-fold r overall,
    and separately within mouse and within rat samples;
  * mouse -> rat (train on all mouse datasets, test on all rat datasets) and
    rat -> mouse;
  * benchmarks on the same folds: size-matched random gene sets (N_RAND draws)
    and a clock trained on the union of all module genes.

Outputs: rev02_species_cv.csv, rev02_benchmarks.csv, rev02_oof_predictions.pkl
"""
import os
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import pearsonr
from sklearn.base import clone
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import rev_common as rc

TARGETS = {"Mortality": "Difference.Hazard.log10",
           "Chrono": "Difference.Chronological_age.Normalized_by_species_max_lifespan"}
N_RAND = int(os.environ.get("REV02_NRAND", "10"))
N_JOBS = int(os.environ.get("REV02_JOBS", "40"))


def base_model(n_jobs=1):
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                         ElasticNetCV(l1_ratio=[.5, .9, 1.0], alphas=np.logspace(-3, 1, 20),
                                      cv=3, max_iter=4000, n_jobs=n_jobs))


def load_training():
    cache = rc.OUT / "_cache_rodent_expr.pkl"
    ann = pd.read_excel(rc.DATA / "rodent" / "Data_annotation_relative_rodents.xlsx").set_index("Sample")
    if cache.exists():
        expr = pd.read_pickle(cache)
    else:
        expr = pd.read_csv(rc.DATA / "rodent" / "Expression_data_relative_rodents_Scaled.csv", index_col=0)
        if not set(ann.index[:5]).issubset(set(expr.columns)):
            expr = expr.T
        expr.index = expr.index.astype(str)
        expr.to_pickle(cache)
    common = [s for s in ann.index if s in expr.columns]
    return ann.loc[common], expr[common]


def safe_r(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return pearsonr(a[ok], b[ok])[0] if ok.sum() > 3 and np.std(a[ok]) > 0 else np.nan


def run_task(kind, outcome, module, genes, X_all, y_all, groups, species, seed=None):
    """X_all: samples x genes (numpy) for `genes`; returns dict of metrics."""
    ok = np.isfinite(y_all)
    X, y, g, sp = X_all[ok], y_all[ok], groups[ok], species[ok]
    oof = np.full(len(y), np.nan)
    for tr, te in GroupKFold(n_splits=5).split(X, y, g):
        oof[te] = clone(base_model()).fit(X[tr], y[tr]).predict(X[te])
    rec = dict(kind=kind, outcome=outcome, module=module, n_genes=X.shape[1], seed=seed,
               r_oof_all=safe_r(oof, y), r_oof_mouse=safe_r(oof[sp == "Mouse"], y[sp == "Mouse"]),
               r_oof_rat=safe_r(oof[sp == "Rat"], y[sp == "Rat"]))
    if kind == "module":
        m, r = sp == "Mouse", sp == "Rat"
        p_rat = clone(base_model()).fit(X[m], y[m]).predict(X[r])
        p_mouse = clone(base_model()).fit(X[r], y[r]).predict(X[m])
        rec["r_mouse_to_rat"] = safe_r(p_rat, y[r])
        rec["r_rat_to_mouse"] = safe_r(p_mouse, y[m])
        oof_rat = np.full(r.sum(), np.nan)
        Xr, yr, gr = X[r], y[r], g[r]
        for tr, te in GroupKFold(n_splits=len(np.unique(gr))).split(Xr, yr, gr):
            oof_rat[te] = clone(base_model()).fit(Xr[tr], yr[tr]).predict(Xr[te])
        rec["r_within_rat_lodo"] = safe_r(oof_rat, yr)
        rec["oof"] = oof
    return rec


def main():
    t0 = time.time()
    ann, expr = load_training()
    mem = pd.read_csv(rc.DATA / "supp" / "module_membership_rodent.csv")
    mem["entrez"] = mem["entrez"].astype(str)
    mod_genes = {f"{m}|{a}": [x for x in g.entrez if x in expr.index]
                 for (m, a), g in mem.groupby(["module", "annotation"])}
    lab = rc.make_labels(list(mod_genes))
    groups = ann["Source"].values
    species = ann["Species"].values
    universe = np.array(expr.index)
    E = expr.T  # samples x genes
    print(f"loaded {E.shape} in {time.time()-t0:.0f}s; species={pd.Series(species).value_counts().to_dict()}")

    rng = np.random.default_rng(rc.SEED)
    tasks = []
    for outcome, col in TARGETS.items():
        y = pd.to_numeric(ann[col], errors="coerce").values
        for key, genes in mod_genes.items():
            tasks.append(("module", outcome, lab[key], genes, y, None))
            for s in range(N_RAND):
                rg = list(rng.choice(universe, size=len(genes), replace=False))
                tasks.append(("random", outcome, lab[key], rg, y, s))
        allg = sorted({g for v in mod_genes.values() for g in v})
        tasks.append(("union_modules", outcome, "All module genes", allg, y, None))

    print(f"{len(tasks)} tasks on {N_JOBS} workers ...")
    res = Parallel(n_jobs=N_JOBS, verbose=5)(
        delayed(run_task)(k, o, m, g, E[g].values, y, groups, species, s)
        for (k, o, m, g, y, s) in tasks)

    oof = {(r["outcome"], r["module"]): r.pop("oof") for r in res if "oof" in r}
    pd.to_pickle({"oof": oof, "species": species, "groups": groups,
                  "y": {o: pd.to_numeric(ann[c], errors="coerce").values for o, c in TARGETS.items()}},
                 rc.OUT / "rev02_oof_predictions.pkl")
    R = pd.DataFrame(res)
    mods = R[R.kind == "module"].drop(columns=["seed"])
    rand = (R[R.kind == "random"].groupby(["outcome", "module"])
            .agg(r_random_mean=("r_oof_all", "mean"), r_random_sd=("r_oof_all", "std"),
                 r_random_max=("r_oof_all", "max")).reset_index())
    mods = mods.merge(rand, on=["outcome", "module"], how="left")
    mods["exceeds_all_random"] = mods.r_oof_all > mods.r_random_max
    rc.save(mods, "rev02_species_cv.csv")
    rc.save(R[R.kind != "module"].drop(columns=[c for c in R.columns if c.startswith("r_mouse_to") or
                                                c.startswith("r_rat_to") or c == "r_within_rat_lodo"],
                                       errors="ignore"), "rev02_benchmarks.csv")
    print(f"done in {(time.time()-t0)/60:.1f} min")
    print(mods.groupby("outcome")[["r_oof_all", "r_oof_mouse", "r_oof_rat", "r_mouse_to_rat",
                                   "r_rat_to_mouse", "r_random_mean"]].median().round(3))


if __name__ == "__main__":
    main()
