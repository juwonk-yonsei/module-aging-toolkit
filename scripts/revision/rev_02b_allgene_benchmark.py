"""All-gene (composite-equivalent) elastic net on the same dataset-grouped folds (Section 3.1).

Same recipe as the module clocks, trained on every gene of the released rodent training matrix
(18,286 genes), GroupKFold(5) by dataset. Benchmark for interpreting module-clock accuracy.

Output: rev02b_allgene_benchmark.csv
"""
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.base import clone
from sklearn.model_selection import GroupKFold

import rev_common as rc
from rev_02_species_cv import TARGETS, base_model, load_training, safe_r


def fit_fold(X, y, tr, te):
    return te, clone(base_model(n_jobs=4)).fit(X[tr], y[tr]).predict(X[te])


def main():
    ann, expr = load_training()
    X_all = expr.T.values
    groups, species = ann["Source"].astype(str).values, ann["Species"].astype(str).values
    rows = []
    for oc, col in TARGETS.items():
        y_all = pd.to_numeric(ann[col], errors="coerce").values
        ok = np.isfinite(y_all)
        X, y, g, sp = X_all[ok], y_all[ok], groups[ok], species[ok]
        folds = list(GroupKFold(n_splits=5).split(X, y, g))
        res = Parallel(n_jobs=5)(delayed(fit_fold)(X, y, tr, te) for tr, te in folds)
        oof = np.full(len(y), np.nan)
        for te, p in res:
            oof[te] = p
        rows.append(dict(kind="all_genes", outcome=oc, n_genes=X.shape[1], r_oof_all=safe_r(oof, y),
                         r_oof_mouse=safe_r(oof[sp == "Mouse"], y[sp == "Mouse"]),
                         r_oof_rat=safe_r(oof[sp == "Rat"], y[sp == "Rat"])))
        print(rows[-1], flush=True)
    rc.save(pd.DataFrame(rows), "rev02b_allgene_benchmark.csv")


if __name__ == "__main__":
    main()
