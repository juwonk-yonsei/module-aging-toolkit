"""Attention-MIL, supplementary re-analysis (Supplementary Methods S4 and Results,
Supplementary Table S16).

Instances: each donor x macrophage annotation (Macrophage, Macrophage_Alveolar) pseudobulk
is split into mini-pseudobulks of ~60 randomly assigned cells (scripts/mil_extract.py);
instance features = 23 module scores; bag = donor; label = IPF vs Control.

For 10 seeds x 5-fold donor-stratified CV: out-of-fold donor-level AUC for attention
pooling, mean pooling and an L2 logistic regression on donor-mean module scores (same folds).
Donor-bootstrap 95% CI of the seed-averaged out-of-fold AUC.

Outputs: rev13_mil_auc.csv, rev13_mil_auc_summary.csv
"""
import os
import sys

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

import rev_common as rc

sys.path.insert(0, str(rc.TOOLKIT / "scripts"))
import mil_model as M  # noqa: E402
from mil_score import MODULES  # noqa: E402

N_SEED = 10
N_BOOT = 2000
FINE = ["Macrophage", "Macrophage_Alveolar"]


def predict(model, bags):
    X, m = M.pad_bags(bags)
    model.eval()
    with torch.no_grad():
        return torch.sigmoid(model(X, m)[0]).cpu().numpy()


def main():
    sc = pd.read_csv(rc.RES / "mil" / "inst_scores_adams_bs60.csv", index_col=0)
    sc = sc[sc.celltype.isin(FINE) & sc.Disease_Identity.isin(["IPF", "Control"])].copy()
    donors = sc.Subject_Identity.values
    uniq = pd.Index(sorted(sc.Subject_Identity.unique()))
    y = (sc.groupby("Subject_Identity").Disease_Identity.first().reindex(uniq) == "IPF").astype(int).values
    X = sc[MODULES].values.astype(np.float32)
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    rows_of = {d: np.where(donors == d)[0] for d in uniq}
    dmean = np.stack([X[rows_of[d]].mean(0) for d in uniq])
    print(f"instances={len(sc)}, donors={len(uniq)} (IPF {y.sum()}, Control {len(y) - y.sum()}); "
          f"instances/donor median {np.median([len(v) for v in rows_of.values()]):.0f}; "
          f"cells/instance median {sc.n_cells.median() if 'n_cells' in sc else np.nan}")
    recs, oof_all = [], {k: [] for k in ["attention", "mean", "logistic"]}
    for seed in range(N_SEED):
        oof = {k: np.full(len(uniq), np.nan) for k in oof_all}
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        for tr, va in skf.split(np.arange(len(uniq)), y):
            bags_tr = [X[rows_of[d]] for d in uniq[tr]]
            bags_va = [X[rows_of[d]] for d in uniq[va]]
            for pool in ["attention", "mean"]:
                model, _ = M.train_eval(bags_tr, y[tr], bags_tr, y[tr], in_dim=len(MODULES), pool=pool,
                                        seed=seed, max_epochs=200, patience=25)
                oof[pool][va] = predict(model, bags_va)
            lr = LogisticRegression(C=1.0, max_iter=2000).fit(dmean[tr], y[tr])
            oof["logistic"][va] = lr.predict_proba(dmean[va])[:, 1]
        for k, v in oof.items():
            recs.append(dict(seed=seed, model=k, auc=roc_auc_score(y, v)))
            oof_all[k].append(v)
        print(f"seed {seed}: " + ", ".join(f"{k} {roc_auc_score(y, v):.3f}" for k, v in oof.items()))
    R = pd.DataFrame(recs)
    rc.save(R, "rev13_mil_auc.csv")
    rng = np.random.default_rng(rc.SEED)
    rows = []
    for k, vs in oof_all.items():
        p = np.mean(vs, 0)
        boots = []
        for _ in range(N_BOOT):
            ii = np.r_[rng.choice(np.where(y == 1)[0], y.sum()), rng.choice(np.where(y == 0)[0], (1 - y).sum())]
            boots.append(roc_auc_score(y[ii], p[ii]))
        a = R[R.model == k].auc
        rows.append(dict(model=k, auc_seed_mean=a.mean(), auc_seed_sd=a.std(), auc_seed_min=a.min(),
                         auc_seed_max=a.max(), auc_of_mean_prediction=roc_auc_score(y, p),
                         ci_lo=np.percentile(boots, 2.5), ci_hi=np.percentile(boots, 97.5),
                         n_donors=len(uniq), n_folds=5, n_seeds=N_SEED))
    S = pd.DataFrame(rows)
    rc.save(S, "rev13_mil_auc_summary.csv")
    print(S.round(3).to_string())


if __name__ == "__main__":
    main()
