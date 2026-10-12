"""Refit the genes of the serialized red mortality module with the published pipeline.

Same fit as scripts/v3/r5/r5_scrambled.py fit_one: median imputation, centering,
SelectKBest(all), ElasticNet at the serialized alpha and l1_ratio, full rodent
training matrix. Output is one row of counts and correlations; no gene names or
coefficients are written.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "v3"))
sys.path.insert(0, str(REPO / "scripts" / "v3" / "r5"))
from module_aging.official_clocks import load_clocks
from rodent_reanalysis import TARGETS, load_matrix
from r5_scrambled import fit_one

OUT = REPO / "results" / "v19"
MODULE = "red"


def main():
    clock = load_clocks(outcome="Mortality")[MODULE]
    ann, genes, _ = load_matrix()
    gene_index = {g: i for i, g in enumerate(genes)}
    y = pd.to_numeric(ann[TARGETS["Mortality"]], errors="coerce").to_numpy(float)
    path = str(REPO / "results" / "v3" / "rodent" / "_expr_samples_by_genes.npy")
    gns = [str(g) for g in clock["genes"]]
    med, mean, coef, _ = fit_one(gns, clock["alpha"], clock["l1_ratio"], y, path, gene_index)
    rho = float(stats.spearmanr(coef, clock["coef"])[0])
    row = dict(module=MODULE, n_genes=len(gns), spearman_rho_coefficients=rho,
               max_abs_diff_training_mean=float(np.max(np.abs(mean - clock["mean"]))),
               max_abs_diff_training_median=float(np.max(np.abs(med - clock["median"]))),
               alpha=float(clock["alpha"]), l1_ratio=float(clock["l1_ratio"]))
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(OUT / "retrain_check_v19.csv", index=False)
    print(pd.DataFrame([row]).to_string(index=False))


if __name__ == "__main__":
    main()
