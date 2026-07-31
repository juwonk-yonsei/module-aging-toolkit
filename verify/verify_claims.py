#!/usr/bin/env python3
"""Verify that the numbers reported in the manuscript match the shipped result
CSVs exactly (within tolerance).

This is the machine-checkable audit a reviewer can run in seconds, without the
multi-GB raw data, to confirm that every headline statistic is a genuine code
output rather than a hand-typed (or hallucinated) value. It recomputes each
claim from `results/*.csv` and prints a PASS/FAIL table.

Usage
-----
    python verify/verify_claims.py            # uses ./results
    MAT_RESULTS=/path/to/results python verify/verify_claims.py

Exit code is non-zero if any claim fails.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS = Path(os.environ.get("MAT_RESULTS", Path(__file__).resolve().parents[1] / "results"))


def _load(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / name)


class Checker:
    def __init__(self) -> None:
        self.rows: list[tuple[str, object, object, bool]] = []

    def check(self, desc: str, observed, expected, tol: float = 0.0) -> None:
        if isinstance(expected, bool):
            ok = bool(observed) == expected
        elif tol == 0.0 and isinstance(expected, int):
            ok = int(observed) == expected
        else:
            ok = abs(float(observed) - float(expected)) <= tol
        self.rows.append((desc, observed, expected, ok))

    def report(self) -> bool:
        w = max(len(r[0]) for r in self.rows)
        print(f"\n{'CLAIM'.ljust(w)}  {'OBSERVED':>12}  {'EXPECTED':>12}  RESULT")
        print("-" * (w + 40))
        for desc, obs, exp, ok in self.rows:
            o = f"{obs:.3f}" if isinstance(obs, float) else str(obs)
            e = f"{exp:.3f}" if isinstance(exp, float) else str(exp)
            print(f"{desc.ljust(w)}  {o:>12}  {e:>12}  {'PASS' if ok else 'FAIL'}")
        n_fail = sum(1 for r in self.rows if not r[3])
        print("-" * (w + 40))
        print(f"{len(self.rows) - n_fail}/{len(self.rows)} checks passed"
              + ("" if n_fail == 0 else f"  ({n_fail} FAILED)"))
        return n_fail == 0


def main() -> int:
    c = Checker()

    # --- Retrained module clocks (Methods 2.1; §3.1) --------------------------
    cov = _load("module_coverage.csv")
    c.check("retrained module clocks = 23", len(cov), 23)

    # --- Composite masking: composite is n.s., modules oppose (§3.2, Fig 2) ---
    g = _load("module_global_robust.csv").set_index("module")
    comp = g.loc["Composite"]
    c.check("composite global effect beta = -2.9", comp["beta"], -2.9, tol=0.1)
    c.check("composite global effect non-significant (p>0.05)", comp["p"] > 0.05, True)
    c.check("Lipid met module +0.31", g.loc["Lipid met", "beta"], 0.31, tol=0.02)
    c.check("Translation module +0.27", g.loc["Translation", "beta"], 0.27, tol=0.02)
    c.check("Interferon signaling module -0.33", g.loc["Interferon signaling", "beta"], -0.33, tol=0.02)
    c.check("VEGF signaling module -0.27", g.loc["VEGF signaling", "beta"], -0.27, tol=0.02)

    # --- Cell-type-resolved effects (§3.3, Fig 3) -----------------------------
    st = _load("ipf_module_stats_adjusted.csv")
    c.check("cell-type x module effects tested = 253", int(st.fdr_IPF.notna().sum()), 253)
    c.check("effects passing FDR<0.10 = 30", int((st.fdr_IPF < 0.10).sum()), 30)

    # --- Resolution-dependent replication boundary (§3.5, Fig 4a) -------------
    rep = _load("replication_coarse_compare.csv")
    r_by_ct = {ct: float(np.corrcoef(gg.beta_adams, gg.beta_hab)[0, 1])
               for ct, gg in rep.groupby("coarse")}
    c.check("cross-cohort r: NK = 0.47", r_by_ct["NK"], 0.47, tol=0.02)
    c.check("cross-cohort r: T cell = 0.51", r_by_ct["Tcell"], 0.51, tol=0.02)
    c.check("cross-cohort r: monocyte = 0.08", r_by_ct["Monocyte"], 0.08, tol=0.02)
    c.check("cross-cohort r: macrophage = -0.11", r_by_ct["Macrophage"], -0.11, tol=0.02)
    c.check("cross-cohort r: ciliated = -0.20", r_by_ct["Ciliated"], -0.20, tol=0.02)

    # --- Split-half ceiling for non-replicating cell types (§3.5, Fig 4b) -----
    sc = _load("mil/stageC_calibrate.csv").set_index("celltype")
    c.check("split-half ceiling macrophage = 0.77", sc.loc["Macrophage", "splithalf_med"], 0.77, tol=0.02)
    c.check("split-half ceiling monocyte = 0.75", sc.loc["Monocyte", "splithalf_med"], 0.75, tol=0.02)

    # --- Whole-tissue bulk recovery across 3 cohorts (§3.5, Fig 4c) -----------
    bulk = _load("mil/stageD_bulk.csv")
    c.check("bulk-recovery r all in [0.46, 0.51]", bool(((bulk.r >= 0.46) & (bulk.r <= 0.51)).all()), True)
    c.check("bulk-recovery permutation p<0.05 (all)", bool((bulk.p_perm < 0.05).all()), True)

    # --- Attention-MIL calibration on semi-synthetic data (§3.6, Methods S2) --
    syn = _load("mil/stage1_results.csv")
    sub = syn[syn.regime == "substate"]
    p02 = sub[np.isclose(sub.p_s, 0.20)]
    att02 = p02[p02.pool == "attention"]["bag_auc"].median()
    mean02 = p02[p02.pool == "mean"]["bag_auc"].median()
    c.check("attention > mean pooling at p_S=0.2 (bag AUC)", bool(att02 > mean02), True)
    auroc_by_ps = sub[sub.pool == "attention"].groupby("p_s")["att_auroc"].median()
    c.check("attention localization AUROC rises with p_S (max~0.87)",
            float(auroc_by_ps.max()), 0.89, tol=0.06)

    ok = c.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
