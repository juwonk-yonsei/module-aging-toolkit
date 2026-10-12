"""Rerun the sub-state permutation of scripts/v3/r5/r5_concordance.py with its stored seed.

Counts the permutations kept per sub-state and compares p with the stored table.
Writes results/v19/substate_rerun_v19.csv only.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "v3" / "r5"))
from r5_concordance import B, OUT, SEED, effects


def main():
    s = pd.read_parquet(OUT / "scores_r5.parquet",
                        columns=["cohort", "variant", "clock", "celltype", "unit", "donor", "disease", "score"])
    s = s[s.clock.ne("Composite") & ~s.clock.astype(str).str.startswith("part:")]
    rng = np.random.default_rng(SEED)
    rows = []
    for variant, kind in [("substate", "marker sub-state"), ("randsplit_001", "random split")]:
        A = s[(s.variant == variant) & (s.cohort == "GSE136831")]
        H = s[(s.variant == variant) & (s.cohort == "GSE135893")]
        eA, eH = effects(A), effects(H)
        labA = A.groupby("donor").disease.first()
        labH = H.groupby("donor").disease.first()
        nulls = [(effects(A, pd.Series(rng.permutation(labA.values), index=labA.index)),
                  effects(H, pd.Series(rng.permutation(labH.values), index=labH.index))) for _ in range(B)]
        for key in sorted(set(eA) & set(eH)):
            r = float(np.corrcoef(eA[key], eH[key])[0, 1])
            nr = np.array([np.corrcoef(a[key], h[key])[0, 1] for a, h in nulls if key in a and key in h])
            k = int(np.sum(np.abs(nr) >= abs(r)))
            rows.append(dict(split=kind, celltype=key[0], substate=key[1], r=r, k=k,
                             B_valid=len(nr), p_two_sided=(1 + k) / (1 + len(nr))))
    C = pd.DataFrame(rows)
    stored = pd.read_csv(OUT / "substate_concordance_r5.csv")
    C = C.merge(stored.rename(columns={"p_two_sided": "p_stored", "B": "B_stored", "r": "r_stored"}),
                on=["split", "celltype", "substate"], how="left")
    C["p_matches"] = np.isclose(C.p_two_sided, C.p_stored, rtol=0, atol=1e-12)
    out = REPO / "results" / "v19"
    out.mkdir(parents=True, exist_ok=True)
    C.to_csv(out / "substate_rerun_v19.csv", index=False)
    print(C.to_string(index=False))


if __name__ == "__main__":
    main()
