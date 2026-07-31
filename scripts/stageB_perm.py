"""Rigorous null for Option B: is the matched-sub-state cross-cohort concordance
better than chance? Uses cached stageB pseudobulks (no re-streaming).

Disease-permutation null: score every (donor x ct x substate) pseudobulk ONCE
(the IPF-vs-Control contrast is invariant to the control-subtraction offset, so
scoring is done once), then permute IPF/Control donor labels within each cohort
(preserving counts) and recompute the cross-cohort concordance NPERM times.

Reports observed r, empirical p, and null mean+/-sd for coarse and each matched
sub-state of Macrophage / Monocyte.
"""
import os
import sys
from pathlib import Path
import numpy as np, pandas as pd
REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
from mil_score import score_matrix, MODULES
from stageB_substate import MIN_CELLS, conc
MIL = REPO / "results" / "mil"
NPERM = 2000


def load(tag):
    return pd.read_pickle(MIL / f"stageB_pb_{tag}.pkl"), pd.read_csv(MIL / f"stageB_gm_{tag}.csv")


def scored(pbdf, gm, gmap, ct, ss=None):
    """Return DataFrame indexed by donor with 23 module scores + disease, or None."""
    sub = gm[(gm.n_cells >= MIN_CELLS) & (gm.celltype == ct)]
    if ss is not None:
        sub = sub[sub.substate == ss]
    n = sub.Disease_Identity.value_counts()
    if n.get("IPF", 0) < 3 or n.get("Control", 0) < 3:
        return None
    md = sub.set_index(sub.group.astype(str))
    pp = tp.preprocess(pbdf[md.index.tolist()], md, species="human", gene_mapping_type=gmap,
                       control_group_column="Disease_Identity", control_group_label="Control")
    S = score_matrix(pp["scaled_diff"])[MODULES].copy()
    S["donor"] = md.loc[S.index, "donor"].values
    S["dis"] = md.loc[S.index, "Disease_Identity"].values
    return S


def eff_from_labels(S, dis_by_donor):
    d = S["donor"].map(dis_by_donor)
    return S[d == "IPF"][MODULES].mean().values - S[d == "Control"][MODULES].mean().values


def main():
    Apb, Agm = load("adams"); Hpb, Hgm = load("hab")
    rng = np.random.default_rng(0)
    rows = []
    for ct in ["Macrophage", "Monocyte"]:
        for name, ss in [("coarse", None), ("monoHi", "monoHi"), ("resHi", "resHi")]:
            SA = scored(Apb, Agm, "Ensembl", ct, ss)
            SH = scored(Hpb, Hgm, "Gene.Symbol", ct, ss)
            if SA is None or SH is None:
                continue
            dA = SA.drop_duplicates("donor").set_index("donor")["dis"].to_dict()
            dH = SH.drop_duplicates("donor").set_index("donor")["dis"].to_dict()
            r_obs = conc(eff_from_labels(SA, dA), eff_from_labels(SH, dH))
            null = np.empty(NPERM)
            keysA, valsA = np.array(list(dA)), np.array(list(dA.values()))
            keysH, valsH = np.array(list(dH)), np.array(list(dH.values()))
            for i in range(NPERM):
                pA = dict(zip(keysA, rng.permutation(valsA)))
                pH = dict(zip(keysH, rng.permutation(valsH)))
                null[i] = conc(eff_from_labels(SA, pA), eff_from_labels(SH, pH))
            null = null[np.isfinite(null)]
            p = (1 + np.sum(null >= r_obs)) / (1 + len(null))
            rows.append({"celltype": ct, "level": name, "r_obs": r_obs,
                         "null_mean": float(null.mean()), "null_sd": float(null.std()),
                         "p_perm": p, "n_IPF_A": sum(v == "IPF" for v in dA.values()),
                         "n_Ctrl_A": sum(v == "Control" for v in dA.values())})
            print(f"  {ct:11s} {name:7s}  r_obs={r_obs:+.2f}   null={null.mean():+.2f}+/-{null.std():.2f}   p={p:.3f}")
    R = pd.DataFrame(rows)
    R.to_csv(MIL / "stageB_perm.csv", index=False)
    print(f"\nsaved: {MIL/'stageB_perm.csv'}")


if __name__ == "__main__":
    main()
