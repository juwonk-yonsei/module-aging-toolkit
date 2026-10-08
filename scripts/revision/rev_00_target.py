"""Construction and scale of the mortality training target (Supplementary Methods S1).

Difference.Hazard.log10 in the released rodent annotation is checked against
    Expected_Hazard.log10(sample) - Expected_Hazard.log10(reference group),
where the reference group of a sample is the set of samples of the same dataset (Source), tissue and sex whose
target is 0 and whose age difference is 0. Reported: fraction reproduced exactly, the association with the
chronological age difference among control animals (overall, per month in mice and per 0.1 of maximum lifespan),
and the target of lifespan-modifying interventions compared at the reference age.

Outputs: rev00_target_construction.csv, rev00_target_by_intervention.csv
"""
import numpy as np
import pandas as pd

import rev_common as rc

TOL = 1e-6


def main():
    a = pd.read_excel(rc.DATA / "rodent" / "Data_annotation_relative_rodents.xlsx")
    tgt, exp_h = "Difference.Hazard.log10", "Expected_Hazard.log10"
    d_age, d_norm = "Difference.Chronological_age.days", "Difference.Chronological_age.Normalized_by_species_max_lifespan"
    key = ["Source", "Tissue", "Sex"]

    ref = a[(a[tgt].abs() < TOL) & (a[d_age].abs() < TOL)]
    ref_h = ref.groupby(key)[exp_h].agg(["mean", "std", "size"]).rename(columns={"mean": "ref_h"})
    m = a.merge(ref_h[["ref_h"]], left_on=key, right_index=True, how="left")
    m["recon"] = m[exp_h] - m["ref_h"]
    has_ref = m["ref_h"].notna()
    exact = has_ref & ((m["recon"] - m[tgt]).abs() < 1e-4)

    ctrl = m[m["Intervention.type"] == "Control"].dropna(subset=[tgt, d_age, d_norm])
    r_ctrl = np.corrcoef(ctrl[tgt], ctrl[d_age])[0, 1]
    mouse = ctrl[ctrl.Species == "Mouse"]
    slope_month = np.polyfit(mouse[d_age] / 30.4375, mouse[tgt], 1)[0]
    slope_norm = np.polyfit(ctrl[d_norm] / 0.1, ctrl[tgt], 1)[0]

    rows = [
        ("samples", len(a)), ("mouse samples", int((a.Species == "Mouse").sum())),
        ("rat samples", int((a.Species == "Rat").sum())), ("datasets (Source)", a.Source.nunique()),
        ("samples with a reference group (Source x Tissue x Sex)", int(has_ref.sum())),
        ("target reproduced exactly (|diff| < 1e-4)", int(exact.sum())),
        ("fraction reproduced exactly", float(exact.mean())),
        ("fraction reproduced exactly among samples with a reference group", float(exact[has_ref].mean())),
        ("control samples with target and age difference", len(ctrl)),
        ("Pearson r, target vs age difference, controls", float(r_ctrl)),
        ("slope, log10 hazard per month of age difference, mouse controls", float(slope_month)),
        ("slope, log10 hazard per 0.1 of species maximum lifespan, controls", float(slope_norm)),
        ("target range (min)", float(a[tgt].min())), ("target range (max)", float(a[tgt].max())),
        ("target median |value|", float(a[tgt].abs().median())),
    ]
    rc.save(pd.DataFrame(rows, columns=["quantity", "value"]), "rev00_target_construction.csv")

    same_age = m[(m[d_age].abs() < TOL) & (m["Intervention.type"] != "Control")]
    by_int = (same_age.groupby("Intervention.type")[tgt]
              .agg(n="size", mean="mean", median="median", sd="std").reset_index())
    by_int["hazard_ratio_at_mean"] = 10 ** by_int["mean"]
    rc.save(by_int, "rev00_target_by_intervention.csv")
    print(pd.DataFrame(rows, columns=["quantity", "value"]).to_string(index=False))
    print(by_int.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
