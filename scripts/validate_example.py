"""Validate the Python tAge preprocessing port on the bundled example data.

Example = Klotho-KO vs WT (Kidney, Skeletal muscle) mouse Ensembl counts.
Expectation (paper Fig.3f): Klotho-KO shows increased MORTALITY tAge vs WT,
while chronological tAge should be similar (age-matched animals).
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "python"))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))

import tage_prep as tp
from tage_predict import predict_tAge

EXT = Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "extdata"
MODELS = REPO / "data" / "clock_models"

counts = pd.read_csv(EXT / "Exprs_example.csv", index_col=0)
meta = pd.read_csv(EXT / "Metadata_example.csv", index_col=0)
print(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples")
print("meta cols:", list(meta.columns))
print(meta.groupby(["Tissue", "Genotype"]).size(), "\n")

clocks = {
    "Chronoage_mouse": MODELS / "EN_Chronoage_Mouse_Multitissue_scaleddiff.pkl",
    "Mortality_rodent": MODELS / "EN_Mortality_Rodents_Multitissue_scaleddiff.pkl",
}

rows = []
for tissue in meta["Tissue"].unique():
    m_t = meta[meta["Tissue"] == tissue]
    c_t = counts[m_t.index]
    pp = tp.preprocess(
        c_t, m_t, species="mouse", gene_mapping_type="Ensembl",
        control_group_column="Genotype", control_group_label="WT",
    )
    X = pp["scaled_diff"]  # samples x gene_list
    ann = m_t.copy()
    for name, path in clocks.items():
        res = predict_tAge(str(path), X, ann, species="mouse", prefix=f"{name}_")
        ann = res
    ann["Tissue"] = tissue
    rows.append(ann)

allann = pd.concat(rows)
cols = [c for c in allann.columns if c.endswith("tAge")]
print("=== per-sample tAge ===")
print(allann[["Tissue", "Genotype"] + cols].round(2).to_string())

print("\n=== group means (by Tissue x Genotype) ===")
summ = allann.groupby(["Tissue", "Genotype"])[cols].mean().round(2)
print(summ.to_string())

print("\n=== KO - WT (mortality) per tissue ===")
for tissue in allann["Tissue"].unique():
    sub = allann[allann["Tissue"] == tissue]
    mcol = [c for c in cols if c.startswith("Mortality")][0]
    ko = sub[sub["Genotype"].str.contains("KO")][mcol].mean()
    wt = sub[sub["Genotype"] == "WT"][mcol].mean()
    flag = "OK (KO>WT)" if ko > wt else "UNEXPECTED"
    print(f"  {tissue:16s} KO={ko:7.2f}  WT={wt:7.2f}  diff={ko-wt:+7.2f}  {flag}")
