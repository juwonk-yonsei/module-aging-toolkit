"""Validate linear module clocks on the Klotho-KO example (expect KO pro-mortality shifts)."""
import os
import sys
from pathlib import Path
import numpy as np, pandas as pd
pd.set_option("display.width", 200); pd.set_option("display.max_columns", 40)
REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import tage_prep as tp
import module_clocks as mc

EXT = Path(os.environ.get("TAGE_DIR", str(REPO / "third_party" / "tAge"))) / "inst" / "extdata"
counts = pd.read_csv(EXT / "Exprs_example.csv", index_col=0)
meta = pd.read_csv(EXT / "Metadata_example.csv", index_col=0)

for tissue in meta["Tissue"].unique():
    m_t = meta[meta["Tissue"] == tissue]
    c_t = counts[m_t.index]
    pp = tp.preprocess(c_t, m_t, species="mouse", gene_mapping_type="Ensembl",
                       control_group_column="Genotype", control_group_label="WT")
    S = mc.module_scores(pp["scaled_diff"], outcome="Mortality")
    S["Genotype"] = m_t["Genotype"].values
    print(f"\n===== {tissue}: mean module score (Mortality) KO vs WT =====")
    g = S.groupby("Genotype").mean(numeric_only=True).T
    g["KO_minus_WT"] = g.get("Klotho KO", 0) - g.get("WT", 0)
    print(g.sort_values("KO_minus_WT", ascending=False).round(3).to_string())
