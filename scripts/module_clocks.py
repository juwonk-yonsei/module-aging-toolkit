"""Linear module clocks from published multi-species coefficients (Supp Table 5C).

Rationale: module-specific clocks are NOT released as .pkl. Supp Table 5C provides
their elastic-net coefficients (mouse Entrez) per module for chrono & mortality.
For GROUP CONTRASTS (disease vs control) the intercept and any internal centering
cancel, so a relative module score  s_m = sum_{g in module m} coef_{g} * x_g
(on the scaled_diff input) is a valid relative module tAge for comparisons.

We validate this on the Klotho example before trusting it (KO should show
pro-mortality module shifts, cf. paper Fig.3h).
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
SUPP = REPO / "data" / "supp"

_MOD_DICT = {
    "turquoise": "Inflammation", "pink": "Mito translation/OxPhos", "orange": "Chromatin",
    "green": "Cell cycle/DNA repl", "blue": "Muscle/Cytoskel/Glycolysis",
    "darkgreen": "Adaptive imm/T cell", "darkred": "mRNA splicing", "brown4": "ECM/EMT",
    "white": "OxPhos/Heme", "darkmagenta": "Interferon", "sienna3": "VEGF",
    "darkslateblue": "ER/UPR", "plum1": "Protein folding", "ivory": "Fatty acid/Peroxisome",
}


def load_module_coefs(outcome="Mortality"):
    """Return dict module -> Series(index=entrez str, value=coef), excluding 'All module genes'."""
    df = pd.read_csv(SUPP / "module_coef_multispecies.csv")
    df = df[(df.outcome == outcome) & (df.module != "All module genes")]
    df["entrez"] = df["entrez"].astype(str)
    return {m: g.set_index("entrez")["coef"] for m, g in df.groupby("module")}


def module_scores(scaled_diff_df, outcome="Mortality"):
    """scaled_diff_df: samples x gene_list (mouse Entrez, str cols). Returns samples x modules."""
    coefs = load_module_coefs(outcome)
    X = scaled_diff_df.copy()
    X.columns = X.columns.map(str)
    X = X.fillna(0.0)  # centered input: missing => at control median => 0 contribution
    out = {}
    for mod, cser in coefs.items():
        genes = [g for g in cser.index if g in X.columns]
        if not genes:
            continue
        out[mod] = X[genes].values @ cser.loc[genes].values
    S = pd.DataFrame(out, index=scaled_diff_df.index)
    S = S.rename(columns={m: f"{m} ({_MOD_DICT.get(m, m)})" for m in S.columns})
    return S
