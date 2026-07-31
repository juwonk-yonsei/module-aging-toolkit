"""Pull healthy-lung scRNA (one CELLxGENE dataset) and build donor x coarse-celltype
pseudobulks (raw counts) efficiently: ONE get_anndata for a stratified donor subset,
then sparse group-sum. For the module-transfer-vs-age test at single-cell modality."""
import os
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
import scipy.sparse as sp
import cellxgene_census

DATASET = "9f222629-9e39-47d0-b83f-e08d610c7479"
MAX_DONORS = 45
MIN_CELLS = 50
OUT = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1]))) / (
    "data/hlca")
OUT.mkdir(parents=True, exist_ok=True)

def log(*a): print(*a, flush=True)

def coarse(ct):
    s = str(ct).lower()
    if "natural killer" in s: return "NK"
    if "cd8" in s or "cd4" in s or "t cell" in s or "thymocyte" in s or "regulatory t" in s: return "T"
    if "b cell" in s or "plasma" in s: return "B"
    if "macrophage" in s: return "Macrophage"
    if "monocyte" in s: return "Monocyte"
    if "dendritic" in s: return "DC"
    if "mast" in s: return "Mast"
    if "ciliated" in s: return "Ciliated"
    if "type ii pneumocyte" in s or "type 2" in s: return "ATII"
    if "type i pneumocyte" in s or "type 1" in s: return "ATI"
    if "club" in s or "secretory" in s or "goblet" in s or "basal" in s: return "Airway_epi"
    if "endothelial" in s or "capillary" in s or "vein" in s or "lymphatic" in s: return "Endothelial"
    if "fibroblast" in s or "smooth muscle" in s or "pericyte" in s or "stromal" in s or "mesothelial" in s: return "Stromal"
    return None

def parse_age(s):
    m = re.search(r"(\d+)-year-old", str(s)); return int(m.group(1)) if m else None

base = (f"dataset_id == '{DATASET}' and tissue_general == 'lung' "
        "and disease == 'normal' and is_primary_data == True")

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    log("querying obs ...")
    obs = cellxgene_census.get_obs(census, "homo_sapiens", value_filter=base,
                                   column_names=["donor_id", "development_stage", "sex"])
    obs["age"] = obs["development_stage"].map(parse_age)
    obs = obs.dropna(subset=["age"])
    dinfo = obs[obs.age >= 20].drop_duplicates("donor_id").sort_values("age")
    if len(dinfo) > MAX_DONORS:                       # stratify across age
        idx = np.linspace(0, len(dinfo) - 1, MAX_DONORS).round().astype(int)
        dinfo = dinfo.iloc[np.unique(idx)]
    donors = dinfo.donor_id.tolist()
    log(f"selected {len(donors)} adult donors, ages {dinfo.age.min()}-{dinfo.age.max()}")

    dlist = ", ".join(f"'{d}'" for d in donors)
    log("pulling anndata (single call) ...")
    ad = cellxgene_census.get_anndata(
        census, "homo_sapiens", measurement_name="RNA", X_name="raw",
        obs_value_filter=f"{base} and donor_id in [{dlist}]",
        column_names={"obs": ["donor_id", "development_stage", "sex", "cell_type"],
                      "var": ["feature_id"]},
    )
    log(f"pulled {ad.n_obs} cells x {ad.n_vars} genes")

ad.obs["coarse"] = ad.obs["cell_type"].map(coarse)
ad.obs["age"] = ad.obs["development_stage"].map(parse_age)
keep = ad.obs["coarse"].notna().values
ad = ad[keep]
grp = (ad.obs["donor_id"].astype(str) + "||" + ad.obs["coarse"].astype(str)).values
uniq, inv = np.unique(grp, return_inverse=True)
log(f"{len(uniq)} donor x celltype groups")

# sparse group-sum: (ngroups x ncells) one-hot @ (ncells x ngenes)
onehot = sp.csr_matrix((np.ones(len(inv)), (inv, np.arange(len(inv)))),
                       shape=(len(uniq), ad.n_obs))
S = onehot @ ad.X                                    # ngroups x ngenes
S = np.asarray(S.todense())

genes = pd.Index(ad.var["feature_id"].values).str.split(".").str[0]
pbdf = pd.DataFrame(S.T, index=genes, columns=uniq).groupby(level=0).sum()

# metadata
ncell = pd.Series(np.bincount(inv), index=uniq)
meta = pd.DataFrame({"key": uniq})
meta["donor"] = [k.split("||")[0] for k in uniq]
meta["celltype"] = [k.split("||")[1] for k in uniq]
meta["n_cells"] = ncell.values
dmap_age = dict(zip(ad.obs["donor_id"].astype(str), ad.obs["age"]))
dmap_sex = dict(zip(ad.obs["donor_id"].astype(str), ad.obs["sex"].astype(str)))
meta["age"] = meta.donor.map(dmap_age)
meta["sex"] = meta.donor.map(dmap_sex)
meta = meta[meta.n_cells >= MIN_CELLS]
pbdf = pbdf[meta.key.tolist()]

pbdf.to_pickle(OUT / "pseudobulk_counts.pkl")
meta.to_csv(OUT / "pseudobulk_meta.csv", index=False)
log(f"saved: {pbdf.shape} genes x pseudobulks; donors={meta.donor.nunique()}")
log(meta.celltype.value_counts().to_string())
