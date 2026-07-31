"""Extract GTEx lung bulk RNA-seq raw counts + donor age/sex from the full GCT.

Streams the ~900MB gzipped gene_reads GCT, keeping only Lung RNASEQ sample columns,
and writes a compact genes(Ensembl, unversioned) x samples count matrix + metadata.
"""
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parents[1] / "data" / "gtex"
GCT = BASE / "gene_reads.gct.gz"
SAMP = BASE / "SampleAttributes.txt"
PHEN = BASE / "SubjectPhenotypes.txt"

# 1. lung RNASEQ sample IDs
sa = pd.read_csv(SAMP, sep="\t", dtype=str)
lung = set(sa.loc[(sa.SMTSD == "Lung") & (sa.SMAFRZE == "RNASEQ"), "SAMPID"])
print(f"lung RNASEQ samples in attributes: {len(lung)}")

# 2. header of GCT (line 3) to find which columns to keep
with pd.io.common.get_handle(GCT, "r", compression="gzip") as h:
    f = h.handle
    f.readline(); f.readline()               # skip #1.2 and dims
    header = f.readline().rstrip("\n").split("\t")
keep_cols = ["Name"] + [c for c in header if c in lung]
print(f"columns kept (Name + lung): {len(keep_cols)-1} lung columns found in GCT")

# 3. stream-read only needed columns
df = pd.read_csv(GCT, sep="\t", skiprows=2, usecols=keep_cols, dtype={"Name": str})
df["Name"] = df["Name"].str.split(".").str[0]          # strip Ensembl version
df = df.groupby("Name").sum()                           # collapse dup unversioned
print("lung count matrix:", df.shape)

# 4. metadata: map SAMPID -> SUBJID -> AGE/SEX
ph = pd.read_csv(PHEN, sep="\t", dtype=str).set_index("SUBJID")
mid = {"20-29": 25, "30-39": 35, "40-49": 45, "50-59": 55, "60-69": 65, "70-79": 75}
meta = pd.DataFrame({"SAMPID": df.columns})
meta["SUBJID"] = meta.SAMPID.str.split("-").str[:2].str.join("-")
meta["AGE_bracket"] = meta.SUBJID.map(ph["AGE"])
meta["age_mid"] = meta.AGE_bracket.map(mid)
meta["SEX"] = meta.SUBJID.map(ph["SEX"])               # 1=male, 2=female
meta["sexM"] = (meta.SEX == "1").astype(int)
meta = meta.set_index("SAMPID")
print(meta.AGE_bracket.value_counts().sort_index())

df.to_pickle(BASE / "lung_counts.pkl")
meta.to_csv(BASE / "lung_meta.csv")
print("saved lung_counts.pkl and lung_meta.csv")
