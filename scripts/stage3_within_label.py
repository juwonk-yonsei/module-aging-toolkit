"""Within-label control for the attention-MIL localization (Fig. 6 of the original
submission; Supplementary Fig. S5d of the revision): does attention separate the
mono-derived/resident axis WITHIN a single fine annotation, or only by recovering the
pre-existing Macrophage vs Macrophage_Alveolar annotation boundary?

Pooled separation (all macrophage-coarse instances) vs within-label separation
(computed inside celltype=='Macrophage' and inside 'Macrophage_Alveolar'),
averaged over seeds. Also reports how well attention alone recovers the
annotation label (AUROC). Annotation-recovery signature: pooled >> within-label,
and AUROC >> 0.5.

Outputs: results/mil/stage3_within_label.csv (+ _raw.csv)
"""
import os
os.environ.setdefault("CUDA_DEVICE_ORDER", "PCI_BUS_ID")
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "2")
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

REPO = Path(os.environ.get("MAT_ROOT", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "third_party"))
import mil_model as M
from mil_score import MODULES

MIL = REPO / "results" / "mil"
MONO_DERIVED = {"SPP1": "ENSG00000118785", "FCN1": "ENSG00000085265",
                "VCAN": "ENSG00000038427", "CTHRC1": "ENSG00000164932",
                "INHBA": "ENSG00000122641", "MERTK": "ENSG00000153208"}
RESIDENT = {"FABP4": "ENSG00000170323", "MARCO": "ENSG00000019169",
            "MRC1": "ENSG00000260314"}
FINE = ["Macrophage", "Macrophage_Alveolar"]   # the two pooled annotations
N_SEED = 10
MIN_N = 30


def load():
    sc = pd.read_csv(MIL / "inst_scores_adams_bs60.csv", index_col=0)
    sc = sc[sc.celltype.isin(FINE) & sc.Disease_Identity.isin(["IPF", "Control"])].copy()
    sc.index = sc.index.astype(str)
    counts = pd.read_pickle(MIL / "inst_counts_adams_bs60.pkl")
    counts.columns = counts.columns.map(str)
    counts = counts[sc.index]
    lib = counts.sum(0).values + 1e-9
    mk = {}
    for name, ens in {**MONO_DERIVED, **RESIDENT}.items():
        if ens in counts.index:
            mk[name] = np.log1p(counts.loc[ens].values / lib * 1e6)
    return sc, pd.DataFrame(mk, index=sc.index)


def oof_attention(sc, seed):
    donors = sc.Subject_Identity.values
    uniq = pd.Index(sorted(sc.Subject_Identity.unique()))
    dlabel = sc.groupby("Subject_Identity").Disease_Identity.first().reindex(uniq)
    y = (dlabel == "IPF").astype(float).values
    X = sc[MODULES].values.astype(np.float32)
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    rows_of = {d: np.where(donors == d)[0] for d in uniq}
    w = np.full(len(sc), np.nan)
    k = min(5, int(y.sum()), int((1 - y).sum()))
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    for tr, va in skf.split(np.arange(len(uniq)), y):
        bags_tr = [X[rows_of[d]] for d in uniq[tr]]
        model, _ = M.train_eval(bags_tr, y[tr], bags_tr, y[tr], in_dim=len(MODULES),
                                pool="attention", seed=seed, max_epochs=200, patience=25)
        for d in uniq[va]:
            r = rows_of[d]
            w[r] = M.bag_attention(model, X[r]) * len(r)
    return w


def sep(att, mk, m, mono, res):
    if m.sum() < MIN_N:
        return np.nan
    rm = [spearmanr(att[m], mk[c].values[m])[0] for c in mono if c in mk]
    rr = [spearmanr(att[m], mk[c].values[m])[0] for c in res if c in mk]
    return np.nanmean(rm) - np.nanmean(rr)


def main():
    sc, mk = load()
    print(f"macrophage-coarse instances={len(sc)}  markers={list(mk.columns)}")
    print("celltype counts:\n", sc.celltype.value_counts().to_string())
    mono, res = list(MONO_DERIVED), list(RESIDENT)
    is_mac = (sc.celltype.values == "Macrophage")   # 1 = mono-derived annotation
    rows = []
    for seed in range(N_SEED):
        att = oof_attention(sc, seed)
        for grp in ["IPF", "Control"]:
            g = sc.Disease_Identity.values == grp
            rec = {"seed": seed, "group": grp, "n": int(g.sum())}
            rec["sep_pooled"] = sep(att, mk, g, mono, res)
            lab = is_mac[g].astype(int)
            if g.sum() >= MIN_N and 0 < lab.sum() < lab.size:
                a = roc_auc_score(lab, att[g])
                rec["auroc_annotation"] = max(a, 1 - a)
            else:
                rec["auroc_annotation"] = np.nan
            for label in FINE:
                m = g & (sc.celltype.values == label)
                rec[f"sep_within_{label}"] = sep(att, mk, m, mono, res)
                rec[f"n_{label}"] = int(m.sum())
            rows.append(rec)
        print(f"seed {seed} done")
    R = pd.DataFrame(rows)
    R.to_csv(MIL / "stage3_within_label_raw.csv", index=False)
    num = R.drop(columns=["seed"]).groupby("group").agg(["mean", "std"])
    num.to_csv(MIL / "stage3_within_label.csv")
    print("\n===== summary (mean, std over seeds) =====")
    print(num.round(3).to_string())
    print(f"\nsaved: {MIL/'stage3_within_label.csv'} (+ _raw.csv)")


if __name__ == "__main__":
    main()
