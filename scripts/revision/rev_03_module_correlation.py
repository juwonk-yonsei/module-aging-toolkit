"""Module redundancy (Section 3.3, Supplementary Fig. S1, Table S6): gene overlap, score correlation in the rodent training
data and in the IPF pseudobulk data, and the effective number of independent modules.

Outputs: rev03_gene_overlap.csv, rev03_corr_training.csv, rev03_corr_ipf.csv,
         rev03_meff.csv, figs/rev03_module_correlation.png
"""
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

import rev_common as rc


def meff(C):
    """Effective number of tests: Nyholt (2004) and Li & Ji (2005)."""
    lam = np.clip(np.linalg.eigvalsh(np.asarray(C, float)), 0, None)
    M = len(lam)
    nyholt = 1 + (M - 1) * (1 - np.var(lam, ddof=1) / M)
    liji = float(np.sum((lam >= 1) + (lam - np.floor(lam))))
    return nyholt, liji


def main():
    clocks = rc.load_module_clocks()
    mem = pd.read_csv(rc.DATA / "supp" / "module_membership_rodent.csv")
    dup = mem.entrez.duplicated(keep=False).sum()
    comp, cg = rc.load_composite()
    *_, coef, _ = rc.linear_parts(comp)
    nz = {g for g, c in zip(cg, coef) if c != 0}
    rows = []
    for lab, d in clocks.items():
        genes = [str(g) for g in d["genes"]]
        *_, mc, _ = rc.linear_parts(d["pipeline"])
        rows.append(dict(module=lab, n_genes=len(genes), n_nonzero_module_coef=int((mc != 0).sum()),
                         in_composite_features=sum(g in set(cg) for g in genes),
                         nonzero_in_composite=sum(g in nz for g in genes)))
    O = pd.DataFrame(rows)
    O.attrs["duplicated_genes"] = int(dup)
    rc.save(O, "rev03_gene_overlap.csv")
    print(f"genes assigned to >1 module: {dup}; composite non-zero genes: {len(nz)}; "
          f"of which in modules: {sum(O.nonzero_in_composite)}")

    # training data: in-sample module predictions on all rodent samples
    from rev_02_species_cv import load_training
    ann, expr = load_training()
    E = expr.T
    E.columns = E.columns.astype(str)
    P = pd.DataFrame({lab: d["pipeline"].predict(E.reindex(columns=[str(g) for g in d["genes"]]).values)
                      for lab, d in clocks.items()}, index=E.index)
    Ct = P.corr(method="spearman")
    rc.save(Ct, "rev03_corr_training.csv", index=True)

    # IPF: per-sample scores, centred within cell type
    S = rc.adams_scores()
    mods = rc.module_names(S)
    Sc = S[mods].sub(S.groupby("celltype")[mods].transform("mean"))
    Ci = Sc.corr(method="spearman")
    rc.save(Ci, "rev03_corr_ipf.csv", index=True)

    rows = []
    for name, C in [("rodent training (in-sample predictions)", Ct), ("IPF pseudobulk (centred within cell type)", Ci)]:
        ny, lj = meff(C.values)
        off = C.values[np.triu_indices(len(C), 1)]
        rows.append(dict(dataset=name, n_modules=len(C), meff_nyholt=ny, meff_liji=lj,
                         median_abs_r=np.median(np.abs(off)), max_abs_r=np.max(np.abs(off)),
                         n_pairs_abs_r_gt_0_7=int((np.abs(off) > 0.7).sum())))
    M = pd.DataFrame(rows)
    rc.save(M, "rev03_meff.csv")
    print(M.round(2).to_string())

    plt = rc.setup_mpl()
    order = leaves_list(linkage(squareform(1 - np.abs(Ct.values), checks=False), "average"))
    labs = Ct.index[order]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), constrained_layout=True)
    for ax, C, title in [(axes[0], Ct, "Rodent training data"), (axes[1], Ci, "IPF pseudobulk (within cell type)")]:
        im = ax.imshow(C.loc[labs, labs].values, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(labs)))
        ax.set_yticks(range(len(labs)))
        ax.set_xticklabels(labs, rotation=90, fontsize=5)
        ax.set_yticklabels(labs if ax is axes[0] else [], fontsize=5)
        ax.set_title(title)
    fig.colorbar(im, ax=axes, shrink=0.6, label="Spearman ρ")
    fig.savefig(rc.FIG / "rev03_module_correlation.png")


if __name__ == "__main__":
    main()
