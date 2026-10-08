"""Confounding and technical sensitivity analyses (Sections 2.9 and 3.5, Supplementary Table S12).

For every analysis the full 11-cell-type x 23-module grid is re-estimated
(score ~ IPF + age + sex [+ extra covariates], IPF and Control donors) and compared with
the main analysis; the 30 main FDR<0.10 effects are tracked explicitly.

  age_overlap      Controls restricted to the IPF age range (>= 54 y)
  tech_covariates  + log cell number + log library size + detected-gene fraction
  depth_group_matched  the deeper disease group binomially thinned so that group-median
                   library sizes match (10 seeds; median beta, fraction of seeds p < 0.05)
  cells_group_matched  the disease group with more cells per donor subsampled (60-cell
                   mini-pseudobulks) so that group-median cell numbers match (20 seeds)
  depth_equalized / cells_equalized   stress tests: every pseudobulk thinned to the
                   cell-type minimum depth / one mini-pseudobulk per donor, each with an
                   information-loss-matched control (*_loss_control) that removes the same
                   average information but preserves the group differences
  min_cells_30 / min_cells_100   alternative inclusion thresholds
  cpm              library-size scaling instead of RLE

Outputs: rev09_sensitivity_long.csv, rev09_sensitivity_summary.csv, rev09_hits_sensitivity.csv
"""
import numpy as np
import pandas as pd
from scipy import stats

import rev_common as rc

N_THIN = 10
N_INST = 20


def grid_effects(S, extra=()):
    rows = []
    mods = [m for m in rc.module_names(S) if m not in ("detected_frac", "log_cells", "log_lib")]
    for ct, g in S.groupby("celltype"):
        g = g[g.Disease_Identity.isin(["IPF", "Control"])].dropna(subset=["age", "sexM", *extra])
        nI, nC = (g.Disease_Identity == "IPF").sum(), (g.Disease_Identity == "Control").sum()
        if nI < 3 or nC < 3:
            continue
        X = np.column_stack([np.ones(len(g)), (g.Disease_Identity == "IPF").astype(float), g.age, g.sexM]
                            + [g[e].values for e in extra])
        Y = g[mods].values.astype(float)
        XtX_inv = np.linalg.pinv(X.T @ X)
        Bt = XtX_inv @ X.T @ Y
        E = Y - X @ Bt
        dof = len(g) - X.shape[1]
        if dof < 2:
            continue
        se = np.sqrt(XtX_inv[1, 1] * (E ** 2).sum(0) / dof)
        p = 2 * stats.t.sf(np.abs(Bt[1] / se), dof)
        for k, m in enumerate(mods):
            rows.append(dict(celltype=ct, module=m, beta=Bt[1, k], se=se[k], p=p[k], n_IPF=nI, n_Control=nC))
    R = pd.DataFrame(rows)
    R["fdr"] = rc.bh(R.p)
    return R


def score_counts(pb, meta, ct_list, clocks, comp, normalization="RLE"):
    recs = []
    for ct in ct_list:
        sub = meta[meta.Manuscript_Identity == ct]
        n = sub.Disease_Identity.value_counts()
        if n.get("IPF", 0) < 3 or n.get("Control", 0) < 3:
            continue
        X, md = rc.preprocess_ct(pb, meta, ct, normalization=normalization)
        S = rc.score_modules(X, clocks, comp).join(
            md[["Subject_Identity", "Disease_Identity", "age", "sexM", "n_cells", "lib_size"]])
        S["celltype"] = ct
        recs.append(S)
    return pd.concat(recs)


def thin(pb, meta, ct_list, rng, matched_loss=False):
    """Binomial thinning to the cell-type minimum library size. matched_loss=True thins
    every pseudobulk by the same fraction (the median of the equalizing fractions), which
    removes as much information on average but preserves the depth differences."""
    cols = meta[meta.Manuscript_Identity.isin(ct_list)].group.tolist()
    out = pb[cols].copy()
    for ct in ct_list:
        g = meta[meta.Manuscript_Identity == ct].group.tolist()
        lib = pb[g].sum(0)
        frac = lib.min() / lib
        for c in g:
            f = float(np.median(frac)) if matched_loss else float(frac[c])
            out[c] = rng.binomial(pb[c].values.astype(np.int64), f)
    return out


def thin_group_matched(pb, meta, ct_list, rng):
    """Remove the systematic IPF-vs-Control depth difference: thin every pseudobulk of the
    deeper group by the ratio of group-median library sizes (the other group untouched)."""
    cols = meta[meta.Manuscript_Identity.isin(ct_list)].group.tolist()
    out = pb[cols].copy()
    for ct in ct_list:
        m = meta[(meta.Manuscript_Identity == ct)]
        lib = pb[m.group].sum(0)
        med_i = lib[m[m.Disease_Identity == "IPF"].group].median()
        med_c = lib[m[m.Disease_Identity == "Control"].group].median()
        deeper, r = ("IPF", med_c / med_i) if med_i > med_c else ("Control", med_i / med_c)
        for c in m[m.Disease_Identity == deeper].group:
            out[c] = rng.binomial(pb[c].values.astype(np.int64), r)
    return out


def pick_instances_group_matched(im, rng):
    """Remove the systematic IPF-vs-Control cell-number difference: in the group with more
    cells per donor, keep the fraction of each donor's mini-pseudobulks equal to the ratio
    of group-median cell numbers."""
    picks = []
    for ct, g in im.groupby("celltype"):
        per = g.groupby(["Subject_Identity", "Disease_Identity"]).n_cells.sum().reset_index()
        med_i = per[per.Disease_Identity == "IPF"].n_cells.median()
        med_c = per[per.Disease_Identity == "Control"].n_cells.median()
        deeper, r = ("IPF", med_c / med_i) if med_i > med_c else ("Control", med_i / med_c)
        for (sid, dis), gg in g.groupby(["Subject_Identity", "Disease_Identity"]):
            k = max(1, int(round(r * len(gg)))) if dis == deeper else len(gg)
            picks.append(gg.sample(k, random_state=int(rng.integers(1e9))))
    return pd.concat(picks)


def pick_instances(im, rng, matched_loss=False):
    """One mini-pseudobulk per donor x cell type; matched_loss=True instead keeps the same
    fraction of every donor's mini-pseudobulks (median of 1/n), preserving cell-number
    differences at a similar average information loss."""
    picks = []
    for ct, g in im.groupby("celltype"):
        n_inst = g.groupby("Subject_Identity").size()
        f = float(np.median(1.0 / n_inst))
        for sid, gg in g.groupby("Subject_Identity"):
            k = max(1, int(round(f * len(gg)))) if matched_loss else 1
            picks.append(gg.sample(k, random_state=int(rng.integers(1e9))))
    return pd.concat(picks)


def main():
    rng = np.random.default_rng(rc.SEED)
    clocks, comp = rc.load_module_clocks(), rc.load_composite()
    pb, meta = rc.load_adams()
    cts = rc.testable_celltypes(meta)
    main_S = rc.adams_scores(pb=pb, meta=meta)
    det = (pb[main_S.index] > 0).mean(0)
    main_S["detected_frac"] = det.reindex(main_S.index).values
    main_S["log_cells"] = np.log(main_S.n_cells)
    main_S["log_lib"] = np.log(main_S.lib_size)
    M = grid_effects(main_S)
    hits = M[M.fdr < 0.10][["celltype", "module", "beta"]].rename(columns={"beta": "beta_main"})
    print(f"main analysis: {len(M)} tests, {len(hits)} FDR<0.10 hits")
    results = {"main": M}

    # age overlap
    S = main_S[(main_S.Disease_Identity != "Control") | (main_S.age >= 54)]
    results["age_overlap"] = grid_effects(S)
    # technical covariates
    results["tech_covariates"] = grid_effects(main_S, extra=("log_cells", "log_lib", "detected_frac"))
    # thresholds and normalization
    _, meta30 = rc.load_adams(min_cells=30)
    _, meta100 = rc.load_adams(min_cells=100)
    results["min_cells_30"] = grid_effects(score_counts(pb, meta30, rc.testable_celltypes(meta30), clocks, comp))
    results["min_cells_100"] = grid_effects(score_counts(pb, meta100, rc.testable_celltypes(meta100), clocks, comp))
    results["cpm"] = grid_effects(score_counts(pb, meta, cts, clocks, comp, normalization="CPM"))

    # depth equalization (binomial thinning) and its information-loss-matched control
    for name, fn in [("depth_group_matched", lambda: thin_group_matched(pb, meta, cts, rng)),
                     ("depth_equalized", lambda: thin(pb, meta, cts, rng)),
                     ("depth_loss_control", lambda: thin(pb, meta, cts, rng, matched_loss=True))]:
        reps = []
        for s in range(N_THIN):
            reps.append(grid_effects(score_counts(fn(), meta, cts, clocks, comp)).assign(rep=s))
        results[name] = pd.concat(reps)

    # cell-number equalization (one ~60-cell mini-pseudobulk per donor x cell type)
    inst = pd.read_pickle(rc.RES / "mil" / "inst_counts_adams_bs60.pkl")
    im = pd.read_csv(rc.RES / "mil" / "inst_meta_adams_bs60.csv")
    im["instance"] = im["instance"].astype(str)
    cov = pd.read_csv(rc.DATA / "ipf" / "subject_covariates.csv").set_index("Subject_Identity")
    im = im.merge(cov[["age", "Sex"]], left_on="Subject_Identity", right_index=True, how="left")
    im["sexM"] = (im.Sex == "M").astype(float)
    im = im[im.celltype.isin(cts)]
    for name, picker in [("cells_group_matched", lambda: pick_instances_group_matched(im, rng)),
                         ("cells_equalized", lambda: pick_instances(im, rng)),
                         ("cells_loss_control", lambda: pick_instances(im, rng, matched_loss=True))]:
        reps = []
        for s in range(N_INST):
            pick = picker()
            cnt = inst[pick.instance].T.groupby(
                (pick.celltype + "||" + pick.Subject_Identity).values).sum().T
            mm = (pick.assign(group=pick.celltype + "||" + pick.Subject_Identity)
                  .groupby("group").agg(Manuscript_Identity=("celltype", "first"),
                                        Subject_Identity=("Subject_Identity", "first"),
                                        Disease_Identity=("Disease_Identity", "first"),
                                        age=("age", "first"), sexM=("sexM", "first"),
                                        n_cells=("n_cells", "sum")).reset_index())
            mm["lib_size"] = cnt[mm.group].sum(0).values
            reps.append(grid_effects(score_counts(cnt, mm, cts, clocks, comp)).assign(rep=s))
        results[name] = pd.concat(reps)

    # ---- assemble
    long = pd.concat([r.assign(analysis=k) for k, r in results.items()], ignore_index=True)
    rc.save(long, "rev09_sensitivity_long.csv")
    rows, hit_rows = [], []
    for k, R in results.items():
        if "rep" in R:
            agg = R.groupby(["celltype", "module"]).agg(beta=("beta", "median"), p=("p", "median"),
                                                        frac_p05=("p", lambda x: np.mean(x < .05))).reset_index()
            agg["fdr"] = rc.bh(agg.p)
        else:
            agg = R.assign(frac_p05=(R.p < .05).astype(float))
        j = M.merge(agg, on=["celltype", "module"], suffixes=("_main", ""))
        h = hits.merge(agg, on=["celltype", "module"], how="left")
        rows.append(dict(analysis=k, n_celltypes=agg.celltype.nunique(), n_tests=len(agg),
                         fdr10_hits=int((agg.fdr < .1).sum()),
                         beta_spearman_vs_main=stats.spearmanr(j.beta_main, j.beta)[0],
                         main_hits_available=int(h.beta.notna().sum()),
                         main_hits_same_sign=int((np.sign(h.beta) == np.sign(h.beta_main)).sum()),
                         main_hits_p05=int((h.p < .05).sum()),
                         main_hits_median_beta_ratio=float(np.nanmedian(h.beta / h.beta_main))))
        hit_rows.append(h.assign(analysis=k))
    Sm = pd.DataFrame(rows)
    rc.save(Sm, "rev09_sensitivity_summary.csv")
    rc.save(pd.concat(hit_rows), "rev09_hits_sensitivity.csv")
    print(Sm.round(3).to_string())


if __name__ == "__main__":
    main()
