"""Formal composite-masking analysis (Sections 2.8 and 3.4, Supplementary Table S11).

The composite mortality clock is linear after imputation/standardization, so each sample's
composite score is exactly  intercept + sum_k c_k  where c_k is the summed contribution of
the composite's own features that belong to module k (23 modules) plus one remainder
component for features outside every module. For the IPF effect (covariate-adjusted OLS
beta, same design for every component) this gives  Delta_composite = sum_k Delta_k.

Masking statistics (all on the composite's raw scale, Delta log10 hazard):
  S   = sum_k |Delta_k|            total component-level disease signal
  N   = |Delta_composite|          net composite signal
  CM  = S - N                      signal cancelled inside the composite
  MI  = CM / S                     fraction cancelled (0 = none, 1 = complete)
Inference: Freedman-Lane permutation (B = 2000, same permutation for every component) for
S, N and MI; stratified donor bootstrap (B = 2000) for 95% CIs. Standardized effects
(beta / residual SD) are reported for the composite and for every module clock.

Outputs: rev07_masking_celltype.csv, rev07_masking_components.csv, rev07_masking_global.csv,
         figs/rev07_masking.png
"""
import numpy as np
import pandas as pd

import rev_common as rc

B = 2000


def composite_components(X, comp_pipe, comp_genes, membership):
    med, mu, sc, coef, b0 = rc.linear_parts(comp_pipe)
    Xc = X.reindex(columns=comp_genes).values
    Xc = np.where(np.isnan(Xc), med[None, :], Xc)
    contrib = ((Xc - mu[None, :]) / sc[None, :]) * coef[None, :]
    lab = np.array([membership.get(g, "Non-module genes") for g in comp_genes])
    comps = {k: contrib[:, lab == k].sum(1) for k in pd.unique(lab)}
    C = pd.DataFrame(comps, index=X.index)
    total = C.sum(1) + b0
    return C, total


def design(md):
    return np.column_stack([np.ones(len(md)), (md.Disease_Identity == "IPF").astype(float),
                            md.age.values, md.sexM.values])


def betas(Y, Xd):
    return (np.linalg.pinv(Xd.T @ Xd) @ Xd.T @ Y)[1]


def mask_stats(d):
    S = np.abs(d).sum(-1)
    N = np.abs(d.sum(-1))
    return S, N, (S - N) / np.where(S > 0, S, np.nan)


def analyse(C, md, rng):
    Y = C.values
    Xd = design(md)
    d_obs = betas(Y, Xd)
    S, N, MI = mask_stats(d_obs)
    # Freedman-Lane: residuals of the reduced (covariates-only) model, permuted jointly
    Xr = Xd[:, [0, 2, 3]]
    H = Xr @ np.linalg.pinv(Xr.T @ Xr) @ Xr.T
    fit, res = H @ Y, Y - H @ Y
    P = np.array([betas(fit + res[rng.permutation(len(Y))], Xd) for _ in range(B)])
    Sp, Np, MIp = mask_stats(P)
    # stratified donor bootstrap
    idx_i = np.where(md.Disease_Identity.values == "IPF")[0]
    idx_c = np.where(md.Disease_Identity.values != "IPF")[0]
    boots = []
    for _ in range(B):
        ii = np.concatenate([rng.choice(idx_i, len(idx_i)), rng.choice(idx_c, len(idx_c))])
        try:
            boots.append(betas(Y[ii], Xd[ii]))
        except np.linalg.LinAlgError:
            continue
    Bt = np.array(boots)
    Sb, Nb, MIb = mask_stats(Bt)
    resid = Y - Xd @ (np.linalg.pinv(Xd.T @ Xd) @ Xd.T @ Y)
    rsd = resid.std(0, ddof=Xd.shape[1])
    out = dict(S=S, N=N, MI=MI, delta_composite=d_obs.sum(),
               p_S=(1 + np.sum(Sp >= S)) / (B + 1), p_N=(1 + np.sum(Np >= N)) / (B + 1),
               p_MI=(1 + np.sum(MIp >= MI)) / (B + 1),
               S_null_median=np.median(Sp), MI_null_median=np.nanmedian(MIp),
               S_lo=np.percentile(Sb, 2.5), S_hi=np.percentile(Sb, 97.5),
               N_lo=np.percentile(Nb, 2.5), N_hi=np.percentile(Nb, 97.5),
               MI_lo=np.nanpercentile(MIb, 2.5), MI_hi=np.nanpercentile(MIb, 97.5),
               comp_lo=np.percentile(Bt.sum(1), 2.5), comp_hi=np.percentile(Bt.sum(1), 97.5))
    comp_rows = pd.DataFrame(dict(component=C.columns, delta=d_obs,
                                  lo=np.percentile(Bt, 2.5, 0), hi=np.percentile(Bt, 97.5, 0),
                                  p_perm=[(1 + np.sum(np.abs(P[:, k]) >= abs(d_obs[k]))) / (B + 1)
                                          for k in range(len(d_obs))],
                                  std_effect=d_obs / rsd))
    return out, comp_rows


def main():
    rng = np.random.default_rng(rc.SEED)
    pb, meta = rc.load_adams()
    comp_pipe, comp_genes = rc.load_composite()
    mem = pd.read_csv(rc.DATA / "supp" / "module_membership_rodent.csv")
    lab = rc.make_labels([f"{m}|{a}" for m, a in mem[["module", "annotation"]].drop_duplicates().values])
    membership = {str(e): lab[f"{m}|{a}"] for e, m, a in mem[["entrez", "module", "annotation"]].values}
    clocks = rc.load_module_clocks()

    ct_rows, comp_tabs, std_rows = [], [], []
    allC, allmd = [], []
    for ct in rc.testable_celltypes(meta):
        X, md = rc.preprocess_ct(pb, meta, ct, diseases=("IPF", "Control"))
        md = md.loc[X.index]
        C, total = composite_components(X, comp_pipe, comp_genes, membership)
        assert np.allclose(total, comp_pipe.predict(X.reindex(columns=comp_genes).values))
        out, rows = analyse(C, md, rng)
        top_pos = rows[rows.component != "Non-module genes"].nlargest(1, "delta").iloc[0]
        top_neg = rows[rows.component != "Non-module genes"].nsmallest(1, "delta").iloc[0]
        ct_rows.append(dict(celltype=ct, n_IPF=int((md.Disease_Identity == "IPF").sum()),
                            n_Control=int((md.Disease_Identity == "Control").sum()), **out,
                            delta_nonmodule=float(rows.loc[rows.component == "Non-module genes", "delta"].iloc[0]),
                            top_positive=f"{top_pos.component} ({top_pos.delta:+.3f})",
                            top_negative=f"{top_neg.component} ({top_neg.delta:+.3f})"))
        comp_tabs.append(rows.assign(celltype=ct))
        # standardized effects: composite vs separately trained module clocks
        S = rc.score_modules(X, clocks, (comp_pipe, comp_genes))
        Xd = design(md)
        bet = np.linalg.pinv(Xd.T @ Xd) @ Xd.T @ S.values
        rsd = (S.values - Xd @ bet).std(0, ddof=Xd.shape[1])
        for k, c in enumerate(S.columns):
            std_rows.append(dict(celltype=ct, score=c, beta=bet[1, k], std_effect=bet[1, k] / rsd[k]))
        allC.append(C.assign(celltype=ct))
        allmd.append(md.assign(celltype=ct))

    T = pd.DataFrame(ct_rows)
    T["fdr_S"] = rc.bh(T.p_S)
    T["fdr_N"] = rc.bh(T.p_N)
    rc.save(T, "rev07_masking_celltype.csv")
    Cm = pd.concat(comp_tabs)
    rc.save(Cm, "rev07_masking_components.csv")
    Sd = pd.DataFrame(std_rows)
    rc.save(Sd, "rev07_standardized_effects.csv")

    # ---- global: cell-type fixed effects, donor-level permutation / bootstrap
    C = pd.concat(allC).fillna(0.0)
    md = pd.concat(allmd)
    comps = [c for c in C.columns if c != "celltype"]
    cts = sorted(md.celltype.unique())
    D = pd.get_dummies(md.celltype)[cts[1:]].astype(float).values

    def gdesign(dis):
        return np.column_stack([np.ones(len(md)), dis, md.age.values, md.sexM.values, D])

    Y = C[comps].values
    dis = (md.Disease_Identity == "IPF").astype(float).values
    d_obs = betas(Y, gdesign(dis))
    S0, N0, MI0 = mask_stats(d_obs)
    donors = md.Subject_Identity.values
    ud = pd.Series(dis, index=donors).groupby(level=0).first()
    P = []
    for _ in range(B):
        perm = pd.Series(rng.permutation(ud.values), index=ud.index)
        P.append(betas(Y, gdesign(perm.reindex(donors).values)))
    P = np.array(P)
    Sp, Np, MIp = mask_stats(P)
    di, dc = ud.index[ud == 1].values, ud.index[ud == 0].values
    rows_by_donor = pd.Series(np.arange(len(md))).groupby(donors).apply(list)
    Bt = []
    for _ in range(B):
        pick = np.concatenate([rng.choice(di, len(di)), rng.choice(dc, len(dc))])
        ii = np.concatenate([rows_by_donor[p] for p in pick])
        Xb = gdesign(dis)[ii]
        Bt.append(betas(Y[ii], Xb))
    Bt = np.array(Bt)
    Sb, Nb, MIb = mask_stats(Bt)
    G = pd.DataFrame(dict(component=comps, delta=d_obs, lo=np.percentile(Bt, 2.5, 0),
                          hi=np.percentile(Bt, 97.5, 0),
                          p_perm=[(1 + np.sum(np.abs(P[:, k]) >= abs(d_obs[k]))) / (B + 1)
                                  for k in range(len(comps))]))
    G["fdr"] = rc.bh(G.p_perm)
    summary = dict(component="__summary__", S=S0, N=N0, MI=MI0, delta_composite=d_obs.sum(),
                   comp_lo=np.percentile(Bt.sum(1), 2.5), comp_hi=np.percentile(Bt.sum(1), 97.5),
                   p_S=(1 + np.sum(Sp >= S0)) / (B + 1), p_N=(1 + np.sum(Np >= N0)) / (B + 1),
                   p_MI=(1 + np.sum(MIp >= MI0)) / (B + 1), S_null_median=np.median(Sp),
                   MI_null_median=np.nanmedian(MIp), S_lo=np.percentile(Sb, 2.5),
                   S_hi=np.percentile(Sb, 97.5), MI_lo=np.nanpercentile(MIb, 2.5),
                   MI_hi=np.nanpercentile(MIb, 97.5))
    G = pd.concat([G, pd.DataFrame([summary])], ignore_index=True)
    rc.save(G, "rev07_masking_global.csv")
    print(T[["celltype", "delta_composite", "comp_lo", "comp_hi", "p_N", "S", "S_lo", "S_hi", "p_S",
             "MI", "MI_lo", "MI_hi", "MI_null_median", "p_MI", "top_positive", "top_negative"]].round(3).to_string())
    print(pd.DataFrame([summary]).round(4).T.to_string())

    # ---- figure
    plt = rc.setup_mpl()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)
    ax = axes[0]
    T2 = T.sort_values("S")
    y = np.arange(len(T2))
    ax.hlines(y, T2.S_lo, T2.S_hi, color="#0072B2", lw=1)
    ax.plot(T2.S, y, "o", color="#0072B2", ms=4, label="Σ|Δ component| (S)")
    ax.plot(T2.S_null_median, y, "|", color="#0072B2", ms=8, alpha=.5, label="S, permutation median")
    ax.hlines(y + .25, T2.N_lo, T2.N_hi, color="#D55E00", lw=1)
    ax.plot(T2.N, y + .25, "s", color="#D55E00", ms=3.5, label="|Δ composite| (N)")
    ax.set_yticks(y + .12)
    ax.set_yticklabels(T2.celltype)
    ax.set_xlabel("IPF effect, Δ log10 hazard (composite scale)")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("a  Component-level vs net composite signal", loc="left")
    ax = axes[1]
    g = G[G.component != "__summary__"].sort_values("delta")
    cols = ["#999999" if c == "Non-module genes" else ("#D55E00" if v > 0 else "#0072B2")
            for c, v in zip(g.component, g.delta)]
    ax.barh(range(len(g)), g.delta, xerr=[g.delta - g.lo, g.hi - g.delta], color=cols,
            error_kw=dict(lw=0.5))
    ax.set_yticks(range(len(g)))
    ax.set_yticklabels(g.component, fontsize=5.5)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("Global IPF effect on composite component (Δ log10 hazard)")
    ax.set_title(f"b  Composite decomposition (net Δ = {d_obs.sum():+.3f})", loc="left")
    fig.savefig(rc.FIG / "rev07_masking.png")


if __name__ == "__main__":
    main()
