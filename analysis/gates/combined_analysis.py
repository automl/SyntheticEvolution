#!/usr/bin/env python3
"""Combined analysis: SHS 3D quality is GATED by how well RNAformer predicts the 2D structure,
crossed with every covariate (family, length, ions, topology, protein family, era, subtype).

Master frame (per exp_db_id):
  3D:  RNA_TM for no-MSA / natural-MSA / RNAformer-SHS  (from merged_long.csv)
  2D:  rf2d_pred     = TRUE RNAformer 2D-prediction accuracy vs GT_dssr (mcc) -- the MEDIATOR
       shs_2dfrom3d  = 2D read (DSSR) from the SHS 3D model vs GT-2D
       msa_2dfrom3d  = 2D read from the natural-MSA 3D model (1v1 only)
  covariates: family / length / ions / ligands / pH / temp / topology / protein family /
              era / orphan / subtype (rna-only vs 1v1 vs multichain)

2D metrics cover the 91 rna-only + 112 1v1 complexes (203 examples); the multichain set has
no RNAformer-2D metric and is kept only for the covariate context.

Writes tables (*.csv) here and a battery of figures to ../figures/combined/.
"""
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
NE = HERE.parent
FIGDIR = NE / "figures" / "combined"; FIGDIR.mkdir(parents=True, exist_ok=True)

# reuse deep_dive_full.load() for the RNA_TM + covariate frame
_spec = importlib.util.spec_from_file_location("ddf", HERE / "deep_dive_full.py")
ddf = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(ddf)
COLORS = ddf.COLORS; LABEL = ddf.LABEL; METHODS = ddf.METHODS

# prim_family is upper-cased (tok() upper-cases), so keys MUST be upper-case to match.
FAM_COLORS = {"SAM": "#d62728", "PURINE": "#ff7f0e", "THF": "#9467bd", "TRNA": "#2ca02c",
              "HDV_RIBOZYME": "#8c564b"}


def load_2d():
    def up(d):
        d["pdbid"] = d["pdbid"].astype(str).str.upper(); return d
    sc = up(pd.read_csv(NE / "RNA_Monomer_2D_self_consistency.csv"))
    d3 = up(pd.read_csv(NE / "RNA_Monomer_2D_and_3D_metricsN100.csv"))
    rp = up(pd.read_csv(NE / "rp1v1_2D_full_with_pairs_af3_rnaformerN100_dssr.csv"))

    def pick(df, method, name):
        return df[df.method == method][["pdbid", "mcc"]].rename(columns={"mcc": name}).set_index("pdbid")[name]

    # rf2d_pred = TRUE RNAformer 2D-prediction accuracy vs GT_dssr (seed pred scored against GT),
    # built by build_rnaformer_2d_gt.py. NOT self_consistency (seed-vs-realised-3D). See that file.
    rft = pd.read_csv(HERE / "rnaformer_2d_vs_gt.csv")
    rft["pdbid"] = rft.exp_db_id.astype(str).str.upper()
    rf = rft.set_index("pdbid")["rf2d_mcc"].rename("rf2d_pred")
    # keep the old self-consistency mcc too, for transparency / comparison
    sc_series = pd.concat([pick(sc, "rnaformerN100_self_consistency", "rf2d_selfcons"),
                           pick(rp, "rnaformerN100_self_consistency", "rf2d_selfcons")])
    sc_series = sc_series[~sc_series.index.duplicated()]
    shs = pd.concat([pick(d3, "rnaformerN100_dssr", "shs_2dfrom3d"),
                     pick(rp, "rnaformerN100_dssr", "shs_2dfrom3d")])
    msa = pick(rp, "alphafold_dssr", "msa_2dfrom3d")
    # other 2D predictors' self-consistency (rna-only) for a predictor comparison
    preds = sc.pivot_table(index="pdbid", columns="method", values="mcc")
    preds = preds.rename(columns=lambda c: c.replace("_self_consistency", "_2dpred"))
    return rf[~rf.index.duplicated()], shs[~shs.index.duplicated()], msa[~msa.index.duplicated()], preds, sc_series


def load_master():
    w = ddf.load()                                   # index=exp_db_id, RNA_TM wide + covariates
    rf, shs, msa, preds, sc_series = load_2d()
    w = w.join(rf).join(shs).join(msa).join(preds).join(sc_series)
    w["has2d"] = w["rf2d_pred"].notna()
    w["len_bin"] = pd.cut(w.rna_len, [0, 40, 70, 120, 10000], labels=["≤40", "41-70", "71-120", ">120"])
    w["pred_bin"] = pd.cut(w.rf2d_pred, [-.01, .5, .8, .95, 1.01],
                           labels=["<0.5", "0.5-0.8", "0.8-0.95", ">0.95"])
    w.reset_index().to_csv(HERE / "master_combined.csv", index=False)
    return w


# ---------- small plotting helpers ----------
def grouped_bar(ax, piv, order, title, ylab, xlab=""):
    piv = piv.reindex(order)
    x = np.arange(len(piv)); wdt = 0.26
    for i, m in enumerate(METHODS):
        ax.bar(x + (i - 1) * wdt, piv[m], wdt, label=LABEL[m], color=COLORS[m])
    ax.set_xticks(x); ax.set_xticklabels(piv.index, fontsize=8)
    ax.set_title(title, fontsize=10); ax.set_ylabel(ylab); ax.set_xlabel(xlab)
    ax.set_ylim(0, 1); ax.legend(fontsize=7, frameon=False)


def means_by(w, by, order=None):
    g = w.groupby(by, observed=True)
    piv = g[METHODS].mean()
    piv["n"] = g.size()
    piv["rf2d_pred"] = g["rf2d_pred"].mean()
    piv["d_shs_nomsa"] = g["d_shs_nomsa"].mean()
    return piv.reindex(order) if order else piv


def save(fig, name):
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(FIGDIR / f"{name}.{e}", bbox_inches="tight")
    plt.close(fig)


def main():
    w = load_master()
    w2 = w[w.has2d].copy()
    print(f"[combined] {len(w)} examples; {len(w2)} with RNAformer 2D metric "
          f"({(w2.group=='rna-only').sum()} rna-only, {(w2.group=='1v1').sum()} 1v1)")

    # ============ headline correlations ============
    print("\n=== correlations (Spearman, n=203 with 2D) ===")
    corr = {}
    for a, b in [("rnaformerN100", "rf2d_pred"), ("d_shs_nomsa", "rf2d_pred"),
                 ("shs_2dfrom3d", "rf2d_pred"), ("rnaformerN100", "shs_2dfrom3d"),
                 ("rf2d_pred", "rna_len"), ("rf2d_pred", "n_stems")]:
        s = w2[[a, b]].dropna()
        rho = s[a].corr(s[b], method="spearman")
        corr[f"{a}~{b}"] = rho
        print(f"  {a:14s} vs {b:12s}: rho={rho:+.2f} (n={len(s)})")
    pd.Series(corr).to_csv(HERE / "combined_correlations.csv")
    # control for baseline: within the MSA-dependent regime (poor no-MSA), does better 2D-pred
    # give a bigger SHS gain? (the family plot is confounded by 'room to improve')
    print("\n=== gain gated by 2D-pred, WITHIN the hard regime (no-MSA RNA_TM < 0.7) ===")
    for name, sub in [("rna-only", w2[(w2.group == 'rna-only') & (w2.noMSA < 0.7)]),
                      ("rna-only+1v1", w2[w2.noMSA < 0.7])]:
        s = sub[["rf2d_pred", "d_shs_nomsa"]].dropna()
        if len(s) > 8:
            print(f"  {name:14s} n={len(s):3d} | gain~2Dpred rho="
                  f"{s.rf2d_pred.corr(s.d_shs_nomsa, method='spearman'):+.2f} | mean gain "
                  f"{sub.d_shs_nomsa.mean():+.3f}")

    # ============ 1) gating: 3D by RNAformer-2D-pred bin (rna-only / 1v1 / all) ============
    order = ["<0.5", "0.5-0.8", "0.8-0.95", ">0.95"]
    rows = []
    for grp, sub in [("rna-only", w2[w2.group == "rna-only"]),
                     ("1v1", w2[w2.group == "1v1"]), ("all-2D", w2)]:
        piv = means_by(sub, "pred_bin", order)
        for b in order:
            if b in piv.index and piv.loc[b, "n"] >= 3:
                r = piv.loc[b]
                rows.append({"subset": grp, "pred_bin": b, "n": int(r.n),
                             "noMSA": r.noMSA, "MSA": r.alphafold, "SHS": r.rnaformerN100,
                             "d_shs_nomsa": r.d_shs_nomsa})
    gate = pd.DataFrame(rows); gate.to_csv(HERE / "combined_gating_by_2dpred.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    for ax, grp in zip(axes, ["rna-only", "1v1", "all-2D"]):
        sub = w2 if grp == "all-2D" else w2[w2.group == grp]
        grouped_bar(ax, means_by(sub, "pred_bin", order)[METHODS], order,
                    f"{grp}: 3D RNA_TM by RNAformer 2D-pred (n={len(sub)})", "mean RNA_TM",
                    "RNAformer 2D-pred mcc")
    fig.suptitle("SHS 3D accuracy is gated by RNAformer 2D-prediction quality "
                 "(good 2D -> SHS recovers the MSA benefit; near-perfect -> beats MSA)")
    save(fig, "fig_gating_by_2dpred")

    # ============ 2) headline scatter: SHS 3D-TM vs RNAformer 2D-pred, families highlighted ====
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    base = w2[~w2.prim_family.isin(FAM_COLORS)]
    ax.scatter(base.rf2d_pred, base.rnaformerN100, s=18, c="#cccccc", label="other", alpha=0.7)
    for fam, col in FAM_COLORS.items():
        s = w2[w2.prim_family == fam]
        ax.scatter(s.rf2d_pred, s.rnaformerN100, s=34, c=col, label=f"{fam} (n={len(s)})",
                   edgecolor="k", linewidth=0.3)
    ax.set_xlabel("RNAformer 2D-prediction accuracy vs GT (mcc)")
    ax.set_ylabel("SHS 3D RNA_TM")
    ax.set_title("SHS 3D accuracy tracks RNAformer 2D-prediction accuracy (rho~0.5)")
    ax.legend(fontsize=8, frameon=False)
    save(fig, "fig_scatter_3d_vs_2dpred")

    # ============ 3) family mediation: per-family 2D-pred vs SHS gain ============
    frows = []
    for fam, g in w2[w2.prim_family != "(no family)"].groupby("prim_family"):
        if len(g) < 3:
            continue
        frows.append({"family": fam, "n": len(g), "rf2d_pred": g.rf2d_pred.mean(),
                      "noMSA": g.noMSA.mean(), "MSA": g.alphafold.mean(), "SHS": g.rnaformerN100.mean(),
                      "d_shs_nomsa": g.d_shs_nomsa.mean(), "d_shs_msa": g.d_shs_msa.mean(),
                      "shs_2dfrom3d": g.shs_2dfrom3d.mean()})
    fam = pd.DataFrame(frows).sort_values("d_shs_nomsa", ascending=False)
    fam.to_csv(HERE / "combined_family.csv", index=False)
    print("\n=== per-family (n>=3, with 2D): RNAformer-2D-pred vs SHS 3D gain ===")
    print(fam[["family", "n", "rf2d_pred", "noMSA", "MSA", "SHS", "d_shs_nomsa"]].round(3).to_string(index=False))
    fig, ax = plt.subplots(figsize=(7, 5.2))
    ax.axhline(0, color="k", lw=0.6)
    ax.scatter(fam.rf2d_pred, fam.d_shs_nomsa, s=40 + 8 * fam.n,
               c=[FAM_COLORS.get(f, "#4477aa") for f in fam.family], edgecolor="k", linewidth=0.4)
    for _, r in fam.iterrows():
        ax.annotate(r.family, (r.rf2d_pred, r.d_shs_nomsa), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("RNAformer 2D-prediction quality (family mean mcc)")
    ax.set_ylabel("SHS 3D gain over no-MSA (Δ RNA_TM)")
    ax.set_title("A large SHS 3D gain needs BOTH accurate RNAformer 2D and a poor no-MSA baseline\n"
                 "(SAM has both; tRNA/rRNA are predicted well but already easy → no room)")
    save(fig, "fig_family_mediation")

    # ============ 4) 2D-pred quality by covariate (family / length / subtype) ============
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    fa = fam.sort_values("rf2d_pred", ascending=False)
    axes[0].bar(range(len(fa)), fa.rf2d_pred, color=[FAM_COLORS.get(f, "#4477aa") for f in fa.family])
    axes[0].set_xticks(range(len(fa))); axes[0].set_xticklabels(fa.family, rotation=40, ha="right", fontsize=7)
    axes[0].set_ylabel("RNAformer 2D-pred mcc"); axes[0].set_title("2D-pred quality by RNA family"); axes[0].set_ylim(0, 1)
    lp = w2.groupby("len_bin", observed=True)["rf2d_pred"].mean().reindex(["≤40", "41-70", "71-120", ">120"])
    axes[1].bar(range(len(lp)), lp.values, color="#4477aa")
    axes[1].set_xticks(range(len(lp))); axes[1].set_xticklabels(lp.index, fontsize=8)
    axes[1].set_ylabel("RNAformer 2D-pred mcc"); axes[1].set_title("2D-pred quality by RNA length"); axes[1].set_ylim(0, 1)
    gp = w2.groupby("group", observed=True)["rf2d_pred"].mean()
    axes[2].bar(range(len(gp)), gp.values, color="#4477aa")
    axes[2].set_xticks(range(len(gp))); axes[2].set_xticklabels(gp.index, fontsize=8)
    axes[2].set_ylabel("RNAformer 2D-pred mcc"); axes[2].set_title("2D-pred quality by subtype"); axes[2].set_ylim(0, 1)
    save(fig, "fig_2dpred_by_covariate")

    # ============ 5) Δ(SHS-noMSA) vs 2D-pred, colored by length bin ============
    fig, ax = plt.subplots(figsize=(7.2, 5))
    ax.axhline(0, color="k", lw=0.6)
    lcolors = {"≤40": "#08519c", "41-70": "#3182bd", "71-120": "#f16913", ">120": "#a50f15"}
    for lb, col in lcolors.items():
        s = w2[w2.len_bin == lb]
        ax.scatter(s.rf2d_pred, s.d_shs_nomsa, s=22, c=col, alpha=0.75, label=f"{lb} (n={len(s)})")
    ax.set_xlabel("RNAformer 2D-prediction quality (mcc)")
    ax.set_ylabel("SHS 3D gain over no-MSA (Δ RNA_TM)")
    ax.set_title("SHS helps most when 2D-pred is good — across all lengths")
    ax.legend(fontsize=8, frameon=False, title="RNA length")
    save(fig, "fig_delta_vs_2dpred_by_length")

    # ============ 6) 2D fidelity: 2D-from-3D vs 2D-pred (does 3D keep the topology) ============
    fig, ax = plt.subplots(figsize=(6.4, 5))
    for fam_, col in FAM_COLORS.items():
        s = w2[w2.prim_family == fam_]
        ax.scatter(s.rf2d_pred, s.shs_2dfrom3d, s=30, c=col, label=fam_, edgecolor="k", linewidth=0.3)
    other = w2[~w2.prim_family.isin(FAM_COLORS)]
    ax.scatter(other.rf2d_pred, other.shs_2dfrom3d, s=16, c="#cccccc", alpha=0.6, zorder=0)
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("RNAformer 2D-prediction mcc")
    ax.set_ylabel("2D read from SHS 3D model (mcc)")
    ax.set_title("The SHS 3D model faithfully carries RNAformer's predicted topology (2D-from-3D vs 2D-pred rho~0.5)")
    ax.legend(fontsize=8, frameon=False)
    save(fig, "fig_2dfrom3d_vs_2dpred")

    # ============ 7) predictor comparison (rna-only self-consistency) ============
    predcols = [c for c in w.columns if c.endswith("_2dpred")]
    if predcols:
        pm = w[predcols].mean().sort_values(ascending=False)
        pm.to_csv(HERE / "combined_predictor_selfconsistency.csv")
        print("\n=== 2D self-consistency by predictor (rna-only mean mcc) ===")
        print(pm.round(3).to_string())

    # length-combined table
    lb_tab = means_by(w2, "len_bin", ["≤40", "41-70", "71-120", ">120"])
    lb_tab.to_csv(HERE / "combined_length.csv")
    print(f"\nDone. tables in {HERE}; figures in {FIGDIR}")


if __name__ == "__main__":
    main()
