#!/usr/bin/env python3
"""Ground-truth-free trust filter for RNAformer-SHS.

self_consistency = inject SHS seed -> predict 3D -> read 2D back -> compare to seed. It needs NO
ground truth, so it is computable at inference time as a confidence score. When AF3 is self-
consistent with the injected topology (the model 'commits' to it), the SHS prediction matches
natural MSA; when it is not, defer to MSA / no-MSA. This is the deployable complement to the
mechanistic rf2d-vs-GT gate (which needs the answer).

rna-only (n=91). Bars = no-MSA / MSA / SHS RNA_TM by RNAformer self-consistency mcc bin.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
NE = HERE.parent
RAW = NE / "iris_eval_neurips_2026_04_05_2026"
FIG = NE / "figures" / "trust"; FIG.mkdir(parents=True, exist_ok=True)
try:
    import sys; sys.path.insert(0, str(HERE)); import figstyle
    figstyle.write_matplotlibrc(str(HERE / "matplotlibrc")); matplotlib.rc_file(str(HERE / "matplotlibrc"))
except Exception:
    matplotlib.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

MCOL = {"noMSA": "#9E9E9E", "alphafold": "#1f77b4", "rnaformerN100": "#d62728"}
MLAB = {"noMSA": "no-MSA", "alphafold": "natural MSA", "rnaformerN100": "RNAformer-SHS"}
METHODS = ["noMSA", "alphafold", "rnaformerN100"]


def tm(m):
    p = RAW / f"{m}_rna.csv"
    try:
        d = pd.read_csv(p)
        if d.shape[1] == 1:
            raise ValueError
    except Exception:
        d = pd.read_csv(p, sep=";")
    d["id"] = d.exp_db_id.astype(str).str.upper()
    return d.set_index("id")["RNA_TM"].astype(float).rename(m)


def star(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"


def main():
    sc = pd.read_csv(NE / "RNA_Monomer_2D_self_consistency.csv"); sc["id"] = sc.pdbid.str.upper()
    scm = sc[sc.method == "rnaformerN100_self_consistency"].set_index("id")["mcc"].rename("selfcons")
    w = pd.concat([scm] + [tm(m) for m in METHODS], axis=1).dropna()

    edges = [-.01, .8, .9, .95, 1.01]
    labs = ["<0.8", "0.8-0.9", "0.9-0.95", "≥0.95"]
    w["b"] = pd.cut(w.selfcons, edges, labels=labs)
    piv = w.groupby("b", observed=True)[METHODS].mean().reindex(labs)
    cnt = w.groupby("b", observed=True).size().reindex(labs)

    print("=== rna-only: SHS vs MSA / no-MSA gated by GT-free self-consistency mcc ===")
    rows = []
    for b in labs:
        g = w[w.b == b]
        pM = wilcoxon(g.rnaformerN100, g.alphafold).pvalue if len(g) >= 6 else np.nan
        pN = wilcoxon(g.rnaformerN100, g.noMSA).pvalue if len(g) >= 6 else np.nan
        rows.append({"bin": b, "n": len(g), "noMSA": g.noMSA.mean(), "MSA": g.alphafold.mean(),
                     "SHS": g.rnaformerN100.mean(), "d_SHS_MSA": (g.rnaformerN100 - g.alphafold).mean(),
                     "p_MSA": pM, "d_SHS_noMSA": (g.rnaformerN100 - g.noMSA).mean(), "p_noMSA": pN})
        print(f"  {b:9s} n={len(g):2d} | noMSA {g.noMSA.mean():.3f} MSA {g.alphafold.mean():.3f} "
              f"SHS {g.rnaformerN100.mean():.3f} | ΔSHS-MSA {rows[-1]['d_SHS_MSA']:+.3f}({star(pM) if pd.notna(pM) else '-'}) "
              f"ΔSHS-noMSA {rows[-1]['d_SHS_noMSA']:+.3f}({star(pN) if pd.notna(pN) else '-'})")
    pd.DataFrame(rows).to_csv(HERE / "selfcons_trust.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    x = np.arange(len(labs)); wd = 0.26
    for i, m in enumerate(METHODS):
        ax.bar(x + (i - 1) * wd, piv[m], wd, label=MLAB[m], color=MCOL[m], zorder=3)
    for xi, b in enumerate(labs):
        g = w[w.b == b]
        if len(g) >= 6:
            pM = wilcoxon(g.rnaformerN100, g.alphafold).pvalue
            pN = wilcoxon(g.rnaformerN100, g.noMSA).pvalue
            top = max(piv.loc[b])
            ax.text(xi, top + 0.03, f"vs MSA:{star(pM)}\nvs noMSA:{star(pN)}", ha="center", fontsize=7)
    # mark the trustworthy regime
    ax.axvspan(1.5, 3.5, color="#2ca02c", alpha=0.06, zorder=0)
    ax.text(2.5, 1.0, "trustworthy SHS (self-consistent)\nSHS ≈ MSA, >> no-MSA",
            ha="center", va="top", fontsize=8, color="#2a7a2a")
    ax.set_xticks(x); ax.set_xticklabels([f"{b}\n(n={int(cnt[b])})" for b in labs])
    ax.set_xlabel("RNAformer SHS self-consistency mcc  (GT-free: seed vs 2D-read-from-3D)")
    ax.set_ylabel("mean RNA_TM"); ax.set_ylim(0, 1.08)
    ax.set_title("A GT-free trust filter for SHS: when AF3 commits to the injected topology\n"
                 "(high self-consistency), RNAformer-SHS reaches natural-MSA accuracy")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(FIG / f"fig_selfcons_trust_gate.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nfigure: {FIG / 'fig_selfcons_trust_gate.png'}")


if __name__ == "__main__":
    main()
