#!/usr/bin/env python3
"""Stratified analyses for the RNAformer-SHS downstream (AF3) evaluation.

Two questions:
  (a) failure analysis: where does RNAformer-SHS degrade vs natural MSA / no-MSA?
  (b) stratification:    performance split by era (pre/post cutoff), subtype
                         (complex vs rna-only), and family-present vs absent.

Data:
  - raw per-method eval CSVs  <method>_{protein_rna,rna}.csv  (the exact files the
    preprocessing uses; comma- or semicolon-delimited),
  - metadata_468.csv          (per exp_db_id: era, subtype, family, topology, biochem),
  - the two pdb_*_from_alphafold_eval.csv + the 4 nar_dbs/*.db  (to build a merged
    RNA-family label and has_family flag; also written back into metadata_468.csv).

Methods contrasted:
  no-MSA        = noMSA
  natural MSA   = alphafold
  RNAformer-SHS = rnaformerN100     (NOT the *SHSandMSA version)

Metrics (kept separate for rna-only vs complexes):
  complex      -> Complex_TM (primary), Dockq (secondary)
  single_chain -> RNA_TM (primary),     RNA_LDDT_allAtom (secondary)

Outputs (in this script's directory): tables (*.csv) + figures (*.png/*.pdf).
Run inside an env with pandas + matplotlib, e.g. `conda run -n af3_eval python gate_analyses.py`.
"""
import ast
import re
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
NE = HERE.parent                                   # neurips_shs_2026/
RAW = NE / "iris_eval_neurips_2026_04_05_2026"
META = NE / "metadata_468.csv"
DBS = NE / "nar_dbs"
OUT = HERE                       # tables land next to the script
FIGDIR = NE / "figures"          # dedicated figures folder within the neurips dir
FIGDIR.mkdir(exist_ok=True)

# reuse the NAR figure style if importable, else fall back to sensible defaults
try:
    sys.path.insert(0, str(HERE))
    import figstyle  # noqa
    figstyle.write_matplotlibrc(str(HERE / "matplotlibrc"))
    matplotlib.rc_file(str(HERE / "matplotlibrc"))
    print("[style] using figstyle.py")
except Exception as exc:  # pragma: no cover
    print(f"[style] figstyle unavailable ({exc}); using local defaults")
    matplotlib.rcParams.update({"font.size": 10, "axes.spines.top": False,
                                "axes.spines.right": False, "figure.dpi": 120})

METHODS = {"noMSA": "no-MSA", "alphafold": "natural MSA", "rnaformerN100": "RNAformer-SHS"}
NOMSA, MSA, SHS = "noMSA", "alphafold", "rnaformerN100"
PRIMARY = {"complex": "Complex_TM", "single_chain": "RNA_TM"}
SECONDARY = {"complex": "Dockq", "single_chain": "RNA_LDDT_allAtom"}
COLORS = {"noMSA": "#9E9E9E", "alphafold": "#1f77b4", "rnaformerN100": "#d62728"}


# --------------------------------------------------------------------------------------
# loaders
# --------------------------------------------------------------------------------------
def read_eval_csv(method, subtype):
    """Load <method>_{protein_rna|rna}.csv, sniffing comma vs semicolon (like the preproc)."""
    fn = f"{method}_{'protein_rna' if subtype == 'complex' else 'rna'}.csv"
    path = RAW / fn
    try:
        df = pd.read_csv(path)
        if df.shape[1] == 1:            # single fat column => wrong delimiter
            raise ValueError("one column")
    except Exception:
        df = pd.read_csv(path, sep=";")
    df["exp_db_id"] = df["exp_db_id"].astype(str).str.upper().str.strip()
    return df


def num(s):
    """Coerce a metric column to float, unwrapping the occasional '[x, ...]' list-string."""
    def one(v):
        if isinstance(v, (int, float)):
            return float(v)
        if isinstance(v, str):
            v = v.strip()
            if v.startswith("["):
                try:
                    x = ast.literal_eval(v)
                    return float(x[0]) if isinstance(x, list) and x else np.nan
                except Exception:
                    return np.nan
            try:
                return float(v)
            except ValueError:
                return np.nan
        return np.nan
    return s.map(one)


def family_tokens(v):
    toks = [t.strip() for t in re.split(r'[\[\]",]', str(v or "")) if t.strip()]
    return {t for t in toks if t.lower() != "none"}


# --------------------------------------------------------------------------------------
# 1. enrich metadata with a merged RNA family (union of the CSV + DB annotations)
# --------------------------------------------------------------------------------------
def enrich_metadata():
    md = pd.read_csv(META, dtype=str, keep_default_na=False)
    md["exp_db_id"] = md["exp_db_id"].str.upper().str.strip()

    fam = {i: set() for i in md["exp_db_id"]}
    # from the two exported CSVs
    for f in ["pdb_protein_rna_from_alphafold_eval.csv", "pdb_rna_rna_from_alphafold_eval.csv"]:
        d = pd.read_csv(NE / f, dtype=str, keep_default_na=False)
        for _, r in d.iterrows():
            pid = str(r.get("PDBId", "")).upper().strip()
            if pid in fam:
                fam[pid] |= family_tokens(r.get("RNAFamily"))
    # from the 4 NAR DBs
    for db, tab in [("a2021_protein_rna.db", "exp_protein_rna"), ("b2021_protein_rna.db", "exp_protein_rna"),
                    ("a2021_rna_only.db", "exp_rna_rna"), ("b2021_rna_only.db", "exp_rna_rna")]:
        p = DBS / db
        if not p.exists():
            continue
        c = sqlite3.connect(p)
        for pid, rf in c.execute(f"SELECT PDBId, RNAFamily FROM {tab}"):
            pid = str(pid).upper().strip()
            if pid in fam:
                fam[pid] |= family_tokens(rf)
        c.close()

    md["RNAFamily_merged"] = md["exp_db_id"].map(lambda i: ";".join(sorted(fam[i])))
    md["has_family"] = md["RNAFamily_merged"].map(lambda s: bool(s))
    md.to_csv(META, index=False)
    print(f"[metadata] merged family: {md['has_family'].sum()}/{len(md)} have a family "
          f"(was 163 CSV / 155 DB); wrote {META.name}")
    return md


# --------------------------------------------------------------------------------------
# 2. merged long frame: one row per (exp_db_id, subtype, method)
# --------------------------------------------------------------------------------------
def build_frame(md):
    md_small = md[["exp_db_id", "subtype", "era", "has_family", "RNAFamily_merged",
                   "RNALength", "RNA_Pseudoknots"]].copy()

    def max_len(v):
        try:
            x = ast.literal_eval(str(v))
            return float(np.max(x)) if isinstance(x, list) and x else float(x)
        except Exception:
            try:
                return float(str(v).split(",")[0])
            except Exception:
                return np.nan
    md_small["rna_len"] = md_small["RNALength"].map(max_len)
    md_small["has_pseudoknot"] = md_small["RNA_Pseudoknots"].map(
        lambda v: bool(family_tokens(v)) or str(v).strip() not in ("", "[]", "None", "nan"))

    # orphan / non-orphan = the paper's MSA-availability split (from the processed csvs).
    # This is the natural "when should synthetic topology be trusted?" axis (af3_rna_MSA is
    # empty in the raw eval CSVs, so we use the curated orphan split instead).
    ss = NE / "2026_neurips_code" / "results" / "csvs_neurips"
    orphan, nonorphan = set(), set()
    if ss.exists():
        for f in ss.glob("*_orphan.csv"):
            if f.name.endswith("_non_orphan.csv"):
                continue
            try:
                orphan |= set(pd.read_csv(f, dtype=str)["exp_db_id"].str.upper())
            except Exception:
                pass
        for f in ss.glob("*_non_orphan.csv"):
            try:
                nonorphan |= set(pd.read_csv(f, dtype=str)["exp_db_id"].str.upper())
            except Exception:
                pass
    md_small["is_orphan"] = md_small["exp_db_id"].str.upper().map(
        lambda i: True if i in orphan else (False if i in nonorphan else np.nan))
    print(f"[orphan] {md_small['is_orphan'].sum():.0f} orphan / "
          f"{(md_small['is_orphan'] == False).sum():.0f} non-orphan (of {len(md_small)})")

    rows = []
    for method in METHODS:
        for sub in ("complex", "single_chain"):
            d = read_eval_csv(method, sub)
            keep = {"exp_db_id": d["exp_db_id"]}
            for m in {PRIMARY[sub], SECONDARY[sub], "RNA_TM", "Complex_TM", "Dockq", "RNA_LDDT_allAtom"}:
                if m in d.columns:
                    keep[m] = num(d[m])
            t = pd.DataFrame(keep)
            t["method"] = method
            t["subtype"] = sub
            rows.append(t)
    long = pd.concat(rows, ignore_index=True)
    # restrict to the curated 468 + attach covariates
    long = long.merge(md_small[["exp_db_id", "subtype", "era", "has_family",
                                 "RNAFamily_merged", "rna_len", "has_pseudoknot", "is_orphan"]],
                      on=["exp_db_id", "subtype"], how="inner")
    long["method_label"] = long["method"].map(METHODS)
    long.to_csv(OUT / "merged_long.csv", index=False)
    print(f"[frame] merged_long.csv: {len(long)} rows "
          f"({long['exp_db_id'].nunique()} ids x {long['method'].nunique()} methods)")
    return long


# --------------------------------------------------------------------------------------
# 3. failure analysis (R2)
# --------------------------------------------------------------------------------------
def wide_primary(long, sub):
    d = long[long["subtype"] == sub]
    metric = PRIMARY[sub]
    w = d.pivot_table(index="exp_db_id", columns="method", values=metric, aggfunc="first")
    cov = d.drop_duplicates("exp_db_id").set_index("exp_db_id")[
        ["era", "has_family", "RNAFamily_merged", "rna_len", "has_pseudoknot", "is_orphan"]]
    w = w.join(cov)
    w = w.dropna(subset=[NOMSA, MSA, SHS])          # need all three to compare
    w["d_shs_vs_msa"] = w[SHS] - w[MSA]
    w["d_shs_vs_nomsa"] = w[SHS] - w[NOMSA]
    w["metric"] = metric
    return w


def failure_analysis(long):
    summary = []
    for sub in ("complex", "single_chain"):
        w = wide_primary(long, sub)
        metric = PRIMARY[sub]
        n = len(w)
        worse_msa = w[w["d_shs_vs_msa"] < 0]
        worse_nomsa = w[w["d_shs_vs_nomsa"] < 0]
        print(f"\n[R2 {sub}] {metric}: n={n} | SHS<MSA: {len(worse_msa)} "
              f"({100*len(worse_msa)/n:.0f}%), mean Δ={w['d_shs_vs_msa'].mean():+.3f} | "
              f"SHS<noMSA: {len(worse_nomsa)} ({100*len(worse_nomsa)/n:.0f}%), "
              f"mean Δ={w['d_shs_vs_nomsa'].mean():+.3f}")
        # ranked failure tables
        cols = [NOMSA, MSA, SHS, "d_shs_vs_msa", "d_shs_vs_nomsa", "era",
                "has_family", "RNAFamily_merged", "rna_len", "has_pseudoknot", "is_orphan"]
        w.sort_values("d_shs_vs_msa")[cols].to_csv(OUT / f"failure_shs_vs_msa_{sub}.csv")
        w.sort_values("d_shs_vs_nomsa")[cols].to_csv(OUT / f"failure_shs_vs_nomsa_{sub}.csv")
        summary.append((sub, metric, n, len(worse_msa), len(worse_nomsa),
                        w["d_shs_vs_msa"].mean(), w["d_shs_vs_nomsa"].mean()))
    pd.DataFrame(summary, columns=["subtype", "metric", "n", "n_SHS<MSA", "n_SHS<noMSA",
                                   "meanΔ_SHS-MSA", "meanΔ_SHS-noMSA"]).to_csv(
        OUT / "failure_summary.csv", index=False)

    # ---- plots ----
    subs = ("complex", "single_chain")
    # (a) Δ distributions
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, sub in zip(axes, subs):
        w = wide_primary(long, sub)
        data = [w["d_shs_vs_msa"].values, w["d_shs_vs_nomsa"].values]
        bp = ax.boxplot(data, tick_labels=["SHS − MSA", "SHS − noMSA"], showmeans=True,
                        widths=0.5, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("#d6272822")
        for i, arr in enumerate(data, 1):
            ax.scatter(np.random.default_rng(0).normal(i, 0.05, len(arr)), arr,
                       s=8, alpha=0.35, color=COLORS[SHS])
        ax.axhline(0, color="k", lw=0.8, ls="--")
        ax.set_title(f"{sub}  ({PRIMARY[sub]}, n={len(w)})")
        ax.set_ylabel("Δ (RNAformer-SHS − baseline)")
    fig.suptitle("Where RNAformer-SHS degrades vs natural MSA / no-MSA (below 0 = worse)")
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(FIGDIR / f"fig_failure_delta.{e}", bbox_inches="tight")
    plt.close(fig)

    # (b) scatter SHS vs MSA, colored by era
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4))
    for ax, sub in zip(axes, subs):
        w = wide_primary(long, sub)
        for era, col in [("pre", "#1f77b4"), ("post", "#d62728")]:
            s = w[w["era"] == era]
            ax.scatter(s[MSA], s[SHS], s=16, alpha=0.6, color=col, label=f"{era} (n={len(s)})")
        lim = [min(w[MSA].min(), w[SHS].min()), max(w[MSA].max(), w[SHS].max())]
        ax.plot(lim, lim, "k--", lw=0.8)
        ax.set_xlabel(f"natural MSA  {PRIMARY[sub]}")
        ax.set_ylabel(f"RNAformer-SHS  {PRIMARY[sub]}")
        ax.set_title(f"{sub}: below diagonal = SHS worse")
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(FIGDIR / f"fig_scatter_shs_vs_msa.{e}", bbox_inches="tight")
    plt.close(fig)

    # (c) performance by method x MSA-availability (orphan) -> "when is SHS worth it?"
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4))
    for ax, sub in zip(axes, subs):
        d = long[(long["subtype"] == sub) & long["is_orphan"].notna()]
        d = d.assign(grp=np.where(d["is_orphan"], "orphan\n(no natural MSA)", "non-orphan\n(has MSA)"))
        piv = d.pivot_table(index="grp", columns="method", values=PRIMARY[sub], aggfunc="mean") \
               .reindex(["orphan\n(no natural MSA)", "non-orphan\n(has MSA)"])
        _grouped_bar(ax, piv, f"{sub}: {PRIMARY[sub]} by MSA availability", f"mean {PRIMARY[sub]}")
    fig.suptitle("When is RNAformer-SHS worth it? (orphan = no natural MSA available)")
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(FIGDIR / f"fig_orphan_msa_availability.{e}", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------
# 4. stratification (R1)
# --------------------------------------------------------------------------------------
def _grouped_bar(ax, piv, title, ylabel):
    methods = [m for m in METHODS if m in piv.columns]
    groups = list(piv.index)
    x = np.arange(len(groups)); w = 0.8 / max(len(methods), 1)
    for i, m in enumerate(methods):
        ax.bar(x + i * w, piv[m].values, w, label=METHODS[m], color=COLORS[m])
    ax.set_xticks(x + w * (len(methods) - 1) / 2)
    ax.set_xticklabels(groups)
    ax.set_title(title); ax.set_ylabel(ylabel)
    ax.legend(fontsize=8, frameon=False)


def stratification(long):
    recs = []
    for sub in ("complex", "single_chain"):
        d = long[long["subtype"] == sub]
        metric = PRIMARY[sub]
        for keys in (["era"], ["has_family"], ["era", "has_family"]):
            g = d.groupby(keys + ["method"])[metric].agg(["mean", "median", "count"]).reset_index()
            g["subtype"] = sub; g["metric"] = metric; g["strata"] = "+".join(keys)
            recs.append(g)
    strata = pd.concat(recs, ignore_index=True)
    strata.to_csv(OUT / "summary_strata.csv", index=False)

    # figure: mean primary metric by method x era, and by method x has_family (2x2)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for j, sub in enumerate(("complex", "single_chain")):
        d = long[long["subtype"] == sub]; metric = PRIMARY[sub]
        piv_era = d.pivot_table(index="era", columns="method", values=metric, aggfunc="mean")
        piv_era = piv_era.reindex(["pre", "post"])
        _grouped_bar(axes[j, 0], piv_era, f"{sub}: {metric} by era (pre/post cutoff)", f"mean {metric}")
        piv_fam = d.assign(fam=np.where(d["has_family"], "family", "no-family")) \
            .pivot_table(index="fam", columns="method", values=metric, aggfunc="mean") \
            .reindex(["family", "no-family"])
        _grouped_bar(axes[j, 1], piv_fam, f"{sub}: {metric} by family presence", f"mean {metric}")
    fig.suptitle("Stratified performance: era (generalization) and family presence")
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(FIGDIR / f"fig_strata.{e}", bbox_inches="tight")
    plt.close(fig)
    # console summary
    print("\n[R1] mean primary metric by method x era:")
    for sub in ("complex", "single_chain"):
        d = long[long["subtype"] == sub]
        print(f"  {sub} ({PRIMARY[sub]}):")
        print(d.pivot_table(index="era", columns="method_label", values=PRIMARY[sub],
                            aggfunc="mean").reindex(["pre", "post"]).round(3).to_string())


def main():
    md = enrich_metadata()
    long = build_frame(md)
    failure_analysis(long)
    stratification(long)
    print(f"\nDone. Tables + figures in {OUT}")


if __name__ == "__main__":
    main()
