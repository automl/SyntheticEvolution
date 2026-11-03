#!/usr/bin/env python3
"""(A) Post-cutoff INTERFACE-metric significance for RNA-protein complexes, and
(B) pre-vs-post FAMILY leakage: do all methods predict post-cutoff RNAs better when the
family was already present pre-cutoff (AF3 training era)?
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import wilcoxon, mannwhitneyu
import re

from paired_stats import subset, paired   # reuse loaders/stat helper

HERE = Path(__file__).resolve().parent
NE = HERE.parent


# ---------------- (A) interface significance ----------------
def interface_significance():
    RP = [("RP full", "RNA-Protein", "", ""),
          ("RP with-homologs", "RNA-Protein", "", "non_orphan"),
          ("RP post-cutoff", "RNA-Protein", "a2021", ""),
          ("RP post-cutoff ORPHAN", "RNA-Protein", "a2021", "orphan")]
    IFACE = ["Dockq", "ICS", "QS_global", "iLDDT", "Complex_LDDT_allAtom", "Complex_TM"]
    rows = []
    print("=== (A) RNA-protein INTERFACE metrics: SHS vs MSA (paired Wilcoxon) ===")
    for metric in IFACE:
        print(f"\n-- {metric} --")
        for label, level, era, split in RP:
            j = subset(level, era, split, metric)
            if j is None or len(j) < 5:
                continue
            dM, pM, loM, hiM = paired(j.SHS.values, j.MSA.values); nM = len(j)
            dN, pN, loN, hiN = paired(j.SHS.values, j.noMSA.values)
            star = "*" if (pd.notna(pM) and pM < 0.05) else " "
            print(f"  {label:24s} n={nM:3d} | MSA {j.MSA.mean():.3f} SHS {j.SHS.mean():.3f} "
                  f"noMSA {j.noMSA.mean():.3f} | Δ(SHS-MSA) {dM:+.3f} p={pM:.3f}{star} "
                  f"CI[{loM:+.3f},{hiM:+.3f}] | Δ(SHS-noMSA) {dN:+.3f} p={pN:.3f}")
            rows.append({"metric": metric, "subset": label, "n": nM, "MSA": j.MSA.mean(),
                         "SHS": j.SHS.mean(), "noMSA": j.noMSA.mean(), "d_SHS_MSA": dM,
                         "p_SHS_MSA": pM, "d_SHS_noMSA": dN, "p_SHS_noMSA": pN})
    pd.DataFrame(rows).to_csv(HERE / "interface_significance.csv", index=False)


# ---------------- (B) pre/post family leakage ----------------
def famtoks(v):
    return {t.strip().upper() for t in re.split(r"[;\[\]\",]", str(v)) if t.strip() and t.strip().lower() != "none"}


def family_leakage():
    md = pd.read_csv(NE / "metadata_468.csv", dtype=str, keep_default_na=False)
    md["exp_db_id"] = md.exp_db_id.str.upper()
    ml = pd.read_csv(HERE / "merged_long.csv")
    tm = ml.pivot_table(index="exp_db_id", columns="method", values="RNA_TM", aggfunc="first")
    d = md.set_index("exp_db_id").join(tm).dropna(subset=["noMSA", "alphafold", "rnaformerN100"])
    d["fams"] = d.RNAFamily_merged.map(famtoks)

    pre_fams = set().union(*d[d.era == "pre"].fams) if (d.era == "pre").any() else set()

    def status(r):
        if not r.fams:
            return "no-family"
        return "seen-in-pre" if (r.fams & pre_fams) else "novel-in-post"
    d["fstatus"] = d.apply(status, axis=1)

    print("\n\n=== (B) POST-cutoff performance by family status (RNA_TM) ===")
    print(f"pre-cutoff families: {len(pre_fams)} distinct")
    post = d[d.era == "post"]
    for sub in ["single_chain", "complex", "ALL"]:
        p = post if sub == "ALL" else post[post.subtype == sub]
        print(f"\n-- post-cutoff {sub} (n={len(p)}) --")
        g = p.groupby("fstatus").agg(n=("noMSA", "size"), noMSA=("noMSA", "mean"),
                                     MSA=("alphafold", "mean"), SHS=("rnaformerN100", "mean"))
        g["d_SHS_MSA"] = g.SHS - g.MSA
        print(g.round(3).to_string())
        # leakage test: MSA on seen-in-pre vs (novel + no-family), Mann-Whitney
        seen = p[p.fstatus == "seen-in-pre"]
        unseen = p[p.fstatus != "seen-in-pre"]
        if len(seen) >= 5 and len(unseen) >= 5:
            for m, lab in [("alphafold", "MSA"), ("rnaformerN100", "SHS"), ("noMSA", "noMSA")]:
                u = mannwhitneyu(seen[m], unseen[m], alternative="greater")
                print(f"    leakage[{lab}]: seen {seen[m].mean():.3f} vs unseen {unseen[m].mean():.3f} "
                      f"(MannWhitney p={u.pvalue:.3f})")
    d.reset_index()[["exp_db_id", "subtype", "era", "fstatus", "RNAFamily_merged",
                     "noMSA", "alphafold", "rnaformerN100"]].to_csv(HERE / "family_leakage.csv", index=False)


if __name__ == "__main__":
    interface_significance()
    family_leakage()
