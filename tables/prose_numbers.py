#!/usr/bin/env python3
"""Paired statistics quoted in the manuscript text that are not printed in a table (Secs. 6 and 7).

Same per-target values as all tables (tables/make_crc_tables.load_condition: mean over 4 AF3 seeds x 5 samples;
RhoFold+: one deterministic prediction). Paired Wilcoxon signed-rank p and 5000x bootstrap 95% CI of the mean difference
(numpy default_rng(0), fresh per comparison).
  python tables/prose_numbers.py      -> results/tables/prose_numbers.md
"""
import glob, importlib.util, statistics as st
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
_s = importlib.util.spec_from_file_location("extra", ROOT / "tables" / "make_extra_tables.py")
E = importlib.util.module_from_spec(_s); _s.loader.exec_module(E)
T = E.T
ROWS = []


def compare(section, what, a, b, ids, metric="RNA_TM"):
    x = np.array([a[i][metric] for i in ids]); y = np.array([b[i][metric] for i in ids]); d = x - y
    rng = np.random.default_rng(0); bs = [rng.choice(d, len(d)).mean() for _ in range(5000)]
    ROWS.append([section, what, len(ids), f"{x.mean():.3f}", f"{y.mean():.3f}", f"{d.mean():+.3f}",
                 f"[{np.percentile(bs, 2.5):+.3f}, {np.percentile(bs, 97.5):+.3f}]", f"{wilcoxon(x, y).pvalue:.1e}",
                 f"{int((d > 0).sum())}/{len(d)}"])


def main():
    L = lambda tag: T.load_condition(tag, "single", 20)[0]
    af3, emp, rnd, shf, rf = L("af3"), L("empty"), L("randomSHS"), L("shuffleSHS"), L("rnaformer")
    ids = sorted(set(E.split_ids("rna_monomers_ids")) & set(af3) & set(emp) & set(rnd) & set(shf) & set(rf))
    S = "6.1 controls"
    compare(S, "SHS (RNAformer) vs no alignment", rf, emp, ids)
    compare(S, "shuffled SHS vs no alignment", shf, emp, ids)
    compare(S, "random SHS vs no alignment", rnd, emp, ids)
    compare(S, "SHS (RNAformer) vs natural MSA", rf, af3, ids)
    src = {"RNAformer": rf, "RNAfold": L("rnafold"), "SPOT-RNA": L("spotrna")}
    rho = {"RNAformer": T.load_rhofold("rnaformerN100Rhofold"), "RNAfold": T.load_rhofold("rnafoldN100Rhofold"),
           "SPOT-RNA": T.load_rhofold("spotrnaN100Rhofold")}
    oids = sorted((set(ids) & set.intersection(*map(set, src.values())) & set.intersection(*map(set, rho.values()))) - {"9cxf"})
    for p, q in (("RNAformer", "RNAfold"), ("RNAformer", "SPOT-RNA"), ("SPOT-RNA", "RNAfold")):
        compare("6.2 ordering (AF3)", f"SHS ({p}) vs SHS ({q})", src[p], src[q], oids)
        compare("6.2 ordering (RhoFold+)", f"SHS ({p}) vs SHS ({q})", rho[p], rho[q], oids)
    sc = E.selfcons_means()
    hi = [i for i in ids if sc.get(i, -1) >= 0.9]; lo = [i for i in ids if i in sc and sc[i] < 0.9]
    compare("6.2 self-consistency >= 0.9", "SHS vs natural MSA", rf, af3, hi)
    compare("6.2 self-consistency >= 0.9", "SHS vs no alignment", rf, emp, hi)
    compare("6.2 self-consistency < 0.9", "SHS vs no alignment", rf, emp, lo)
    b0, dn, gt = L("rnafold"), L("denoise_rnafold"), L("dssr")
    did = sorted(set(ids) & set(b0) & set(dn) & set(gt))
    compare("6.3 denoising (AF3)", "denoised vs SHS (RNAfold)", dn, b0, did)
    compare("6.3 denoising (AF3)", "reference topology vs SHS (RNAfold)", gt, b0, did)
    r0, rd, rg = (T.load_rhofold(t) for t in ("rnafoldN100Rhofold", "rnafoldN100DenoiseRhofold", "dssrN100Rhofold"))
    rid = sorted(set(ids) & set(r0) & set(rd) & set(rg))
    compare("6.3 denoising (RhoFold+)", "denoised vs SHS (RNAfold)", rd, r0, rid)
    compare("6.3 denoising (RhoFold+)", "reference topology vs SHS (RNAfold)", rg, r0, rid)
    c = pd.read_csv(ROOT / "results" / "analysis" / "construction" / "construct91_N100_summary.csv")
    tm = c.pivot(index="pid", columns="step", values="mean_tm"); rc = c.pivot(index="pid", columns="step", values="recall")
    tmd = {s: {i: {"RNA_TM": v} for i, v in tm[s].items()} for s in ("wc", "all")}
    compare("6.3 construction (91 RNAs, N=100)", f"full reference vs Watson-Crick only (mean recall {rc['wc'].mean():.2f} -> "
            f"{rc['all'].mean():.2f})", tmd["all"], tmd["wc"], sorted(tm.index))
    vals = {}
    for lab, prefixes in (("msa", ("af3", "af3msa")), ("ref", ("dssr",)), ("nomsa", ("nomsa",)), ("shuffle", ("shuffle",)),
                          ("random", ("random",))):
        acc = {}
        for p in prefixes:
            for db in glob.glob(f"{E.DB_LONG}/{p}_s*_m*.db"):
                for r in (dict(x) for x in E._rows(db)):
                    t = T.transform(r)
                    if t.get("RNA_TM") is not None: acc.setdefault(E.idof(r["FileName"]), []).append(t["RNA_TM"])
        vals[lab] = {i: {"RNA_TM": st.mean(v)} for i, v in acc.items()}
    lids = sorted(set.intersection(*[set(v) for v in vals.values()]))
    for lab, name in (("ref", "reference-topology SHS"), ("shuffle", "shuffled SHS"), ("random", "random SHS")):
        compare("6.4 RNAs > 200 nt", f"{name} vs no alignment", vals[lab], vals["nomsa"], lids)
    for split, idset, metrics, sec in (("protein_rna_a2021_ids", "complex", T.COMPLEX, "7.1 post-cutoff complexes"),
                                       ("protein_rna_a2021_orphan_ids", "complex", T.COMPLEX, "7.1 post-cutoff orphan complexes"),
                                       ("rna_monomers_a2021_orphan_ids", "single", T.SINGLE, "7.2 post-cutoff orphan single chains"),
                                       ("rna_monomers_orphan_ids", "single", T.SINGLE, "7.2 single chains without homologs")):
        Sx, Mx = T.load_condition("rnaformer", idset, 20)[0], T.load_condition("af3", idset, 20)[0]
        for m in metrics:
            sid = sorted(i for i in E.split_ids(split) if i in Sx and i in Mx and Sx[i].get(m) is not None and Mx[i].get(m) is not None)
            compare(sec, f"SHS (RNAformer) vs natural MSA, {m}", Sx, Mx, sid, m)
    out = ROOT / "results" / "tables" / "prose_numbers.md"
    head = ["Section", "Comparison (A vs B)", "n", "mean A", "mean B", "A - B", "95% CI", "Wilcoxon p", "A > B"]
    out.write_text("# Paired statistics quoted in the manuscript text\n\nGenerated by tables/prose_numbers.py. "
                   "Per-target means over 4 AF3 seeds x 5 samples (RhoFold+: one prediction).\n\n"
                   + "| " + " | ".join(head) + " |\n|" + "---|" * len(head) + "\n"
                   + "".join("| " + " | ".join(map(str, r)) + " |\n" for r in ROWS))
    print(out.read_text())


if __name__ == "__main__":
    main()
