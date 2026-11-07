#!/usr/bin/env python3
"""LaTeX tables that make_crc_tables.py does not produce (results/tables/<name>.{tex,csv}):

  gates_single.tex      tab:gates          single-chain TM by 2D MCC of the encoded RNAformer topology and by self-consistency
  gates_complex1v1.tex  tab:gates_complex  1v1-complex DockQ by 2D MCC of the encoded RNAformer topology
  seed_robustness.tex   tab:seed_robustness per-seed condition means, 4x5 mean, s.d. across seeds (from seed_robustness.csv)
  long_rna.tex          tab:long_rna       25 RNAs of 202-746 nt, five inputs (metric DBs in $DB_LONG)
  orphan_tests.tex      tab:orphan_tests   paired SHS vs natural-MSA tests on post-cutoff orphans (orphan_parity.csv)
  order_rna_single_chain.tex   tab:order_rna_single_chain   2D-quality ordering, AF3 (4x5) and RhoFold+ (single run), n=90
  rf_rna_single_chain_other_splits.tex   tab:rf_rna_single_chain_other_splits   RhoFold+ SHS sources on further splits
All tables are also appended to results/tables/all_tables.tex (section "extra tables"), the file to copy from.

Same protocol and loaders as make_crc_tables.py (per-target mean over 4 seeds x 5 samples; mean +- s.d. over targets).
Run after make_crc_tables.py, analysis/gates/build_rnaformer_2d_gt.py and analysis/benchmark/orphan_parity.py.
"""
import csv, glob, importlib.util, sqlite3, statistics as st, sys
from pathlib import Path
import numpy as np, pandas as pd
import math
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from paths import DB_LONG, MC_CSV, MD_CSV, OUT_ROOT, idof, selfcons_means, split_ids  # noqa: E402
_s = importlib.util.spec_from_file_location("crc_tables", ROOT / "tables" / "make_crc_tables.py")
T = importlib.util.module_from_spec(_s); _s.loader.exec_module(T)
OUT = ROOT / "results" / "tables"
PROTO = r"Mean $\pm$ s.d.\ over targets of the per-target mean over 4 AF3 seeds $\times$ 5 samples"
rng = np.random.default_rng(0)


TEX = []


def write(name, label, caption, header, rows, colspec, body=None):
    body = body or "\n".join(" & ".join(r) + r" \\" for r in rows)
    wide = len(colspec) >= 7
    tex = (f"\\begin{{table}}[t]\n\\centering\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
           + ("\\resizebox{\\textwidth}{!}{%\n" if wide else "") + f"\\begin{{tabular}}{{{colspec}}}\n"
           f"\\toprule\n{' & '.join(header)} \\\\\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}\n"
           + ("}\n" if wide else "") + "\\end{table}\n")
    (OUT / f"{name}.tex").write_text(tex); TEX.append(tex)
    pd.DataFrame(rows, columns=header).to_csv(OUT / f"{name}.csv", index=False)
    print(f"  {name}.tex ({len(rows)} rows)")


def sgn(x): return f"${x:+.3f}$"


def gates():
    mc = pd.read_csv(MC_CSV); mc["id"] = mc.exp_db_id.str.lower()
    for idset, group, metric, name, label, ylab in (
            ("single", "rna-only", "RNA_TM", "gates_single", "tab:gates", "TM-score"),
            ("complex", "1v1", "Dockq", "gates_complex1v1", "tab:gates_complex", "DockQ")):
        w = mc[mc.group == group].drop_duplicates("id").set_index("id")
        for tag, c in (("empty", "N"), ("rnaformer", "S"), ("af3", "M")):
            d, _ = T.load_condition(tag, idset, 20); w[c] = pd.Series({i: d[i][metric] for i in w.index if i in d})
        w = w.dropna(subset=["N", "S", "M"])
        subsets = [(f"2D MCC ${lab}$", w[(w.rf2d_pred >= lo) & (w.rf2d_pred < hi)]) for lab, lo, hi in
                   (("<0.5", -1.01, .5), ("0.5\\text{--}0.8", .5, .8), ("0.8\\text{--}0.95", .8, .95), (">0.95", .95, 1.01))]
        if idset == "single":
            sc = pd.Series(selfcons_means()); w["sc"] = sc
            subsets += [("Self-consistency $\\ge 0.9$", w[w.sc >= 0.9]), ("Self-consistency $<0.9$", w[w.sc < 0.9])]
        rows = []
        for lab, g in subsets:
            d = (g.S - g.N).values
            p = wilcoxon(g.S, g.N).pvalue if len(g) > 5 and np.any(d != 0) else float("nan")
            rows.append([lab, str(len(g)), f"{g.N.mean():.3f}", f"{g.S.mean():.3f}", f"{g.M.mean():.3f}", sgn(d.mean()),
                         "--" if np.isnan(p) else f"{p:.1e}"])
        what = "single-chain RNAs" if idset == "single" else "one-RNA-one-protein complexes"
        cap = (f"{ylab} of {what} ($n{{=}}{len(w)}$) by the 2D MCC of the encoded RNAformer topology against the DSSR reference "
               f"pairs (all base pairs)" + (" and by self-consistency (MCC between the encoded topology and the base pairs of the "
               "SHS-guided prediction; needs no ground truth)" if idset == "single" else "") +
               f". {PROTO}; $\\Delta$ = SHS $-$ no alignment, paired Wilcoxon $p$.")
        write(name, label, cap, ["Subset", "$n$", "No alignment", "SHS", "Natural MSA", "$\\Delta$", "$p$"], rows, "lcccccc")


def seed_robustness():
    d = pd.read_csv(OUT / "seed_robustness.csv")
    d = d[~((d.condition == "MSA removed baseline") & (d.idset == "complex"))]   # not used in the paper
    rows = [[T.DISPLAY.get(r.condition, r.condition), "single" if r.idset == "single" else "complex", "TM" if r.metric == "RNA_TM" else "DockQ", str(r.n)]
            + [f"{getattr(r, c):.3f}" for c in ("seed1_mean", "seed0_mean", "seed42_mean", "seed137_mean", "mean_4x5")]
            + [f"{r.sd_across_seed_means:.4f}", f"{r.mean_per_target_sd_across_samples:.3f}"] for r in d.itertuples()]
    cap = ("Seed robustness. Condition means per AF3 seed (each the mean over targets of the per-target mean over 5 samples), "
           "the 4$\\times$5 mean used throughout, the s.d.\\ of the four seed means, and the mean per-target s.d.\\ across the "
           "5 samples of a seed.")
    write("seed_robustness_table", "tab:seed_robustness", cap,
          ["Input", "Set", "Metric", "$n$", "Seed 1", "Seed 0", "Seed 42", "Seed 137", "4$\\times$5", "s.d.\\ seeds", "s.d.\\ samples"],
          rows, "lllcccccccc")


def long_rna():
    conds = [("Natural MSA", ("af3", "af3msa")), ("SHS (reference topology)", ("dssr",)), ("No alignment", ("nomsa",)),
             ("Shuffled SHS", ("shuffle",)), ("Random SHS", ("random",))]
    vals = {}
    for lab, prefixes in conds:
        acc = {}
        for p in prefixes:
            for db in glob.glob(f"{DB_LONG}/{p}_s*_m*.db"):
                for r in (dict(x) for x in _rows(db)):
                    t = T.transform(r); i = idof(r["FileName"])
                    for m in ("RNA_TM", "RNA_RMSD", "RNA_LDDT_allAtom"):
                        if t.get(m) is not None: acc.setdefault(i, {}).setdefault(m, []).append(t[m])
        vals[lab] = {i: {m: st.mean(v) for m, v in d.items()} for i, d in acc.items()}
    ids = sorted(set.intersection(*[set(v) for v in vals.values()]))
    rows = []
    for lab, _ in conds:
        rows.append([lab] + [f"{st.mean(vals[lab][i][m] for i in ids):.3f} $\\pm$ {st.stdev(vals[lab][i][m] for i in ids):.3f}"
                             for m in ("RNA_TM", "RNA_RMSD", "RNA_LDDT_allAtom")])
    tm = {lab: st.mean(vals[lab][i]["RNA_TM"] for i in ids) for lab, _ in conds}
    gap = (tm["SHS (reference topology)"] - tm["No alignment"]) / (tm["Natural MSA"] - tm["No alignment"])
    cap = (f"RNAs longer than 200 nt ($n{{=}}{len(ids)}$, 202--746\\,nt). SHS are generated from the reference (DSSR) topology; "
           f"shuffled and random SHS have the same depth. {PROTO}. Reference-topology SHS recover {gap * 100:.0f}\\% of the "
           f"TM-score gap between no alignment and natural MSA. RMSD on P atoms; TM-score on C3$'$ atoms; all-atom LDDT. "
           + _runtime_note())
    write("long_rna", "tab:long_rna", cap, ["AlphaFold~3 input", "TM-score ($\\uparrow$)", "RMSD (\\AA, $\\downarrow$)",
                                            "LDDT ($\\uparrow$)"], rows, "lccc")


def _runtime_note():
    from paths import TIMING_CSV
    r = list(csv.DictReader(open(TIMING_CSV)))
    frac = [float(x["orig_msa_s"]) / float(x["orig_total_s"]) for x in r]
    msa = [float(x["orig_msa_s"]) for x in r]; gen = [float(x["shs_gen_s"]) for x in r]
    return (f"Homolog search takes {st.mean(frac) * 100:.0f}\\% of AF3's end-to-end runtime on these RNAs "
            f"(mean {st.mean(msa):.0f}\\,s), SHS generation {st.mean(gen):.0f}\\,s.")


def _rows(db):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True); con.row_factory = sqlite3.Row
    return con.execute("SELECT * FROM pred_rna_rna").fetchall()


def orphan_tests():
    """paired SHS vs natural MSA on post-cutoff orphans, from the same per-target values as all tables."""
    boot = np.random.default_rng(0); rows = []
    for idset, split, metric, lab in (("single", "rna_monomers_a2021_orphan_ids", "RNA_TM", "Single-chain RNAs (TM-score)"),
                                      ("complex", "protein_rna_a2021_orphan_ids", "Dockq", "RNA-protein complexes (DockQ)")):
        M, _ = T.load_condition("af3", idset, 20); S, _ = T.load_condition("rnaformer", idset, 20)
        ids = sorted(i for i in split_ids(split) if i in M and i in S and M[i].get(metric) is not None and S[i].get(metric) is not None)
        m_ = np.array([M[i][metric] for i in ids]); s_ = np.array([S[i][metric] for i in ids]); d = s_ - m_
        bs = [boot.choice(d, len(d)).mean() for _ in range(5000)]
        rows.append([lab, str(len(ids)), f"{m_.mean():.3f}", f"{s_.mean():.3f}", sgn(d.mean()),
                     f"$[{np.percentile(bs, 2.5):+.3f}, {np.percentile(bs, 97.5):+.3f}]$", f"{wilcoxon(s_, m_).pvalue:.2f}"])
    cap = ("Post-cutoff orphans: RNAformer-seeded SHS vs.\\ natural MSA. Paired differences of per-target means "
           "(4 AF3 seeds $\\times$ 5 samples), 5000$\\times$ bootstrap 95\\% CI and Wilcoxon signed-rank $p$. Neither difference "
           "is significant; SHS and natural MSAs are on par in this regime.")
    write("orphan_tests", "tab:orphan_tests", cap, ["Subset", "$n$", "Natural MSA", "SHS", "$\\Delta$ (SHS $-$ MSA)",
                                                    "95\\% CI", "$p$"], rows, "lcccccc")

SRC = [("rnafold", "rnafoldN100Rhofold", "SHS (RNAfold)"), ("spotrna", "spotrnaN100Rhofold", "SHS (SPOT-RNA)"),
       ("rnaformer", "rnaformerN100Rhofold", "SHS (RNAformer)")]
M3 = ("RNA_TM", "RNA_RMSD", "RNA_LDDT_allAtom")
HEAD3 = ["TM-score ($\\uparrow$)", "RMSD (\\AA, $\\downarrow$)", "LDDT ($\\uparrow$)"]


def _block(data, ids, labels, with_sd=True):
    """rows of mean (+- s.d.) over ids; bold = best per metric within the block."""
    stats = [[(st.mean(d[i][m] for i in ids), st.stdev(d[i][m] for i in ids)) for m in M3] for d in data]
    best = [(min if m == "RNA_RMSD" else max)(s[k][0] for s in stats) for k, m in enumerate(M3)]
    out = []
    for lab, s in zip(labels, stats):
        cells = []
        for k, (mu, sd) in enumerate(s):
            c = f"{mu:.3f} $\\pm$ {sd:.3f}" if with_sd else f"{mu:.3f}"
            cells.append(f"\\textbf{{{c}}}" if abs(mu - best[k]) < 5e-4 else c)
        out.append([lab] + cells)
    return out


def order_table():
    af3 = [T.load_condition(t, "single", 20)[0] for t, _, _ in SRC]; rho = [T.load_rhofold(r) for _, r, _ in SRC]
    ids = sorted(i for i in split_ids("rna_monomers_ids") - {"9cxf"}
                 if all(i in d and all(d[i].get(m) is not None for m in M3) for d in af3 + rho))
    labels = [l for _, _, l in SRC]
    a, r = _block(af3, ids, labels), _block(rho, ids, labels)
    body = "\n".join(("\\multirow{3}{*}{" + p + "}" if k == 0 else "") + " & " + " & ".join(row) + r" \\"
                      for p, blk in (("AlphaFold~3", a), ("RhoFold+", r)) for k, row in enumerate(blk))
    body = body.replace("\n\\multirow{3}{*}{RhoFold+}", "\n\\midrule\n\\multirow{3}{*}{RhoFold+}")
    cap = (f"Cross-architecture 2D-quality ordering on single-chain RNAs ($n{{=}}{len(ids)}$; 9CXF excluded because the RhoFold+ "
           "relaxation failed for SPOT-RNA-seeded SHS). AlphaFold~3: mean $\\pm$ s.d.\\ over targets of the per-target mean over "
           "4 seeds $\\times$ 5 samples; RhoFold+: one deterministic prediction per target. TM-score on C3$'$ atoms "
           "(US-align); RMSD on P atoms; LDDT is all-atom RNA LDDT.")
    write("order_rna_single_chain", "tab:order_rna_single_chain", cap, ["3D predictor", "SHS source"] + HEAD3,
          [[p] + row for p, blk in (("AlphaFold 3", a), ("RhoFold+", r)) for row in blk], "llccc", body=body)


def rf_other_splits():
    rho = [T.load_rhofold(r) for _, r, _ in SRC]; labels = [l for _, _, l in SRC]; parts, rows = [], []
    for split, name in (("rna_monomers_b2021_orphan_ids", "Pre-cutoff orphan"), ("rna_monomers_b2021_non_orphan_ids", "Pre-cutoff non-orphan"),
                        ("rna_monomers_a2021_non_orphan_ids", "Post-cutoff non-orphan")):
        ids = sorted(i for i in split_ids(split) if all(i in d and all(d[i].get(m) is not None for m in M3) for d in rho))
        blk = _block(rho, ids, labels, with_sd=False)
        parts.append("\n".join(("\\multirow{3}{*}{" + f"{name} (n={len(ids)})" + "}" if k == 0 else "") + " & " + " & ".join(row)
                               + r" \\" for k, row in enumerate(blk)))
        rows += [[f"{name} (n={len(ids)})"] + row for row in blk]
    cap = ("RhoFold+ with SHS from different 2D sources on further single-chain splits (one deterministic prediction per target; "
           "mean over targets; bold = best per split).")
    write("rf_rna_single_chain_other_splits", "tab:rf_rna_single_chain_other_splits", cap, ["Split", "RhoFold+ input"] + HEAD3,
          rows, "llccc", body="\n\\midrule\n".join(parts))


STAGES = [("local_wc", "Local WC only"), ("init_wc", "Global WC topology"), ("step1_nc", "+ non-canonical"),
          ("step2_pk", "+ pseudoknots"), ("step3_multi", "+ multi-partner pairs"), ("final_all", "Full reference")]


def construction_table():
    f = ROOT / "results" / "dssr2d" / "construction_2d.tsv.gz"
    d = pd.read_csv(f, sep="\t")
    per_n = d.groupby(["target", "stage", "N"]).mcc_vs_reference.agg(["mean", "count"]).reset_index()
    if not (per_n["count"] == 20).all(): sys.exit("construction_2d: not all (target, stage, N) have 20 predictions")
    nN = per_n.N.nunique(); rows = []
    for t in ("4P8Z", "6WJS"):
        g = per_n[per_n.target == t]
        rows.append([t] + [f"{g[g.stage == s]['mean'].mean():.3f} $\\pm$ {g[g.stage == s]['mean'].std():.3f}" for s, _ in STAGES])
    Ns = ", ".join(str(int(x)) for x in sorted(per_n.N.unique()))
    cap = ("2D accuracy of the AF3 predictions along the construction sweep: MCC between the base pairs DSSR extracts from the "
           "predicted structure and the DSSR reference pairs, for SHS encoding progressively richer topologies. Per SHS depth, "
           f"mean over 4 AF3 seeds $\\times$ 5 samples; values are mean $\\pm$ s.d.\\ across the {nN} depths ($N\\in\\{{{Ns}\\}}$). "
           "Each nucleotide carries at most one encoded partner, so ``+ multi-partner pairs'' re-assigns partners rather than "
           "encoding base triples.")
    write("construction_stage_mcc", "tab:intervention_stage_secondary_structure_summary_all", cap,
          ["Target"] + [n for _, n in STAGES], rows, "l" + "c" * len(STAGES))


C91 = [("wc", "Nested Watson--Crick and G--U pairs"), ("nc", "+ nested non-canonical pairs"), ("pk", "+ pseudoknotted pairs"),
       ("all", "Full reference (incl.\\ multi-partner pairs)")]


def construct91_table():
    """91-RNA construction sweep (SHS N=100): cumulative topology steps built from the DSSR reference pairs."""
    c = pd.read_csv(ROOT / "results" / "analysis" / "construction" / "construct91_N100_summary.csv")
    if not (c.n == 20).all(): sys.exit("construct91_N100_summary: not all (RNA, step) have 20 predictions")
    tm = c.pivot(index="pid", columns="step", values="mean_tm"); rc = c.pivot(index="pid", columns="step", values="recall")
    rows = []
    for s, lab in C91:
        d = (tm[s] - tm["wc"]).values
        rows.append([lab, f"{rc[s].mean():.2f}", f"{tm[s].mean():.3f} $\\pm$ {tm[s].std():.3f}",
                     "--" if s == "wc" else sgn(d.mean()), "--" if s == "wc" else f"{wilcoxon(tm[s], tm['wc']).pvalue:.1e}",
                     "--" if s == "wc" else f"{int((d > 0).sum())}/{len(d)}"])
    cap = (f"Construction sweep on all {len(tm)} single-chain RNAs (SHS depth $N{{=}}100$). Topologies are built cumulatively "
           "from the DSSR reference pairs; recall is the fraction of reference pairs encoded (precision 1 by construction). "
           f"TM-score: {PROTO}. $\\Delta$TM and Wilcoxon $p$ are paired against the Watson--Crick step; the last column counts "
           "RNAs with a higher TM-score than at the Watson--Crick step.")
    write("construct91_steps", "tab:construct91", cap,
          ["Encoded topology", "Recall", "TM-score ($\\uparrow$)", "$\\Delta$TM", "$p$", "Improved"], rows, "lccccc")


BENCH = [("af3", "Natural MSA"), ("rnafold", "SHS (RNAfold)"), ("spotrna", "SHS (SPOT-RNA)"), ("rnaformer", "SHS (RNAformer)")]


def benchmark_orphans_table():
    """Main-text benchmark: natural MSA vs. SHS on post-cutoff orphans, single chains next to complexes."""
    panels = [("single", "rna_monomers_a2021_orphan_ids", T.SINGLE, "Single-chain RNAs"),
              ("complex", "protein_rna_a2021_orphan_ids", T.COMPLEX_RED, "RNA-protein complexes")]
    cols, ns, pmin = [], [], 1.0
    for idset, split, metrics, _ in panels:
        data = [T.load_condition(t, idset, 20)[0] for t, _ in BENCH]
        ids = sorted(i for i in split_ids(split) if all(i in d and all(d[i].get(m) is not None for m in metrics) for d in data))
        ns.append(len(ids))
        for m in metrics:
            v = [[d[i][m] for i in ids] for d in data]
            cols.append((m, [(st.mean(x), st.stdev(x)) for x in v]))
            pmin = min(pmin, wilcoxon(v[3], v[0]).pvalue)               # RNAformer-seeded SHS vs. natural MSA
    if pmin < 0.05: sys.exit("benchmark_orphans: a difference is significant; caption would be wrong")
    rows = []
    for k, (_, lab) in enumerate(BENCH):
        cells = []
        for m, stats in cols:
            means = [x[0] for x in stats]; best = min(means) if m in T.LOWER_BETTER else max(means)
            c = T.fmt(*stats[k]); cells.append(r"\textbf{" + c + "}" if abs(stats[k][0] - best) < 5e-4 else c)
        rows.append([lab] + cells)
    cap = ("AlphaFold~3 with natural MSA vs.\\ SHS on post-cutoff orphans, i.e., targets deposited after the AF3 and RNAformer "
           f"training cutoffs for which the AF3 pipeline retrieves no RNA homologs. {PROTO}; bold = best mean per metric. "
           f"Differences between natural MSA and RNAformer-seeded SHS are not statistically significant (paired Wilcoxon, all "
           f"$p{{\\ge}}{math.floor(pmin * 100) / 100:.2f}$; App.~\\ref{{app:orphan_tests}}). TM-score on C3$'$ atoms; RMSD on P atoms; "
           "LDDT is all-atom RNA LDDT; DockQ scores the RNA-protein interface; cRMSD over all polymer atoms after superposition, "
           "before outlier rejection; cLDDT is all-atom complex LDDT.")
    hdr = ([""] + [rf"\multicolumn{{3}}{{c}}{{{t} ($n{{=}}{n}$)}}" for (_, _, _, t), n in zip(panels, ns)])
    tex = ("\\begin{table}[t]\n\\centering\n\\caption{" + cap + "}\n\\label{tab:benchmark_orphans}\n\\resizebox{\\textwidth}{!}{%\n"
           "\\begin{tabular}{lcccccc}\n\\toprule\n" + " & ".join(hdr) + " \\\\\n\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"
           "AlphaFold~3 input & " + " & ".join(T.HEAD[m] for m, _ in cols) + " \\\\\n\\midrule\n"
           + "\n".join(" & ".join(r) + r" \\" for r in rows) + "\n\\bottomrule\n\\end{tabular}\n}\n\\end{table}\n")
    (OUT / "benchmark_orphans.tex").write_text(tex); TEX.append(tex)
    pd.DataFrame(rows, columns=["AlphaFold 3 input"] + [f"{m} ({'single' if k < 3 else 'complex'})" for k, (m, _) in enumerate(cols)]
                 ).to_csv(OUT / "benchmark_orphans.csv", index=False)
    print(f"  benchmark_orphans.tex (n={ns}, min p={pmin:.3f})")



def _rna_name(desc):
    d = desc.lower()
    for key, name in (("delta", "HDV ribozyme"), ("hdv", "HDV ribozyme"), ("alkylating", "self-alkylating ribozyme"),
                      ("cpeb3", "CPEB3 ribozyme"), ("p4-p6", "P4-P6 ribozyme domain")):
        if key in d: return name
    return desc


def ribozyme_table():
    """Protein-bound ribozymes that AF3 folds poorly without an RNA alignment (RNA TM-score < 0.5): natural MSA vs. SHS."""
    md = pd.read_csv(MD_CSV, sep=None, engine="python"); md["id"] = md.exp_db_id.str.lower()
    rib = md[md.RNADescription.fillna("").str.lower().str.contains("ribozyme")].set_index("id")
    S, M, N = (T.load_condition(t, "complex", 20)[0] for t in ("rnaformer", "af3", "empty"))
    ids = sorted(i for i in rib.index if all(i in d and d[i].get("RNA_TM") is not None for d in (S, M, N)) and N[i]["RNA_TM"] < 0.5)
    s_ = np.array([S[i]["RNA_TM"] for i in ids]); m_ = np.array([M[i]["RNA_TM"] for i in ids]); d = s_ - m_
    boot = np.random.default_rng(0); bs = [boot.choice(d, len(d)).mean() for _ in range(5000)]
    lo, hi, p = np.percentile(bs, 2.5), np.percentile(bs, 97.5), wilcoxon(s_, m_).pvalue
    rows = [[i.upper(), _rna_name(rib.RNADescription[i]), f"{M[i]['RNA_TM']:.3f}", f"{S[i]['RNA_TM']:.3f}"] for i in ids]
    body = "\n".join(" & ".join(r) + r" \\" for r in rows) + "\n\\midrule\n" + \
        f"Mean & & {m_.mean():.3f} & {s_.mean():.3f} \\\\"
    cap = (f"Protein-bound ribozymes for which AF3 without an RNA alignment reaches an RNA TM-score below 0.5 ($n{{=}}{len(ids)}$ of "
           f"{len(rib)} ribozymes in the complex set). RNA TM-score, per-target mean over 4 AF3 seeds $\\times$ 5 samples. "
           f"SHS $-$ natural MSA: ${d.mean():+.3f}$, bootstrap 95\\% CI $[{lo:+.3f}, {hi:+.3f}]$, paired Wilcoxon $p{{=}}{p:.2f}$. "
           "Post hoc subset; reported as a hypothesis, not as evidence.")
    write("ribozymes", "tab:ribozymes", cap, ["PDB ID", "RNA", "Natural MSA", "SHS (RNAformer)"], rows, "llcc", body=body)
    print(f"  ribozymes: n={len(ids)} SHS {s_.mean():.3f} MSA {m_.mean():.3f} d={d.mean():+.3f} [{lo:+.3f},{hi:+.3f}] p={p:.3f}")


def confidence_table():
    import json, os
    from paths import prediction_files, IDS_DIR
    from scipy.stats import spearmanr
    singles = {x.strip().lower() for x in open(IDS_DIR / "single_91.txt") if x.strip()}

    def conf(tag):
        acc = {}
        for f in prediction_files(tag, "single", "summary_confidences.json"):
            i = idof(os.path.basename(os.path.dirname(os.path.dirname(f))))
            if not i or i not in singles: continue
            v = json.load(open(f)).get("ptm")
            if v is not None: acc.setdefault(i, []).append(float(v))
        return {i: st.mean(x) for i, x in acc.items()}
    conds = [("af3", "Natural MSA"), ("rnaformer", "SHS (RNAformer)"), ("empty", "No alignment"),
             ("shuffleSHS", "Shuffled SHS"), ("randomSHS", "Random SHS")]
    C = {t: conf(t) for t, _ in conds}
    ids = sorted(set.intersection(*[set(v) for v in C.values()]))
    rows = []
    for t, lab in conds:
        x = np.array([C[t][i] for i in ids]); base = np.array([C["empty"][i] for i in ids])
        if t == "empty": rows.append([lab, f"{x.mean():.3f} $\\pm$ {x.std(ddof=1):.3f}", "--", "--"]); continue
        rows.append([lab, f"{x.mean():.3f} $\\pm$ {x.std(ddof=1):.3f}", sgn((x - base).mean()), f"{wilcoxon(x, base).pvalue:.1e}"])
    tm, _ = T.load_condition("rnaformer", "single", 20)
    r = spearmanr([C["rnaformer"][i] for i in ids if i in tm], [tm[i]["RNA_TM"] for i in ids if i in tm])
    cap = (f"AF3 confidence (pTM) on single-chain RNAs ($n{{=}}{len(ids)}$) by alignment input. Per-target mean over 4 AF3 seeds "
           f"$\\times$ 5 samples; mean $\\pm$ s.d.\\ over targets; $\\Delta$ and paired Wilcoxon $p$ relative to no alignment. "
           f"With SHS, pTM tracks the TM-score of the prediction (Spearman $\\rho{{=}}{r.correlation:.2f}$).")
    write("confidence_ptm", "tab:confidence", cap, ["AlphaFold~3 input", "pTM", "$\\Delta$ vs.\\ no alignment", "$p$"], rows, "lccc")


def near_floor_table():
    d = pd.read_csv(Path(OUT_ROOT) / "near_floor" / "near_floor.csv")
    rows = [[("all" if r.floor == 0 else f"$\\ge {r.floor:.2f}$"), str(int(r.n)), f"{r.MSA:.3f}", f"{r.SHS:.3f}", sgn(r.d_SHS_MSA),
             f"{r.frac_SHS_ge_MSA:.0f}\\%"] for r in d.itertuples()]
    cap = ("Post-cutoff single-chain orphans restricted to solvable targets (best TM-score of any input above the floor). For "
           "orphans the natural-MSA input contains no homologs. TM-score, per-target mean over 4 AF3 seeds $\\times$ 5 samples; "
           "last column: share of targets where SHS reach at least the natural-MSA TM-score.")
    write("near_floor", "tab:near_floor", cap, ["Floor", "$n$", "Natural MSA", "SHS", "$\\Delta$", "SHS $\\ge$ MSA"], rows, "lccccc")


def novelty_table():
    import json
    seq = pd.read_csv(Path(OUT_ROOT) / "novelty" / "sequence_similarity.csv")
    rows = []
    for sub, lab in (("single_chain", "Single-chain RNAs"), ("complex", "RNA-protein complexes")):
        s8 = seq[(seq.subtype == sub) & (seq.threshold.round(2) == 0.80)].iloc[0]
        s9 = seq[(seq.subtype == sub) & (seq.threshold.round(2) == 0.90)].iloc[0]
        tm = np.array(list(json.load(open(Path(OUT_ROOT) / "novelty" / f"structure_maxTM_{sub}.json")).values()), float)
        rows.append([lab, str(int(s8.n_post)), f"{int(s9.novel)}/{int(s9.n_post)}", f"{int(s8.novel)}/{int(s8.n_post)}",
                     f"{np.median(tm):.2f}", f"{int((tm > 0.5).sum())}/{len(tm)}", f"{int((tm > 0.7).sum())}/{len(tm)}"])
    cap = ("Novelty of the post-cutoff targets relative to the pre-cutoff set. Sequence: targets without a pre-cutoff RNA "
           "sequence at $\\ge$90\\% or $\\ge$80\\% identity (CD-HIT-EST-2D). Structure: maximum TM-score (US-align) of each "
           "post-cutoff reference structure to any pre-cutoff reference structure.")
    write("novelty", "tab:novelty", cap, ["Set", "$n$", "Seq.-novel (90\\%)", "Seq.-novel (80\\%)", "Median max-TM",
                                          "max-TM $>0.5$", "max-TM $>0.7$"], rows, "lcccccc")


def append_to_all_tables():
    f = OUT / "all_tables.tex"; base = f.read_text().split("% ===== extra tables")[0].rstrip("\n")
    f.write_text(base + "\n\n% ===== extra tables (tables/make_extra_tables.py) =====\n" + "\n".join(TEX))


if __name__ == "__main__":
    gates(); seed_robustness(); long_rna(); orphan_tests(); order_table(); rf_other_splits()
    confidence_table(); near_floor_table(); novelty_table()
    if (ROOT / "results" / "dssr2d" / "construction_2d.tsv.gz").exists(): construction_table()
    else: print("  construction_stage_mcc: results/dssr2d/construction_2d.tsv.gz missing (run eval/dssr2d) -> skipped")
    construct91_table()
    benchmark_orphans_table()
    ribozyme_table()
    append_to_all_tables()
