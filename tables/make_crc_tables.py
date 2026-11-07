#!/usr/bin/env python3
"""Regenerate every AF3 results table of the paper from the per-shard metric DBs.

Protocol: each target is predicted with 4 AF3 seeds (1, 0, 42, 137) x 5 diffusion samples; a target's
score is the mean over its 20 predictions; tables report mean +- s.d. over targets. Targets are the paper
split (data/splits/neurips_ids.py) restricted to complete cases (every condition in the table has a value
for every metric). Metric transforms: RMSD lists -> first element; complex metrics fall back to the RNA metric
for single chains.

  python tables/make_crc_tables.py                 # all tables -> results/tables/{<label>.tex,.csv}, all_tables.tex
  python tables/make_crc_tables.py --min-preds 20  # require the full 4x5 protocol per target (default: 20)
RhoFold+ rows: one deterministic run per input (no seeds), scored by the same engine ($DB_RHOFOLD/<tag>_single_s1_m0.db).
Complex RMSD (cRMSD) of DB-backed conditions is read from $CRMSD_DIR (eval/crmsd): all polymer atoms, RMSD before
PyMOL outlier rejection; the metric DBs store the post-rejection value in Complex_RMSD[0].
"""
import argparse, ast, csv, glob, gzip, math, os, sqlite3, statistics as st, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
from paths import CRC, CRMSD_DIR, DB_DIR, DB_RHOFOLD, SEEDS, idof, split_ids  # noqa: E402

SINGLE = ["RNA_TM", "RNA_RMSD", "RNA_LDDT_allAtom"]
COMPLEX = ["Dockq", "ICS", "QS_global", "Complex_TM", "Complex_RMSD", "Complex_LDDT_allAtom"]
COMPLEX_RED = ["Dockq", "Complex_RMSD", "Complex_LDDT_allAtom"]
LOWER_BETTER = {"RNA_RMSD", "Complex_RMSD"}
HEAD = {"RNA_TM": r"TM-score ($\uparrow$)", "RNA_RMSD": r"RMSD (\AA, $\downarrow$)", "RNA_LDDT_allAtom": r"LDDT ($\uparrow$)",
        "Dockq": r"DockQ ($\uparrow$)", "ICS": r"ICS ($\uparrow$)", "QS_global": r"QS-global ($\uparrow$)",
        "Complex_TM": r"cTM ($\uparrow$)", "Complex_RMSD": r"cRMSD (\AA, $\downarrow$)", "Complex_LDDT_allAtom": r"cLDDT ($\uparrow$)"}
LABEL = {r["tag"]: r["paper_label"] for r in csv.DictReader(open(CRC / "eval" / "conditions.tsv"), delimiter="\t")}
MAIN4 = ["af3", "rnafold", "spotrna", "rnaformer"]
# LaTeX display names (CSV outputs keep the conditions.tsv paper_label); per-table overrides below
DISPLAY = {"MSA": "Natural MSA", "MSA removed baseline": "No alignment", "SHS (Random)": "Random SHS",
           "SHS (Shuffled)": "Shuffled SHS", "SHS (GT)": "SHS (reference topology)", "SHS (RNAfold Denoise)": "SHS (RNAfold, denoised)",
           "RhoFold+ native MSA": "Native RhoFold+ MSA"}
for _n in (20, 50, 100, 200, 500, 1000):
    DISPLAY[f"SHS (RNAformerN{_n}) + MSA"] = f"SHS (RNAformer, $N{{=}}{_n}$) + natural MSA"
for _n in (20, 500, 1000, 5000, 10000, 15000, 20000):
    DISPLAY[f"SHS (RNAformer, N={_n})"] = f"SHS (RNAformer, $N{{=}}{_n}$)"
DISPLAY_OVERRIDE = {t: {"SHS (RNAformer)": "SHS (RNAformer, $N{=}100$)"} for t in
                    ("rna_single_chain_scaling", "rna_protein_complexes_scaling", "shs_and_msa_depth_sweep")}
METRIC_NOTE = {"RNA_RMSD": r"RMSD on P atoms", "RNA_TM": r"TM-score on C3$'$ atoms (US-align)", "RNA_LDDT_allAtom": r"LDDT is all-atom RNA LDDT",
               "Dockq": r"DockQ, ICS and QS-global score the RNA-protein interface", "Complex_TM": r"cTM is the complex TM-score",
               "Complex_RMSD": r"cRMSD over all polymer atoms after superposition, before outlier rejection",
               "Complex_LDDT_allAtom": r"cLDDT is all-atom complex LDDT"}


def metric_note(metrics):
    parts = [METRIC_NOTE[m] for m in metrics if m in METRIC_NOTE]
    if "Dockq" in metrics and not ({"ICS", "QS_global"} & set(metrics)):
        parts = [x.replace("DockQ, ICS and QS-global score", "DockQ scores") for x in parts]
    return "; ".join(parts) + "."

# (label, split ids name, idset, metrics, rows, caption); a row is a condition tag or a RhoFold+ row ("rhofold", tag, display)
def rf(tag, name): return ("rhofold", tag, name)
SINGLE_SPLITS = [("rna_single_chain", "rna_monomers_ids", "RNA single chain (n={n})"),
                 ("rna_single_chain_post_cutoff", "rna_monomers_a2021_ids", "RNA single chain post-cutoff (n={n})"),
                 ("rna_single_chain_without_homologs", "rna_monomers_orphan_ids", "RNA single chain without homologs (n={n})"),
                 ("rna_single_chain_with_homologs", "rna_monomers_non_orphan_ids", "RNA single chain with homologs (n={n})"),
                 ("rna_single_chain_pre_cutoff", "rna_monomers_b2021_ids", "RNA single chain pre-cutoff (n={n})"),
                 ("rna_single_chain_pre_cutoff_without_homologs", "rna_monomers_b2021_orphan_ids", "RNA single chain pre-cutoff without homologs (n={n})"),
                 ("rna_single_chain_pre_cutoff_with_homologs", "rna_monomers_b2021_non_orphan_ids", "RNA single chain pre-cutoff with homologs (n={n})"),
                 ("rna_single_chain_post_cutoff_without_homologs", "rna_monomers_a2021_orphan_ids", "RNA single chain post-cutoff without homologs (n={n})"),
                 ("rna_single_chain_post_cutoff_with_homologs", "rna_monomers_a2021_non_orphan_ids", "RNA single chain post-cutoff with homologs (n={n})")]
COMPLEX_SPLITS = [("rna_protein_complexes", "protein_rna_ids", "RNA-protein complexes (n={n})"),
                  ("rna_protein_complexes_without_homologs", "protein_rna_orphan_ids", "RNA-protein complexes without homologs (n={n})"),
                  ("rna_protein_complexes_with_homologs", "protein_rna_non_orphan_ids", "RNA-protein complexes with homologs (n={n})"),
                  ("rna_protein_complexes_pre_cutoff", "protein_rna_b2021_ids", "RNA-protein complexes pre-cutoff (n={n})"),
                  ("rna_protein_complexes_pre_cutoff_without_homologs", "protein_rna_b2021_orphan_ids", "RNA-protein complexes pre-cutoff without homologs (n={n})"),
                  ("rna_protein_complexes_pre_cutoff_with_homologs", "protein_rna_b2021_non_orphan_ids", "RNA-protein complexes pre-cutoff with homologs (n={n})"),
                  ("rna_protein_complexes_post_cutoff", "protein_rna_a2021_ids", "RNA-protein complexes post-cutoff (n={n})"),
                  ("rna_protein_complexes_post_cutoff_without_homologs", "protein_rna_a2021_orphan_ids", "RNA-protein complexes post-cutoff without homologs (n={n})"),
                  ("rna_protein_complexes_post_cutoff_with_homologs", "protein_rna_a2021_non_orphan_ids", "RNA-protein complexes post-cutoff with homologs (n={n})"),
                  ("rna_protein_1v1_complexes", "protein_rna_1v1_ids", "RNA-protein 1v1 complexes (n={n})"),
                  ("rna_protein_multichain_complexes", "protein_rna_multichain_ids", "RNA-protein multichain complexes (n={n})")]
RHO = [rf("rnafoldN100Rhofold", "SHS (RNAfold)"), rf("spotrnaN100Rhofold", "SHS (SPOT-RNA)"), rf("rnaformerN100Rhofold", "SHS (RNAformer)")]
RHO_NATIVE = RHO + [rf("stdmsaRhofold", "RhoFold+ native MSA")]
DEPTHS = ["rnaformerN20", "rnaformer", "rnaformerN500", "rnaformerN1000", "rnaformerN5000", "rnaformerN10000", "rnaformerN15000", "rnaformerN20000"]
TABLES = (
    [(l, s, "single", SINGLE, MAIN4, c) for l, s, c in SINGLE_SPLITS]
    + [(l, s, "complex", COMPLEX, MAIN4, c) for l, s, c in COMPLEX_SPLITS]
    + [("rna_protein_complexes_post_cutoff_without_homologs_reduced", "protein_rna_a2021_orphan_ids", "complex", COMPLEX_RED, MAIN4,
        "RNA-protein complexes post-cutoff without homologs (n={n}); main-text metrics"),
       ("rna_single_chain_shuffle", "rna_monomers_ids", "single", SINGLE, ["empty", "randomSHS", "shuffleSHS", "rnaformer"],
        "Equal-depth negative controls, RNA single chain (n={n})"),
       ("rna_single_chain_controls", "rna_monomers_ids", "single", SINGLE, ["empty", "randomSHS", "shuffleSHS", "rnaformer"],
        "Equal-depth controls on single-chain RNAs (n={n}). All SHS variants use $N{{=}}100$"),
       ("order_rna_single_chain_af3", "rna_monomers_ids", "single", SINGLE, ["rnafold", "spotrna", "rnaformer"],
        "Cross-architecture ordering, AlphaFold 3 part (n={n}; 9CXF excluded as in the RhoFold+ part)"),
       ("denoise_main", "rna_monomers_ids", "single", SINGLE, ["rnafold", "denoise_rnafold", "dssr"],
        "AF3 denoising ladder, RNA single chain (n={n})"),
       ("shs_and_msa_depth_sweep", "rna_monomers_ids", "single", SINGLE,
        ["af3", "rnaformer", "synMSA_N20", "synMSA_N50", "synMSA_N100", "synMSA_N200", "synMSA_N500", "synMSA_N1000"],
        "SHS combined with natural MSA, RNA single chain (n={n})"),
       ("rna_single_chain_scaling", "rna_monomers_ids", "single", SINGLE, ["af3"] + DEPTHS,
        "SHS depth scaling (RNAformer SHS, N=20 to 20000), RNA single chain (n={n})"),
       ("rna_protein_complexes_scaling", "protein_rna_ids", "complex", COMPLEX, ["af3"] + DEPTHS,
        "SHS depth scaling (RNAformer SHS, N=20 to 20000), RNA-protein complexes (n={n})")]
    + [("rf_" + l, s, "single", SINGLE, RHO, "RhoFold+, " + c + "; single run")
       for l, s, c in SINGLE_SPLITS if s in ("rna_monomers_ids", "rna_monomers_orphan_ids", "rna_monomers_non_orphan_ids",
                                              "rna_monomers_b2021_ids", "rna_monomers_a2021_ids", "rna_monomers_a2021_orphan_ids",
                                              "rna_monomers_b2021_orphan_ids", "rna_monomers_b2021_non_orphan_ids", "rna_monomers_a2021_non_orphan_ids")]
    + [("denoise_main_rf", "rna_monomers_ids", "single", SINGLE,
        [rf("rnafoldN100Rhofold", "SHS (RNAfold)"), rf("rnafoldN100DenoiseRhofold", "SHS (RNAfold Denoise)"), rf("dssrN100Rhofold", "SHS (GT)")],
        "RhoFold+ denoising ladder, RNA single chain (n={n}); single run"),
       ("rf_rna_single_chain_native_msa", "rna_monomers_ids", "single", SINGLE, RHO_NATIVE,
        "RhoFold+ SHS variants vs. the native RhoFold+ MSA pipeline, RNA single chain (n={n}); single run")]
)
EXCLUDE = {"order_rna_single_chain_af3": {"9cxf"}}


def num(x):
    if x is None or x == "": return None
    if isinstance(x, (int, float)): v = float(x)
    else:
        try: v = ast.literal_eval(str(x))
        except Exception: return None
        if isinstance(v, (list, tuple)): v = v[0] if v else None
        try: v = float(v)
        except (TypeError, ValueError): return None
    return None if math.isnan(v) else v


def transform(row):
    """Metric preprocessing: first RMSD element; complex metrics fall back to RNA metrics."""
    out = {m: num(row.get(m)) for m in set(SINGLE + COMPLEX)}
    if out["Complex_RMSD"] is None: out["Complex_RMSD"] = out["RNA_RMSD"]
    if out["Complex_TM"] is None: out["Complex_TM"] = out["RNA_TM"]
    if out["RNA_TM"] == 0.0: out["RNA_TM"] = None   # US-align failure sentinel
    return out


def crmsd_values(db):
    """{id: paper cRMSD} for the shard of metric DB `db` (eval/crmsd output); a missing file is an error, so DB and
    paper cRMSD definitions are never mixed in one table."""
    p = Path(CRMSD_DIR) / (Path(db).stem + ".tsv.gz")
    if not p.is_file():
        raise FileNotFoundError(f"{p}: run eval/crmsd/run_crmsd_array.sh for this shard")
    return {r["id"]: float(r["rmsd_all"]) for r in csv.DictReader(gzip.open(p, "rt"), delimiter="\t") if r["status"] == "ok"}


_cache = {}
def load_condition(tag, idset, min_preds):
    """{id: {metric: mean over seeds x samples}} for a DB-backed condition."""
    key = (tag, idset)
    if key in _cache: return _cache[key]
    tbl = "pred_rna_rna" if idset == "single" else "pred_protein_rna"
    acc, npred = {}, {}
    for seed in SEEDS:
        for db in sorted(glob.glob(f"{DB_DIR}/{tag}_{idset}_s{seed}_m*.db")):
            crmsd = crmsd_values(db) if idset == "complex" else None
            con = sqlite3.connect(f"file:{db}?mode=ro", uri=True); con.row_factory = sqlite3.Row
            for r in con.execute(f"SELECT * FROM {tbl}"):
                i = idof(r["FileName"]); t = transform(dict(r))
                if crmsd is not None: t["Complex_RMSD"] = crmsd.get(i)
                npred[i] = npred.get(i, 0) + 1
                for m, v in t.items():
                    if v is not None: acc.setdefault(i, {}).setdefault(m, []).append(v)
    res = {i: {m: st.mean(v) for m, v in d.items()} for i, d in acc.items() if npred.get(i, 0) >= min_preds}
    _cache[key] = (res, npred)
    return _cache[key]


def load_rhofold(tag):
    """{id: metrics} of one RhoFold+ run (single deterministic prediction per target)."""
    con = sqlite3.connect(f"file:{DB_RHOFOLD}/{tag}_single_s1_m0.db?mode=ro", uri=True); con.row_factory = sqlite3.Row
    return {idof(r["FileName"]): transform(dict(r)) for r in con.execute("SELECT * FROM pred_rna_rna")}


def fmt(m, s): return f"{m:.3f} $\\pm$ {s:.3f}"


def build(label, split, idset, metrics, rows, caption, min_preds):
    data, names = [], []
    for r in rows:
        if isinstance(r, tuple): d = load_rhofold(r[1]); names.append(r[2])
        else: d, _ = load_condition(r, idset, min_preds); names.append(LABEL[r])
        data.append(d)
    ids = split_ids(split) - EXCLUDE.get(label, set())
    ids = sorted(i for i in ids if all(i in d and all(d[i].get(m) is not None for m in metrics) for d in data))
    stats = [[(st.mean(d[i][m] for i in ids), st.stdev(d[i][m] for i in ids)) if len(ids) > 1 else (float("nan"),) * 2
              for m in metrics] for d in data]
    best = [min if metrics[k] in LOWER_BETTER else max for k in range(len(metrics))]
    best = [best[k](s[k][0] for s in stats) for k in range(len(metrics))]
    out = Path(CRC / "results" / "tables"); out.mkdir(parents=True, exist_ok=True)
    with open(out / f"{label}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["condition", "n"] + [x for m in metrics for x in (m + "_mean", m + "_sd")])
        for nm, s in zip(names, stats): w.writerow([nm, len(ids)] + [f"{v:.6f}" for ms in s for v in ms])
    cols = "l" + "c" * len(metrics)
    rho = all(isinstance(r, tuple) for r in rows)
    disp = [DISPLAY_OVERRIDE.get(label, {}).get(nm, DISPLAY.get(nm, nm)) for nm in names]
    note = metric_note(metrics)
    tex = [r"\begin{table}[htbp]", r"\centering",
           r"\caption{" + caption.format(n=len(ids)) + (r". Mean $\pm$ s.d. over targets (one deterministic RhoFold+ prediction per target)"
           if all(isinstance(r, tuple) for r in rows) else r". Mean $\pm$ s.d. over targets of the per-target mean over 4 AF3 seeds $\times$ 5 samples")
           + "; bold = best mean per metric. " + note + "}",
           r"\label{tab:" + label + "}"] + ([r"\resizebox{\textwidth}{!}{%"] if len(metrics) >= 6 else []) + [
           r"\begin{tabular}{" + cols + "}", r"\toprule",
           ("RhoFold+ input & " if rho else "AlphaFold~3 input & ") + " & ".join(HEAD[m] for m in metrics) + r" \\", r"\midrule"]
    for nm, s in zip(disp, stats):
        cells = [(r"\textbf{" + fmt(*ms) + "}") if abs(ms[0] - best[k]) < 5e-4 else fmt(*ms) for k, ms in enumerate(s)]
        tex.append(nm + " & " + " & ".join(cells) + r" \\")
    tex += [r"\bottomrule", r"\end{tabular}"] + (["}"] if len(metrics) >= 6 else []) + [r"\end{table}", ""]
    (out / f"{label}.tex").write_text("\n".join(tex))
    return "\n".join(tex), len(ids)


def seed_robustness(min_preds):
    """Per-seed means, the 4x5 mean, the spread across seed means and across samples; single RNA_TM / complex DockQ."""
    rows = []
    for tag in ("af3", "empty", "rnafold", "spotrna", "rnaformer", "dssr"):
        for idset, metric in (("single", "RNA_TM"), ("complex", "Dockq")):
            if tag == "dssr" and idset == "complex": continue
            ids = split_ids("rna_monomers_ids" if idset == "single" else "protein_rna_ids")
            per_seed = {}
            for seed in SEEDS:
                vals = {}
                for db in glob.glob(f"{DB_DIR}/{tag}_{idset}_s{seed}_m*.db"):
                    tbl = "pred_rna_rna" if idset == "single" else "pred_protein_rna"
                    for fn, v in sqlite3.connect(f"file:{db}?mode=ro", uri=True).execute(f"SELECT FileName,{metric} FROM {tbl}"):
                        v = num(v); i = idof(fn)
                        if v is not None and i in ids and not (metric == "RNA_TM" and v == 0): vals.setdefault(i, []).append(v)
                per_seed[seed] = vals
            common = sorted(set.intersection(*[set(v) for v in per_seed.values() if v] or [set()]))
            if not common: continue
            seed_means = {s: st.mean(st.mean(per_seed[s][i]) for i in common) for s in SEEDS if per_seed[s]}
            all_mean = st.mean(st.mean(sum((per_seed[s].get(i, []) for s in SEEDS), [])) for i in common)
            samp_sd = st.mean(st.pstdev(per_seed[s][i]) for s in SEEDS for i in common if len(per_seed[s].get(i, [])) > 1)
            rows.append([LABEL[tag], idset, metric, len(common)] + [seed_means.get(s, float("nan")) for s in SEEDS]
                        + [all_mean, st.pstdev(seed_means.values()), samp_sd])
    out = CRC / "results" / "tables" / "seed_robustness.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["condition", "idset", "metric", "n"] + [f"seed{s}_mean" for s in SEEDS]
                   + ["mean_4x5", "sd_across_seed_means", "mean_per_target_sd_across_samples"])
        w.writerows(rows)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--min-preds", type=int, default=20)
    a = ap.parse_args()
    alltex = []
    for spec in TABLES:
        tex, n = build(*spec, a.min_preds); alltex.append(tex); print(f"{spec[0]:60s} n={n}")
    (CRC / "results" / "tables" / "all_tables.tex").write_text("\n".join(alltex))
    for r in seed_robustness(a.min_preds):
        print("robustness", r[0], r[1], f"4x5={r[8]:.3f} sd_seeds={r[9]:.4f}")
