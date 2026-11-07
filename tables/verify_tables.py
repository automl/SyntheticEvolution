#!/usr/bin/env python3
"""Independent check of every generated table (results/tables/*.tex) and of the copy-paste folder.

The expected numbers are recomputed here with separate code directly from the raw sources (AF3 metric DBs + cRMSD
files, RhoFold+ DBs, long-RNA DBs, AF3 summary_confidences.json, self-consistency files, the DSSR sweep file, analysis
CSVs) and compared with the numbers printed in each .tex cell (same rounding). Also checks that every file in the
copy-paste folder is byte-identical to the repo tables.

  python tables/verify_tables.py --dest <copy-paste dir>     -> <dest>/VERIFY_REPORT.md, exit code 1 on any mismatch
"""
import argparse, ast, csv, glob, gzip, importlib.util, json, math, os, re, sqlite3, statistics as st, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from paths import (CRMSD_DIR, DB_DIR, DB_LONG, DB_RHOFOLD, MC_CSV, MD_CSV, OUT_ROOT, OUTPUTS, SELFCONS_DIR, TIMING_CSV,  # noqa: E402
                   split_ids)
_s = importlib.util.spec_from_file_location("crc_tables", ROOT / "tables" / "make_crc_tables.py")
T = importlib.util.module_from_spec(_s); _s.loader.exec_module(T)   # only TABLES / EXCLUDE (table definitions) are used
TAB = ROOT / "results" / "tables"
SINGLE = ["RNA_TM", "RNA_RMSD", "RNA_LDDT_allAtom"]
REPORT = []


# ---------------------------------------------------------------- independent loaders
def pid(s):
    m = re.search(r"([0-9][A-Za-z0-9]{3})", s or ""); return m.group(1).lower() if m else None


def val(x):
    if x is None or x == "": return None
    if isinstance(x, str):
        try: x = ast.literal_eval(x)
        except Exception: return None
    if isinstance(x, (list, tuple)): x = x[0] if x else None
    try: x = float(x)
    except (TypeError, ValueError): return None
    return None if math.isnan(x) else x


_af3 = {}
def af3(tag, idset):
    """{id: {metric: mean over all 20 predictions}} read straight from the metric DBs (+ cRMSD files for complexes)."""
    if (tag, idset) in _af3: return _af3[(tag, idset)]
    tbl = "pred_rna_rna" if idset == "single" else "pred_protein_rna"; acc, cnt = {}, {}
    for db in sorted(glob.glob(f"{DB_DIR}/{tag}_{idset}_s*_m*.db")):
        cr = {}
        if idset == "complex":
            for r in csv.DictReader(gzip.open(Path(CRMSD_DIR) / (Path(db).stem + ".tsv.gz"), "rt"), delimiter="\t"):
                if r["status"] == "ok": cr[r["id"]] = float(r["rmsd_all"])
        con = sqlite3.connect(f"file:{db}?mode=ro", uri=True); con.row_factory = sqlite3.Row
        for r in con.execute(f"SELECT * FROM {tbl}"):
            r = dict(r); i = pid(r["FileName"]); cnt[i] = cnt.get(i, 0) + 1
            v = {m: val(r.get(m)) for m in SINGLE + T.COMPLEX}
            if idset == "complex": v["Complex_RMSD"] = cr.get(i)
            if v["Complex_RMSD"] is None: v["Complex_RMSD"] = v["RNA_RMSD"]
            if v["Complex_TM"] is None: v["Complex_TM"] = v["RNA_TM"]
            if v["RNA_TM"] == 0.0: v["RNA_TM"] = None
            for m, x in v.items():
                if x is not None: acc.setdefault(i, {}).setdefault(m, []).append(x)
    # correctly rounded mean (statistics.mean): exact ties between inputs must stay exact zeros for Wilcoxon
    _af3[(tag, idset)] = {i: {m: st.mean(x) for m, x in d.items()} for i, d in acc.items() if cnt[i] >= 20}
    return _af3[(tag, idset)]


def rho(tag):
    con = sqlite3.connect(f"file:{DB_RHOFOLD}/{tag}_single_s1_m0.db?mode=ro", uri=True); con.row_factory = sqlite3.Row
    out = {}
    for r in con.execute("SELECT * FROM pred_rna_rna"):
        r = dict(r); v = {m: val(r.get(m)) for m in SINGLE}
        if v["RNA_TM"] == 0.0: v["RNA_TM"] = None
        out[pid(r["FileName"])] = {m: x for m, x in v.items() if x is not None}
    return out


# ---------------------------------------------------------------- tex parsing / comparison
def tex_rows(name):
    t = (TAB / f"{name}.tex").read_text()
    body = t.split("\\midrule", 1)[1].split("\\bottomrule")[0]
    rows = []
    for line in body.replace("\\midrule", "").splitlines():
        line = line.strip()
        if not line.endswith("\\\\"): continue
        cells = [c.strip() for c in line[:-2].split("&")]
        cells = [re.sub(r"\\textbf\{(.*)\}$", r"\1", c) for c in cells]
        rows.append(cells)
    return t, rows


def ms(xs): return f"{st.mean(xs):.3f} $\\pm$ {st.stdev(xs):.3f}"


def check(name, expected, got, note=""):
    bad = [(k, e, g) for k, (e, g) in enumerate(zip(expected, got)) if e != g]
    if len(expected) != len(got): bad.append(("rows", len(expected), len(got)))
    REPORT.append((name, "PASS" if not bad else "FAIL", note if not bad else f"{len(bad)} mismatches, e.g. {bad[:2]}"))


def numeric_cells(rows, start):
    return [r[start:] for r in rows]


# ---------------------------------------------------------------- per-table expectations
def crc_tables():
    for label, split, idset, metrics, rows, caption in T.TABLES:
        data = [rho(r[1]) if isinstance(r, tuple) else af3(r, idset) for r in rows]
        ids = sorted(i for i in split_ids(split) - T.EXCLUDE.get(label, set())
                     if all(i in d and all(d[i].get(m) is not None for m in metrics) for d in data))
        exp = [[ms([d[i][m] for i in ids]) for m in metrics] for d in data]
        _, got = tex_rows(label)
        ok_n = re.search(r"n=" + str(len(ids)) + r"(?!\d)", (TAB / f"{label}.tex").read_text()) is not None
        check(label, exp, numeric_cells(got, 1), f"{len(rows)} rows x {len(metrics)} metrics, n={len(ids)}" + ("" if ok_n else " (n in caption NOT found)"))
        if not ok_n: REPORT[-1] = (label, "FAIL", f"n={len(ids)} not in caption")



def benchmark_orphans():
    """main-text benchmark table = the verified orphan split tables, side by side (cells must be identical)."""
    _, got = tex_rows("benchmark_orphans")
    _, sc = tex_rows("rna_single_chain_post_cutoff_without_homologs")
    _, cx = tex_rows("rna_protein_complexes_post_cutoff_without_homologs")
    cidx = [T.COMPLEX.index(m) for m in T.COMPLEX_RED]
    exp = [a[1:] + [b[1 + k] for k in cidx] for a, b in zip(sc, cx)]
    check("benchmark_orphans", exp, numeric_cells(got, 1), "cells identical to the two verified orphan split tables")


def ribozymes():
    md = pd.read_csv(MD_CSV, sep=None, engine="python"); md["id"] = md.exp_db_id.str.lower()
    rib = set(md[md.RNADescription.fillna("").str.lower().str.contains("ribozyme")].id)
    S, M, N = af3("rnaformer", "complex"), af3("af3", "complex"), af3("empty", "complex")
    ids = sorted(i for i in rib if all(i in d and "RNA_TM" in d[i] for d in (S, M, N)) and N[i]["RNA_TM"] < 0.5)
    exp = [[f"{M[i]['RNA_TM']:.3f}", f"{S[i]['RNA_TM']:.3f}"] for i in ids]
    exp.append([f"{np.mean([M[i]['RNA_TM'] for i in ids]):.3f}", f"{np.mean([S[i]['RNA_TM'] for i in ids]):.3f}"])
    _, got = tex_rows("ribozymes"); check("ribozymes", exp, [r[2:] for r in got], f"n={len(ids)} + mean row")

def gates():
    mc = pd.read_csv(MC_CSV); mc["id"] = mc.exp_db_id.str.lower()
    sc = {}
    for f in glob.glob(f"{SELFCONS_DIR}/rnaformer_single_s*_m*.tsv.gz"):
        for r in csv.DictReader(gzip.open(f, "rt"), delimiter="\t"):
            if r["status"] == "ok": sc.setdefault(r["id"].lower(), []).append(float(r["selfcons_mcc"]))
    sc = {i: sum(v) / len(v) for i, v in sc.items()}
    for idset, group, metric, name in (("single", "rna-only", "RNA_TM", "gates_single"), ("complex", "1v1", "Dockq", "gates_complex1v1")):
        w = mc[mc.group == group].drop_duplicates("id").set_index("id")
        N, S, M = af3("empty", idset), af3("rnaformer", idset), af3("af3", idset)
        ids = [i for i in w.index if i in N and i in S and i in M and metric in N[i] and metric in S[i] and metric in M[i]]
        subsets = [[i for i in ids if lo <= w.rf2d_pred[i] < hi] for lo, hi in ((-1.01, .5), (.5, .8), (.8, .95), (.95, 1.01))]
        if idset == "single": subsets += [[i for i in ids if sc.get(i, -9) >= 0.9], [i for i in ids if i in sc and sc[i] < 0.9]]
        exp = []
        for g in subsets:
            n_, s_, m_ = [N[i][metric] for i in g], [S[i][metric] for i in g], [M[i][metric] for i in g]
            d = np.array(s_) - np.array(n_)
            p = wilcoxon(s_, n_).pvalue if len(g) > 5 and np.any(d != 0) else float("nan")
            exp.append([str(len(g)), f"{np.mean(n_):.3f}", f"{np.mean(s_):.3f}", f"{np.mean(m_):.3f}", f"${d.mean():+.3f}$",
                        "--" if np.isnan(p) else f"{p:.1e}"])
        _, got = tex_rows(name); check(name, exp, numeric_cells(got, 1), f"{len(exp)} subsets")


def seed_robustness():
    _, got = tex_rows("seed_robustness_table"); exp = []
    for tag in ("af3", "empty", "rnafold", "spotrna", "rnaformer", "dssr"):
        for idset, metric in (("single", "RNA_TM"), ("complex", "Dockq")):
            if tag in ("dssr", "empty") and idset == "complex": continue
            ids = split_ids("rna_monomers_ids" if idset == "single" else "protein_rna_ids"); per = {}
            tbl = "pred_rna_rna" if idset == "single" else "pred_protein_rna"
            for seed in (1, 0, 42, 137):
                vals = {}
                for db in glob.glob(f"{DB_DIR}/{tag}_{idset}_s{seed}_m*.db"):
                    for fn, v in sqlite3.connect(f"file:{db}?mode=ro", uri=True).execute(f"SELECT FileName,{metric} FROM {tbl}"):
                        v = val(v); i = pid(fn)
                        if v is not None and i in ids and not (metric == "RNA_TM" and v == 0): vals.setdefault(i, []).append(v)
                per[seed] = vals
            common = sorted(set.intersection(*[set(v) for v in per.values()]))
            sm = [st.mean(st.mean(per[s][i]) for i in common) for s in (1, 0, 42, 137)]
            allm = st.mean(st.mean(sum((per[s][i] for s in per), [])) for i in common)
            sd_s = st.mean(st.pstdev(per[s][i]) for s in per for i in common if len(per[s][i]) > 1)
            exp.append([str(len(common))] + [f"{x:.3f}" for x in sm + [allm]] + [f"{st.pstdev(sm):.4f}", f"{sd_s:.3f}"])
    check("seed_robustness_table", exp, [r[3:] for r in got], "recomputed per seed from the DBs")


def long_rna():
    conds = [("af3", "af3msa"), ("dssr",), ("nomsa",), ("shuffle",), ("random",)]; data = []
    for prefixes in conds:
        acc = {}
        for p in prefixes:
            for db in glob.glob(f"{DB_LONG}/{p}_s*_m*.db"):
                con = sqlite3.connect(f"file:{db}?mode=ro", uri=True); con.row_factory = sqlite3.Row
                for r in con.execute("SELECT * FROM pred_rna_rna"):
                    r = dict(r); i = pid(r["FileName"])
                    for m in SINGLE:
                        x = val(r.get(m))
                        if m == "RNA_TM" and x == 0.0: x = None
                        if x is not None: acc.setdefault(i, {}).setdefault(m, []).append(x)
        data.append({i: {m: st.mean(v) for m, v in d.items()} for i, d in acc.items()})
    ids = sorted(set.intersection(*[set(d) for d in data]))
    exp = [[ms([d[i][m] for i in ids]) for m in SINGLE] for d in data]
    _, got = tex_rows("long_rna"); check("long_rna", exp, numeric_cells(got, 1), f"n={len(ids)}")
    t = (TAB / "long_rna.tex").read_text(); tm = [st.mean(d[i]["RNA_TM"] for i in ids) for d in data]
    gap = f"recover {(tm[1] - tm[2]) / (tm[0] - tm[2]) * 100:.0f}\\%"
    rows = list(csv.DictReader(open(TIMING_CSV)))
    frac = st.mean(float(r["orig_msa_s"]) / float(r["orig_total_s"]) for r in rows)
    ok = gap in t and f"takes {frac * 100:.0f}\\%" in t and f"n{{=}}{len(ids)}$" in t
    REPORT.append(("long_rna caption", "PASS" if ok else "FAIL", "gap recovery, runtime share and n"))


def orphan_tests():
    _, got = tex_rows("orphan_tests"); exp = []; boot = np.random.default_rng(0)
    for idset, split, metric in (("single", "rna_monomers_a2021_orphan_ids", "RNA_TM"), ("complex", "protein_rna_a2021_orphan_ids", "Dockq")):
        M, S = af3("af3", idset), af3("rnaformer", idset)
        ids = sorted(i for i in split_ids(split) if i in M and i in S and metric in M[i] and metric in S[i])
        m_, s_ = [M[i][metric] for i in ids], [S[i][metric] for i in ids]
        d = np.array(s_) - np.array(m_); bs = [boot.choice(d, len(d)).mean() for _ in range(5000)]
        exp.append([str(len(ids)), f"{np.mean(m_):.3f}", f"{np.mean(s_):.3f}", f"${d.mean():+.3f}$",
                    f"$[{np.percentile(bs, 2.5):+.3f}, {np.percentile(bs, 97.5):+.3f}]$", f"{wilcoxon(s_, m_).pvalue:.2f}"])
    check("orphan_tests", exp, numeric_cells(got, 1), "means, diff, bootstrap CI and Wilcoxon p recomputed")


def order_and_rf_splits():
    tags = [("rnafold", "rnafoldN100Rhofold"), ("spotrna", "spotrnaN100Rhofold"), ("rnaformer", "rnaformerN100Rhofold")]
    A = [af3(t, "single") for t, _ in tags]; R = [rho(r) for _, r in tags]
    ids = sorted(i for i in split_ids("rna_monomers_ids") - {"9cxf"} if all(i in d and all(d[i].get(m) is not None for m in SINGLE) for d in A + R))
    exp = [[ms([d[i][m] for i in ids]) for m in SINGLE] for d in A + R]
    _, got = tex_rows("order_rna_single_chain"); check("order_rna_single_chain", exp, [r[2:] for r in got], f"n={len(ids)}")
    exp = []
    for split in ("rna_monomers_b2021_orphan_ids", "rna_monomers_b2021_non_orphan_ids", "rna_monomers_a2021_non_orphan_ids"):
        ids = sorted(i for i in split_ids(split) if all(i in d and all(d[i].get(m) is not None for m in SINGLE) for d in R))
        exp += [[f"{st.mean(d[i][m] for i in ids):.3f}" for m in SINGLE] for d in R]
    _, got = tex_rows("rf_rna_single_chain_other_splits"); check("rf_rna_single_chain_other_splits", exp, [r[2:] for r in got])


def confidence():
    single = {x.strip().lower() for x in open(ROOT / "data" / "ids" / "single_91.txt") if x.strip()}
    def ptm(tag):
        acc = {}
        for f in glob.glob(str(OUTPUTS / "af3" / f"{tag}_single" / "*" / "*" / "seed-*_sample-*" / "summary_confidences.json")):
            i = pid(Path(f).parent.parent.name); sub = Path(f).parent.parent.parent.name
            seed = int(re.search(r"seed-(\d+)_", f).group(1))
            if seed == 1 and sub not in ("seed_1", "seeds_1_0_42_137"): continue   # seed 1 is read from its own folder only
            if i in single:
                v = json.load(open(f)).get("ptm")
                if v is not None: acc.setdefault(i, []).append(float(v))
        if any(len(v) != 20 for v in acc.values()): raise ValueError(f"{tag}: not 20 predictions per target")
        return {i: st.mean(v) for i, v in acc.items()}
    tags = ["af3", "rnaformer", "empty", "shuffleSHS", "randomSHS"]; C = {t: ptm(t) for t in tags}
    ids = sorted(set.intersection(*[set(v) for v in C.values()])); base = np.array([C["empty"][i] for i in ids]); exp = []
    for t in tags:
        x = np.array([C[t][i] for i in ids])
        exp.append([f"{x.mean():.3f} $\\pm$ {x.std(ddof=1):.3f}"] + (["--", "--"] if t == "empty" else
                   [f"${(x - base).mean():+.3f}$", f"{wilcoxon(x, base).pvalue:.1e}"]))
    _, got = tex_rows("confidence_ptm"); check("confidence_ptm", exp, numeric_cells(got, 1), f"n={len(ids)}, read from summary_confidences.json")


def near_floor_novelty():
    d = pd.read_csv(Path(OUT_ROOT) / "near_floor" / "near_floor.csv")
    exp = [[str(int(r.n)), f"{r.MSA:.3f}", f"{r.SHS:.3f}", f"${r.d_SHS_MSA:+.3f}$", f"{r.frac_SHS_ge_MSA:.0f}\\%"] for r in d.itertuples()]
    _, got = tex_rows("near_floor"); check("near_floor", exp, numeric_cells(got, 1), "vs near_floor.csv")
    seq = list(csv.DictReader(open(Path(OUT_ROOT) / "novelty" / "sequence_similarity.csv"))); exp = []
    for sub in ("single_chain", "complex"):
        g = {round(float(r["threshold"]), 2): r for r in seq if r["subtype"] == sub}
        tm = [float(v) for v in json.load(open(Path(OUT_ROOT) / "novelty" / f"structure_maxTM_{sub}.json")).values()]
        exp.append([g[0.8]["n_post"], f"{g[0.9]['novel']}/{g[0.9]['n_post']}", f"{g[0.8]['novel']}/{g[0.8]['n_post']}",
                    f"{st.median(tm):.2f}", f"{sum(x > 0.5 for x in tm)}/{len(tm)}", f"{sum(x > 0.7 for x in tm)}/{len(tm)}"])
    _, got = tex_rows("novelty"); check("novelty", exp, numeric_cells(got, 1), "vs novelty CSV/JSON")


def construction():
    acc = {}
    for r in csv.DictReader(gzip.open(ROOT / "results" / "dssr2d" / "construction_2d.tsv.gz", "rt"), delimiter="\t"):
        acc.setdefault((r["target"], r["stage"], int(r["N"])), []).append(float(r["mcc_vs_reference"]))
    if any(len(v) != 20 for v in acc.values()):
        REPORT.append(("construction_stage_mcc", "FAIL", "not 20 predictions per cell")); return
    exp = []
    for t in ("4P8Z", "6WJS"):
        row = []
        for s in ("local_wc", "init_wc", "step1_nc", "step2_pk", "step3_multi", "final_all"):
            per_n = [st.mean(v) for (tt, ss, n), v in sorted(acc.items()) if tt == t and ss == s]
            row.append(ms(per_n))
        exp.append(row)
    _, got = tex_rows("construction_stage_mcc"); check("construction_stage_mcc", exp, numeric_cells(got, 1), "3360 DSSR runs")


def copies(dest):
    dest = Path(dest); bad = []
    for f in list((dest / "main").glob("*.tex")) + list((dest / "appendix_sections").glob("*.tex")):
        if f.read_bytes() != (TAB / f.name).read_bytes(): bad.append(f.name)
    app = (dest / "appendix_additional_results.tex").read_text()
    n_in = 0
    for name in re.findall(r"\\label\{(tab:[^}]*)\}", app):
        src = [f for f in TAB.glob("*.tex") if f.name != "all_tables.tex" and f"\\label{{{name}}}" in f.read_text()]
        if not src or src[0].read_text() not in app: bad.append(f"appendix:{name}")
        n_in += 1
    REPORT.append(("copy-paste folder", "PASS" if not bad else "FAIL",
                   f"main/ + appendix_sections/ byte-identical, {n_in} appendix tables verbatim" if not bad else f"differs: {bad}"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--dest", required=True); a = ap.parse_args()
    for f in (crc_tables, benchmark_orphans, ribozymes, gates, seed_robustness, long_rna, orphan_tests, order_and_rf_splits, confidence, near_floor_novelty,
              construction, lambda: copies(a.dest)):
        try: f()
        except Exception as e: REPORT.append((getattr(f, "__name__", "copies"), "FAIL", f"{type(e).__name__}: {e}"))
    nfail = sum(r[1] != "PASS" for r in REPORT)
    lines = ["# Table verification report", "", f"{len(REPORT)} checks, {nfail} failed. Every number is recomputed from the raw "
             "sources with code independent of the table generators and compared with the printed .tex cells.", "",
             "| Table | Result | Detail |", "|---|---|---|"] + [f"| {n} | {r} | {d} |" for n, r, d in REPORT]
    Path(a.dest, "VERIFY_REPORT.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    sys.exit(1 if nfail else 0)
