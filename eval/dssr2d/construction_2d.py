#!/usr/bin/env python3
"""2D accuracy of the 3D predictions of the 4P8Z / 6WJS construction sweep (appendix construction-stage table).

For every prediction of af3_case_studies/<T>_intervention_full_sweep (seed 1) and ..._multiseed (seeds 0, 42, 137),
5 samples each, 6 topology stages x 14 SHS depths: DSSR base pairs of the predicted structure vs DSSR base pairs of the
reference structure (data/gt_intervention/<t>_gt.cif); MCC over the full symmetric L x L matrix (same definition as the
submitted table, which used seed 1, sample 0 per depth: `--validate` reproduces those 168 values).

DSSR (x3dna-dssr v2.2.1) is licensed and NOT distributed; set DSSR=/path/to/x3dna-dssr. Outputs (shipped, precomputed):
  python eval/dssr2d/construction_2d.py --task K --ntasks M   -> results/dssr2d/construction/task<K>.tsv.gz
  python eval/dssr2d/construction_2d.py --merge                 -> results/dssr2d/construction_2d.tsv.gz
  python eval/dssr2d/construction_2d.py --validate              -> compares seed-1 sample-0 MCC with the submitted values
"""
import argparse, glob, gzip, json, os, re, subprocess, sys, tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis"))
from paths import OUTPUTS  # noqa: E402
DSSR = os.environ.get("DSSR", os.path.expanduser("~/DSSR/x3dna-dssr"))
CS = OUTPUTS / "af3_case_studies"
OUT = ROOT / "results" / "dssr2d"
NAME = re.compile(r"^(?P<t>[0-9a-z]{4})_custom_rnamsa_n(?P<n>\d+)_.*_intervention_(?P<stage>[a-z0-9_]+)$")


def dssr_pairs(cif, work):
    out, link = os.path.join(work, "out.json"), os.path.join(work, "in.cif")   # DSSR aborts on very long paths
    for f in (out, link):
        if os.path.lexists(f): os.remove(f)
    os.symlink(os.path.realpath(cif), link)
    subprocess.run([DSSR, "-i=in.cif", "--json", "-o=out.json"], cwd=work, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    d = json.load(open(out)); nts = {n["nt_id"]: n for n in d.get("nts", [])}
    L = len(nts)
    P = [(nts[p["nt1"]]["index_chain"] - 1, nts[p["nt2"]]["index_chain"] - 1) for p in d.get("pairs", [])
         if nts[p["nt1"]]["chain_name"] == nts[p["nt2"]]["chain_name"]]
    return P, L


def matrix(P, L):
    m = np.zeros((L, L), bool)
    for i, j in P:
        if 0 <= i < L and 0 <= j < L and i != j: m[i, j] = m[j, i] = True
    return m


def mcc(p, t):
    tp = (p & t).sum(); tn = (~p & ~t).sum(); fp = (p & ~t).sum(); fn = (~p & t).sum()
    den = np.sqrt(float(tp + fp) * float(tp + fn) * float(tn + fp) * float(tn + fn))
    return float((tp * tn - fp * fn) / den) if den else 0.0


def jobs():
    out = []
    for T in ("4P8Z", "6WJS"):
        for sub in (f"{T}_intervention_full_sweep", f"{T}_intervention_full_sweep_multiseed"):
            for d in sorted((CS / sub).iterdir()):
                m = NAME.match(d.name)
                if not m: continue
                for s in sorted(d.glob("seed-*_sample-*")):
                    seed, sample = re.match(r"seed-(\d+)_sample-(\d+)", s.name).groups()
                    out.append((T, m["stage"], int(m["n"]), int(seed), int(sample), str(s / "model.cif")))
    return out


def run(task, ntasks):
    todo = jobs()[task::ntasks]; refs = {}; rows = []
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as w:
        for T, stage, n, seed, sample, cif in todo:
            if T not in refs: refs[T] = dssr_pairs(str(ROOT / "data" / "gt_intervention" / f"{T.lower()}_gt.cif"), w)
            R, L = refs[T]; P, Lp = dssr_pairs(cif, w)
            rows.append((T, stage, n, seed, sample, Lp, len(set(map(tuple, map(sorted, P)))), f"{mcc(matrix(P, L), matrix(R, L)):.6f}"))
    (OUT / "construction").mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / "construction" / f"task{task:03d}.tsv.gz", "wt") as f:
        f.write("target\tstage\tN\tseed\tsample\tn_nt_pred\tn_pairs_pred\tmcc_vs_reference\n")
        for r in rows: f.write("\t".join(map(str, r)) + "\n")
    print(f"task {task}: {len(rows)} predictions")


def merge(ntasks):
    parts = sorted((OUT / "construction").glob("task*.tsv.gz"))
    lines = []
    for p in parts: lines += gzip.open(p, "rt").read().splitlines()[1:]
    n = len(jobs())
    if len(lines) != n: sys.exit(f"merge: {len(lines)} rows, expected {n} ({len(parts)} task files)")
    with gzip.open(OUT / "construction_2d.tsv.gz", "wt") as f:
        f.write("target\tstage\tN\tseed\tsample\tn_nt_pred\tn_pairs_pred\tmcc_vs_reference\n" + "\n".join(sorted(lines)) + "\n")
    print(f"merged {len(lines)} rows -> {OUT / 'construction_2d.tsv.gz'}")


def validate():
    import pandas as pd
    old = pd.read_csv(ROOT / "data" / "gates" / "inputs" / "all_secondary_structure_data_rna_rna.csv")
    old = old[old.experiment.str.contains("full_sweep", na=False)]
    refs = {}; worst = 0.0; k = 0
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as w:
        for r in old.itertuples():
            T = r.pdbid; stage = r.intervention_stage; n = int(r.shs_depth)
            d = [x for x in (CS / f"{T}_intervention_full_sweep").iterdir()
                 if (m := NAME.match(x.name)) and m["stage"] == stage and int(m["n"]) == n][0]
            top = str(d / "seed-1_sample-0" / "model.cif")
            if T not in refs: refs[T] = dssr_pairs(str(ROOT / "data" / "gt_intervention" / f"{T.lower()}_gt.cif"), w)
            R, L = refs[T]; P, _ = dssr_pairs(top, w)
            worst = max(worst, abs(mcc(matrix(P, L), matrix(R, L)) - r.mcc)); k += 1
    print(f"validate: {k} seed-1 sample-0 predictions of the sweep, max |MCC - submitted| = {worst:.2e}")
    if worst > 1e-6: sys.exit("VALIDATION FAILED")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--merge", action="store_true"); ap.add_argument("--validate", action="store_true")
    a = ap.parse_args()
    if a.validate: validate()
    elif a.merge: merge(a.ntasks)
    else: run(a.task, a.ntasks)
