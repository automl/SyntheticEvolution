#!/usr/bin/env python
"""Eval the 91-RNA construction sweep: US-align every predicted cif vs GT,
mean TM per (RNA, step), recall per (RNA, step); pooled recall->TM + per-step means.
Writes construct91_tm.csv (per-prediction) + construct91_summary.csv (per RNA,step)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import json, os, glob, re, subprocess, sys
from multiprocessing import Pool

SCR = CONSTRUCT_DIR
UA = USALIGN
GTDIR = GT_SINGLE
OUTBASE = case_study_dir(CONSTRUCT_RUN)   # CONSTRUCT_N=100 -> construct91_N100_multiseed
PAIRS = f"{CONSTRUCT_DATA}/construct_pairs_91.json"   # frozen output of build_construct_91.py
IDS = str(IDS_DIR/"single_91.txt")
STEPS = ["wc","nc","pk","all"]

def recall_map():
    recs = json.load(open(PAIRS))
    npairs = {}
    for r in recs:
        npairs[r["id"]] = len(r["pairs"])
    out = {}
    for pid in {r["pdb_id"] for r in recs}:
        nall = npairs.get(f"{pid}_all", 0) or 1
        for s in STEPS:
            out[(pid,s)] = round(npairs.get(f"{pid}_{s}",0)/nall, 4)
    return out

def find_cifs(pid, step):
    # output dir name starts with "{pid}_{step}_custom_rnamsa_...construct" (lowercase)
    pat = f"{OUTBASE}/{pid}_{step}_custom_rnamsa_*construct"
    dirs = glob.glob(pat)
    cifs = []
    for d in dirs:
        cifs += glob.glob(f"{d}/seed-*_sample-*/model.cif")
    return cifs

def usalign_tm(args):
    pid, step, cif = args
    gt = f"{GTDIR}/{pid}.cif"
    try:
        out = subprocess.run([UA, cif, gt, "-mol","RNA","-ter","0"],
                             capture_output=True, text=True, timeout=120).stdout
    except Exception:
        return (pid, step, cif, None)
    tm = None
    for line in out.splitlines():
        if "normalized by length of Structure_2" in line:
            m = re.search(r"TM-score=\s*([0-9.]+)", line)
            if m: tm = float(m.group(1)); break
    m2 = re.search(r"seed-(\d+)_sample-(\d+)", cif)
    seed, samp = (m2.group(1), m2.group(2)) if m2 else ("?","?")
    return (pid, step, seed, samp, tm)

def main():
    ids = [x.strip().lower() for x in open(IDS) if x.strip()]
    rc = recall_map()
    jobs = []
    for pid in ids:
        for step in STEPS:
            for cif in find_cifs(pid, step):
                jobs.append((pid, step, cif))
    print(f"US-align jobs: {len(jobs)} (expect ~{len(ids)*4*4*5})")
    if not jobs:
        print("No cifs found yet -- folds not complete."); sys.exit(1)
    with Pool(16) as p:
        res = p.map(usalign_tm, jobs)
    # per-prediction csv
    import csv
    with open(f"{SCR}/{CONSTRUCT_PREFIX}_tm.csv","w",newline="") as f:
        w = csv.writer(f); w.writerow(["pid","step","seed","sample","tm"])
        for r in res: w.writerow(r)
    # aggregate per (pid, step)
    from collections import defaultdict
    agg = defaultdict(list)
    for pid, step, seed, samp, tm in res:
        if tm is not None: agg[(pid,step)].append(tm)
    import statistics as St
    with open(f"{SCR}/{CONSTRUCT_PREFIX}_summary.csv","w",newline="") as f:
        w = csv.writer(f); w.writerow(["pid","step","recall","n","mean_tm","sd_tm"])
        for pid in ids:
            for step in STEPS:
                v = agg.get((pid,step),[])
                if v:
                    w.writerow([pid,step,rc.get((pid,step)),len(v),
                                round(St.mean(v),4), round(St.pstdev(v),4)])
    # per-step pooled means (across RNAs, mean of per-RNA means)
    per_rna = {(pid,step): St.mean(agg[(pid,step)]) for (pid,step) in agg}
    print("\nPer-step aggregate (mean over 91 RNA means):")
    for step in STEPS:
        vals = [per_rna[(pid,step)] for pid in ids if (pid,step) in per_rna]
        recs = [rc[(pid,step)] for pid in ids if (pid,step) in per_rna]
        if vals:
            print(f"  {step}: meanTM={St.mean(vals):.3f}  meanRecall={St.mean(recs):.3f}  nRNA={len(vals)}")
    # pooled spearman recall vs TM (per RNA,step points)
    pts = [(rc[(pid,step)], per_rna[(pid,step)]) for (pid,step) in per_rna]
    try:
        from scipy.stats import spearmanr
        rho,pval = spearmanr([a for a,_ in pts],[b for _,b in pts])
        print(f"\nPooled Spearman(recall, TM) over {len(pts)} points: rho={rho:.3f} p={pval:.2e}")
    except Exception as e:
        print("scipy unavailable:", e)
    print(f"\nwrote {SCR}/{CONSTRUCT_PREFIX}_tm.csv and _summary.csv")

if __name__=="__main__":
    main()
