#!/usr/bin/env python3
"""uTY4-Q4: is the post-cutoff parity 'shared failure'? Restrict the post-cutoff single-chain ORPHAN
split to solvable targets (best of SHS/MSA >= threshold) and re-compare SHS vs MSA (both multiseed;
on orphans MSA == no-MSA by construction). If SHS still matches MSA among solvable targets, the
parity is genuine, not two failures tying. Run in synEvo.
Writes reproduce_rebuttal/results/uTY4/near_floor.csv."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, glob, re, csv, os, statistics as st
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
DB = DB_DIR
SP = f"{OUT_ROOT}/paper_splits"
OUT = f"{OUT_ROOT}/uTY4"; os.makedirs(OUT, exist_ok=True)
def idof(s):
    m = re.search(r"([0-9][A-Za-z0-9]{3})", s or ""); return m.group(1).lower() if m else None
def ms(method):
    acc = {}
    for db in glob.glob(f"{DB}/{method}_single_s*_m*.db"):
        for fn, v in sqlite3.connect(db).execute("SELECT FileName,RNA_TM FROM pred_rna_rna"):
            i = idof(fn)
            try: v = float(v)
            except: continue
            if i and v != 0: acc.setdefault(i, []).append(v)
    return {i: st.mean(x) for i, x in acc.items()}
def pub_ids(f):
    p = f"{SP}/{f}.csv"; d = ";" if ";" in open(p).readline() else ","
    return {str(r.get("exp_db_id", "")).lower() for r in csv.DictReader(open(p), delimiter=d) if r.get("exp_db_id")}
SHS, MSA = ms("rnaformer"), ms("af3")
orph = split_ids("RNA_Monomers_a2021_alphafold_orphan")
ids = sorted(i for i in orph if i in SHS and i in MSA)
shs = np.array([SHS[i] for i in ids]); msa = np.array([MSA[i] for i in ids])
print(f"post-cutoff single-chain orphan, n_total={len(ids)}")
print(f"{'floor best>=':12}{'n':>4}{'MSA(=noMSA)':>13}{'SHS':>8}{'d(SHS-MSA)':>12}{'SHS>=MSA%':>10}")
rows = []
for t in (0.0, 0.30, 0.45, 0.50, 0.60):
    m = np.maximum(shs, msa) >= t; n = int(m.sum())
    if n < 3: print(f"{t:<12.2f}{n:>4}  (too few)"); continue
    d = (shs[m] - msa[m])
    frac = (shs[m] >= msa[m]).mean() * 100
    print(f"{t:<12.2f}{n:>4}{msa[m].mean():>13.3f}{shs[m].mean():>8.3f}{d.mean():>+12.3f}{frac:>9.0f}%")
    rows.append(dict(floor=t, n=n, MSA=msa[m].mean(), SHS=shs[m].mean(), d_SHS_MSA=d.mean(), frac_SHS_ge_MSA=frac))
pd.DataFrame(rows).to_csv(f"{OUT}/near_floor.csv", index=False)
print(f"wrote {OUT}/near_floor.csv")
