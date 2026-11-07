#!/usr/bin/env python3
"""Long-RNA mechanism ladder (25 RNAs, 202-746 nt).
Mean RNA_TM per condition over the multiseed long-run DBs (results/long_run/<cond>_s*_m*.db).
Conditions: natural MSA (af3msa), reference-topology SHS (dssr), no-MSA, shuffle-SHS, random-SHS.
Also reports the fraction of the MSA-over-no-MSA gap recovered by reference-SHS. Run in synEvo."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, glob, re, statistics as st, os
import pandas as pd
DB = DB_LONG
OUT = f"{OUT_ROOT}/long_rna"; os.makedirs(OUT, exist_ok=True)
def idof(s):
    m = re.search(r"([0-9][A-Za-z0-9]{3})", s or ""); return m.group(1).lower() if m else None
def ms(cond):
    acc = {}
    for db in glob.glob(f"{DB}/{cond}_s*_m*.db"):
        for fn, v in sqlite3.connect(db).execute("SELECT FileName,RNA_TM FROM pred_rna_rna"):
            i = idof(fn)
            try: v = float(v)
            except: continue
            if i and v != 0: acc.setdefault(i, []).append(v)
    return {i: st.mean(x) for i, x in acc.items()}
COND = {"natural MSA": "af3msa", "reference-topology SHS": "dssr", "no-MSA": "nomsa",
        "shuffle-SHS": "shuffle", "random-SHS": "random"}
data = {k: ms(v) for k, v in COND.items()}
common = sorted(set.intersection(*[set(d) for d in data.values() if d]))
means = {k: st.mean([d[i] for i in common]) for k, d in data.items()}
print(f"long-RNA ladder, mean RNA_TM over n={len(common)} RNAs (multiseed):")
for k in COND:
    print(f"  {k:26} {means[k]:.3f}")
gap = (means["reference-topology SHS"] - means["no-MSA"]) / (means["natural MSA"] - means["no-MSA"])
print(f"  reference-SHS recovers {gap*100:.0f}% of the natural-MSA over no-MSA gap")
pd.DataFrame([{"condition": k, "RNA_TM": means[k]} for k in COND]).to_csv(f"{OUT}/long_run_ladder.csv", index=False)
