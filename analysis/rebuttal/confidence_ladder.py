import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import json, glob, re, os, statistics as st
import numpy as np
from scipy.stats import spearmanr, wilcoxon
def idof(s):
    m=re.search(r"([0-9][A-Za-z0-9]{3})",s or ""); return m.group(1).lower() if m else None
singles=set(x.strip().lower() for x in open(IDS_DIR/"single_91.txt") if x.strip())
def conf(tag, metric="ptm"):
    acc={}
    for f in prediction_files(tag, "single", "summary_confidences.json"):
        i=idof(os.path.basename(os.path.dirname(os.path.dirname(f))))
        if not i or i not in singles: continue
        try: v=json.load(open(f)).get(metric)
        except: continue
        if v is not None: acc.setdefault(i,[]).append(float(v))
    return {i:st.mean(x) for i,x in acc.items()}
# CRC: every condition pools seeds 1/0/42/137 x 5 samples, each seed read once (prediction_files).
# (The posted rebuttal numbers used the multiseed dirs for SHS/noMSA/MSA and seed-1 dirs for shuffle/random.)
TAGS={"SHS":"rnaformer","noMSA":"empty","MSA":"af3","shuffle":"shuffleSHS","random":"randomSHS"}
C={k:conf(t,"ptm") for k,t in TAGS.items()}
common=set.intersection(*[set(C[k]) for k in C])
print(f"single-chain ids with ptm in ALL methods: n={len(common)}")
for k in ["MSA","SHS","shuffle","random","noMSA"]:
    print(f"  {k:8} mean ptm = {st.mean([C[k][i] for i in common]):.3f}")
print()
base=[C['noMSA'][i] for i in common]
for k in ["SHS","shuffle","random","MSA"]:
    d=[C[k][i]-C['noMSA'][i] for i in common]
    p=wilcoxon([C[k][i] for i in common], base).pvalue
    print(f"  {k}-noMSA ptm: {np.mean(d):+.3f} (p={p:.4f})")
# key: SHS vs equal-depth shuffle/random
for k in ["shuffle","random"]:
    d=[C['SHS'][i]-C[k][i] for i in common]
    print(f"  SHS-{k} ptm (topology vs depth): {np.mean(d):+.3f} (p={wilcoxon([C['SHS'][i] for i in common],[C[k][i] for i in common]).pvalue:.4f})")
# does SHS confidence track accuracy? correlate ptm with RNA_TM
import sqlite3
def tm(method):
    acc={}
    for db in glob.glob(f"{DB_DIR}/{method}_single_s*_m*.db"):
        for fn,v in sqlite3.connect(f"file:{db}?mode=ro",uri=True).execute("SELECT FileName,RNA_TM FROM pred_rna_rna"):
            i=idof(fn)
            try:v=float(v)
            except:continue
            if i and v!=0: acc.setdefault(i,[]).append(v)
    return {i:st.mean(x) for i,x in acc.items()}
TM=tm("rnaformer")
ids=[i for i in common if i in TM]
r=spearmanr([C['SHS'][i] for i in ids],[TM[i] for i in ids])
print(f"\n  corr(SHS ptm, SHS RNA_TM) spearman={r.correlation:+.3f} (p={r.pvalue:.4f}, n={len(ids)})")
