import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import json, glob, re, os, statistics as st
import numpy as np
from scipy.stats import wilcoxon
def idof(s):
    m=re.search(r"([0-9][A-Za-z0-9]{3})",s or ""); return m.group(1).lower() if m else None
singles=set(x.strip().lower() for x in open(IDS_DIR/"single_91.txt") if x.strip())
def scan(tag, seeds=SEEDS):
    pt={}; rk={}
    for f in prediction_files(tag, "single", "summary_confidences.json", seeds):
        i=idof(os.path.basename(os.path.dirname(os.path.dirname(f))))
        if not i or i not in singles: continue
        try: d=json.load(open(f))
        except: continue
        if d.get("ptm") is not None: pt.setdefault(i,[]).append(float(d["ptm"]))
        if d.get("ranking_score") is not None: rk.setdefault(i,[]).append(float(d["ranking_score"]))
    return {i:st.mean(x) for i,x in pt.items()}, {i:st.mean(x) for i,x in rk.items()}
# CRC: all conditions pooled over seeds 1/0/42/137 x 5 samples (posted rebuttal: shuffle/random were seed 1 only).
TAGS={"SHS":"rnaformer","noMSA":"empty","MSA":"af3","shuffle":"shuffleSHS","random":"randomSHS"}
PT={}; RK={}
for k,v in TAGS.items():
    PT[k],RK[k]=scan(v); print(f"scanned {k}: {len(PT[k])} ids", flush=True)
for name,C in [("ranking_score",RK),("ptm",PT)]:
    common=set.intersection(*[set(C[k]) for k in C])
    print(f"\n### {name} (n={len(common)}) ###")
    for k in ["MSA","SHS","noMSA","shuffle","random"]:
        print(f"  {k:8} {st.mean([C[k][i] for i in common]):.3f}")
    print(f"  SHS-noMSA={np.mean([C['SHS'][i]-C['noMSA'][i] for i in common]):+.3f}  SHS-shuffle={np.mean([C['SHS'][i]-C['shuffle'][i] for i in common]):+.3f}  SHS-random={np.mean([C['SHS'][i]-C['random'][i] for i in common]):+.3f}")
# single-seed-matched: seed 1 only for all three
S1pt,_=scan(TAGS["SHS"],seeds=(1,)); PT["shuffle"],_=scan(TAGS["shuffle"],seeds=(1,)); PT["random"],_=scan(TAGS["random"],seeds=(1,))
comm=set(S1pt)&set(PT["shuffle"])&set(PT["random"])
print(f"\n### single-seed-matched SHS(seed-1) vs shuffle/random, ptm (n={len(comm)}) ###")
print(f"  SHS(s1)={st.mean([S1pt[i] for i in comm]):.3f} shuffle={st.mean([PT['shuffle'][i] for i in comm]):.3f} random={st.mean([PT['random'][i] for i in comm]):.3f}")
print(f"  SHS-shuffle={np.mean([S1pt[i]-PT['shuffle'][i] for i in comm]):+.3f} (p={wilcoxon([S1pt[i] for i in comm],[PT['shuffle'][i] for i in comm]).pvalue:.4f})")
print(f"  SHS-random ={np.mean([S1pt[i]-PT['random'][i] for i in comm]):+.3f} (p={wilcoxon([S1pt[i] for i in comm],[PT['random'][i] for i in comm]).pvalue:.4f})")
print("DONE")
