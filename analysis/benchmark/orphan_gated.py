import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, glob, re, csv, statistics as st
import pandas as pd, numpy as np
from scipy.stats import wilcoxon
DB=DB_DIR; SPl=f"{OUT_ROOT}/paper_splits"
MC=MC_CSV
def idof(s):
    m=re.search(r"([0-9][A-Za-z0-9]{3})",s or ""); return m.group(1).lower() if m else None
def ms(method):
    acc={}
    for db in glob.glob(f"{DB}/{method}_single_s*_m*.db"):
        for fn,v in sqlite3.connect(db).execute("SELECT FileName,RNA_TM FROM pred_rna_rna"):
            i=idof(fn)
            try:v=float(v)
            except:continue
            if i and v!=0: acc.setdefault(i,[]).append(v)
    return {i:st.mean(x) for i,x in acc.items()}
def pub_ids(f):
    p=f"{SPl}/{f}.csv"; d=";" if ";" in open(p).readline() else ","
    return {str(r.get("exp_db_id","")).lower() for r in csv.DictReader(open(p),delimiter=d) if r.get("exp_db_id")}
SHS=ms("rnaformer"); NOM=ms("empty")
orph=split_ids("RNA_Monomers_a2021_alphafold_orphan")
w=pd.read_csv(MC); w["idl"]=w.exp_db_id.astype(str).str.lower()
w.loc[w.group=="rna-only","rf2d_selfcons"]=w.loc[w.group=="rna-only","idl"].map(selfcons_means())  # 4 seeds x 5 samples
w=w[w.idl.isin(orph)].copy()
w["S"]=w.idl.map(SHS); w["N"]=w.idl.map(NOM); w=w.dropna(subset=["S","N"])
print(f"post-cutoff single-chain ORPHANS with SHS+noMSA(multiseed): n={len(w)}")
print(f"overall: SHS={w.S.mean():.3f} noMSA={w.N.mean():.3f} d(SHS-noMSA)={ (w.S-w.N).mean():+.3f}")
print()
for gate,lab in [("rf2d_selfcons","self-consistency"),("rf2d_pred","GT-MCC")]:
    ww=w.dropna(subset=[gate])
    print(f"=== ORPHANS gated by {lab} ({gate}) ===")
    for cut in (0.7,0.8,0.9,0.95):
        g=ww[ww[gate]>=cut]
        if len(g)==0: print(f"  >={cut}: n=0"); continue
        d=(g.S-g.N)
        p=wilcoxon(g.S.values,g.N.values).pvalue if len(g)>=6 else float('nan')
        frac=(g.S.values>g.N.values).mean()*100
        print(f"  >={cut}: n={len(g):2d} SHS={g.S.mean():.3f} noMSA={g.N.mean():.3f} d={d.mean():+.3f} (SHS>noMSA in {frac:.0f}%, p={p:.3f})")
    print()
