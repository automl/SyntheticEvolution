import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, glob, re, csv, statistics as st
import pandas as pd, numpy as np
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
    return acc  # keep raw lists to inspect
def pub_ids(f):
    p=f"{SPl}/{f}.csv"; d=";" if ";" in open(p).readline() else ","
    return {str(r.get("exp_db_id","")).lower() for r in csv.DictReader(open(p),delimiter=d) if r.get("exp_db_id")}
SHSr=ms("rnaformer"); NOMr=ms("empty")
orph=sorted(split_ids("RNA_Monomers_a2021_alphafold_orphan"))
w=pd.read_csv(MC); w["idl"]=w.exp_db_id.astype(str).str.lower()
w.loc[w.group=="rna-only","rf2d_selfcons"]=w.loc[w.group=="rna-only","idl"].map(selfcons_means())  # 4 seeds x 5 samples
w=w.set_index("idl")
print(f"orphan ids in split: {len(orph)} -> {orph}")
print(f"\n{'id':6}{'selfcons':>9}{'GTmcc':>7}{'#SHSpred':>9}{'SHS':>7}{'#NOMpred':>9}{'noMSA':>7}{'SHS-noMSA':>11}")
rows=[]
for i in orph:
    if i not in SHSr or i not in NOMr: continue
    sc = w.loc[i,'rf2d_selfcons'] if i in w.index else np.nan
    mc = w.loc[i,'rf2d_pred'] if i in w.index else np.nan
    s=st.mean(SHSr[i]); n=st.mean(NOMr[i])
    rows.append((i,sc,mc,len(SHSr[i]),s,len(NOMr[i]),n,s-n))
for r in sorted(rows, key=lambda x:-(x[1] if x[1]==x[1] else -9)):
    sc = f"{r[1]:.3f}" if r[1]==r[1] else "  NaN"
    mc = f"{r[2]:.3f}" if r[2]==r[2] else "  NaN"
    print(f"{r[0]:6}{sc:>9}{mc:>7}{r[3]:>9}{r[4]:>7.3f}{r[5]:>9}{r[6]:>7.3f}{r[7]:>+11.3f}")
print(f"\ntotal orphans with both SHS+noMSA: {len(rows)}")
# selfcons>=0.9 subset
hi=[r for r in rows if r[1]==r[1] and r[1]>=0.9]
print(f"selfcons>=0.9: n={len(hi)} ; all SHS>noMSA? {all(r[7]>0 for r in hi)} ; mean d={np.mean([r[7] for r in hi]):+.3f}")
