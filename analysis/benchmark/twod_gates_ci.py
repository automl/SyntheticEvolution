import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, glob, re, csv, statistics as st
import pandas as pd, numpy as np
np.random.seed(0)
DB=DB_DIR; MC=MC_CSV
MD=MD_CSV
def idof(s):
    m=re.search(r"([0-9][A-Za-z0-9]{3})",s or ""); return m.group(1).lower() if m else None
def ms(method,idset,col):
    tbl="pred_rna_rna" if idset=="single" else "pred_protein_rna"; acc={}
    for db in glob.glob(f"{DB}/{method}_{idset}_s*_m*.db"):
        for fn,v in sqlite3.connect(db).execute(f"SELECT FileName,{col} FROM {tbl}"):
            i=idof(fn)
            try:v=float(v)
            except:continue
            if i and v!=0: acc.setdefault(i,[]).append(v)
    return {i:st.mean(x) for i,x in acc.items()}
def bootci(d,n=5000):
    d=np.array(d); bs=[np.mean(np.random.choice(d,len(d),replace=True)) for _ in range(n)]
    return np.percentile(bs,2.5),np.percentile(bs,97.5)
S=ms("rnaformer","single","RNA_TM"); M=ms("af3","single","RNA_TM"); N=ms("empty","single","RNA_TM")
w=pd.read_csv(MC); w=w[w.group=="rna-only"].copy(); w["idl"]=w.exp_db_id.astype(str).str.lower(); w=w.dropna(subset=["rf2d_pred"])
g=w[w.rf2d_pred>=0.9]; ids=[i for i in g.idl if i in S and i in N]
d=np.array([S[i] for i in ids])-np.array([N[i] for i in ids])
lo,hi=bootci(d); print(f"GT-gate >=0.9 SHS-noMSA: mean={d.mean():+.3f} 95%CI=[{lo:+.3f},{hi:+.3f}] n={len(ids)}")
# ribozymes
Sc=ms("rnaformer","complex","RNA_TM"); Mc=ms("af3","complex","RNA_TM"); Nc=ms("empty","complex","RNA_TM")
dd=";" if ";" in open(MD).readline() else ","
rib={str(r["exp_db_id"]).lower() for r in csv.DictReader(open(MD),delimiter=dd) if "ribozyme" in (r.get("RNADescription") or "").lower()}
hard=sorted(i for i in rib if i in Sc and i in Mc and i in Nc and Nc[i]<0.5)
a=np.array([Sc[i] for i in hard]); m=np.array([Mc[i] for i in hard]); n=np.array([Nc[i] for i in hard])
lo,hi=bootci(a-m); print(f"ribozyme SHS-MSA: mean={ (a-m).mean():+.3f} 95%CI=[{lo:+.3f},{hi:+.3f}] n={len(hard)}")
lo,hi=bootci(a-n); print(f"ribozyme SHS-noMSA: mean={ (a-n).mean():+.3f} 95%CI=[{lo:+.3f},{hi:+.3f}] n={len(hard)}")
