import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, glob, re, statistics as st
import numpy as np
from scipy.stats import wilcoxon
DB=DB_DIR
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
base=ms("rnafold"); den=ms("denoise_rnafold"); ref=ms("dssr")
ids=sorted(set(base)&set(den)&set(ref))
b=np.array([base[i] for i in ids]); d=np.array([den[i] for i in ids]); r=np.array([ref[i] for i in ids])
print(f"AF3 denoise ladder MULTISEED (n={len(ids)}):")
print(f"  base(rnafold)={b.mean():.3f}  denoise={d.mean():.3f}  reference(dssr)={r.mean():.3f}")
print(f"  denoise - base = {(d-b).mean():+.4f}  (p={wilcoxon(d,b).pvalue:.4f})   [draft AF3: +0.002, p=0.20]")
print(f"  reference - base = {(r-b).mean():+.4f}   [draft AF3 ref: +0.082]")
