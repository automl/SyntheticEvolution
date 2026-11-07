import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, glob, re, json, statistics as st
import numpy as np
from scipy.stats import spearmanr, pearsonr
DB=DB_DIR
def idof(s):
    m=re.search(r"([0-9][A-Za-z0-9]{3})",s or ""); return m.group(1).lower() if m else None
def ms(method,col):
    acc={}
    for db in glob.glob(f"{DB}/{method}_complex_s*_m*.db"):
        for fn,v in sqlite3.connect(db).execute(f"SELECT FileName,{col} FROM pred_protein_rna"):
            i=idof(fn)
            try:v=float(v)
            except:continue
            if i and v!=0: acc.setdefault(i,[]).append(v)
    return {i:st.mean(x) for i,x in acc.items()}
stab=json.load(open(str(CRC/"data"/"gates"/"structure_maxTM_complex.json")))
stab={idof(k):v for k,v in stab.items()}
MSA=ms("af3","Dockq"); SHS=ms("rnaformer","Dockq"); NOM=ms("empty","Dockq")
ids=[i for i in stab if i in MSA and i in SHS and i in NOM]
S=np.array([stab[i] for i in ids]); dM=np.array([MSA[i]-NOM[i] for i in ids]); dS=np.array([SHS[i]-NOM[i] for i in ids])
print(f"n complexes with structTM + all three scores: {len(ids)}")
print(f"corr(structTM, MSA-noMSA):  spearman={spearmanr(S,dM).correlation:+.3f} (p={spearmanr(S,dM).pvalue:.3f})")
print(f"corr(structTM, SHS-noMSA):  spearman={spearmanr(S,dS).correlation:+.3f} (p={spearmanr(S,dS).pvalue:.3f})")
print()
print(f"{'structTM bin':16}{'n':>4}{'noMSA':>8}{'MSA':>8}{'SHS':>8}{'MSA-noMSA':>11}{'SHS-noMSA':>11}")
for lab,lo,hi in [("< 0.45 (novel)",0,0.45),("0.45-0.7",0.45,0.7),("> 0.7 (similar)",0.7,2)]:
    m=[j for j,i in enumerate(ids) if lo<=stab[i]<hi]
    if not m: continue
    nom=np.mean([NOM[ids[j]] for j in m]); msa=np.mean([MSA[ids[j]] for j in m]); shs=np.mean([SHS[ids[j]] for j in m])
    print(f"{lab:16}{len(m):>4}{nom:>8.3f}{msa:>8.3f}{shs:>8.3f}{msa-nom:>+11.3f}{shs-nom:>+11.3f}")
