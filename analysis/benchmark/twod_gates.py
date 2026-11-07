import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, glob, re, csv, statistics as st
import pandas as pd, numpy as np
from scipy.stats import wilcoxon
DB=DB_DIR; MC=MC_CSV
MD=MD_CSV; FAM=FAM_CSV
BINS=[-1.01,.5,.8,.95,1.01]; LAB=["<0.5","0.5-0.8","0.8-0.95",">0.95"]
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
S_s=ms("rnaformer","single","RNA_TM"); M_s=ms("af3","single","RNA_TM"); N_s=ms("empty","single","RNA_TM")
S_c=ms("rnaformer","complex","Dockq"); M_c=ms("af3","complex","Dockq"); N_c=ms("empty","complex","Dockq")
w0=pd.read_csv(MC); w0["idl"]=w0.exp_db_id.astype(str).str.lower()
w0.loc[w0.group=="rna-only","rf2d_selfcons"]=w0.loc[w0.group=="rna-only","idl"].map(selfcons_means())  # 4 seeds x 5 samples
def strat(df,gate,S,M,N,lab):
    w=df.dropna(subset=[gate]).copy()
    w["S"]=w.idl.map(S);w["M"]=w.idl.map(M);w["N"]=w.idl.map(N);w=w.dropna(subset=["S","M","N"]);w["bin"]=pd.cut(w[gate],BINS,labels=LAB)
    print(f"\n## {lab} (rows={len(w)}) ##  bin n | noMSA SHS MSA | SHS-noMSA")
    for b in LAB:
        g=w[w.bin==b]
        if len(g)<2: continue
        print(f"  {b:9} {len(g):>3} | {g.N.mean():.3f} {g.S.mean():.3f} {g.M.mean():.3f} | {(g.S-g.N).mean():+.3f}")
strat(w0[w0.group=="rna-only"],"rf2d_pred",S_s,M_s,N_s,"SINGLE 2D-MCC")
# rf2d_pred = MCC of the encoded RNAformer topology vs DSSR reference (analysis/gates/build_rnaformer_2d_gt.py);
# available for single chains and 1v1 complexes
strat(w0[w0.group=="1v1"],"rf2d_pred",S_c,M_c,N_c,"COMPLEX (1v1) 2D-MCC")
print("\n## SELF-CONSISTENCY (single) ##  cut n | noMSA SHS MSA | SHS-noMSA (p)")
ws=w0[w0.group=="rna-only"].dropna(subset=["rf2d_selfcons"]).copy(); ws["S"]=ws.idl.map(S_s);ws["M"]=ws.idl.map(M_s);ws["N"]=ws.idl.map(N_s)
for cut in (0.9,0.95):
    g=ws[ws.rf2d_selfcons>=cut].dropna(subset=["S","M","N"]); p=wilcoxon(g.S.values,g.N.values).pvalue
    print(f"  >={cut} {len(g):>3} | {g.N.mean():.3f} {g.S.mean():.3f} {g.M.mean():.3f} | {(g.S-g.N).mean():+.3f} (p={p:.4f})")
print("\n## GT-MCC cumulative (single) ##  cut n | noMSA SHS MSA | SHS-noMSA (p)")
wg=w0[w0.group=="rna-only"].dropna(subset=["rf2d_pred"]).copy(); wg["S"]=wg.idl.map(S_s);wg["M"]=wg.idl.map(M_s);wg["N"]=wg.idl.map(N_s)
for cut in (0.9,0.95):
    g=wg[wg.rf2d_pred>=cut].dropna(subset=["S","M","N"]); p=wilcoxon(g.S.values,g.N.values).pvalue
    print(f"  >={cut} {len(g):>3} | {g.N.mean():.3f} {g.S.mean():.3f} {g.M.mean():.3f} | {(g.S-g.N).mean():+.3f} (p={p:.4f})  SHS-MSA={ (g.S-g.M).mean():+.4f}")
print("\n## RIBOZYMES (autonomous, hard=multiseed noMSA<0.5) ##")
S_cr=ms("rnaformer","complex","RNA_TM"); M_cr=ms("af3","complex","RNA_TM"); N_cr=ms("empty","complex","RNA_TM")
d=";" if ";" in open(MD).readline() else ","
rib={str(r["exp_db_id"]).lower() for r in csv.DictReader(open(MD),delimiter=d) if "ribozyme" in (r.get("RNADescription") or "").lower()}
for crit,fl in [("multiseed noMSA<0.5",lambda i:N_cr.get(i,9)<0.5)]:
    hard=sorted(i for i in rib if i in S_cr and i in M_cr and i in N_cr and fl(i))
    a=[S_cr[i] for i in hard]; m=[M_cr[i] for i in hard]; n=[N_cr[i] for i in hard]
    pM=wilcoxon(a,m).pvalue; pN=wilcoxon(a,n).pvalue
    print(f"  {crit}: n={len(hard)} SHS={st.mean(a):.3f} MSA={st.mean(m):.3f} noMSA={st.mean(n):.3f}")
    print(f"    SHS-noMSA={st.mean(a)-st.mean(n):+.3f}(p={pN:.4f})  SHS-MSA={st.mean(a)-st.mean(m):+.3f}(p={pM:.4f})")
print("\n## FAMILY LEAKAGE ##")
fam=pd.read_csv(FAM); fam["id"]=fam.exp_db_id.astype(str).str.lower(); fam["seen"]=fam.fstatus.str.contains("seen-in-pre",na=False)
Nall={**N_s,**N_c}; Mall={**M_s,**M_c}; Sall={**S_s,**S_c}
fam["N"]=fam.id.map(Nall); fam["M"]=fam.id.map(Mall); fam["S"]=fam.id.map(Sall); fam=fam.dropna(subset=["N","M","S"])
for lab,g in [("seen",fam[fam.seen]),("unseen",fam[~fam.seen])]:
    print(f"  {lab:7} n={len(g):3d} | noMSA={g.N.mean():.3f} MSA={g.M.mean():.3f} SHS={g.S.mean():.3f}")
for c in ["N","M","S"]:
    b=fam[fam.seen][c].mean()-fam[~fam.seen][c].mean(); print(f"  boost {c}: {b:+.3f}")
sg=fam[fam.seen]; ug=fam[~fam.seen]
print(f"  MSA-SHS gap: unseen={(ug.M-ug.S).mean():+.3f}(p={wilcoxon(ug.M,ug.S).pvalue:.4f}) seen={(sg.M-sg.S).mean():+.3f}(p={wilcoxon(sg.M,sg.S).pvalue:.4f})")
