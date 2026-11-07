#!/usr/bin/env python3
"""Sample-dimension robustness. For each key comparison, the paired delta at EACH diffusion sample m in {0..4}
(per-target value = mean over the 4 AF3 seeds at that sample) and over all 20 predictions, showing the delta is
stable across samples. Metric RNA_TM (single/ribozyme) / DockQ (complex orphans)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, csv, re, glob, os, numpy as np, statistics as st
from scipy.stats import wilcoxon
DB=DB_DIR; SP=f'{OUT_ROOT}/paper_splits'
MD=MD_CSV; RNG=np.random.default_rng(0)
def idof(s):
    m=re.search(r'([0-9][A-Za-z0-9]{3})',s or ''); return m.group(1).lower() if m else None
def fnum(x):
    try: return float(x)
    except: return None
def ms_mean(method, idset, col, sample=None):
    """per-id mean over seeds; if sample given, only that sample; else all 20."""
    tbl='pred_rna_rna' if idset=='single' else 'pred_protein_rna'
    pat=f'{DB}/{method}_{idset}_s*_m{sample}.db' if sample is not None else f'{DB}/{method}_{idset}_s*_m*.db'
    acc={}
    for db in glob.glob(pat):
        try: c=sqlite3.connect(db)
        except: continue
        for fn,v in c.execute(f'SELECT FileName,{col} FROM {tbl}'):
            i=idof(fn); v=fnum(v)
            if i and v is not None and (col!='RNA_TM' or v!=0): acc.setdefault(i,[]).append(v)
    return {i:st.mean(vs) for i,vs in acc.items() if vs}
def pub_map(fname, metric):
    p=f'{SP}/{fname}.csv'
    if not os.path.exists(p): return {}
    d=';' if ';' in open(p).readline() else ','
    return {str(r.get('exp_db_id','')).lower():fnum(r.get(metric)) for r in csv.DictReader(open(p),delimiter=d) if fnum(r.get(metric)) is not None}
def delta(shs, base, ids):
    keep=sorted(set(shs)&set(base)&set(ids))
    if len(keep)<3: return None
    d=np.array([shs[i] for i in keep])-np.array([base[i] for i in keep])
    p=wilcoxon([shs[i] for i in keep],[base[i] for i in keep]).pvalue if (len(d)>=6 and not np.allclose(d,0)) else float('nan')
    return d.mean(), p, len(keep)

# id subsets
orph_s=split_ids('rna_monomers_a2021_orphan_ids')
orph_c=split_ids('protein_rna_a2021_orphan_ids')
allsingle=split_ids('rna_monomers_ids')
d=';' if ';' in open(MD).readline() else ','
rib={str(r['exp_db_id']).lower() for r in csv.DictReader(open(MD),delimiter=d) if 'ribozyme' in (r.get('RNADescription') or '').lower()}
nom_c=ms_mean('empty','complex','RNA_TM'); hard={i for i in rib if nom_c.get(i,1)<0.5}

CLAIMS=[
 ("orphan single SHS-MSA (RNA_TM)","rnaformer","af3","single","RNA_TM",orph_s,"MSA"),
 ("orphan complex SHS-MSA (Dockq)","rnaformer","af3","complex","Dockq",orph_c,"MSA"),
 ("single gain SHS-noMSA (RNA_TM)","rnaformer",None,"single","RNA_TM",allsingle,"noMSA"),
 ("ribozyme SHS-MSA (RNA_TM)","rnaformer","af3","complex","RNA_TM",hard,"MSA"),
]
nom_s=ms_mean('empty','single','RNA_TM')
print(f"{'claim':34}{'s0':>11}{'s1':>8}{'s2':>8}{'s3':>8}{'s4':>8}{'mean20':>9}{'spread':>8}")
for lab,mth,base_mth,idset,col,ids,basekind in CLAIMS:
    ds=[]
    for m in range(5):
        shs=ms_mean(mth,idset,col,sample=m)
        base = ms_mean(base_mth,idset,col,sample=m) if base_mth else (nom_s if idset=='single' else nom_c)
        r=delta(shs,base,ids); ds.append(r[0] if r else float('nan'))
    shs15=ms_mean(mth,idset,col); base15=ms_mean(base_mth,idset,col) if base_mth else (nom_s if idset=='single' else nom_c)
    r15=delta(shs15,base15,ids)
    spread=(max(ds)-min(ds)) if all(x==x for x in ds) else float('nan')
    print(f"{lab:34}"+"".join(f"{x:>+8.4f}" if i==0 else f"{x:>+8.4f}" for i,x in enumerate(ds[:1]))+"".join(f"{x:>+8.4f}" for x in ds[1:])+f"{r15[0]:>+9.4f}{spread:>8.4f}")
print("\n(each column = paired mean delta at that sample, averaged over 4 seeds; spread=max-min across samples)")
