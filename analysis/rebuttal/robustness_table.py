#!/usr/bin/env python3
"""Rebuttal 'correction & robustness' table. For each benchmark method, compares the
REPORTED single seed (seed-1, sample-0 = published) against the full 4-seed picture
(seed 1 published + our seeds 0/42/137), all at sample-0 (paper convention), plus the
across-sample spread. Shows means change by <=~0.01 and nothing flips.
single metric = RNA_TM (n=91); complex metric = DockQ (n=377)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, csv, re, glob, os, numpy as np, statistics as st
DB=DB_DIR; SP=f'{OUT_ROOT}/paper_splits'
OUT=f'{OUT_ROOT}/seed_robustness_table.csv'
def idof(s):
    m=re.search(r'([0-9][A-Za-z0-9]{3})',s or ''); return m.group(1).lower() if m else None
def fnum(x):
    try: return float(x)
    except: return None
def pub(fname, metric):  # submitted seed-1 sample-0 values (published CSVs), restricted to the paper split
    import gzip
    pfx = fname.split('_', 2)[-1] if fname.startswith('RNA_Monomers_') else fname.split('_', 1)[-1]
    suf = 'rna' if fname.startswith('RNA_Monomers_') else 'protein_rna'
    ids = split_ids('rna_monomers_ids' if suf == 'rna' else 'protein_rna_ids')
    p = CRC/'data'/'published_seed1'/f'{pfx}_{suf}.csv.gz'
    if not p.exists(): return {}
    lines = gzip.open(p, 'rt').read().splitlines(); d = ';' if ';' in lines[0] else ','
    return {str(r.get('exp_db_id','')).lower(): fnum(r.get(metric)) for r in csv.DictReader(lines, delimiter=d)
            if str(r.get('exp_db_id','')).lower() in ids and fnum(r.get(metric)) is not None}
def our_seed_sample(method, idset, col, seed, sample):  # id->val
    tbl='pred_rna_rna' if idset=='single' else 'pred_protein_rna'
    db=f'{DB}/{method}_{idset}_s{seed}_m{sample}.db'
    if not os.path.exists(db): return {}
    c=sqlite3.connect(db); out={}
    for fn,v in c.execute(f'SELECT FileName,{col} FROM {tbl}'):
        i=idof(fn); v=fnum(v)
        if i and v is not None and (col!='RNA_TM' or v!=0): out[i]=v
    return out

PUBPFX={'rnaformer':'rnaformerN100','dssr':'dssrN100','rnafold':'rnafoldN100','spotrna':'spotrnaN100','af3':'alphafold'}
rows=[]
print(f"{'method':10}{'idset':8}{'metric':7}{'seed1(rep)':>11}{'seed0':>8}{'seed42':>8}{'seed137':>9}{'4seed_mean':>11}{'seed_sd':>8}{'samp_sd':>8}{'|Δ|':>7}")
for method,pfx in PUBPFX.items():
    for idset in ['single','complex']:
        if method=='dssr' and idset=='complex': continue
        metric='RNA_TM' if idset=='single' else 'Dockq'
        pubf=('RNA_Monomers_' if idset=='single' else 'RNA-Protein_')+pfx
        p=pub(pubf,metric)
        seedmeans={1: (st.mean(p.values()) if p else None)}
        # per-seed sample-0 means (ours)
        for s in (0,42,137):
            d=our_seed_sample(method,idset,metric,s,0)
            seedmeans[s]=st.mean(d.values()) if d else None
        vals=[seedmeans[s] for s in (1,0,42,137) if seedmeans[s] is not None]
        if len(vals)<4:
            print(f"{method:10}{idset:8}{metric:7}  incomplete ({len(vals)}/4 seeds)"); continue
        m4=st.mean(vals); seed_sd=st.pstdev(vals)
        # across-sample spread: per id std over 5 samples, averaged over seeds+ids
        samp_sds=[]
        for s in (0,42,137):
            perid={}
            for m in range(5):
                for i,v in our_seed_sample(method,idset,metric,s,m).items(): perid.setdefault(i,[]).append(v)
            samp_sds+=[st.pstdev(vs) for vs in perid.values() if len(vs)>=2]
        samp_sd=st.mean(samp_sds) if samp_sds else float('nan')
        dabs=abs(m4-seedmeans[1])
        print(f"{method:10}{idset:8}{metric:7}{seedmeans[1]:>11.3f}{seedmeans[0]:>8.3f}{seedmeans[42]:>8.3f}{seedmeans[137]:>9.3f}{m4:>11.3f}{seed_sd:>8.3f}{samp_sd:>8.3f}{dabs:>7.3f}")
        rows.append([method,idset,metric,seedmeans[1],seedmeans[0],seedmeans[42],seedmeans[137],m4,seed_sd,samp_sd,dabs])
os.makedirs(os.path.dirname(OUT),exist_ok=True)
with open(OUT,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['method','idset','metric','seed1_reported','seed0','seed42','seed137','mean_4seed','seed_std','sample_std','abs_delta_vs_reported']); w.writerows(rows)
print(f"\nwrote {OUT}")
print("seed_sd = std across the 4 seed-means (sample-0); samp_sd = mean per-id std across 5 samples;")
print("|Δ| = |4-seed mean - reported single seed|.  All small => reported results are seed/sample robust.")
