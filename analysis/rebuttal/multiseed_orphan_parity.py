#!/usr/bin/env python3
"""Multiseed redo of the uTY4-Q1 orphan-parity paired tests.
SHS (rnaformer) and MSA (af3) = OUR multiseed ensemble, averaged per exp_db_id over
seeds {0,42,137} x samples {0..4} (from reproduce_rebuttal/results/benchmark_multiseed/*.db).
noMSA = single-seed published (low-variance floor; benchmark noMSA multiseed not folded).
Orphan subsets = paper's a2021 (post-cutoff) orphan IDs. Metric: single=RNA_TM, complex=Dockq.
Compares to the single-seed paired_stats result. Writes tables/multiseed_orphan_parity.csv."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import sqlite3, csv, re, glob, os, numpy as np, statistics as st
from scipy.stats import wilcoxon
DB=DB_DIR
SP=f'{OUT_ROOT}/paper_splits'
OUT=f'{OUT_ROOT}/multiseed_orphan_parity.csv'
RNG=np.random.default_rng(0)

def idof(s):
    m=re.search(r'([0-9][A-Za-z0-9]{3})', s or ''); return m.group(1).lower() if m else None
def fnum(x):
    try: return float(x)
    except: return None

def ms_mean(method, idset, col):
    """per-id list of values across all seed x sample DBs, then mean per id."""
    tbl='pred_rna_rna' if idset=='single' else 'pred_protein_rna'
    acc={}
    for db in glob.glob(f'{DB}/{method}_{idset}_s*_m*.db'):
        try: c=sqlite3.connect(db)
        except: continue
        for fn,v in c.execute(f'SELECT FileName,{col} FROM {tbl}'):
            i=idof(fn); v=fnum(v)
            if i and v is not None and (col!='RNA_TM' or v!=0): acc.setdefault(i,[]).append(v)
    return {i:st.mean(vs) for i,vs in acc.items() if vs}

def pub_ids(fname, col='exp_db_id'):
    p=f'{SP}/{fname}.csv'
    if not os.path.exists(p): return set()
    d=';' if ';' in open(p).readline() else ','
    return {str(r.get(col,'')).lower() for r in csv.DictReader(open(p),delimiter=d) if r.get(col)}
def pub_map(fname, metric):
    p=f'{SP}/{fname}.csv'
    if not os.path.exists(p): return {}
    d=';' if ';' in open(p).readline() else ','
    out={}
    for r in csv.DictReader(open(p),delimiter=d):
        i=str(r.get('exp_db_id','')).lower(); v=fnum(r.get(metric))
        if i and v is not None: out[i]=v
    return out

def paired(a,b):
    d=np.array(a)-np.array(b)
    p=wilcoxon(a,b,zero_method='wilcox').pvalue if (len(d)>=6 and not np.allclose(d,0)) else float('nan')
    bs=[RNG.choice(d,len(d),replace=True).mean() for _ in range(5000)]
    lo,hi=np.percentile(bs,[2.5,97.5])
    return d.mean(),p,lo,hi

# multiseed means
shs_s=ms_mean('rnaformer','single','RNA_TM'); msa_s=ms_mean('af3','single','RNA_TM')
shs_c=ms_mean('rnaformer','complex','Dockq'); msa_c=ms_mean('af3','complex','Dockq')
# single-seed noMSA (floor)
nom_s=ms_mean('empty','single','RNA_TM')   # multiseed no-MSA (was single-seed CSV)
nom_c=ms_mean('empty','complex','Dockq')
# orphan id sets (a2021 post-cutoff)
orph_s=split_ids('RNA_Monomers_a2021_alphafold_orphan')
orph_c=split_ids('RNA-Protein_a2021_alphafold_orphan')

rows=[]
print(f"{'subset':26}{'n':>4}{'MSA':>7}{'SHS':>7}{'SHSms-MSAms':>13}{'p':>9}{'  boot95_CI':>18}")
def run(label, shs, msa, ids, ss_ref):
    keep=sorted(set(shs)&set(msa)&set(ids))
    if len(keep)<3: print(f"  {label}: only {len(keep)} matched"); return
    xs=[shs[i] for i in keep]; ys=[msa[i] for i in keep]
    dm,p,lo,hi=paired(xs,ys)
    print(f"{label:26}{len(keep):>4}{st.mean(ys):>7.3f}{st.mean(xs):>7.3f}{dm:>+13.4f}{p:>9.3f}  [{lo:+.3f},{hi:+.3f}]  (1seed:{ss_ref})")
    rows.append([label,len(keep),st.mean(ys),st.mean(xs),dm,p,lo,hi])
print("=== SHS vs MSA  (BOTH multiseed) on a2021 post-cutoff orphans ===")
run("single ORPHAN (RNA_TM)", shs_s, msa_s, orph_s, "+0.029 p=.24")
run("complex ORPHAN (Dockq)", shs_c, msa_c, orph_c, "+0.004 p=.80")
os.makedirs(os.path.dirname(OUT),exist_ok=True)
with open(OUT,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['subset','n','MSA_mean','SHS_mean','d_SHS-MSA_multiseed','wilcoxon_p','ci_lo','ci_hi']); w.writerows(rows)
print(f"\nwrote {OUT}")
print("Read: CI spanning 0 => parity (the honest uTY4-Q1 answer), now seed-robust.")
