#!/usr/bin/env python
"""Build the 91-RNA construction experiment.
 - seq+pairs from dssr (nat_meth primary, a2021/b2021 fallback)
 - classify GT pairs into cumulative recall sweep: wc -> +nc -> +pk -> +all
 - emit construct_pairs_91.json  (list of {id:"<pid>_<step>", pdb_id, pairs})
 - synthesize a minimal AF3 data.json per RNA (seq carrier for the generator)
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import json, os, ast

RB = RNABENCH_HOME
SCR = CONSTRUCT_DIR
IDS = str(IDS_DIR/"single_91.txt")
GTDIR = GT_SINGLE
DATAOUT = f"{SCR}/construct91_datafiles"
DSSR_PRIMARY = f"{RB}/dssr_secondary_structure_data_nat_meth.json"
DSSR_FALL = [f"{RB}/a2021_secondary_structure_data_dssr.json",
             f"{RB}/b2021_secondary_structure_data_dssr.json"]
STEPS = ["wc","nc","pk","all"]
WC = {("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}

def load_dssr(files):
    m = {}
    for f in files:
        for r in json.load(open(f)):
            pid = r.get("pdb_id","").lower()
            if pid in m: continue
            pr = r.get("pairs")
            if isinstance(pr, str): pr = ast.literal_eval(pr)
            pairs = [(int(a),int(b)) for a,b in pr]
            m[pid] = {"seq": r.get("sequence","").upper().replace("T","U"), "pairs": pairs}
    return m

def crosses(a,b):
    (i,j),(k,l)=a,b
    return (i<k<j<l) or (k<i<l<j)

def classify(seq, pairs):
    P = sorted({(min(a,b),max(a,b)) for a,b in pairs})
    deg={}
    for i,j in P: deg[i]=deg.get(i,0)+1; deg[j]=deg.get(j,0)+1
    R=[p for p in P if deg[p[0]]==1 and deg[p[1]]==1]
    nested=[]; pk=[]
    for p in sorted(R, key=lambda x:(x[0],x[1])):
        (pk if any(crosses(p,q) for q in nested) else nested).append(p)
    def canon(p):
        i,j=p
        return i<len(seq) and j<len(seq) and (seq[i],seq[j]) in WC
    wc=[p for p in nested if canon(p)]; nc=[p for p in nested if not canon(p)]
    S={"wc":sorted(wc),"nc":sorted(wc+nc),"pk":sorted(wc+nc+pk),"all":sorted(P)}
    stats={"n_all":len(P),"n_multi":len(P)-len(R),"n_pk":len(pk),"n_nc":len(nc),"n_wc":len(wc)}
    return S, stats

def synth_datajson(pid, seq):
    return {"dialect":"alphafold3","version":1,"name":pid.upper(),
            "sequences":[{"rna":{"id":"A","sequence":seq,"modifications":[],
                                 "unpairedMsa":f">query\n{seq}"}}],
            "modelSeeds":[1,0,42,137],"bondedAtomPairs":None,"userCCD":None}

def main():
    ids=[x.strip().lower() for x in open(IDS) if x.strip()]
    dssr=load_dssr([DSSR_PRIMARY]+DSSR_FALL)
    os.makedirs(DATAOUT, exist_ok=True)
    out=[]; audit=[]; miss=[]
    for pid in ids:
        rec=dssr.get(pid)
        if rec is None: miss.append((pid,"no_dssr")); continue
        if not os.path.exists(f"{GTDIR}/{pid}.cif"): miss.append((pid,"no_gt")); continue
        seq=rec["seq"]
        if not seq: miss.append((pid,"empty_seq")); continue
        S,st=classify(seq, rec["pairs"])
        nall=st["n_all"] or 1
        rs={k:round(len(v)/nall,3) for k,v in S.items()}
        for step in STEPS:
            out.append({"id":f"{pid}_{step}","pdb_id":pid,"pairs":[[i,j] for i,j in S[step]]})
        json.dump(synth_datajson(pid,seq), open(f"{DATAOUT}/{pid}_data.json","w"), indent=2)
        audit.append((pid,len(seq),st,rs))
    op=f"{SCR}/construct_pairs_91.json"
    json.dump(out, open(op,"w"))
    print(f"usable={len(audit)}/91 missing={len(miss)} -> {len(out)} step-records")
    print("datafiles ->", DATAOUT, " pairs ->", op)
    if miss: print("MISSING:", miss)
    print("\nSample recall sweeps (wc/nc/pk/all as fraction of GT pairs):")
    for pid,L,st,rs in audit[:10]:
        print(f"  {pid} L={L:3d} nGT={st['n_all']:3d} multi={st['n_multi']:2d} pk={st['n_pk']:2d} nc={st['n_nc']:2d} | {rs['wc']:.2f} {rs['nc']:.2f} {rs['pk']:.2f} {rs['all']:.2f}")
    import statistics as St
    wcr=[rs['wc'] for *_,rs in audit]
    print(f"\nwc-start recall: min={min(wcr):.2f} med={St.median(wcr):.2f} max={max(wcr):.2f}")
    print(f"RNAs w/ multiplets: {sum(1 for *_,st,_ in [(a[0],a[1],a[2],a[3]) for a in audit] if st['n_multi']>0)}/{len(audit)}",
          f"| w/ pseudoknots: {sum(1 for a in audit if a[2]['n_pk']>0)}/{len(audit)}")
    # how many RNAs have a genuine sweep (wc<all)?
    movers=sum(1 for a in audit if a[3]['wc']<0.999)
    print(f"RNAs where wc<all (real sweep room): {movers}/{len(audit)}")

if __name__=="__main__":
    main()
