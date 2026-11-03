import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import csv, statistics as St
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
SCR=CONSTRUCT_DIR
STEPS=["wc","nc","pk","all"]
rows=list(csv.DictReader(open(f"{SCR}/{CONSTRUCT_PREFIX}_summary.csv")))
by=defaultdict(dict); rec=defaultdict(dict)
for r in rows:
    by[r["pid"]][r["step"]]=float(r["mean_tm"]); rec[r["pid"]][r["step"]]=float(r["recall"])
pids=[p for p in by if all(s in by[p] for s in STEPS)]
def sem(v): return St.pstdev(v)/ (len(v)**0.5) if len(v)>1 else 0.0

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))
# LEFT: pooled scatter (recall, TM) 364 pts + binned mean
X=[rec[p][s] for p in pids for s in STEPS]; Y=[by[p][s] for p in pids for s in STEPS]
ax1.scatter(X,Y,s=10,alpha=0.25,color="#666666")
bins=np.linspace(0,1,9); bx=[];bm=[];be=[]
for i in range(len(bins)-1):
    v=[y for x,y in zip(X,Y) if bins[i]<=x<bins[i+1]]
    if len(v)>=5: bx.append((bins[i]+bins[i+1])/2); bm.append(St.mean(v)); be.append(sem(v))
ax1.errorbar(bx,bm,yerr=be,marker="o",color="#0072B2",lw=2,capsize=3,label="binned mean")
ax1.set_xlabel("Topology recall vs GT (precision ~ 1)"); ax1.set_ylabel("RNA TM-score to GT")
from scipy.stats import spearmanr
_r=spearmanr(X,Y)
ax1.set_title(f"91 RNAs, all steps pooled (n={len(X)}, SHS N={CONSTRUCT_N})\nSpearman rho={_r.correlation:.2f}, p={_r.pvalue:.0e}"); ax1.legend(); ax1.grid(alpha=0.3)
# RIGHT: within-RNA dTM vs wc baseline, mean +/- SEM per step
mr=[St.mean([rec[p][s] for p in pids]) for s in STEPS]
md=[St.mean([by[p][s]-by[p]["wc"] for p in pids]) for s in STEPS]
ed=[sem([by[p][s]-by[p]["wc"] for p in pids]) for s in STEPS]
ax2.axhline(0,color="grey",lw=0.8,ls="--")
ax2.errorbar(mr,md,yerr=ed,marker="o",ms=9,color="#D55E00",lw=2,capsize=4)
for x,y,s in zip(mr,md,STEPS): ax2.annotate(s,(x,y),textcoords="offset points",xytext=(6,6),fontsize=9)
ax2.set_xlabel("Topology recall vs GT (precision ~ 1)"); ax2.set_ylabel("Within-RNA dTM (vs WC-only step)")
from scipy.stats import wilcoxon
_up=sum(by[p]["all"]>by[p]["wc"] for p in pids); _w=wilcoxon([by[p]["all"] for p in pids],[by[p]["wc"] for p in pids])
ax2.set_title(f"Within-RNA: 3D improves on topology completion\nTM(all)>TM(wc): {_up}/{len(pids)}, Wilcoxon p={_w.pvalue:.0e}"); ax2.grid(alpha=0.3)
fig.suptitle("Construction sweep across 91 single-chain RNAs: AF3 3D quality tracks topology recall (precision fixed ~1)",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.96])
out=f"{OUT_ROOT}/{CONSTRUCT_PREFIX}_recall_TM.png"
fig.savefig(out,dpi=150)
print("plot ->",out)
