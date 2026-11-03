import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import csv, statistics as St
from collections import defaultdict
SCR=CONSTRUCT_DIR
STEPS=["wc","nc","pk","all"]
rows=list(csv.DictReader(open(f"{SCR}/{CONSTRUCT_PREFIX}_summary.csv")))
by=defaultdict(dict); rec=defaultdict(dict)
for r in rows:
    by[r["pid"]][r["step"]]=float(r["mean_tm"]); rec[r["pid"]][r["step"]]=float(r["recall"])
pids=[p for p in by if all(s in by[p] for s in STEPS)]
print(f"RNAs with all 4 steps: {len(pids)}")
from scipy.stats import wilcoxon, spearmanr, binomtest
d=[by[p]["all"]-by[p]["wc"] for p in pids]
pos=sum(1 for x in d if x>0); neg=sum(1 for x in d if x<0)
print(f"\n[all vs wc] mean Delta={St.mean(d):+.3f} median={St.median(d):+.3f} up:{pos}/{len(d)} down:{neg}")
print(f"  Wilcoxon p={wilcoxon(d).pvalue:.2e}  sign p={binomtest(pos,pos+neg,0.5).pvalue:.2e}")
rhos=[]
for p in pids:
    xs=[rec[p][s] for s in STEPS]; ys=[by[p][s] for s in STEPS]
    if len(set(xs))>1:
        r,_=spearmanr(xs,ys)
        if r==r: rhos.append(r)
posr=sum(1 for r in rhos if r>0)
print(f"\n[per-RNA Spearman] median rho={St.median(rhos):+.3f} positive:{posr}/{len(rhos)} sign p={binomtest(posr,len(rhos),0.5).pvalue:.2e}")
print("\n[within-RNA step dTM vs wc baseline]:")
for s in STEPS:
    dd=[by[p][s]-by[p]["wc"] for p in pids]
    print(f"  {s}: dTM={St.mean(dd):+.3f} meanRecall={St.mean([rec[p][s] for p in pids]):.3f}")
sub=[p for p in pids if rec[p]["wc"]<0.6]
d2=[by[p]["all"]-by[p]["wc"] for p in sub]
print(f"\n[wc-recall<0.6, n={len(sub)}] mean d(all-wc)={St.mean(d2):+.3f} median={St.median(d2):+.3f} up:{sum(1 for x in d2 if x>0)}/{len(d2)} Wilcoxon p={wilcoxon(d2).pvalue:.2e}")
