#!/usr/bin/env python
"""4P8Z/6WJS intervention sweep, MULTISEED: pool seed-1 (single-seed dir) + seeds
0/42/137 (multiseed dir) per step over all depths x samples. Per-step mean TM and
SD across the 4 per-seed means -> error bars. Regenerate construction plot."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import glob, re, subprocess, os
from multiprocessing import Pool
from collections import defaultdict
import statistics as St

SCR = CONSTRUCT_DIR
UA = USALIGN

TARGETS = {
 "4P8Z": {"gt": f"{GT_INTERVENTION}/4p8z_gt.cif",
   "single": case_study_dir("4P8Z_intervention_full_sweep"),
   "multi": case_study_dir("4P8Z_intervention_full_sweep_multiseed"),
   "steps": ["local_wc","init_wc","step1_nc","step2_pk","step3_multi","final_all"],
   "recall": [0.27,0.52,0.61,0.70,0.88,1.00]},
 "6WJS": {"gt": f"{GT_INTERVENTION}/6wjs_gt.cif",
   "single": case_study_dir("6WJS_intervention_full_sweep"),
   "multi": case_study_dir("6WJS_intervention_full_sweep_multiseed"),
   "steps": ["local_wc","init_wc","step1_nc","step2_pk","step3_multi","final_all"],
   "recall": [0.42,0.56,0.62,0.76,0.84,1.00]},
}

def tm(args):
    cif, gt = args
    try:
        out = subprocess.run([UA, cif, gt, "-mol","RNA","-ter","0"],
                             capture_output=True, text=True, timeout=120).stdout
    except Exception:
        return (cif, None)
    for line in out.splitlines():
        if "normalized by length of Structure_2" in line:
            m = re.search(r"TM-score=\s*([0-9.]+)", line)
            if m: return (cif, float(m.group(1)))
    return (cif, None)

def collect(target):
    cfg = TARGETS[target]; gt = cfg["gt"]
    jobs = []
    meta = {}
    for step in cfg["steps"]:
        for base in (cfg["single"], cfg["multi"]):
            for cif in glob.glob(f"{base}/*intervention_{step}/seed-*_sample-*/model.cif"):
                jobs.append((cif, gt)); meta[cif] = step
    with Pool(16) as p:
        res = p.map(tm, jobs)
    # group by (step, seed)
    byseedstep = defaultdict(list)
    for cif, t in res:
        if t is None: continue
        step = meta[cif]
        sd = re.search(r"seed-(\d+)_sample", cif)
        seed = sd.group(1) if sd else "?"
        byseedstep[(step, seed)].append(t)
    return byseedstep

def main():
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2,5.2))
    colors = {"4P8Z":"#0072B2","6WJS":"#D55E00"}
    summary_lines = []
    for target in ("4P8Z","6WJS"):
        cfg = TARGETS[target]
        bss = collect(target)
        xs, ys, es = [], [], []
        print(f"\n=== {target} (multiseed, 4 seeds) ===")
        for step, rec in zip(cfg["steps"], cfg["recall"]):
            seedmeans = [St.mean(bss[(step,s)]) for s in ("1","0","42","137") if bss.get((step,s))]
            npred = sum(len(bss[(step,s)]) for s in ("1","0","42","137") if bss.get((step,s)))
            if not seedmeans: continue
            m = St.mean(seedmeans); sd = St.pstdev(seedmeans) if len(seedmeans)>1 else 0.0
            xs.append(rec); ys.append(m); es.append(sd)
            print(f"  {step:12s} recall={rec:.2f} meanTM={m:.3f} +/-{sd:.3f} (seeds={len(seedmeans)}, n={npred})")
            summary_lines.append(f"{target},{step},{rec},{m:.4f},{sd:.4f},{len(seedmeans)},{npred}")
        ax.errorbar(xs, ys, yerr=es, marker="o", ms=7, lw=2, capsize=4,
                    color=colors[target], label=target)
        for x,y,st in zip(xs,ys,cfg["steps"]):
            ax.annotate(st.replace("_","\n"), (x,y), textcoords="offset points",
                        xytext=(6,-4), fontsize=6.5, color=colors[target], alpha=0.8)
    ax.set_xlabel("Topology recall vs ground truth (precision ~ 1 throughout)")
    ax.set_ylabel("RNA TM-score to ground truth")
    ax.set_title("AF3 3D quality tracks topology recall, not precision\n"
                 "(multiseed: mean +/- SD across 4 AF3 seeds)")
    ax.grid(alpha=0.3); ax.legend()
    ax.set_xlim(0,1.05)
    out = f"{OUT_ROOT}/construction_recall_TM_multiseed.png"
    fig.tight_layout(); fig.savefig(out, dpi=150)
    with open(f"{SCR}/construction_ms_summary.csv","w") as f:
        f.write("target,step,recall,mean_tm,sd_tm,nseed,npred\n")
        f.write("\n".join(summary_lines)+"\n")
    print(f"\nplot -> {out}")

if __name__=="__main__":
    main()
