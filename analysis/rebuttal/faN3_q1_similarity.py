#!/usr/bin/env python3
"""faN3-Q1 post-vs-pre-cutoff SIMILARITY (leakage) quantification, reproducible for the revision.

For each subtype (single_chain / complex) and each POST-cutoff RNA, how similar is it to its
NEAREST PRE-cutoff neighbour?
  - sequence: cd-hit-est-2d, fraction of post with a pre neighbour at >=80% / >=90% identity.
  - structure: US-align max TM of each post to any pre GT structure (median, frac>0.5, frac>0.7).
Splits (pre=b2021, post=a2021) come from family_leakage.csv `era`+`subtype`, verified to match the
paper's a2021 split. Sequences come from the AF3 input JSONs (inputs bundle), structures from the reference
mmCIFs. Requires: cd-hit-est-2d (RnaBench env) and US-align (eval/build_usalign.sh).

Usage:
  conda activate RnaBench   # for cd-hit-est-2d
  python faN3_q1_similarity.py --mode sequence            # both subtypes, fast
  python faN3_q1_similarity.py --mode structure --subtype single_chain
  python faN3_q1_similarity.py --mode structure --subtype complex   # ~32k US-align pairs, slow
Outputs csv/json under results/analysis/faN3_q1/.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import argparse, json, glob, os, re, subprocess, statistics as st
import pandas as pd

FAM = FAM_CSV
GT = {"single_chain": GT_SINGLE, "complex": GT_COMPLEX}
USALIGN = USALIGN
OUT = f"{OUT_ROOT}/faN3_q1"; os.makedirs(OUT, exist_ok=True)


def splits(subtype):
    f = pd.read_csv(FAM); f = f[f.subtype == subtype]
    pre = sorted(set(f[f.era == "pre"].exp_db_id.str.lower()))
    post = sorted(set(f[f.era == "post"].exp_db_id.str.lower()))
    return pre, post


def rna_seq(idlow):
    fmap = rna_seq._map if hasattr(rna_seq, "_map") else None
    if fmap is None:
        fmap = {**input_jsons("single"), **input_jsons("complex")}
        rna_seq._map = fmap
    p = fmap.get(idlow.upper())
    if not p:
        return None
    d = json.load(open(p))
    for s in d.get("sequences", []):
        r = s.get("rna") or s.get("RNA")
        if r and r.get("sequence"):
            return r["sequence"]
    return None


def sequence_mode():
    rows = []
    for sub in ("single_chain", "complex"):
        pre, post = splits(sub)
        for name, ids in (("pre", pre), ("post", post)):
            recs = [f">{i}\n{rna_seq(i)}" for i in ids if rna_seq(i)]
            open(f"{OUT}/{sub}_{name}.fasta", "w").write("\n".join(recs) + "\n")
        for c, n in (("0.90", 8), ("0.80", 7)):
            subprocess.run(["cd-hit-est-2d", "-i", f"{OUT}/{sub}_pre.fasta", "-i2", f"{OUT}/{sub}_post.fasta",
                            "-o", f"{OUT}/{sub}_novel_{c}", "-c", c, "-n", str(n), "-M", "8000", "-T", "4"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            novel = sum(1 for l in open(f"{OUT}/{sub}_novel_{c}") if l.startswith(">"))
            tot = sum(1 for l in open(f"{OUT}/{sub}_post.fasta") if l.startswith(">"))
            hit = tot - novel
            print(f"  {sub:12} seq-id>={c}: {hit}/{tot} post have a pre neighbour ({novel}/{tot} sequence-novel)")
            rows.append({"subtype": sub, "threshold": c, "n_post": tot, "with_pre_neighbour": hit, "novel": novel})
    pd.DataFrame(rows).to_csv(f"{OUT}/sequence_similarity.csv", index=False)
    print(f"  wrote {OUT}/sequence_similarity.csv")


def cif(gtdir, i):
    g = glob.glob(f"{gtdir}/{i}*.cif")
    return g[0] if g else None


def structure_mode(subtype):
    gtdir = GT[subtype]; pre, post = splits(subtype)
    res = {}
    for po in post:
        pc = cif(gtdir, po)
        if not pc:
            continue
        best = 0.0
        for pr in pre:
            rc = cif(gtdir, pr)
            if not rc:
                continue
            try:
                out = subprocess.run([USALIGN, pc, rc], capture_output=True, text=True, timeout=120).stdout
            except Exception:
                continue
            tms = [float(x) for x in re.findall(r"TM-score=\s*([0-9.]+)", out)]
            if tms:
                best = max(best, max(tms))
        res[po] = best
        print(f"{po}\t{best:.3f}", flush=True)
    json.dump(res, open(f"{OUT}/structure_maxTM_{subtype}.json", "w"))
    v = list(res.values())
    if v:
        print(f"# {subtype}: post max-TM to nearest pre  n={len(v)} median={st.median(v):.2f} "
              f">0.5={sum(x>0.5 for x in v)}/{len(v)} >0.7={sum(x>0.7 for x in v)}/{len(v)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["sequence", "structure"], required=True)
    ap.add_argument("--subtype", choices=["single_chain", "complex"], default="single_chain")
    a = ap.parse_args()
    if a.mode == "sequence":
        sequence_mode()
    else:
        structure_mode(a.subtype)
