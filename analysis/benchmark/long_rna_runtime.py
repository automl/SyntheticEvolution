#!/usr/bin/env python3
"""Long-RNA runtime: AF3 homolog search vs SHS generation (data/long_rna/pipeline_timing_comparison.csv).
 - homolog search as a fraction of the original AF3 end-to-end runtime (~87%),
 - SHS generation vs MSA search speedup (~73x on the search step),
 - the two longest RNAs (544, 746 nt) that found zero homologs despite long search."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import csv, statistics as st
rows = list(csv.DictReader(open(TIMING_CSV)))
frac = [float(r["orig_msa_s"]) / float(r["orig_total_s"]) for r in rows if float(r["orig_total_s"]) > 0]
msa = [float(r["orig_msa_s"]) for r in rows]
gen = [float(r["shs_gen_s"]) for r in rows]
print(f"n={len(rows)} long RNAs")
print(f"homolog search as fraction of original end-to-end runtime: mean={st.mean(frac)*100:.1f}%  median={st.median(frac)*100:.1f}%")
print(f"MSA search mean={st.mean(msa):.1f}s ; SHS generation mean={st.mean(gen):.1f}s ; ratio={st.mean(msa)/st.mean(gen):.0f}x")
zero = [r for r in rows if int(float(r.get("msa_depth", 0) or 0)) == 0]
print(f"RNAs with zero homologs (msa_depth=0): {len(zero)}; longest-search among them:")
for r in sorted(zero, key=lambda r: -float(r["orig_msa_s"]))[:4]:
    print(f"  {r['id']}: msa_depth=0 search={float(r['orig_msa_s']):.0f}s")
