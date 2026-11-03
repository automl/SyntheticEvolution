#!/usr/bin/env python3
"""Write eval shard files from eval/conditions.tsv.

One shard = (tag, AF3 output dir, seed, sample, idset, class); 5 diffusion samples per seed (20 predictions per
target). With the unzipped outputs bundle the dirs are absolute bundle paths (af3/<tag>_<idset>/seed_1 | seeds_*);
otherwise seed 1 comes from `seed1_dir`, seeds 0/42/137 from `multiseed_dir` (relative to $AF_OUTPUT).

  python eval/make_shards.py --out eval/shards/all.tsv                       # every shard
  python eval/make_shards.py --out X.tsv --seeds 1 --tags rnaformer empty    # subset
  python eval/make_shards.py --out X.tsv --missing                           # only shards without a DB
"""
import argparse, csv, os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from paths import bundle_layout, prediction_root  # noqa: E402
ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True, type=Path)
ap.add_argument("--seeds", nargs="+", type=int, default=[1, 0, 42, 137])
ap.add_argument("--tags", nargs="+")
ap.add_argument("--idsets", nargs="+")
ap.add_argument("--missing", action="store_true", help="skip shards whose DB already exists in $EVAL_DB_DIR (default results/db)")
a = ap.parse_args()
db_dir = Path(os.environ.get("EVAL_DB_DIR", ROOT / "results" / "db"))   # where run_eval_array.sh writes

rows = []
for c in csv.DictReader(open(ROOT / "eval" / "conditions.tsv"), delimiter="\t"):
    if a.tags and c["tag"] not in a.tags: continue
    if a.idsets and c["idset"] not in a.idsets: continue
    for seed in a.seeds:
        d = prediction_root(c["tag"], c["idset"], seed) if bundle_layout() else c["seed1_dir"] if seed == 1 else c["multiseed_dir"]
        for sample in range(5):
            name = f'{c["tag"]}_{c["idset"]}_s{seed}_m{sample}'
            if a.missing and (db_dir / f"{name}.db").exists(): continue
            rows.append([c["tag"], d, seed, sample, c["idset"], c["class"]])
a.out.parent.mkdir(parents=True, exist_ok=True)
with open(a.out, "w", newline="") as f:
    csv.writer(f, delimiter="\t", lineterminator="\n").writerows(rows)
print(f"{len(rows)} shards -> {a.out}")
