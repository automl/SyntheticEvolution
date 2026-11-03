#!/usr/bin/env python3
"""Report DB coverage per condition: expected 4 seeds x 5 samples, rows (targets) per DB vs. the target set.
Writes results/completeness.tsv; prints conditions that are not complete."""
import csv, os, sqlite3, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
from paths import CRC, DB_DIR, SEEDS  # noqa: E402

N = {"single": 91, "complex": 377}
rows, bad = [], 0
for c in csv.DictReader(open(CRC / "eval" / "conditions.tsv"), delimiter="\t"):
    tbl = "pred_rna_rna" if c["idset"] == "single" else "pred_protein_rna"
    for seed in SEEDS:
        for m in range(5):
            name = f'{c["tag"]}_{c["idset"]}_s{seed}_m{m}'; db = Path(DB_DIR) / f"{name}.db"
            n = sqlite3.connect(f"file:{db}?mode=ro", uri=True).execute(f"select count(*) from {tbl}").fetchone()[0] if db.exists() else 0
            rows.append([c["tag"], c["idset"], seed, m, n, N[c["idset"]]])
            if n < N[c["idset"]]:
                bad += 1; print(f"incomplete: {name}: {n}/{N[c['idset']]} targets")
with open(CRC / "results" / "completeness.tsv", "w", newline="") as f:
    csv.writer(f, delimiter="\t", lineterminator="\n").writerows([["tag", "idset", "seed", "sample", "targets", "expected"]] + rows)
print(f"{len(rows)} shards checked, {bad} incomplete (dssr single is 89/91 by design: 7EFG, 7VFT never folded)")
