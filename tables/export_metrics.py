#!/usr/bin/env python3
"""Flatten all metric DBs into one small per-prediction table (git-tracked):
results/metrics_per_prediction.csv.gz  -- tag, idset, target, seed, sample, metrics (+ AF3 confidences)."""
import csv, glob, gzip, io, os, re, sqlite3, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
from paths import CRC, DB_DIR, DB_LONG, DB_RHOFOLD, idof  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_crc_tables import COMPLEX, SINGLE, num  # noqa: E402

EXTRA = ["iLDDT", "af3_rna_pTM", "af3_rna_ipTM", "af3_ranking_score", "af3_rna_pLDDT_avg"]
COLS = SINGLE + [c for c in COMPLEX if c not in SINGLE] + EXTRA
out = CRC / "results" / "metrics_per_prediction.csv.gz"
n = 0
with open(out, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz, \
        io.TextIOWrapper(gz, encoding="utf-8", newline="") as f:  # mtime=0: byte-identical reruns
    w = csv.writer(f); w.writerow(["tag", "idset", "target", "seed", "sample"] + COLS)
    for d, kind in ((DB_DIR, None), (DB_RHOFOLD, None), (DB_LONG, "long")):
        for db in sorted(glob.glob(f"{d}/*.db")):
            m = re.match(r"(.+)_(single|complex)_s(\d+)_m(\d)\.db$", os.path.basename(db)) if kind is None else \
                re.match(r"(.+)_s(\d+)_m(\d)\.db$", os.path.basename(db))
            if not m: continue
            tag, idset, seed, sample = (m.group(1), m.group(2), m.group(3), m.group(4)) if kind is None else \
                                       (m.group(1), "long", m.group(2), m.group(3))
            tbl = "pred_protein_rna" if idset == "complex" else "pred_rna_rna"
            con = sqlite3.connect(f"file:{db}?mode=ro", uri=True); con.row_factory = sqlite3.Row
            cols = {r[1] for r in con.execute(f"pragma table_info({tbl})")}
            for r in con.execute(f"SELECT * FROM {tbl}"):
                w.writerow([tag, idset, idof(r["FileName"]), seed, sample] +
                           [("" if (v := num(r[c]) if c in cols else None) is None else f"{v:.5g}") for c in COLS]); n += 1
print(f"{n} predictions -> {out} ({out.stat().st_size/1e6:.1f} MB)")
