#!/usr/bin/env python3
"""Select AF3 input JSONs for a target set and (optionally) re-seed them.

Inputs are located as <src_root>/[config80/]<idx>/<PDBID>...json (one JSON per subdir; the first
four characters of the file name are the PDB id).

  list mode   (no --seeds): write --list with the matching subdirs relative to src_root (no copies).
  reseed mode (--seeds S..): write <dst_root>/<PDBID>/<same file name> with ONLY `modelSeeds`
              replaced, and --list with <PDBID>.
"""
import argparse, json, sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--src_root", required=True, type=Path)
ap.add_argument("--ids", required=True, type=Path, help="one PDB id per line")
ap.add_argument("--list", required=True, type=Path, help="output list file for fold_array.sh")
ap.add_argument("--seeds", nargs="+", type=int)
ap.add_argument("--dst_root", type=Path)
a = ap.parse_args()

ids = [l.strip().upper() for l in a.ids.read_text().split() if l.strip()]
found = {}
for j in sorted(a.src_root.glob("*/*.json")) + sorted(a.src_root.glob("config80/*/*.json")):
    pid = j.name[:4].upper()
    if pid in found and found[pid] != j:
        sys.exit(f"duplicate input for {pid}: {found[pid]} {j}")
    found[pid] = j
missing = [i for i in ids if i not in found]
if missing:
    sys.exit(f"{len(missing)} ids without input under {a.src_root}: {missing[:10]}")

lines = []
for pid in ids:
    src = found[pid]
    if a.seeds:
        d = json.loads(src.read_text())
        d["modelSeeds"] = a.seeds
        out = a.dst_root / pid / src.name
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(d, indent=2))
        lines.append(pid)
    else:
        lines.append(str(src.parent.relative_to(a.src_root)))
a.list.parent.mkdir(parents=True, exist_ok=True)
a.list.write_text("\n".join(lines) + "\n")
print(f"{len(lines)} inputs -> {a.list}" + (f" (re-seeded {a.seeds} into {a.dst_root})" if a.seeds else ""))
