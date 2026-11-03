#!/usr/bin/env python3
"""
Rewrite every JSON in a directory of numbered subfolders with a new set of seeds.

Reads each *.json inside the numbered subdirectories of --in_dir, replaces ONLY the
top-level `modelSeeds` value with the provided list of seeds, and writes the result
into --out_dir, mirroring the subfolder structure and file names exactly.

Nothing else in the JSON is touched: key order, formatting (indent=2), and every
other field are preserved as-is.

Usage:
    python set_seeds.py --in_dir config80 --seeds 42 137 0 --out_dir config80_multiseed
    python set_seeds.py --in_dir config80 --seeds 42,137,0 --out_dir config80_multiseed
"""

import argparse
import json
import os
from pathlib import Path


def long_path(p):
    """Return a path string usable past Windows' 260-char MAX_PATH limit.

    On Windows, prefix the absolute path with the extended-length marker
    ``\\?\``; on other platforms return the path unchanged. The very long
    rnaformer file names push full paths over MAX_PATH, so open()/mkdir()
    fail with FileNotFoundError without this.
    """
    ap = os.path.abspath(str(p))
    if os.name == "nt" and not ap.startswith("\\\\?\\"):
        return "\\\\?\\" + ap
    return ap


def parse_seeds(values):
    """Accept either space-separated (--seeds 42 137 0) or comma-separated
    (--seeds 42,137,0) seed lists, or a mix."""
    seeds = []
    for v in values:
        for part in v.split(","):
            part = part.strip()
            if part:
                seeds.append(int(part))
    if not seeds:
        raise argparse.ArgumentTypeError("no seeds provided")
    return seeds


def main():
    ap = argparse.ArgumentParser(description="Set modelSeeds on every JSON in numbered subfolders.")
    ap.add_argument("--in_dir", required=True, help="Input directory with numbered subfolders of JSONs.")
    ap.add_argument("--seeds", required=True, nargs="+", help="Seed list, e.g. --seeds 42 137 0")
    ap.add_argument("--out_dir", required=True, help="Output directory to write the updated tree into.")
    args = ap.parse_args()

    seeds = parse_seeds(args.seeds)
    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)

    if not in_dir.is_dir():
        raise SystemExit(f"--in_dir not found: {in_dir}")

    json_paths = sorted(in_dir.rglob("*.json"))
    if not json_paths:
        raise SystemExit(f"No .json files found under {in_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)

    for src in json_paths:
        with open(long_path(src), encoding="utf-8") as fh:
            data = json.load(fh)

        # Change ONLY the seeds. Everything else is left untouched.
        data["modelSeeds"] = seeds

        rel = src.relative_to(in_dir)
        dst = out_dir / rel
        os.makedirs(long_path(dst.parent), exist_ok=True)
        with open(long_path(dst), "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2)

    print(f"Updated {len(json_paths)} JSON files with seeds {seeds}")
    print(f"Written to {out_dir}")


if __name__ == "__main__":
    main()
