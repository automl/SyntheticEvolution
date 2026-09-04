#!/usr/bin/env python3
"""Prepare per-target trRosetta custom SS files from base-pair CSV annotations.

The generated files are numeric LxL matrices suitable for:
    predict.py -ss <file> -ss_fmt spot_prob

Important:
- trRosetta's predict.py applies `ss += ss.T` for `spot_prob` inputs.
- This script therefore writes only the upper triangle values to avoid doubling.
"""

from __future__ import annotations

import argparse
import ast
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create custom SS matrices for trRosetta from RNA_Monomers_base_pairs.csv.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--pairs-csv",
        default="./RNA_Monomers_base_pairs.csv",
        help="CSV path containing columns like pdb_id, method, pairs.",
    )
    parser.add_argument(
        "--fasta-dir",
        default="./prepared_shs_inputs/fasta",
        help="Directory with FASTA files used by the runner. Used to validate lengths.",
    )
    parser.add_argument(
        "--method",
        required=True,
        help="Method to extract from CSV (e.g. rnaformer_pred, rnafold_pred, spotrna_pred).",
    )
    parser.add_argument(
        "--pairs-column",
        default="pairs",
        choices=["pairs", "wc_pairs", "wobble_pairs", "nc_pairs", "lone_pairs", "multiplets"],
        help="Which pair column to convert.",
    )
    parser.add_argument(
        "--output-dir",
        default="./prepared_custom_ss/spot_prob",
        help="Directory to write per-target spot_prob matrix files.",
    )
    parser.add_argument(
        "--output-suffix",
        default=".ssprob.txt",
        help="Output filename suffix. Final name is <PDB_ID><suffix>.",
    )
    parser.add_argument(
        "--ids-file",
        default=None,
        help="Optional file with one PDB ID per line. Limits export to these IDs.",
    )
    parser.add_argument(
        "--index-base",
        type=int,
        default=0,
        choices=[0, 1],
        help="Pair index base in CSV.",
    )
    parser.add_argument(
        "--pair-value",
        type=float,
        default=1.0,
        help="Value used for each paired position in the upper triangle.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args()


def read_fasta_sequence(fasta_path: Path) -> str:
    seq_lines: list[str] = []
    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
    return "".join(seq_lines).upper().replace("T", "U")


def load_id_filter(ids_file: str | None) -> set[str] | None:
    if not ids_file:
        return None
    ids: set[str] = set()
    with open(ids_file, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip().upper()
            if line:
                ids.add(line)
    return ids


def find_fasta_for_id(fasta_dir: Path, pdb_id: str) -> Path | None:
    for ext in (".fasta", ".fa", ".fna"):
        candidate = fasta_dir / f"{pdb_id}{ext}"
        if candidate.is_file():
            return candidate
    return None


def parse_pair_list(pair_text: str) -> list[tuple[int, int]]:
    pair_text = (pair_text or "").strip()
    if not pair_text:
        return []
    parsed = ast.literal_eval(pair_text)
    if not isinstance(parsed, list):
        raise ValueError(f"Pairs field is not a list: {pair_text[:80]}")

    pairs: list[tuple[int, int]] = []
    for item in parsed:
        if not (isinstance(item, (tuple, list)) and len(item) == 2):
            raise ValueError(f"Invalid pair item: {item}")
        i, j = item
        if not isinstance(i, int) or not isinstance(j, int):
            raise ValueError(f"Pair indices must be integers: {item}")
        pairs.append((i, j))
    return pairs


def build_upper_tri_matrix(length: int, pairs: list[tuple[int, int]], index_base: int, value: float) -> np.ndarray:
    mat = np.zeros((length, length), dtype=np.float32)

    for i_raw, j_raw in pairs:
        i = i_raw - index_base
        j = j_raw - index_base
        if i < 0 or j < 0:
            raise ValueError(f"Negative index after base conversion: ({i_raw}, {j_raw})")
        if i >= length or j >= length:
            raise ValueError(f"Pair index out of range for length {length}: ({i_raw}, {j_raw})")
        if i == j:
            continue
        if i > j:
            i, j = j, i
        mat[i, j] = max(mat[i, j], value)

    return mat


def main() -> None:
    args = parse_args()

    pairs_csv = Path(args.pairs_csv).resolve()
    fasta_dir = Path(args.fasta_dir).resolve()
    out_dir = Path(args.output_dir).resolve()

    if not pairs_csv.is_file():
        raise FileNotFoundError(f"Pairs CSV not found: {pairs_csv}")
    if not fasta_dir.is_dir():
        raise FileNotFoundError(f"FASTA dir not found: {fasta_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)

    id_filter = load_id_filter(args.ids_file)

    method_rows: dict[str, dict[str, str]] = {}
    duplicates: defaultdict[str, int] = defaultdict(int)

    with open(pairs_csv, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        required_cols = {"pdb_id", "method", args.pairs_column}
        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV missing required columns: {sorted(missing)}")

        for row in reader:
            method = (row.get("method") or "").strip()
            if method != args.method:
                continue
            pdb_id = (row.get("pdb_id") or "").strip().upper()
            if not pdb_id:
                continue
            if id_filter is not None and pdb_id not in id_filter:
                continue
            if pdb_id in method_rows:
                duplicates[pdb_id] += 1
                continue
            method_rows[pdb_id] = row

    if not method_rows:
        raise ValueError(f"No rows found for method='{args.method}' with current filters.")

    exported = 0
    skipped = 0

    manifest_path = out_dir / f"manifest_{args.method}_{args.pairs_column}.tsv"
    manifest_rows: list[dict[str, str]] = []

    for pdb_id in sorted(method_rows):
        row = method_rows[pdb_id]
        fasta_path = find_fasta_for_id(fasta_dir, pdb_id)
        if fasta_path is None:
            manifest_rows.append(
                {
                    "pdb_id": pdb_id,
                    "status": "skipped",
                    "reason": "missing_fasta",
                    "length": "0",
                    "n_pairs": "0",
                    "output_file": "",
                }
            )
            skipped += 1
            continue

        sequence = read_fasta_sequence(fasta_path)
        length = len(sequence)

        try:
            pairs = parse_pair_list(row.get(args.pairs_column, ""))
            mat = build_upper_tri_matrix(length, pairs, args.index_base, args.pair_value)
        except Exception as exc:  # noqa: BLE001
            manifest_rows.append(
                {
                    "pdb_id": pdb_id,
                    "status": "skipped",
                    "reason": f"parse_error:{exc}",
                    "length": str(length),
                    "n_pairs": "0",
                    "output_file": "",
                }
            )
            skipped += 1
            continue

        out_file = out_dir / f"{pdb_id}{args.output_suffix}"
        if out_file.exists() and not args.overwrite:
            manifest_rows.append(
                {
                    "pdb_id": pdb_id,
                    "status": "skipped",
                    "reason": "exists_use_overwrite",
                    "length": str(length),
                    "n_pairs": str(len(pairs)),
                    "output_file": str(out_file),
                }
            )
            skipped += 1
            continue

        np.savetxt(out_file, mat, fmt="%.3f")

        manifest_rows.append(
            {
                "pdb_id": pdb_id,
                "status": "exported",
                "reason": "ok",
                "length": str(length),
                "n_pairs": str(len(pairs)),
                "output_file": str(out_file),
            }
        )
        exported += 1

    with open(manifest_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["pdb_id", "status", "reason", "length", "n_pairs", "output_file"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Pairs CSV     : {pairs_csv}")
    print(f"Method        : {args.method}")
    print(f"Pairs column  : {args.pairs_column}")
    print(f"FASTA dir     : {fasta_dir}")
    print(f"Output dir    : {out_dir}")
    print(f"Exported      : {exported}")
    print(f"Skipped       : {skipped}")
    print(f"Manifest      : {manifest_path}")
    if duplicates:
        print(f"Duplicates ignored for IDs: {sorted(duplicates)}")


if __name__ == "__main__":
    main()
