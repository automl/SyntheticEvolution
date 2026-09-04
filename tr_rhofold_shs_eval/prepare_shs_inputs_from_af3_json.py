#!/usr/bin/env python3
"""Prepare TRRosetta/RhoFold inputs from AlphaFold3 SHS JSON files.

This script extracts, for each AF3 JSON data file:
- RNA sequence -> FASTA
- synthetic unpaired MSA -> A3M

By default, only single-RNA examples are exported because TRRosettaRNA and
RhoFold are RNA-only predictors.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Tuple

# Repo root: this folder sits inside the SyntheticEvolution checkout.
REPO_ROOT = Path(__file__).resolve().parent.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract SHS FASTA/A3M inputs from AF3 JSON files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--af3-json-dir",
        default=str(REPO_ROOT / "data" / "datafiles" / "alphafold3"),
        help="Directory containing *_data.json files.",
    )
    parser.add_argument(
        "--out-dir",
        default="./prepared_shs_inputs",
        help="Output directory (fasta/, a3m/, metadata/).",
    )
    parser.add_argument(
        "--include-non-rna-only",
        action="store_true",
        help="Include entries that are not single-RNA-only.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing FASTA/A3M files.",
    )
    return parser.parse_args()


def normalize_msa(msa_text: str, seq_id: str, sequence: str) -> str:
    msa_text = (msa_text or "").strip()
    if msa_text:
        if not msa_text.endswith("\n"):
            msa_text += "\n"
        return msa_text

    # Fallback if no SHS MSA is embedded in the JSON.
    return f">{seq_id}\n{sequence}\n"


def classify_json(data: Dict) -> Tuple[bool, str]:
    sequences: List[Dict] = data.get("sequences", [])
    n_rna = sum(1 for item in sequences if "rna" in item)
    n_protein = sum(1 for item in sequences if "protein" in item)
    n_dna = sum(1 for item in sequences if "dna" in item)

    is_single_rna_only = (n_rna == 1 and n_protein == 0 and n_dna == 0)
    reason = f"rna={n_rna},protein={n_protein},dna={n_dna}"
    return is_single_rna_only, reason


def extract_rna_payload(data: Dict) -> Dict | None:
    for item in data.get("sequences", []):
        if "rna" in item:
            return item["rna"]
    return None


def main() -> None:
    args = parse_args()

    af3_json_dir = Path(args.af3_json_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    fasta_dir = out_dir / "fasta"
    a3m_dir = out_dir / "a3m"
    meta_dir = out_dir / "metadata"

    fasta_dir.mkdir(parents=True, exist_ok=True)
    a3m_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    json_files = sorted(af3_json_dir.glob("*.json"))
    if not json_files:
        raise FileNotFoundError(f"No *.json files found in: {af3_json_dir}")

    manifest_path = meta_dir / "manifest.tsv"
    ids_path = meta_dir / "exported_ids.txt"

    exported = 0
    skipped = 0
    rows: List[Dict[str, str]] = []
    exported_ids: List[str] = []

    for json_path in json_files:
        pdb_id = json_path.stem[:4].upper()
        status = "skipped"
        reason = ""
        seq_len = "0"
        msa_lines = "0"

        try:
            with open(json_path, "r", encoding="utf-8") as fh:
                data = json.load(fh)

            is_single_rna_only, class_reason = classify_json(data)
            if not is_single_rna_only and not args.include_non_rna_only:
                reason = f"filtered_non_single_rna ({class_reason})"
                skipped += 1
            else:
                rna = extract_rna_payload(data)
                if rna is None:
                    reason = "no_rna_payload"
                    skipped += 1
                else:
                    sequence = (rna.get("sequence") or "").strip().upper().replace("T", "U")
                    msa_text = normalize_msa(rna.get("unpairedMsa", ""), pdb_id, sequence)

                    if not sequence:
                        reason = "empty_rna_sequence"
                        skipped += 1
                    else:
                        fasta_path = fasta_dir / f"{pdb_id}.fasta"
                        a3m_path = a3m_dir / f"{pdb_id}.a3m"

                        if (fasta_path.exists() or a3m_path.exists()) and not args.overwrite:
                            reason = "exists_use_overwrite"
                            skipped += 1
                        else:
                            with open(fasta_path, "w", encoding="utf-8") as f_fasta:
                                f_fasta.write(f">{pdb_id}\n{sequence}\n")

                            with open(a3m_path, "w", encoding="utf-8") as f_a3m:
                                f_a3m.write(msa_text)

                            status = "exported"
                            reason = class_reason
                            seq_len = str(len(sequence))
                            msa_lines = str(len([ln for ln in msa_text.splitlines() if ln.startswith(">")]))
                            exported += 1
                            exported_ids.append(pdb_id)

        except Exception as exc:  # noqa: BLE001
            skipped += 1
            reason = f"error:{exc}"

        rows.append(
            {
                "pdb_id": pdb_id,
                "status": status,
                "reason": reason,
                "seq_len": seq_len,
                "msa_nseq": msa_lines,
                "json_file": str(json_path),
            }
        )

    with open(manifest_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["pdb_id", "status", "reason", "seq_len", "msa_nseq", "json_file"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    with open(ids_path, "w", encoding="utf-8") as fh:
        for pdb_id in exported_ids:
            fh.write(f"{pdb_id}\n")

    print(f"Input JSON dir : {af3_json_dir}")
    print(f"Output dir     : {out_dir}")
    print(f"Exported       : {exported}")
    print(f"Skipped        : {skipped}")
    print(f"Manifest       : {manifest_path}")
    print(f"IDs list       : {ids_path}")


if __name__ == "__main__":
    main()
