#!/usr/bin/env python3
"""Convert TRRosetta/RhoFold prediction PDBs into evaluator-compatible CIF files.

The SyntheticEvolution evaluator expects CIF files and later renames them to the
AlphaFold-style format fold_<PDB>_s1_model_0.cif. This script only needs output
filenames starting with the 4-letter PDB code, and injects minimal _entity_poly
metadata required by parseAFcif2DB.py.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

try:
    from Bio.PDB import MMCIFIO, PDBParser
except ImportError as exc:  # pragma: no cover - runtime dependency guard
    raise SystemExit(
        "Missing dependency: Bio.PDB (Biopython). Install with: pip install biopython"
    ) from exc


PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL",
}

RNA_RESIDUES = {
    "A", "C", "G", "U", "I",
    "RA", "RC", "RG", "RU", "RI",
    "ADE", "CYT", "GUA", "URA",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert predicted PDB outputs to CIF for SyntheticEvolution evaluation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--algorithm", required=True, choices=["trrosetta", "rhofold"])
    parser.add_argument("--input-pred-dir", required=True, help="Prediction root directory.")
    parser.add_argument("--output-cif-dir", required=True, help="Output directory for CIF files.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def clean_pdb_id(name: str) -> str | None:
    candidate = name[:4].upper()
    if re.fullmatch(r"[A-Z0-9]{4}", candidate):
        return candidate
    return None


def select_source_pdb(algorithm: str, target_dir: Path) -> Path | None:
    if algorithm == "trrosetta":
        preferred = sorted((target_dir / "structures").glob("*_model1.pdb"))
        if preferred:
            return preferred[0]
        fallback = sorted(target_dir.glob("**/*.pdb"))
        return fallback[0] if fallback else None

    # rhofold
    relaxed = sorted(target_dir.glob("relaxed_*_model.pdb"))
    if relaxed:
        return relaxed[-1]
    unrelaxed = target_dir / "unrelaxed_model.pdb"
    if unrelaxed.is_file():
        return unrelaxed
    fallback = sorted(target_dir.glob("**/*.pdb"))
    return fallback[0] if fallback else None


def classify_chain_type(chain) -> str:
    prot = 0
    rna = 0

    for residue in chain:
        if residue.id[0] != " ":
            continue
        resn = residue.get_resname().strip().upper()
        if resn in PROTEIN_RESIDUES:
            prot += 1
        if resn in RNA_RESIDUES:
            rna += 1

    if prot > rna:
        return "polypeptide(L)"
    return "polyribonucleotide"


def format_cif_token(token: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_.+-]+", token):
        return token
    return f"'{token}'"


def append_entity_poly_if_missing(cif_path: Path, chain_type_rows: list[tuple[str, str]]) -> None:
    content = cif_path.read_text(encoding="utf-8")
    if "_entity_poly.type" in content and "_entity_poly.pdbx_strand_id" in content:
        return

    lines = ["", "loop_", "_entity_poly.type", "_entity_poly.pdbx_strand_id"]
    for entity_type, chain_id in chain_type_rows:
        lines.append(f"{format_cif_token(entity_type)} {format_cif_token(chain_id)}")
    lines.append("")

    cif_path.write_text(content + "\n".join(lines), encoding="utf-8")


def convert_one(source_pdb: Path, out_cif: Path, pdb_id: str) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(source_pdb))

    mmcif_io = MMCIFIO()
    mmcif_io.set_structure(structure)
    mmcif_io.save(str(out_cif))

    chain_rows: list[tuple[str, str]] = []
    model = next(structure.get_models())
    for idx, chain in enumerate(model.get_chains(), start=1):
        chain_id = (chain.id or "").strip() or f"X{idx}"
        chain_rows.append((classify_chain_type(chain), chain_id))

    if not chain_rows:
        chain_rows.append(("polyribonucleotide", "A"))

    append_entity_poly_if_missing(out_cif, chain_rows)


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_pred_dir).resolve()
    output_dir = Path(args.output_cif_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input prediction dir not found: {input_dir}")

    converted = 0
    skipped = 0

    for target_dir in sorted(p for p in input_dir.iterdir() if p.is_dir()):
        pdb_id = clean_pdb_id(target_dir.name)
        if pdb_id is None:
            skipped += 1
            continue

        source_pdb = select_source_pdb(args.algorithm, target_dir)
        if source_pdb is None:
            print(f"[skip] {target_dir.name}: no source PDB found")
            skipped += 1
            continue

        out_cif = output_dir / f"{pdb_id}_{args.algorithm}.cif"
        if out_cif.exists() and not args.overwrite:
            print(f"[skip] {pdb_id}: output exists ({out_cif.name})")
            skipped += 1
            continue

        try:
            convert_one(source_pdb, out_cif, pdb_id)
            converted += 1
            print(f"[ok]   {pdb_id}: {source_pdb} -> {out_cif}")
        except Exception as exc:  # noqa: BLE001
            skipped += 1
            print(f"[fail] {pdb_id}: {exc}")

    print(f"Converted: {converted}")
    print(f"Skipped:   {skipped}")
    print(f"Output:    {output_dir}")


if __name__ == "__main__":
    main()
