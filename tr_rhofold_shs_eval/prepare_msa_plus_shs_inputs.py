#!/usr/bin/env python3
"""Build combined "natural MSA + SHS" A3M inputs for RhoFold.

For each target, concatenates:
  1) the natural BLASTN-derived MSA A3M (from prepared_db_msa_inputs/a3m), and
  2) the SHS A3M produced by an upstream method (e.g. rnaformer N100),
into a single A3M, dropping the duplicate query row from the SHS file.
The query sequence and aligned-column count must match between the two sources.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stack natural-DB MSA + SHS A3M into combined inputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--db-msa-dir", required=True,
                        help="Folder with fasta/ and a3m/ from BLASTN MSA prep.")
    parser.add_argument("--shs-dir", required=True,
                        help="Folder with fasta/ and a3m/ from SHS prep.")
    parser.add_argument("--output-dir", required=True,
                        help="Output folder; fasta/ and a3m/ are created underneath.")
    return parser.parse_args()


def read_first_record(a3m_path: Path) -> tuple[str, str]:
    header: str | None = None
    seq_chunks: list[str] = []
    with open(a3m_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    break
                header = line
            else:
                seq_chunks.append(line.strip())
    if header is None:
        raise ValueError(f"No header in {a3m_path}")
    return header, "".join(seq_chunks)


def strip_first_record(a3m_path: Path) -> str:
    """Return the A3M content with its first record (query) removed."""
    out_lines: list[str] = []
    seen_first_header = False
    past_first_record = False
    with open(a3m_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if not seen_first_header:
                    seen_first_header = True
                    continue
                past_first_record = True
                out_lines.append(line)
            else:
                if seen_first_header and not past_first_record:
                    continue
                out_lines.append(line)
    return "\n".join(out_lines) + ("\n" if out_lines else "")


def aligned_len(seq: str) -> int:
    return sum(1 for c in seq if not c.islower())


def main() -> None:
    args = parse_args()

    db_root = Path(args.db_msa_dir).resolve()
    shs_root = Path(args.shs_dir).resolve()
    out_root = Path(args.output_dir).resolve()

    db_a3m_dir = db_root / "a3m"
    shs_a3m_dir = shs_root / "a3m"
    shs_fasta_dir = shs_root / "fasta"

    for p in (db_a3m_dir, shs_a3m_dir, shs_fasta_dir):
        if not p.is_dir():
            raise FileNotFoundError(f"Missing input dir: {p}")

    out_fasta = out_root / "fasta"
    out_a3m = out_root / "a3m"
    out_fasta.mkdir(parents=True, exist_ok=True)
    out_a3m.mkdir(parents=True, exist_ok=True)

    fastas = sorted(list(shs_fasta_dir.glob("*.fa")) + list(shs_fasta_dir.glob("*.fasta")))
    if not fastas:
        raise FileNotFoundError(f"No FASTA in {shs_fasta_dir}")

    n_ok = 0
    n_skip = 0
    n_only_shs = 0
    mismatches: list[str] = []

    for fa in fastas:
        seq_id = fa.stem.upper()
        shs_a3m = shs_a3m_dir / f"{seq_id}.a3m"
        if not shs_a3m.is_file():
            print(f"[{seq_id}] missing SHS a3m, skip", file=sys.stderr)
            n_skip += 1
            continue

        db_a3m = db_a3m_dir / f"{seq_id}.a3m"
        shs_header, shs_seq = read_first_record(shs_a3m)

        merged_lines: list[str] = []

        if db_a3m.is_file():
            db_header, db_seq = read_first_record(db_a3m)
            if db_seq.upper() != shs_seq.upper():
                mismatches.append(seq_id)
                print(f"[{seq_id}] query mismatch between DB and SHS a3m", file=sys.stderr)
                n_skip += 1
                continue
            db_aligned = aligned_len(db_seq)
            shs_aligned = aligned_len(shs_seq)
            if db_aligned != shs_aligned:
                mismatches.append(seq_id)
                print(f"[{seq_id}] aligned-col count mismatch: db={db_aligned} shs={shs_aligned}", file=sys.stderr)
                n_skip += 1
                continue
            with open(db_a3m, "r", encoding="utf-8") as fh:
                db_text = fh.read()
            if not db_text.endswith("\n"):
                db_text += "\n"
            merged_lines.append(db_text)
            tail = strip_first_record(shs_a3m)
            if tail:
                merged_lines.append(tail)
        else:
            with open(shs_a3m, "r", encoding="utf-8") as fh:
                merged_lines.append(fh.read())
            n_only_shs += 1

        out_path = out_a3m / f"{seq_id}.a3m"
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write("".join(merged_lines))

        shutil.copyfile(fa, out_fasta / fa.name)
        n_ok += 1

    print(
        f"Built {n_ok} merged A3Ms (only-SHS fallback for {n_only_shs}); "
        f"skipped {n_skip}; mismatches={len(mismatches)}"
    )
    if mismatches:
        print("Mismatched targets: " + ", ".join(mismatches))


if __name__ == "__main__":
    main()
