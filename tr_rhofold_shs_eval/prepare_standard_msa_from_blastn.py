#!/usr/bin/env python3
"""Generate standard-database A3M files from FASTA using RhoFold BLASTN helper."""

from __future__ import annotations

import argparse
import logging
import string
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate standard A3M alignments with RhoFold BLASTN pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fasta-dir", required=True, help="Directory with input FASTA files.")
    parser.add_argument("--output-a3m-dir", required=True, help="Output directory for generated .a3m files.")
    parser.add_argument("--rhofold-root", default="../../RhoFold", help="Path to RhoFold repository root.")
    parser.add_argument(
        "--database-dpath",
        default=None,
        help="Database dir containing rnacentral.fasta and nt. Defaults to <rhofold-root>/database.",
    )
    parser.add_argument(
        "--binary-dpath",
        default=None,
        help="Binary dir containing blastn and helper perl scripts. Defaults to <rhofold-root>/rhofold/data/bin.",
    )
    parser.add_argument("--n-cpu", type=int, default=4, help="Threads passed to BLASTN.")
    parser.add_argument("--ids-file", default=None, help="Optional text file with one ID per line.")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument(
        "--sanitize-a3m-for-trrosetta",
        action="store_true",
        default=True,
        help="Filter generated A3M to rows compatible with trRosettaRNA parser.",
    )
    parser.add_argument("--log-file", default=None)
    return parser.parse_args()


def setup_logging(log_file: str | None) -> logging.Logger:
    logger = logging.getLogger("prepare_standard_msa")
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.handlers = [stream_handler]

    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def read_id_filter(ids_file: str | None) -> set[str] | None:
    if not ids_file:
        return None
    ids = set()
    with open(ids_file, "r", encoding="utf-8") as fh:
        for line in fh:
            token = line.strip().upper()
            if token:
                ids.add(token)
    return ids


def read_query_length(fasta_path: Path) -> int:
    length = 0
    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            length += len(line)
    return length


def _trrosetta_clean_seq(seq: str) -> str:
    table = str.maketrans(dict.fromkeys(string.ascii_lowercase))
    return (
        seq.rstrip()
        .replace("W", "A")
        .replace("R", "A")
        .replace("Y", "C")
        .replace("E", "A")
        .replace("I", "A")
        .replace("P", "G")
        .replace("T", "U")
        .translate(table)
    )


def sanitize_a3m_for_trrosetta(a3m_path: Path, query_length: int, logger: logging.Logger) -> None:
    rows: list[tuple[str, str]] = []
    header: str | None = None
    seq_chunks: list[str] = []

    with open(a3m_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    rows.append((header, "".join(seq_chunks)))
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if header is not None:
        rows.append((header, "".join(seq_chunks)))

    if not rows:
        raise ValueError(f"Empty A3M file: {a3m_path}")

    kept: list[tuple[str, str]] = []
    for h, s in rows:
        cleaned = _trrosetta_clean_seq(s)
        if len(cleaned) == query_length:
            kept.append((h, s))

    if not kept:
        # Keep query row as a last resort to avoid hard failure.
        kept = [rows[0]]
        logger.warning(
            "A3M sanitize removed all rows for %s; keeping only query row.",
            a3m_path.name,
        )

    with open(a3m_path, "w", encoding="utf-8") as fh:
        for h, s in kept:
            fh.write(f"{h}\n{s}\n")

    logger.info(
        "Sanitized %s for trRosetta: kept %d/%d rows (query_len=%d)",
        a3m_path.name,
        len(kept),
        len(rows),
        query_length,
    )


def main() -> None:
    args = parse_args()
    logger = setup_logging(args.log_file)

    fasta_dir = Path(args.fasta_dir).resolve()
    out_dir = Path(args.output_a3m_dir).resolve()
    rhofold_root = Path(args.rhofold_root).resolve()

    if not fasta_dir.is_dir():
        raise FileNotFoundError(f"FASTA dir not found: {fasta_dir}")
    if not rhofold_root.is_dir():
        raise FileNotFoundError(f"RhoFold root not found: {rhofold_root}")

    sys.path.insert(0, str(rhofold_root))
    from rhofold.data.balstn import BLASTN  # pylint: disable=import-outside-toplevel

    database_dpath = Path(args.database_dpath).resolve() if args.database_dpath else (rhofold_root / "database")
    binary_dpath = Path(args.binary_dpath).resolve() if args.binary_dpath else (rhofold_root / "rhofold" / "data" / "bin")

    db_rnacentral = database_dpath / "rnacentral.fasta"
    db_nt = database_dpath / "nt"
    blastn_bin = binary_dpath / "blastn"
    parse_perl = binary_dpath / "parse_blastn_local.pl"
    reformat_perl = binary_dpath / "reformat.pl"

    missing: list[Path] = []
    for req in (blastn_bin, parse_perl, reformat_perl):
        if not req.exists():
            missing.append(req)
    if missing:
        msg = "\n".join(str(p) for p in missing)
        raise FileNotFoundError(f"Missing required BLASTN/MSA resources:\n{msg}")

    databases: list[str] = []
    if db_rnacentral.exists():
        databases.append(str(db_rnacentral))
    if db_nt.exists():
        databases.append(str(db_nt))
    if not databases:
        raise FileNotFoundError(
            "No usable BLAST database found. Expected at least one of:\n"
            f"{db_rnacentral}\n{db_nt}"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    id_filter = read_id_filter(args.ids_file)

    fasta_files = sorted(
        list(fasta_dir.glob("*.fa"))
        + list(fasta_dir.glob("*.fasta"))
        + list(fasta_dir.glob("*.fna"))
    )
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in: {fasta_dir}")

    blast = BLASTN(
        binary_dpath=str(binary_dpath),
        databases=databases,
        n_cpu=args.n_cpu,
    )

    ok = 0
    fail = 0
    skip = 0

    logger.info("Found %d FASTA files.", len(fasta_files))
    logger.info("Using database_dpath: %s", database_dpath)
    logger.info("Using binary_dpath: %s", binary_dpath)
    logger.info("Using BLAST databases: %s", ", ".join(databases))

    for fasta_path in fasta_files:
        seq_id = fasta_path.stem.upper()
        if id_filter is not None and seq_id not in id_filter:
            skip += 1
            continue

        out_a3m = out_dir / f"{seq_id}.a3m"
        if args.skip_existing and out_a3m.is_file() and out_a3m.stat().st_size > 0:
            logger.info("[%s] Existing A3M found, skipping.", seq_id)
            skip += 1
            continue

        try:
            blast.query(str(fasta_path), str(out_a3m), logger)
            if out_a3m.is_file() and out_a3m.stat().st_size > 0:
                if args.sanitize_a3m_for_trrosetta:
                    query_len = read_query_length(fasta_path)
                    sanitize_a3m_for_trrosetta(out_a3m, query_len, logger)
                logger.info("[%s] A3M ready: %s", seq_id, out_a3m)
                ok += 1
            else:
                logger.error("[%s] BLAST finished but A3M missing/empty: %s", seq_id, out_a3m)
                fail += 1
        except Exception as exc:  # noqa: BLE001
            logger.error("[%s] Failed generating standard MSA: %s", seq_id, exc)
            fail += 1

    logger.info("Summary: success=%d failed=%d skipped=%d", ok, fail, skip)
    if fail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
