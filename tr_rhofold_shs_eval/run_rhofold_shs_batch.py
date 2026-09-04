#!/usr/bin/env python3
"""Run RhoFold on FASTA inputs with configurable MSA sourcing."""

from __future__ import annotations

import argparse
import logging
import os
import shlex
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch RhoFold predictions from prepared inputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fasta-dir", default="./prepared_shs_inputs/fasta")
    parser.add_argument("--a3m-dir", default="./prepared_shs_inputs/a3m")
    parser.add_argument(
        "--msa-mode",
        default="provided",
        choices=["provided", "auto_db", "single_seq"],
        help="MSA source mode.",
    )
    parser.add_argument("--output-dir", default="./predictions_rhofold")
    parser.add_argument(
        "--rhofold-bin",
        default="../../RhoFold/inference.py",
        help="RhoFold script path or command.",
    )
    parser.add_argument("--checkpoint", default=None)
    parser.add_argument(
        "--database-dpath",
        default=None,
        help="RhoFold DB dir for --msa-mode auto_db (contains rnacentral.fasta and nt).",
    )
    parser.add_argument(
        "--binary-dpath",
        default=None,
        help="RhoFold binary dir for --msa-mode auto_db (contains blastn + helper perl scripts).",
    )
    parser.add_argument("--device", default="cuda", choices=["cuda", "cpu"])
    parser.add_argument("--relax-steps", default="1000")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--ids-file", default=None, help="Optional file with one PDB ID per line.")
    parser.add_argument("--log-file", default=None)
    parser.add_argument(
        "--rhofold-args",
        nargs=argparse.REMAINDER,
        default=[],
        help="Extra args passed verbatim to RhoFold.",
    )
    return parser.parse_args()


def setup_logging(log_file: str | None) -> logging.Logger:
    logger = logging.getLogger("rhofold_shs")
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
            line = line.strip().upper()
            if line:
                ids.add(line)
    return ids


def resolve_command(rhofold_bin: str) -> tuple[list[str], str | None]:
    tokens = shlex.split(rhofold_bin)
    if not tokens:
        raise ValueError("--rhofold-bin resolved to an empty command")

    first = tokens[0]
    looks_like_script = first.endswith(".py") or os.sep in first or Path(first).is_file()

    if len(tokens) == 1 and looks_like_script:
        script = str(Path(first).expanduser().resolve())
        return [sys.executable, script], str(Path(script).parent)

    cwd = None
    if len(tokens) >= 2:
        first_name = Path(tokens[0]).name.lower()
        second = tokens[1]
        if first_name.startswith("python") and (second.endswith(".py") or Path(second).is_file()):
            cwd = str(Path(second).expanduser().resolve().parent)
    return tokens, cwd


def get_primary_pdb(pred_dir: Path) -> Path | None:
    relaxed = sorted(pred_dir.glob("relaxed_*_model.pdb"))
    if relaxed:
        return relaxed[-1]
    unrelaxed = pred_dir / "unrelaxed_model.pdb"
    if unrelaxed.is_file():
        return unrelaxed
    return None


def main() -> None:
    args = parse_args()
    logger = setup_logging(args.log_file)

    fasta_dir = Path(args.fasta_dir).resolve()
    a3m_dir = Path(args.a3m_dir).resolve() if args.a3m_dir else None
    out_root = Path(args.output_dir).resolve()

    if not fasta_dir.is_dir():
        raise FileNotFoundError(f"FASTA dir not found: {fasta_dir}")
    if args.msa_mode == "provided":
        if a3m_dir is None or not a3m_dir.is_dir():
            raise FileNotFoundError(f"A3M dir not found: {a3m_dir}")

    out_root.mkdir(parents=True, exist_ok=True)
    id_filter = read_id_filter(args.ids_file)
    cmd_prefix, cmd_cwd = resolve_command(args.rhofold_bin)

    fasta_files = sorted(list(fasta_dir.glob("*.fa")) + list(fasta_dir.glob("*.fasta")) + list(fasta_dir.glob("*.fna")))
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in: {fasta_dir}")

    ok = 0
    fail = 0
    skip = 0

    logger.info("Found %d FASTA files.", len(fasta_files))

    for fasta_path in fasta_files:
        seq_id = fasta_path.stem.upper()
        if id_filter is not None and seq_id not in id_filter:
            skip += 1
            continue

        seq_out = out_root / seq_id
        seq_out.mkdir(parents=True, exist_ok=True)

        if args.skip_existing and get_primary_pdb(seq_out) is not None:
            logger.info("[%s] Existing prediction found, skipping.", seq_id)
            skip += 1
            continue

        cmd = list(cmd_prefix)
        cmd += [
            "--input_fas",
            str(fasta_path),
            "--output_dir",
            str(seq_out),
            "--device",
            args.device,
            "--relax_steps",
            str(args.relax_steps),
        ]

        if args.msa_mode == "provided":
            assert a3m_dir is not None
            a3m_path = a3m_dir / f"{seq_id}.a3m"
            if not a3m_path.is_file():
                logger.warning("[%s] Missing A3M file: %s", seq_id, a3m_path)
                fail += 1
                continue
            cmd += ["--input_a3m", str(a3m_path)]
        elif args.msa_mode == "single_seq":
            cmd += ["--single_seq_pred", "True"]
        else:
            if args.database_dpath:
                cmd += ["--database_dpath", str(Path(args.database_dpath).expanduser().resolve())]
            if args.binary_dpath:
                cmd += ["--binary_dpath", str(Path(args.binary_dpath).expanduser().resolve())]

        if args.checkpoint:
            cmd += ["--ckpt", str(Path(args.checkpoint).expanduser().resolve())]

        cmd += list(args.rhofold_args)

        logger.info("[%s] Running: %s", seq_id, " ".join(cmd))
        rc = subprocess.run(cmd, cwd=cmd_cwd, check=False).returncode
        primary_pdb = get_primary_pdb(seq_out)
        if primary_pdb is None:
            logger.error("[%s] Prediction failed (no output PDB).", seq_id)
            fail += 1
            continue
        if rc != 0:
            unrelaxed = seq_out / "unrelaxed_model.pdb"
            if primary_pdb == unrelaxed:
                logger.warning("[%s] Amber relaxation failed; using unrelaxed model: %s", seq_id, unrelaxed)
            else:
                logger.error("[%s] Prediction failed (rc=%d).", seq_id, rc)
                fail += 1
                continue

        logger.info("[%s] Done.", seq_id)
        ok += 1

    logger.info("Summary: success=%d failed=%d skipped=%d", ok, fail, skip)
    if fail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
