#!/usr/bin/env python3
"""Run trRosettaRNA on FASTA inputs with either provided or single-seq MSA."""

from __future__ import annotations

import argparse
import logging
import re
import subprocess
import sys
import textwrap
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch trRosettaRNA predictions from prepared inputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fasta-dir", default="./prepared_shs_inputs/fasta")
    parser.add_argument("--a3m-dir", default="./prepared_shs_inputs/a3m")
    parser.add_argument(
        "--msa-mode",
        default="provided",
        choices=["provided", "single_seq"],
        help="MSA source: use provided A3M files or single-sequence fallback.",
    )
    parser.add_argument("--output-dir", default="./predictions_trrosetta")
    parser.add_argument("--trrosetta-root", default="../../trRosettaRNA_v1.1")
    parser.add_argument("--gpu", type=int, default=-1, help="GPU index. -1 means CPU-only.")
    parser.add_argument("--n-cpu", type=int, default=8)
    parser.add_argument("--n-decoys", type=int, default=5)
    parser.add_argument("--skip-spotrna", action="store_true")
    parser.add_argument(
        "--ss-dir",
        default=None,
        help="Optional directory containing one SS file per target: <ID><ss-suffix>.",
    )
    parser.add_argument(
        "--ss-fmt",
        default="spot_prob",
        choices=["spot_prob", "dot_bracket", "ct"],
        help="Format passed to trRosetta when --ss-dir is provided.",
    )
    parser.add_argument(
        "--ss-suffix",
        default=".ssprob.txt",
        help="Suffix for SS files in --ss-dir. Final lookup is <ID><ss-suffix>.",
    )
    parser.add_argument(
        "--require-ss",
        action="store_true",
        help="Fail targets with missing SS files when --ss-dir is used.",
    )
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--ids-file", default=None, help="Optional file with one PDB ID per line.")
    parser.add_argument("--log-file", default=None)
    return parser.parse_args()


def setup_logging(log_file: str | None) -> logging.Logger:
    logger = logging.getLogger("trrosetta_shs")
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


def parse_fasta_sequence(fasta_path: Path) -> str:
    seq_lines = []
    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
    seq = "".join(seq_lines).upper().replace("T", "U")
    bad = set(seq) - set("ACGUN-")
    if bad:
        raise ValueError(f"Invalid characters in sequence for {fasta_path.name}: {sorted(bad)}")
    return seq


def stream_cmd(cmd: list[str], cwd: Path, logger: logging.Logger) -> int:
    logger.info("Running: %s", " ".join(cmd))
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        logger.info(line.rstrip())
    return proc.wait()


def make_dot_bracket(seq: str) -> str:
    return "." * len(seq)


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


def safe_id(name: str) -> str:
    return re.sub(r"[^\w\-]", "_", name)


def write_single_seq_a3m(path: Path, seq_id: str, sequence: str) -> None:
    path.write_text(f">{seq_id}\n{sequence}\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    logger = setup_logging(args.log_file)

    fasta_dir = Path(args.fasta_dir).resolve()
    a3m_dir = Path(args.a3m_dir).resolve() if args.a3m_dir else None
    ss_dir = Path(args.ss_dir).resolve() if args.ss_dir else None
    out_root = Path(args.output_dir).resolve()
    trrosetta_root = Path(args.trrosetta_root).resolve()

    predict_py = trrosetta_root / "predict.py"
    fold_py = trrosetta_root / "fold.py"
    params_dir = trrosetta_root / "params" / "model_1"

    if not fasta_dir.is_dir():
        raise FileNotFoundError(f"FASTA dir not found: {fasta_dir}")
    if args.msa_mode == "provided":
        if a3m_dir is None or not a3m_dir.is_dir():
            raise FileNotFoundError(f"A3M dir not found: {a3m_dir}")
    if ss_dir is not None and not ss_dir.is_dir():
        raise FileNotFoundError(f"SS dir not found: {ss_dir}")
    if not predict_py.is_file() or not fold_py.is_file() or not params_dir.is_dir():
        raise FileNotFoundError("trRosettaRNA installation layout not found. Check --trrosetta-root.")

    out_root.mkdir(parents=True, exist_ok=True)

    id_filter = read_id_filter(args.ids_file)

    fasta_files = sorted(list(fasta_dir.glob("*.fa")) + list(fasta_dir.glob("*.fasta")))
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in: {fasta_dir}")

    ok = 0
    fail = 0
    skip = 0

    logger.info("Found %d FASTA files.", len(fasta_files))

    for fasta_path in fasta_files:
        seq_id = safe_id(fasta_path.stem.upper())
        if id_filter is not None and seq_id not in id_filter:
            skip += 1
            continue

        seq_out = out_root / seq_id
        struct_dir = seq_out / "structures"
        geom_dir = seq_out / "geometry"
        seq_out.mkdir(parents=True, exist_ok=True)
        struct_dir.mkdir(parents=True, exist_ok=True)
        geom_dir.mkdir(parents=True, exist_ok=True)

        out_pdb = struct_dir / f"{seq_id}_model1.pdb"
        npz_path = geom_dir / f"{seq_id}.npz"

        if args.skip_existing and out_pdb.is_file():
            logger.info("[%s] Existing PDB found, skipping.", seq_id)
            skip += 1
            continue

        seq = parse_fasta_sequence(fasta_path)

        # Keep a copy of used inputs under each prediction folder for traceability.
        local_fasta = seq_out / f"{seq_id}.fasta"
        local_a3m = seq_out / f"{seq_id}.a3m"
        local_fasta.write_text(f">{seq_id}\n" + "\n".join(textwrap.wrap(seq, 60)) + "\n", encoding="utf-8")

        if args.msa_mode == "provided":
            assert a3m_dir is not None
            a3m_path = a3m_dir / f"{seq_id}.a3m"
            if not a3m_path.is_file():
                logger.warning("[%s] Missing A3M file: %s", seq_id, a3m_path)
                fail += 1
                continue
            local_a3m.write_text(a3m_path.read_text(encoding="utf-8"), encoding="utf-8")
        else:
            write_single_seq_a3m(local_a3m, seq_id, seq)

        ss_args: list[str] = []
        if ss_dir is not None:
            ss_path = ss_dir / f"{seq_id}{args.ss_suffix}"
            if ss_path.is_file():
                ss_args = ["-ss", str(ss_path), "-ss_fmt", args.ss_fmt]
            elif args.require_ss:
                logger.warning("[%s] Missing SS file: %s", seq_id, ss_path)
                fail += 1
                continue
            else:
                logger.info("[%s] No SS file found at %s, falling back to SPOT-RNA.", seq_id, ss_path)
        elif args.skip_spotrna:
            dbn_path = seq_out / f"{seq_id}.dbn"
            dbn_path.write_text(make_dot_bracket(seq) + "\n", encoding="utf-8")
            ss_args = ["-ss", str(dbn_path), "-ss_fmt", "dot_bracket"]

        geom_cmd = [
            sys.executable,
            str(predict_py),
            "-i",
            str(local_a3m),
            "-o",
            str(npz_path),
            "-mdir",
            str(params_dir),
            "-cpu",
            str(args.n_cpu),
        ] + ss_args
        if args.gpu >= 0:
            geom_cmd.extend(["-gpu", str(args.gpu)])

        rc = stream_cmd(geom_cmd, trrosetta_root, logger)
        if rc != 0 or not npz_path.is_file():
            logger.error("[%s] Geometry prediction failed.", seq_id)
            fail += 1
            continue

        fold_cmd = [
            sys.executable,
            str(fold_py),
            "-npz",
            str(npz_path),
            "-fa",
            str(local_fasta),
            "-out",
            str(out_pdb),
            "-nm",
            str(args.n_decoys),
            "-cpu",
            str(args.n_cpu),
        ]

        rc = stream_cmd(fold_cmd, trrosetta_root, logger)
        if rc != 0 or not out_pdb.is_file():
            logger.error("[%s] Folding failed.", seq_id)
            fail += 1
            continue

        logger.info("[%s] Done.", seq_id)
        ok += 1

    logger.info("Summary: success=%d failed=%d skipped=%d", ok, fail, skip)
    if fail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
