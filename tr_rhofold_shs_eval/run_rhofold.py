#!/usr/bin/env python3
"""
RhoFold RNA structure prediction pipeline.
Processes a directory of FASTA files and runs RhoFold on each sequence.

Usage:
    python run_rhofold.py --input_dir <fasta_dir> --output_dir <out_dir> [options]
"""

import argparse
import logging
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

# ─── Logging ────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)


# ─── Helpers ────────────────────────────────────────────────────────────────

def parse_fasta(fasta_path: Path) -> list[tuple[str, str]]:
    """Return list of (header, sequence) tuples from a FASTA file."""
    records = []
    header, seq_parts = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:].split()[0]  # use first token as ID
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if header is not None:
            records.append((header, "".join(seq_parts)))
    return records


def write_single_fasta(header: str, sequence: str, path: Path) -> None:
    with open(path, "w") as fh:
        fh.write(f">{header}\n{sequence}\n")


def run_rhofold(
    input_fasta: Path,
    output_dir: Path,
    rhofold_bin: str,   # now should be "inference" script path
    checkpoint: str | None,
    device: str,
    extra_args: list[str],
) -> bool:
    """
    Call RhoFold on a single FASTA file.
    Returns True on success, False on failure.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build command in a flexible way:
    # - script path:     /path/to/inference.py      -> python /path/to/inference.py ...
    # - CLI executable:  rhofold                    -> rhofold ...
    # - full command:    python -m rhofold.infer    -> python -m rhofold.infer ...
    rhofold_tokens = shlex.split(rhofold_bin)
    if not rhofold_tokens:
        raise ValueError("--rhofold_bin resolved to an empty command")

    first = rhofold_tokens[0]
    looks_like_script = (
        first.endswith(".py")
        or os.sep in first
        or Path(first).is_file()
    )

    cmd_cwd = None
    if len(rhofold_tokens) == 1 and looks_like_script:
        log.warning("Interpreting --rhofold_bin as script path: %s", first)
        script_path = str(Path(first).expanduser().resolve())
        cmd = [sys.executable, script_path]
        cmd_cwd = str(Path(script_path).parent)
    else:
        log.warning("Interpreting --rhofold_bin as command: %s", rhofold_bin)
        cmd = rhofold_tokens
        if len(rhofold_tokens) >= 2:
            first_lower = Path(rhofold_tokens[0]).name.lower()
            second = rhofold_tokens[1]
            if first_lower.startswith("python") and (second.endswith(".py") or Path(second).is_file()):
                cmd_cwd = str(Path(second).expanduser().resolve().parent)

    cmd += ["--input_fas", str(input_fasta),
            "--output_dir", str(output_dir)]

    if checkpoint:
        cmd += ["--ckpt", str(Path(checkpoint).expanduser().resolve())]

    cmd += ["--device", device]

    # add any extra arguments
    cmd += extra_args

    log.info("Running: %s", " ".join(cmd))
    if cmd_cwd:
        log.info("Working directory: %s", cmd_cwd)
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cmd_cwd)
    elapsed = time.time() - t0

    if result.returncode != 0:
        log.error("RhoFold FAILED for %s (%.1fs)\nSTDERR:\n%s", input_fasta.name, elapsed, result.stderr)
        return False

    log.info("Finished %s in %.1fs", input_fasta.name, elapsed)
    if result.stdout:
        log.debug("STDOUT:\n%s", result.stdout)
    return True


# ─── Main ───────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="RhoFold batch prediction pipeline")

    parser.add_argument("--input_dir", required=True,
                        help="Directory containing input FASTA files (*.fa, *.fasta, *.fna)")
    parser.add_argument("--output_dir", required=True,
                        help="Root output directory; one sub-folder is created per sequence")
    parser.add_argument("--rhofold_bin", default="rhofold",
                        help="RhoFold command or script path (e.g. 'rhofold' or '/path/to/inference.py')")
    parser.add_argument("--checkpoint", default=None,
                        help="Path to RhoFold model checkpoint (.pt)")
    parser.add_argument("--device", default="cuda",
                        choices=["cuda", "cpu"],
                        help="Compute device (default: cuda)")
    parser.add_argument("--split_multi", action="store_true",
                        help="If a FASTA has multiple records, run each as a separate job")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip sequences whose output directory already contains a .pdb file")
    parser.add_argument("--log_file", default=None,
                        help="Optional path to write a log file")
    parser.add_argument("--rhofold_args", nargs=argparse.REMAINDER, default=[],
                        help="Extra arguments passed verbatim to RhoFold (append after '--')")

    args = parser.parse_args()

    # ── File handler for log ──────────────────────────────────────────────
    if args.log_file:
        fh = logging.FileHandler(args.log_file)
        fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        log.addHandler(fh)

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.is_dir():
        log.error("--input_dir does not exist: %s", input_dir)
        sys.exit(1)

    # ── Collect FASTA files ───────────────────────────────────────────────
    fasta_files = sorted(
        f for f in input_dir.iterdir()
        if f.suffix.lower() in {".fa", ".fasta", ".fna"}
    )

    if not fasta_files:
        log.error("No FASTA files found in %s", input_dir)
        sys.exit(1)

    log.info("Found %d FASTA file(s) in %s", len(fasta_files), input_dir)

    # ── Process each file ─────────────────────────────────────────────────
    success_count = 0
    fail_count = 0
    skip_count = 0

    tmp_dir = output_dir / "_tmp_fastas"
    tmp_dir.mkdir(exist_ok=True)

    for fasta_path in fasta_files:
        records = parse_fasta(fasta_path)
        if not records:
            log.warning("No records found in %s — skipping", fasta_path.name)
            continue

        # Decide which records to process
        if args.split_multi and len(records) > 1:
            work_items = records  # each record becomes its own job
        else:
            # Treat whole file as one unit; use first header for naming
            work_items = [(records[0][0], None)]  # sentinel: None = use original file

        for header, sequence in work_items:
            safe_name = header.replace("/", "_").replace(" ", "_")
            seq_out_dir = output_dir / safe_name

            # ── Skip if output already exists ──────────────────────────
            if args.skip_existing:
                existing_pdbs = list(seq_out_dir.glob("*.pdb"))
                if existing_pdbs:
                    log.info("Skipping %s (output already exists)", safe_name)
                    skip_count += 1
                    continue

            # ── Prepare input FASTA ────────────────────────────────────
            if sequence is None:
                # Use original file directly
                run_fasta = fasta_path
            else:
                run_fasta = tmp_dir / f"{safe_name}.fasta"
                write_single_fasta(header, sequence, run_fasta)

            # ── Run RhoFold ────────────────────────────────────────────
            ok = run_rhofold(
                input_fasta=run_fasta,
                output_dir=seq_out_dir,
                rhofold_bin=args.rhofold_bin,
                checkpoint=args.checkpoint,
                device=args.device,
                extra_args=args.rhofold_args,
            )

            if ok:
                success_count += 1
            else:
                fail_count += 1

    # ── Summary ───────────────────────────────────────────────────────────
    log.info(
        "Pipeline complete — success: %d | failed: %d | skipped: %d",
        success_count, fail_count, skip_count,
    )

    if fail_count > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
