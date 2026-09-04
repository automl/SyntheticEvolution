#!/usr/bin/env python3
"""
run_trrosettarna.py
====================
Pipeline wrapper for trRosettaRNA v1.1 (Zenodo 10616895).

For each FASTA in --input_dir, this script:
  1. Validates and reformats the sequence.
  2. Generates a single-sequence A3M file (no external MSA search needed).
  3. Runs geometry prediction:  predict.py -i seq.a3m -o seq.npz [-gpu N|-cpu N]
  4. Runs 3D structure folding: fold.py   -npz seq.npz -fa seq.fasta -out model.pdb
  5. Collects output PDBs into --output_dir/<seq_id>/.

Note on MSA: without real homologs the geometry prediction quality is reduced,
but single-sequence mode still works and is fine for testing / small RNAs.
If you have a real A3M from Infernal/BLASTN, place it at:
  <input_dir>/<seq_id>.a3m   — it will be used automatically instead.

Usage (standalone, outside SLURM):
    python run_trrosettarna.py \
        --input_dir  data/ \
        --output_dir trrosettarna_predictions/ \
        --trrosetta_root /abs/path/to/trRosettaRNA_v1.1 \
        --gpu 0

Write bundled example FASTAs:
    python run_trrosettarna.py --write_examples --input_dir data/
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import textwrap
import time
from pathlib import Path
from typing import Optional

# =============================================================================
#  Example RNA sequences (used when --write_examples is passed)
#  These are small, well-characterised RNAs good for a quick smoke-test.
# =============================================================================

EXAMPLE_SEQUENCES = {
    # Adenine riboswitch aptamer domain (add A-riboswitch, ~70 nt)
    "add_riboswitch": (
        "UUAUAUAAUUCUUGGAGGCUUCGGGCCAAGGACUAAACUCUUUGAUUCUUGGAGGCUCUUGGAGGCUUUGUUUU"
    ),
    # P4-P6 domain of Tetrahymena group I intron (a classic RNA folding benchmark)
    "p4p6_domain": (
        "GGAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGG"
        "GUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAU"
        "GGAUGCAGUUCAAAACCAAACCGUUC"
    ),
    # Hammerhead ribozyme (minimal motif, ~40 nt) — fast to fold
    "hammerhead_ribozyme": (
        "CUGAUGAGUCCGUGAGGACGAAACAGC"
        "AAAAAACGUGUCUCUCAAUGAUGAGCCC"
    ),
}

# =============================================================================
#  Helpers
# =============================================================================


def setup_logging(log_file: Optional[str], level: int = logging.INFO) -> logging.Logger:
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
    )
    return logging.getLogger(__name__)


def clean_rna_sequence(seq: str) -> str:
    """Upper-case, T→U, strip whitespace/numbers. Raises ValueError on bad chars."""
    seq = re.sub(r"[\s\d]", "", seq).upper().replace("T", "U")
    bad = set(seq) - set("ACGUN")
    if bad:
        raise ValueError(f"Invalid characters in sequence: {bad}")
    return seq


def parse_fasta(fasta_path: Path) -> list:
    """Return list of (header, sequence) tuples from a FASTA file."""
    records = []
    header, parts = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header = line[1:].split()[0]  # use first token as ID
                parts = []
            else:
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def safe_id(name: str) -> str:
    """Sanitise sequence name for use as a directory/file name."""
    return re.sub(r"[^\w\-]", "_", name)


def run(cmd: list, logger: logging.Logger, cwd: Optional[str] = None,
        env: Optional[dict] = None) -> int:
    """Run a subprocess, stream stdout/stderr to logger, return exit code."""
    logger.info("Running: %s", " ".join(str(c) for c in cmd))
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=cwd,
        env=env or os.environ.copy(),
    )
    for line in proc.stdout:
        logger.info(line.rstrip())
    proc.wait()
    return proc.returncode


# =============================================================================
#  Core prediction logic
# =============================================================================

def make_single_seq_a3m(seq_id: str, sequence: str, a3m_path: Path) -> None:
    """
    Write a minimal single-sequence A3M file for trRosettaRNA.
    Format: just a FASTA-style header + sequence, no size header line.
    """
    with open(a3m_path, "w") as f:
        f.write(f">{seq_id}\n")
        f.write(sequence + "\n")


def make_dot_bracket(sequence: str) -> str:
    """Generate a flat (all-unpaired) dot-bracket string. Used to bypass SPOT-RNA."""
    return "." * len(sequence)


# =============================================================================
#  Core prediction logic
# =============================================================================

def predict_one(
    seq_id: str,
    sequence: str,
    out_dir: Path,
    trrosetta_root: Path,
    gpu: int,
    n_cpu: int,
    n_decoys: int,
    skip_spotrna: bool,
    logger: logging.Logger,
) -> bool:
    """
    Run the full trRosettaRNA v1.1 pipeline for a single RNA sequence.

    Step 1 — geometry prediction (transformer network):
        predict.py -i <seq.a3m> -o <seq.npz> -mdir <params/model_1> -gpu <N>

    Step 2 — 3D folding (PyRosetta energy minimisation):
        fold.py -npz <seq.npz> -fa <seq.fasta> -out <model_1.pdb> --cpu <N> -nm <N>
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    predict_py = trrosetta_root / "predict.py"
    fold_py = trrosetta_root / "fold.py"
    params_dir = trrosetta_root / "params" / "model_1"

    # ── Write FASTA ────────────────────────────────────────────────────────────
    input_fasta = out_dir / f"{seq_id}.fasta"
    with open(input_fasta, "w") as f:
        f.write(f">{seq_id}\n")
        f.write("\n".join(textwrap.wrap(sequence, 60)) + "\n")

    # ── Prepare A3M (use pre-existing one if available next to the input FASTA) ─
    # Convention: place <seq_id>.a3m in the same input_dir for real MSA support.
    geom_dir = out_dir / "geometry"
    geom_dir.mkdir(exist_ok=True)
    a3m_path = out_dir / f"{seq_id}.a3m"
    if not a3m_path.exists():
        logger.info("[%s] No A3M found — generating single-sequence A3M.", seq_id)
        make_single_seq_a3m(seq_id, sequence, a3m_path)
    else:
        logger.info("[%s] Using existing A3M: %s", seq_id, a3m_path)

    # ── Prepare secondary structure (optional SPOT-RNA bypass) ────────────────
    ss_args = []
    if skip_spotrna:
        ss_file = out_dir / f"{seq_id}.dbn"
        if not ss_file.exists():
            logger.info("[%s] Writing flat dot-bracket SS (SPOT-RNA bypassed).", seq_id)
            ss_file.write_text(make_dot_bracket(sequence) + "\n")
        ss_args = ["-ss", str(ss_file.resolve()), "-ss_fmt", "dot_bracket"]

    # ── Step 1: geometry prediction ────────────────────────────────────────────
    geom_npz = geom_dir / f"{seq_id}.npz"
    logger.info("[%s] Step 1/2 — geometry prediction", seq_id)
    t0 = time.time()

    geom_cmd = [
        sys.executable, str(predict_py),
        "-i", str(a3m_path.resolve()),
        "-o", str(geom_npz.resolve()),
        "-mdir", str(params_dir),
        "-cpu", str(n_cpu),
    ] + ss_args
    if gpu >= 0:
        geom_cmd += ["-gpu", str(gpu)]

    rc = run(geom_cmd, logger, cwd=str(trrosetta_root))
    if rc != 0:
        logger.error("[%s] Geometry prediction FAILED (exit %d)", seq_id, rc)
        return False
    logger.info("[%s] Geometry done in %.1f s", seq_id, time.time() - t0)

    if not geom_npz.exists():
        logger.error("[%s] Expected NPZ not found: %s", seq_id, geom_npz)
        return False

    # ── Step 2: 3D structure folding ───────────────────────────────────────────
    struct_dir = out_dir / "structures"
    struct_dir.mkdir(exist_ok=True)
    out_pdb = struct_dir / f"{seq_id}_model1.pdb"

    logger.info("[%s] Step 2/2 — 3D folding (%d decoys, %d CPUs)", seq_id, n_decoys, n_cpu)
    t1 = time.time()

    fold_cmd = [
        sys.executable, str(fold_py),
        "-npz", str(geom_npz.resolve()),
        "-fa", str(input_fasta.resolve()),
        "-out", str(out_pdb.resolve()),
        "-nm", str(n_decoys),
        "-cpu", str(n_cpu),
    ]
    rc = run(fold_cmd, logger, cwd=str(trrosetta_root))
    if rc != 0:
        logger.error("[%s] Folding FAILED (exit %d)", seq_id, rc)
        return False
    logger.info("[%s] Folding done in %.1f s", seq_id, time.time() - t1)

    # ── Collect final PDBs ─────────────────────────────────────────────────────
    pdbs = list(struct_dir.glob("*.pdb"))
    if pdbs:
        logger.info("[%s] ✓ %d PDB(s) generated:", seq_id, len(pdbs))
        for p in sorted(pdbs):
            logger.info("     %s", p)
    else:
        logger.warning("[%s] No PDB files found after folding.", seq_id)

    return True


# =============================================================================
#  Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Batch trRosettaRNA v1.1 predictions on a directory of FASTA files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing .fa / .fasta input files.")
    parser.add_argument("--output_dir", required=True,
                        help="Root output directory. One sub-folder per sequence.")
    parser.add_argument("--trrosetta_root", required=True,
                        help="Absolute path to unzipped trRosettaRNA_v1.1/ directory.")
    parser.add_argument("--gpu", type=int, default=-1,
                        help="GPU index for geometry prediction (-1 = CPU only).")
    parser.add_argument("--n_cpu", type=int, default=8,
                        help="CPU cores for geometry prediction and Rosetta folding.")
    parser.add_argument("--n_decoys", type=int, default=5,
                        help="Number of Rosetta decoys to generate per sequence.")
    parser.add_argument("--skip_spotrna", action="store_true",
                        help="Bypass SPOT-RNA by using a flat dot-bracket SS. "
                             "Faster but lower geometry quality.")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip sequences that already have a PDB output.")
    parser.add_argument("--split_multi", action="store_true",
                        help="Process each record in multi-sequence FASTAs separately.")
    parser.add_argument("--log_file", default=None,
                        help="Optional path for a persistent log file.")
    parser.add_argument("--write_examples", action="store_true",
                        help="Write example FASTA files to --input_dir and exit.")
    args = parser.parse_args()

    logger = setup_logging(args.log_file)

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    trrosetta_root = Path(args.trrosetta_root)

    # ── Write example FASTAs and exit ─────────────────────────────────────────
    if args.write_examples:
        input_dir.mkdir(parents=True, exist_ok=True)
        for name, seq in EXAMPLE_SEQUENCES.items():
            fasta = input_dir / f"{name}.fasta"
            with open(fasta, "w") as f:
                f.write(f">{name}\n")
                f.write("\n".join(textwrap.wrap(seq, 60)) + "\n")
            logger.info("Written: %s  (%d nt)", fasta, len(seq))
        logger.info("Done. Now run without --write_examples to predict.")
        return

    # ── Validate paths ─────────────────────────────────────────────────────────
    if not input_dir.is_dir():
        logger.error("INPUT_DIR does not exist: %s", input_dir)
        sys.exit(1)
    if not trrosetta_root.is_dir():
        logger.error("trRosettaRNA root not found: %s", trrosetta_root)
        sys.exit(1)
    predict_py = trrosetta_root / "predict.py"
    if not predict_py.exists():
        logger.error("predict.py not found in trRosettaRNA root: %s", predict_py)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Collect FASTA files ────────────────────────────────────────────────────
    fasta_files = sorted(
        list(input_dir.glob("*.fa")) + list(input_dir.glob("*.fasta"))
    )
    if not fasta_files:
        logger.error("No .fa / .fasta files found in %s", input_dir)
        sys.exit(1)
    logger.info("Found %d FASTA file(s) in %s", len(fasta_files), input_dir)

    # ── Build job list ─────────────────────────────────────────────────────────
    jobs: list = []     # (seq_id, sequence)
    for fasta_path in fasta_files:
        try:
            records = parse_fasta(fasta_path)
        except Exception as e:
            logger.warning("Could not parse %s: %s — skipping", fasta_path.name, e)
            continue

        if len(records) > 1 and not args.split_multi:
            logger.warning(
                "%s contains %d records but --split_multi not set; "
                "using only first record.", fasta_path.name, len(records)
            )
            records = [records[0]]

        for header, raw_seq in records:
            try:
                seq = clean_rna_sequence(raw_seq)
            except ValueError as e:
                logger.warning("Skipping %s / %s: %s", fasta_path.name, header, e)
                continue
            if len(seq) == 0:
                logger.warning("Empty sequence in %s / %s — skipping", fasta_path.name, header)
                continue
            jobs.append((safe_id(header), seq))

    if not jobs:
        logger.error("No valid sequences to process. Exiting.")
        sys.exit(1)

    logger.info("Sequences to predict: %d", len(jobs))

    # ── Run predictions ────────────────────────────────────────────────────────
    successes, skipped, failures = 0, 0, 0

    for seq_id, sequence in jobs:
        seq_out = output_dir / seq_id

        # Skip if already done
        if args.skip_existing and list(seq_out.glob("**/*.pdb")):
            logger.info("[%s] Already has PDB output — skipping.", seq_id)
            skipped += 1
            continue

        logger.info("=" * 60)
        logger.info("Processing: %s  (%d nt)", seq_id, len(sequence))
        logger.info("=" * 60)

        ok = predict_one(
            seq_id=seq_id,
            sequence=sequence,
            out_dir=seq_out,
            trrosetta_root=trrosetta_root,
            gpu=args.gpu,
            n_cpu=args.n_cpu,
            n_decoys=args.n_decoys,
            skip_spotrna=args.skip_spotrna,
            logger=logger,
        )
        if ok:
            successes += 1
        else:
            failures += 1

    # ── Summary ───────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("DONE — success: %d | skipped: %d | failed: %d",
                successes, skipped, failures)
    logger.info("=" * 60)

    if failures > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
