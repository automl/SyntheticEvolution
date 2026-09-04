#!/usr/bin/env python3
"""Run SyntheticEvolution evaluation for converted TRRosetta/RhoFold CIF predictions."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# This folder lives inside the SyntheticEvolution checkout, so the repo root is
# its parent. Anchoring on the file, not the cwd, keeps the defaults valid when
# the script is called from elsewhere.
REPO_ROOT = Path(__file__).resolve().parent.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate converted predictions via SyntheticEvolution/run_evaluation.sh",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--algorithm",
        required=True,
        help="Output label used in copied CSV names (e.g. trrosetta, trrosetta_stdmsa).",
    )
    parser.add_argument("--pred-cif-dir", required=True, help="Directory with converted prediction CIFs.")
    parser.add_argument("--synthevolution-root", default=str(REPO_ROOT),
                        help="SyntheticEvolution repository root (this folder's parent).")
    parser.add_argument("--gt-dir", default=None, help="Defaults to <synthevolution-root>/evaluation/predictions/gt")
    parser.add_argument("--eval-type", default="rna_rna", choices=["rna_rna", "protein_rna"])
    parser.add_argument("--output-dir", default="./evaluation_results")
    return parser.parse_args()


def parse_float_like(value: str) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None

    # Handles values like "[1.23, 1.10]" by taking the first numeric token.
    match = re.search(r"-?\d+(?:\.\d+)?", text)
    if not match:
        return None
    return float(match.group(0))


def summarize_csv(csv_path: Path, out_summary_path: Path) -> None:
    if not csv_path.is_file():
        return

    numeric_cols = [
        "Complex_RMSD",
        "RNA_RMSD",
        "Protein_RMSD",
        "Complex_TM",
        "RNA_TM",
        "Protein_TM",
        "Complex_LDDT",
        "RNA_LDDT",
        "Protein_LDDT",
    ]

    sums = {k: 0.0 for k in numeric_cols}
    counts = {k: 0 for k in numeric_cols}
    n_rows = 0

    with open(csv_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            n_rows += 1
            for col in numeric_cols:
                if col not in row:
                    continue
                val = parse_float_like(row[col])
                if val is None:
                    continue
                sums[col] += val
                counts[col] += 1

    lines = [f"rows={n_rows}"]
    for col in numeric_cols:
        if counts[col] > 0:
            lines.append(f"{col}_mean={sums[col] / counts[col]:.4f}")

    out_summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def sanitize_algorithm_label(label: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", label.strip())
    return safe or "algorithm"


def main() -> None:
    args = parse_args()
    algo_label = sanitize_algorithm_label(args.algorithm)

    synthetic_root = Path(args.synthevolution_root).resolve()
    pred_cif_dir = Path(args.pred_cif_dir).resolve()
    gt_dir = Path(args.gt_dir).resolve() if args.gt_dir else (synthetic_root / "evaluation" / "predictions" / "gt")
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    run_eval_sh = synthetic_root / "run_evaluation.sh"
    results_dir = synthetic_root / "results"

    if not run_eval_sh.is_file():
        raise FileNotFoundError(f"run_evaluation.sh not found: {run_eval_sh}")
    if not pred_cif_dir.is_dir():
        raise FileNotFoundError(f"Prediction CIF dir not found: {pred_cif_dir}")
    if not gt_dir.is_dir():
        raise FileNotFoundError(f"GT dir not found: {gt_dir}")

    cmd = [
        "bash",
        str(run_eval_sh),
        str(pred_cif_dir),
        str(gt_dir),
        args.eval_type,
    ]

    print("Running evaluation command:")
    print(" ".join(cmd))
    env = os.environ.copy()
    python_bin = str(Path(sys.executable).resolve().parent)
    env["PATH"] = f"{python_bin}:{env.get('PATH', '')}"
    subprocess.run(cmd, cwd=str(synthetic_root), check=True, env=env)

    pred_csv_name = f"pred_{args.eval_type}.csv"
    exp_csv_name = f"exp_{args.eval_type}.csv"

    src_pred_csv = results_dir / pred_csv_name
    src_exp_csv = results_dir / exp_csv_name

    dst_pred_csv = output_dir / f"{algo_label}_{pred_csv_name}"
    dst_exp_csv = output_dir / f"{algo_label}_{exp_csv_name}"

    if src_pred_csv.is_file():
        shutil.copy2(src_pred_csv, dst_pred_csv)
        print(f"Copied: {dst_pred_csv}")
    else:
        print(f"Warning: missing {src_pred_csv}")

    if src_exp_csv.is_file():
        shutil.copy2(src_exp_csv, dst_exp_csv)
        print(f"Copied: {dst_exp_csv}")
    else:
        print(f"Warning: missing {src_exp_csv}")

    summary_path = output_dir / f"{algo_label}_{args.eval_type}_summary.txt"
    summarize_csv(dst_pred_csv, summary_path)
    if summary_path.is_file():
        print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
