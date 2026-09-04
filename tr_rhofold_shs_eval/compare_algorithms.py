#!/usr/bin/env python3
"""Aggregate and compare RNA-RNA evaluation results across algorithms.

Inputs:
- TRRosetta eval CSV (typically trrosetta_pred_rna_rna.csv)
- RhoFold eval CSV (typically rhofold_pred_rna_rna.csv)
- AlphaFold baseline CSV (default: ../results/csvs/All_alphafold.csv)

Outputs:
- merged_per_target.csv
- summary_means.csv
- summary_overlap_counts.txt
- mean_metrics_barplot.png (if matplotlib is available)
"""

from __future__ import annotations

import argparse
import csv
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Dict, Iterable, List


# Repo root: this folder sits inside the SyntheticEvolution checkout.
REPO_ROOT = Path(__file__).resolve().parent.parent

METRICS = ["Complex_RMSD", "RNA_RMSD", "Complex_TM", "RNA_TM", "Complex_LDDT", "RNA_LDDT"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare TRRosetta/RhoFold evaluation output against AlphaFold baseline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--trrosetta-csv", required=True)
    parser.add_argument("--rhofold-csv", required=True)
    parser.add_argument(
        "--alphafold-csv",
        default=str(REPO_ROOT / "results" / "csvs" / "All_alphafold.csv"),
        help="AlphaFold baseline CSV with exp_db_id and metric columns.",
    )
    parser.add_argument("--output-dir", default="./comparison")
    return parser.parse_args()


def parse_float_like(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    match = re.search(r"-?\d+(?:\.\d+)?", text)
    if not match:
        return None
    return float(match.group(0))


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return list(reader)


def detect_id(row: Dict[str, str]) -> str | None:
    for key in ("exp_db_id", "PDBId", "pdb_id", "id"):
        value = row.get(key)
        if value:
            token = str(value).strip().upper()
            if re.fullmatch(r"[A-Z0-9]{4}", token):
                return token

    file_name = row.get("FileName") or row.get("file_name") or ""
    if file_name:
        stem = Path(file_name).stem.upper()
        m = re.search(r"([A-Z0-9]{4})", stem)
        if m:
            return m.group(1)
    return None


def to_metric_map(rows: Iterable[Dict[str, str]]) -> Dict[str, Dict[str, float | None]]:
    out: Dict[str, Dict[str, float | None]] = {}
    for row in rows:
        pid = detect_id(row)
        if not pid:
            continue
        values = {}
        for metric in METRICS:
            values[metric] = parse_float_like(row.get(metric))
        out[pid] = values
    return out


def mean(values: List[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def write_merged_table(
    out_path: Path,
    common_ids: List[str],
    tr_map: Dict[str, Dict[str, float | None]],
    rho_map: Dict[str, Dict[str, float | None]],
    af_map: Dict[str, Dict[str, float | None]],
) -> None:
    fieldnames = ["pdb_id"]
    for prefix in ("trrosetta", "rhofold", "alphafold"):
        for metric in METRICS:
            fieldnames.append(f"{prefix}_{metric}")
    for metric in METRICS:
        fieldnames.append(f"delta_trrosetta_vs_alphafold_{metric}")
        fieldnames.append(f"delta_rhofold_vs_alphafold_{metric}")

    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for pid in common_ids:
            row: Dict[str, str | float | None] = {"pdb_id": pid}
            for metric in METRICS:
                tr_val = tr_map[pid].get(metric)
                rho_val = rho_map[pid].get(metric)
                af_val = af_map[pid].get(metric)

                row[f"trrosetta_{metric}"] = tr_val
                row[f"rhofold_{metric}"] = rho_val
                row[f"alphafold_{metric}"] = af_val

                row[f"delta_trrosetta_vs_alphafold_{metric}"] = (
                    None if tr_val is None or af_val is None else tr_val - af_val
                )
                row[f"delta_rhofold_vs_alphafold_{metric}"] = (
                    None if rho_val is None or af_val is None else rho_val - af_val
                )

            writer.writerow(row)


def write_summary(
    out_csv: Path,
    common_ids: List[str],
    tr_map: Dict[str, Dict[str, float | None]],
    rho_map: Dict[str, Dict[str, float | None]],
    af_map: Dict[str, Dict[str, float | None]],
) -> Dict[str, Dict[str, float | None]]:
    summary: Dict[str, Dict[str, float | None]] = {
        "trrosetta": {},
        "rhofold": {},
        "alphafold": {},
    }

    for metric in METRICS:
        tr_vals = [v for pid in common_ids if (v := tr_map[pid][metric]) is not None]
        rho_vals = [v for pid in common_ids if (v := rho_map[pid][metric]) is not None]
        af_vals = [v for pid in common_ids if (v := af_map[pid][metric]) is not None]

        summary["trrosetta"][metric] = mean(tr_vals)
        summary["rhofold"][metric] = mean(rho_vals)
        summary["alphafold"][metric] = mean(af_vals)

    fieldnames = ["algorithm"] + METRICS
    with open(out_csv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for algo in ("trrosetta", "rhofold", "alphafold"):
            row: Dict[str, str | float | None] = {"algorithm": algo}
            row.update(summary[algo])
            writer.writerow(row)

    return summary


def write_overlap_counts(path: Path, tr_ids: set[str], rho_ids: set[str], af_ids: set[str], common: set[str]) -> None:
    lines = [
        f"trrosetta_ids={len(tr_ids)}",
        f"rhofold_ids={len(rho_ids)}",
        f"alphafold_ids={len(af_ids)}",
        f"common_ids={len(common)}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def maybe_plot(summary: Dict[str, Dict[str, float | None]], out_path: Path) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return

    metrics_to_plot = [m for m in METRICS if summary["alphafold"].get(m) is not None]
    if not metrics_to_plot:
        return

    algos = ["trrosetta", "rhofold", "alphafold"]
    x = list(range(len(metrics_to_plot)))
    width = 0.25

    fig, ax = plt.subplots(figsize=(12, 5))
    for i, algo in enumerate(algos):
        vals = [summary[algo].get(m) for m in metrics_to_plot]
        vals = [float(v) if v is not None else float("nan") for v in vals]
        offsets = [xi + (i - 1) * width for xi in x]
        ax.bar(offsets, vals, width=width, label=algo)

    ax.set_xticks(x)
    ax.set_xticklabels(metrics_to_plot, rotation=30, ha="right")
    ax.set_title("Mean Metrics On Common ID Set")
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()

    tr_csv = Path(args.trrosetta_csv).resolve()
    rho_csv = Path(args.rhofold_csv).resolve()
    af_csv = Path(args.alphafold_csv).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    for required in (tr_csv, rho_csv, af_csv):
        if not required.is_file():
            raise FileNotFoundError(f"Missing input CSV: {required}")

    tr_map = to_metric_map(read_csv_rows(tr_csv))
    rho_map = to_metric_map(read_csv_rows(rho_csv))
    af_map = to_metric_map(read_csv_rows(af_csv))

    tr_ids = set(tr_map.keys())
    rho_ids = set(rho_map.keys())
    af_ids = set(af_map.keys())
    common = tr_ids & rho_ids & af_ids

    common_ids = sorted(common)

    write_merged_table(out_dir / "merged_per_target.csv", common_ids, tr_map, rho_map, af_map)
    summary = write_summary(out_dir / "summary_means.csv", common_ids, tr_map, rho_map, af_map)
    write_overlap_counts(out_dir / "summary_overlap_counts.txt", tr_ids, rho_ids, af_ids, common)
    maybe_plot(summary, out_dir / "mean_metrics_barplot.png")

    print(f"TRRosetta IDs: {len(tr_ids)}")
    print(f"RhoFold IDs:   {len(rho_ids)}")
    print(f"AlphaFold IDs: {len(af_ids)}")
    print(f"Common IDs:    {len(common_ids)}")
    print(f"Output dir:    {out_dir}")


if __name__ == "__main__":
    main()
