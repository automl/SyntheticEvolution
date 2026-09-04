#!/usr/bin/env python3
"""Compare two evaluation CSVs for the same algorithm (e.g. SHS vs standard MSA)."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


METRICS = ["Complex_RMSD", "RNA_RMSD", "Complex_TM", "RNA_TM", "Complex_LDDT", "RNA_LDDT"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare two evaluation CSVs on their common target IDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--left-csv", required=True, help="First CSV, e.g. SHS run.")
    parser.add_argument("--right-csv", required=True, help="Second CSV, e.g. standard MSA run.")
    parser.add_argument("--left-label", default="shs")
    parser.add_argument("--right-label", default="stdmsa")
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def parse_float_like(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    m = re.search(r"-?\d+(?:\.\d+)?", text)
    if not m:
        return None
    return float(m.group(0))


def detect_id(row: dict[str, str]) -> str | None:
    for key in ("exp_db_id", "PDBId", "pdb_id", "id"):
        value = row.get(key)
        if value:
            token = value.strip().upper()
            if re.fullmatch(r"[A-Z0-9]{4}", token):
                return token

    file_name = row.get("FileName") or row.get("file_name") or ""
    if file_name:
        stem = Path(file_name).stem.upper()
        m = re.search(r"([A-Z0-9]{4})", stem)
        if m:
            return m.group(1)
    return None


def read_metric_map(path: Path) -> dict[str, dict[str, float | None]]:
    out: dict[str, dict[str, float | None]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            pid = detect_id(row)
            if not pid:
                continue
            out[pid] = {metric: parse_float_like(row.get(metric)) for metric in METRICS}
    return out


def mean(values: list[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def maybe_plot(summary_rows: list[dict[str, str | float | None]], out_path: Path) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return

    labels = [str(summary_rows[0]["label"]), str(summary_rows[1]["label"])]
    metrics = METRICS
    x = list(range(len(metrics)))
    width = 0.35

    fig, ax = plt.subplots(figsize=(11, 5))
    for idx, row in enumerate(summary_rows):
        vals = [row[m] for m in metrics]
        vals = [float(v) if v is not None else float("nan") for v in vals]
        offsets = [xi + (idx - 0.5) * width for xi in x]
        ax.bar(offsets, vals, width=width, label=labels[idx])

    ax.set_xticks(x)
    ax.set_xticklabels(metrics, rotation=30, ha="right")
    ax.set_title("Mean Metrics On Common ID Set")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()

    left_csv = Path(args.left_csv).resolve()
    right_csv = Path(args.right_csv).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    for required in (left_csv, right_csv):
        if not required.is_file():
            raise FileNotFoundError(f"Missing input CSV: {required}")

    left_map = read_metric_map(left_csv)
    right_map = read_metric_map(right_csv)

    left_ids = set(left_map.keys())
    right_ids = set(right_map.keys())
    common_ids = sorted(left_ids & right_ids)

    merged_fields = ["pdb_id"]
    for metric in METRICS:
        merged_fields += [f"{args.left_label}_{metric}", f"{args.right_label}_{metric}", f"delta_{args.left_label}_minus_{args.right_label}_{metric}"]

    with open(out_dir / "merged_per_target.csv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=merged_fields)
        writer.writeheader()
        for pid in common_ids:
            row: dict[str, str | float | None] = {"pdb_id": pid}
            for metric in METRICS:
                lv = left_map[pid].get(metric)
                rv = right_map[pid].get(metric)
                row[f"{args.left_label}_{metric}"] = lv
                row[f"{args.right_label}_{metric}"] = rv
                row[f"delta_{args.left_label}_minus_{args.right_label}_{metric}"] = None if lv is None or rv is None else lv - rv
            writer.writerow(row)

    summary_rows: list[dict[str, str | float | None]] = []
    for label, metric_map in ((args.left_label, left_map), (args.right_label, right_map)):
        row: dict[str, str | float | None] = {"label": label}
        for metric in METRICS:
            vals = [v for pid in common_ids if (v := metric_map[pid][metric]) is not None]
            row[metric] = mean(vals)
        summary_rows.append(row)

    with open(out_dir / "summary_means.csv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["label", *METRICS])
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    overlap_lines = [
        f"{args.left_label}_ids={len(left_ids)}",
        f"{args.right_label}_ids={len(right_ids)}",
        f"common_ids={len(common_ids)}",
    ]
    (out_dir / "summary_overlap_counts.txt").write_text("\n".join(overlap_lines) + "\n", encoding="utf-8")

    maybe_plot(summary_rows, out_dir / "mean_metrics_barplot.png")

    print(f"{args.left_label} IDs:  {len(left_ids)}")
    print(f"{args.right_label} IDs: {len(right_ids)}")
    print(f"Common IDs:       {len(common_ids)}")
    print(f"Output dir:       {out_dir}")


if __name__ == "__main__":
    main()
