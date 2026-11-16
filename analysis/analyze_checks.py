#!/usr/bin/env python3
"""Analyze an exported check run: do the best HPO configs hold up?

Reads the check run's ``results/`` and the HPO export it checks (for the HPO
losses and the training RNAs). Writes into ``analysis/`` next to ``results/``:

- ``trials.csv``: per trial, the loss on all, train and held-out RNAs.
- ``summary.csv``: per config (and baseline), HPO loss, reproduction, and the
  mean and standard deviation over the new seeds on train and held-out RNAs.
- ``stability.csv``: per config and RNA, the mean pairwise F1 between the
  seeds' predicted pairs (1 = the same prediction for every seed; RNAs whose
  scoring failed are left out), and the mean and standard deviation of the
  RNA's loss (a failed RNA counts with loss 1.0).
- ``predictions.txt``: per RNA the target and every trial's prediction in
  dot-bracket notation, with its loss.
- ``seed_spread.png``, ``per_rna_loss.png`` (where a config gains over the
  baseline) and ``stability.png``.
"""

import argparse
import csv
import json
import logging
import math
import statistics
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from analysis.analyze_results import _configure_logging, prepare_output_directory
from analysis.modules.dotbracket import pairs_to_dotbracket

LOGGER = logging.getLogger("analyze_checks")
SERIES = "#2a78d6"
INK = "#0b0b0b"
MUTED = "#52514e"
# Categorical slots in fixed order; the first three stay apart for color-blind readers.
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


def _mean(values):
    return statistics.fmean(values) if values else math.nan


def _sd(values):
    return statistics.stdev(values) if len(values) > 1 else math.nan


def _f1(a: set, b: set) -> float:
    return 1.0 if not a and not b else 2 * len(a & b) / (len(a) + len(b))


def load(results: Path, hpo_results: Path) -> tuple[list[dict], dict[str, float], dict[str, str], dict[str, str]]:
    """Trials with their results, HPO losses by config id, split and sequence by RNA id."""
    manifest = json.loads((results / "manifest.json").read_text())
    dataset = json.loads((results / "run_metadata.json").read_text())["dataset"]
    hpo_rows = {row["id"]: row["sequence"] for row in json.loads(
        (hpo_results / "run_metadata.json").read_text())["dataset"]}
    split = {}
    for row in dataset:
        if row["id"] in hpo_rows and hpo_rows[row["id"]] != row["sequence"]:
            raise ValueError(f"RNA {row['id']} has another sequence in the HPO dataset")
        split[row["id"]] = "train" if row["id"] in hpo_rows else "heldout"

    with (hpo_results / "summary" / "full.csv").open(newline="") as stream:
        hpo_losses = {row["id"]: float(row["objective_to_minimize"])
                      for row in csv.DictReader(stream) if row["state"] == "success"}

    trials = []
    for trial in manifest["trials"]:
        path = results / "trials" / trial["name"] / "trial_result.json"
        if not path.is_file():
            LOGGER.warning("Trial %s has no result (failed or not run); skipped", trial["name"])
            continue
        result = json.loads(path.read_text())
        base = trial.get("base")
        group = f"c{base}" if base is not None else "no MSA" if trial.get("check") == "baseline" else trial["name"]
        trials.append(dict(trial, group=group, check=trial.get("check", "seed"), result=result))
    return trials, hpo_losses, split, {row["id"]: row["sequence"] for row in dataset}


def trial_rows(trials: list[dict], split: dict[str, str]) -> list[dict]:
    rows = []
    for trial in trials:
        scores = trial["result"]["scores"]
        losses = {part: [s["loss"] for s in scores if split[s["id"]] == part] for part in ("train", "heldout")}
        rows.append(dict(
            name=trial["name"], group=trial["group"], check=trial["check"], base=trial.get("base"),
            shs_seed=trial["result"].get("shs_seed"), af3_seed=trial["result"].get("af3_seed"),
            loss=trial["result"]["loss"], loss_train=_mean(losses["train"]), loss_heldout=_mean(losses["heldout"]),
        ))
    return rows


def summary_rows(rows: list[dict], hpo_losses: dict[str, float]) -> list[dict]:
    summary = []
    for group in dict.fromkeys(row["group"] for row in rows):
        members = [row for row in rows if row["group"] == group]
        seeds = [row for row in members if row["check"] != "reproduction"]
        reproduction = [row["loss_train"] for row in members if row["check"] == "reproduction"]
        base = members[0]["base"]
        hpo = hpo_losses.get(str(base), math.nan) if base is not None else math.nan
        summary.append(dict(
            group=group, hpo_loss=hpo,
            reproduction_train=reproduction[0] if reproduction else math.nan,
            reproduction_diff=reproduction[0] - hpo if reproduction else math.nan,
            seeds=len(seeds),
            train_mean=_mean([r["loss_train"] for r in seeds]), train_sd=_sd([r["loss_train"] for r in seeds]),
            heldout_mean=_mean([r["loss_heldout"] for r in seeds if not math.isnan(r["loss_heldout"])]),
            heldout_sd=_sd([r["loss_heldout"] for r in seeds if not math.isnan(r["loss_heldout"])]),
        ))
    return summary


def stability_rows(trials: list[dict], split: dict[str, str]) -> list[dict]:
    rows = []
    for group in dict.fromkeys(trial["group"] for trial in trials):
        seeds = [t for t in trials if t["group"] == group and t["check"] != "reproduction"]
        by_rna = {}
        for trial in seeds:
            for score in trial["result"]["scores"]:
                # A failed RNA (scoring error) has no prediction to compare.
                pairs = None if "error" in score else {tuple(pair) for pair in score["predicted_pairs"]}
                by_rna.setdefault(score["id"], []).append((pairs, score["loss"]))
        for rna_id, entries in by_rna.items():
            predictions = [pairs for pairs, _ in entries if pairs is not None]
            losses = [loss for _, loss in entries]
            rows.append(dict(
                group=group, rna_id=rna_id, split=split[rna_id], seeds=len(entries),
                pairwise_f1=_mean([_f1(a, b) for a, b in combinations(predictions, 2)]),
                loss_mean=_mean(losses), loss_sd=_sd(losses),
            ))
    return rows


def write_predictions(trials: list[dict], split: dict[str, str], sequences: dict[str, str], path: Path) -> None:
    """Target and every trial's prediction per RNA as dot-bracket lines.

    Pairs that dot-bracket cannot show (crossing, or a base in two pairs) are
    left out and counted at the end of the line.
    """
    width = max(len(t["name"]) for t in trials)
    lines = []
    for rna_id, sequence in sequences.items():
        lines.append(f"{rna_id} ({split[rna_id]})")
        lines.append(f"  {'sequence':<{width}}  {sequence}")
        target = None
        for trial in trials:
            score = next(s for s in trial["result"]["scores"] if s["id"] == rna_id)
            if target is None:
                structure, skipped = pairs_to_dotbracket(score["target_pairs"], len(sequence))
                lines.append(f"  {'target':<{width}}  {structure}" + (f"  [{len(skipped)} pairs not shown]" if skipped else ""))
                target = structure
            structure, skipped = pairs_to_dotbracket(score["predicted_pairs"], len(sequence))
            note = f"  [{len(skipped)} pairs not shown]" if skipped else ""
            if "error" in score:
                structure, note = f"{'(scoring failed)':<{len(sequence)}}", f"  [{score['error'][:80]}]"
            lines.append(f"  {trial['name']:<{width}}  {structure}  loss {score['loss']:.3f}{note}")
        lines.append("")
    path.write_text("\n".join(lines))


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _style(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    ax.grid(axis="y", color="#e5e4e0", linewidth=0.8)
    ax.set_axisbelow(True)


def plot_seed_spread(rows: list[dict], summary: list[dict], path: Path) -> None:
    """Loss per config: new seeds (dots), their mean, reproduction and HPO loss."""
    groups = [row["group"] for row in summary]
    parts = [("train", "Train RNAs")] + (
        [("heldout", "Held-out RNAs")] if any(not math.isnan(r["loss_heldout"]) for r in rows) else [])
    fig, axes = plt.subplots(1, len(parts), figsize=(1.5 + 1.3 * len(groups) * len(parts), 4),
                             sharey=True, squeeze=False)
    for ax, (part, title) in zip(axes[0], parts):
        _style(ax)
        for x, group in enumerate(groups):
            seeds = [r[f"loss_{part}"] for r in rows if r["group"] == group and r["check"] != "reproduction"]
            # Spread the seeds a little, so equal losses stay visible.
            offsets = [0.2 * (i / (len(seeds) - 1) - 0.5) if len(seeds) > 1 else 0 for i in range(len(seeds))]
            ax.scatter([x + o for o in offsets], seeds, s=36, color=SERIES, edgecolor="white", linewidth=1,
                       zorder=3, label="new seeds")
            if seeds:
                ax.hlines(_mean(seeds), x - 0.25, x + 0.25, color=INK, linewidth=2, zorder=4, label="seed mean")
            if part == "train":
                row = summary[x]
                if not math.isnan(row["hpo_loss"]):
                    ax.hlines(row["hpo_loss"], x - 0.35, x + 0.35, color=MUTED, linewidth=1.5,
                              linestyles="dashed", zorder=4, label="HPO loss")
                ax.scatter([x + 0.3], [row["reproduction_train"]], marker="o", s=70, facecolor="none",
                           edgecolor=INK, linewidth=1.5, zorder=5, label="reproduction (HPO seeds)")
        ax.set_xticks(range(len(groups)), groups)
        ax.set_title(title, color=INK, loc="left")
    axes[0][0].set_ylabel("loss (mean over RNAs)", color=INK)
    handles, labels = axes[0][0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    fig.legend(unique.values(), unique.keys(), loc="upper center", ncol=len(unique), frameon=False,
               bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_per_rna(trials: list[dict], split: dict[str, str], path: Path) -> None:
    """Loss per RNA: each config's seed mean (dot) and range (line); baseline as ring.

    Shows on which RNAs a config gains over the baseline, and where it depends on the seed.
    """
    seeds = [t for t in trials if t["check"] != "reproduction"]
    groups = list(dict.fromkeys(t["group"] for t in seeds))
    rna_ids = list(split)
    losses = {}
    for trial in seeds:
        for score in trial["result"]["scores"]:
            losses.setdefault((trial["group"], score["id"]), []).append(score["loss"])
    configs = [g for g in groups if g != "no MSA"]
    step = 0.6 / max(len(groups) - 1, 1)
    fig, ax = plt.subplots(figsize=(7, 1.2 + 0.45 * len(rna_ids)))
    _style(ax)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color="#e5e4e0", linewidth=0.8)
    for index, group in enumerate(groups):
        offset = -0.3 + index * step if len(groups) > 1 else 0
        baseline = group == "no MSA"
        color = MUTED if baseline else CATEGORICAL[configs.index(group) % len(CATEGORICAL)]
        for y, rna_id in enumerate(rna_ids):
            values = losses.get((group, rna_id))
            if not values:
                continue
            ax.hlines(y + offset, min(values), max(values), color=color, linewidth=1.5, zorder=3)
            ax.scatter([_mean(values)], [y + offset], s=40, zorder=4, linewidth=1.5,
                       facecolor="none" if baseline else color, edgecolor=color if baseline else "white",
                       label=group if y == 0 else None)
    heldout = [y for y, rna_id in enumerate(rna_ids) if split[rna_id] == "heldout"]
    if heldout:
        ax.axhline(heldout[0] - 0.5, color=MUTED, linewidth=0.8)
    ax.set_yticks(range(len(rna_ids)), [f"{r} (held-out)" if split[r] == "heldout" else r for r in rna_ids])
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("loss (dot: mean over seeds, line: seed range)", color=INK)
    ax.set_title("Loss per RNA", color=INK, loc="left")
    fig.legend(loc="upper center", ncol=len(groups), frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.5 / (1.2 + 0.45 * len(rna_ids))))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_stability(rows: list[dict], path: Path) -> None:
    """Heatmap: mean pairwise F1 of the seeds' predictions per RNA and config."""
    groups = list(dict.fromkeys(row["group"] for row in rows))
    rna_ids = list(dict.fromkeys(row["rna_id"] for row in rows))
    value = {(row["group"], row["rna_id"]): row["pairwise_f1"] for row in rows}
    split = {row["rna_id"]: row["split"] for row in rows}
    grid = [[value.get((group, rna_id), math.nan) for group in groups] for rna_id in rna_ids]
    fig, ax = plt.subplots(figsize=(2.5 + 0.9 * len(groups), 1.2 + 0.32 * len(rna_ids)))
    image = ax.imshow(grid, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(groups)), groups)
    ax.set_yticks(range(len(rna_ids)), [f"{r} (held-out)" if split[r] == "heldout" else r for r in rna_ids])
    ax.tick_params(colors=MUTED, labelcolor=INK, length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Prediction agreement between seeds (pairwise F1)", color=INK, loc="left")
    fig.colorbar(image, ax=ax, fraction=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def analyze(results: Path, hpo_results: Path, output: Path) -> list[dict]:
    trials, hpo_losses, split, sequences = load(results, hpo_results)
    if not trials:
        raise ValueError(f"No finished trials in {results}")
    rows = trial_rows(trials, split)
    summary = summary_rows(rows, hpo_losses)
    stability = stability_rows(trials, split)
    write_csv(rows, output / "trials.csv")
    write_csv(summary, output / "summary.csv")
    write_csv(stability, output / "stability.csv")
    write_predictions(trials, split, sequences, output / "predictions.txt")
    plot_seed_spread(rows, summary, output / "seed_spread.png")
    plot_per_rna(trials, split, output / "per_rna_loss.png")
    plot_stability(stability, output / "stability.png")
    for row in summary:
        LOGGER.info(
            "%s: HPO %.4f | reproduction %+.4f | train %.4f ± %.4f | held-out %.4f ± %.4f (%d seeds)",
            row["group"], row["hpo_loss"], row["reproduction_diff"], row["train_mean"], row["train_sd"],
            row["heldout_mean"], row["heldout_sd"], row["seeds"],
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", type=Path, help="Exported check run (checks/results)")
    parser.add_argument("--hpo-results", type=Path,
                        help="Exported HPO run that was checked (default: results/ two levels up)")
    parser.add_argument("--output", type=Path, help="Output directory (default: sibling analysis directory)")
    args = parser.parse_args()
    _configure_logging()
    results = args.results.expanduser().resolve()
    hpo_results = (args.hpo_results or results.parent.parent / "results").expanduser().resolve()
    output = (args.output or results.parent / "analysis").expanduser().resolve()
    prepare_output_directory(output)
    analyze(results, hpo_results, output)
    LOGGER.info("Check analysis written to %s", output)


if __name__ == "__main__":
    main()
