"""Plot the full optimization history from NePS's ``summary/full.csv``."""

import csv
import logging
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis.results_loader import ResultsContext

LOGGER = logging.getLogger("analysis.optimization_history")


def _float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return math.nan


def _plot_loss_over_time(success: list[dict], crashed: list[dict], start: float, path: Path) -> None:
    hours = lambda row: (_float(row["time_end"]) - start) / 3600
    success = sorted(success, key=hours)
    x = [hours(row) for row in success]
    y = [_float(row["objective_to_minimize"]) for row in success]
    best = [min(y[: i + 1]) for i in range(len(y))]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.scatter(x, y, s=10, alpha=0.5, color="tab:blue", label=f"success ({len(success)})")
    ax.step(x, best, where="post", color="black", label="best so far")
    crashed_x = [hours(row) for row in crashed if not math.isnan(hours(row))]
    if crashed_x:
        ax.plot(crashed_x, [1.0] * len(crashed_x), "|", color="tab:red", markersize=12,
                transform=ax.get_xaxis_transform(), label=f"crashed ({len(crashed)})")
    ax.set_xlabel("hours since run start (trial end)")
    ax.set_ylabel("loss")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _plot_loss_histogram(losses: list[float], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(losses, bins=40, color="tab:blue")
    ax.axvline(min(losses), color="black", linestyle="--", label=f"best {min(losses):.3f}")
    ax.set_xlabel("loss")
    ax.set_ylabel("successful trials")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _plot_parameters(rows: list[dict], path: Path) -> None:
    columns = [
        column for column in rows[0]
        if column.startswith("config.")
        and len({row[column] for row in rows}) > 1
        and not any(math.isnan(_float(row[column])) for row in rows)
    ]
    if not columns:
        LOGGER.warning("No varying numeric parameters in full.csv; skipping parameter plot")
        return
    ncols = min(3, len(columns))
    nrows = math.ceil(len(columns) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.2 * nrows), squeeze=False)
    for ax, column in zip(axes.flat, columns):
        for state, color in (("success", "tab:blue"), ("crashed", "tab:red")):
            subset = [row for row in rows if row["state"] == state]
            # Crashed trials have no real loss; draw them at the top edge.
            y = [_float(row["objective_to_minimize"]) if state == "success" else 1.0 for row in subset]
            ax.scatter([_float(row[column]) for row in subset], y, s=8, alpha=0.5, color=color,
                       label=state, transform=ax.transData if state == "success" else ax.get_xaxis_transform())
        ax.set_title(column.removeprefix("config."), fontsize=9)
        ax.set_ylabel("loss")
    for ax in axes.flat[len(columns):]:
        ax.set_visible(False)
    axes.flat[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run(context: ResultsContext, output_dir: Path) -> None:
    """Plot loss over time, the loss distribution and each parameter against loss."""
    full_csv = context.results_dir / "summary" / "full.csv"
    if not full_csv.is_file():
        LOGGER.warning("Skipping optimization history: %s is missing (re-export the run)", full_csv)
        return
    with full_csv.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    success = [row for row in rows if row["state"] == "success"]
    crashed = [row for row in rows if row["state"] == "crashed"]
    if not success:
        raise ValueError(f"No successful trials in {full_csv}")
    start = min(_float(row["time_started"]) for row in rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    _plot_loss_over_time(success, crashed, start, output_dir / "loss_over_time.png")
    _plot_loss_histogram([_float(row["objective_to_minimize"]) for row in success],
                         output_dir / "loss_histogram.png")
    _plot_parameters(success + crashed, output_dir / "parameters_vs_loss.png")
    LOGGER.info("Plotted %d trials (%d success, %d crashed, %d other)", len(rows),
                len(success), len(crashed), len(rows) - len(success) - len(crashed))
