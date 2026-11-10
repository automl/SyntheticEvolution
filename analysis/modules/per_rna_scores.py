"""Plot the per-RNA loss of every trajectory configuration."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis.results_loader import ResultsContext


def _losses(context: ResultsContext) -> dict[str, list[float]]:
    """Per-RNA loss for each trajectory config, keyed by RNA chain ID."""
    losses: dict[str, list[float]] = {}
    for index, config_id in enumerate(context.config_ids):
        for score in context.trial_results[config_id]["scores"]:
            losses.setdefault(score["id"], [float("nan")] * len(context.config_ids))[index] = score["loss"]
    return losses


def run(context: ResultsContext, output_dir: Path) -> None:
    """Write a heatmap and a line plot of per-RNA loss along the trajectory."""
    losses = _losses(context)
    rna_ids = list(losses)
    labels = [str(config_id) for config_id in context.config_ids]
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(1 + 0.55 * len(labels), 1 + 0.4 * len(rna_ids)))
    image = ax.imshow([losses[rna_id] for rna_id in rna_ids], cmap="viridis_r", vmin=0, vmax=1, aspect="auto")
    for row, rna_id in enumerate(rna_ids):
        for column, value in enumerate(losses[rna_id]):
            ax.text(column, row, f"{value:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if value > 0.5 else "black")
    ax.set_xticks(range(len(labels)), labels, fontsize=8)
    ax.set_yticks(range(len(rna_ids)), rna_ids, fontsize=8)
    ax.set_xlabel("trajectory config ID")
    fig.colorbar(image, ax=ax, label="loss")
    fig.tight_layout()
    fig.savefig(output_dir / "per_rna_heatmap.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(labels)), 4.5))
    for rna_id in rna_ids:
        ax.plot(range(len(labels)), losses[rna_id], marker="o", markersize=3, label=rna_id)
    mean = [context.trial_results[config_id]["loss"] for config_id in context.config_ids]
    ax.plot(range(len(labels)), mean, color="black", linewidth=2.5, label="trial loss")
    ax.set_xticks(range(len(labels)), labels, fontsize=8)
    ax.set_xlabel("trajectory config ID")
    ax.set_ylabel("loss")
    ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(output_dir / "per_rna_lines.png", dpi=150)
    plt.close(fig)
