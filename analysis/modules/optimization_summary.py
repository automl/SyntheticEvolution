"""Write a small summary of the exported optimization trajectory."""

import json
from pathlib import Path
from typing import Any

from analysis.results_loader import ResultsContext


def run(context: ResultsContext, output_dir: Path) -> None:
    """Write trajectory losses and the NePS best-configuration text."""
    trajectory = [
        {
            "config_id": config_id,
            "loss": context.trial_results[config_id].get("loss"),
        }
        for config_id in context.config_ids
    ]
    report: dict[str, Any] = {
        "config_ids": context.config_ids,
        "trajectory": trajectory,
        "best_config": context.best_config_text,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "optimization_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
