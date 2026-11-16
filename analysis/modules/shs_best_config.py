"""Generate static SHS plots for selected exported configurations."""

import importlib.util
import json
import logging
import re
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

from analysis.results_loader import ResultsContext

LOGGER = logging.getLogger("analysis.shs_best_config")


ROOT = Path(__file__).resolve().parents[2]
FOLDERS = {"rna": "RNA", "protein": "Protein"}


def _load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load SHS module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _best_config_id(best_config_text: str) -> int:
    match = re.search(r"^Config ID:\s*(\d+)\s*$", best_config_text, re.MULTILINE)
    if match is None:
        raise ValueError("Could not find the best Config ID in best_config.txt")
    return int(match.group(1))


def _selected_config_ids(context: ResultsContext) -> list[int]:
    config = context.analysis_config.get("shs", {})
    if not isinstance(config, dict):
        raise ValueError("Analysis config field 'shs' must be an object")
    top_n = config.get("top_n", 1)
    if not isinstance(top_n, int) or isinstance(top_n, bool) or top_n < 1:
        raise ValueError("Analysis config field 'shs.top_n' must be a positive integer")
    if top_n == 1:
        return [_best_config_id(context.best_config_text)]
    return sorted(
        context.config_ids,
        key=lambda config_id: context.trial_results[config_id].get("loss", float("inf")),
    )[:top_n]


def _shs_modules(chain_type: str) -> tuple[Any, Any]:
    """Load the pair_map and SHS analyzer for one chain type."""
    folder = FOLDERS[chain_type]
    pair_map = _load_module(f"pair_map_{chain_type}", ROOT / f"SHS-Generator-{folder}" / "pair_map.py")
    # Both analyzers run `from pair_map import ...`; give each its own type's module.
    sys.modules["pair_map"] = pair_map
    analyzer = _load_module(
        f"shs_analyzer_{chain_type}", ROOT / "analysis" / f"SHS-Analyzer-{folder}" / "shs_analyzer.py"
    )
    return pair_map, analyzer


def _plot_config(config_id: int, context: ResultsContext, output_dir: Path, modules: dict[str, Any]) -> int:
    config_dir = context.results_dir / "configs" / f"config_{config_id}" / "shs"
    task_dirs = sorted(path for path in config_dir.glob("rna_*") if path.is_dir())
    if not task_dirs:
        raise FileNotFoundError(f"No exported SHS tasks found for Config ID {config_id}")

    dataset = context.run_metadata["dataset"]
    plot_count = 0
    for task_dir in task_dirs:
        # Task rna_<index> holds dataset row <index>, for every chain type.
        chain_type = dataset[int(task_dir.name.split("_")[1])].get("type", "rna")
        if chain_type not in modules:
            modules[chain_type] = _shs_modules(chain_type)
        pair_map_module, analyzer = modules[chain_type]
        af3_input_path = task_dir / "af3_input.json"
        request_path = task_dir / "shs_generation_request.json"
        af3_input = json.loads(af3_input_path.read_text())
        request = json.loads(request_path.read_text())
        pair_map = pair_map_module.build_pair_map(
            request["sequence"], request["pairs"]
        )
        data = analyzer.MsaData.from_text(
            json.dumps(af3_input),
            pairs=pair_map,
        )
        _write_msa_views(data, task_output := output_dir / f"config_{config_id}" / task_dir.name)
        features = analyzer.Features(data)
        task_output.mkdir(parents=True, exist_ok=True)
        # Covariance needs at least two sample rows.
        if len(data.aligned) < 2:
            LOGGER.warning(
                "Skipping covariance plots for Config ID %d, %s: MSA has only %d sample",
                config_id, task_dir.name, len(data.aligned),
            )
        else:
            analyzer.plot_covariance(
                data.seq, features.covariance, show=False,
                out=str(task_output / "covariance.png"),
            )
            analyzer.plot_covariance_classification(
                data.seq, features.covariance, pairs=data.pairs, show=False,
                out=str(task_output / "covariance_classification.png"),
            )
            plot_count += 2
        analyzer.plot_deletions_and_insertion(
            features.col_deletion_rate,
            features.col_insertion_mean,
            pairs=data.pairs,
            show=False,
            out=str(task_output / "deletion_insertion.png"),
        )
        plt.close("all")
        plot_count += 1
    return plot_count


def _write_msa_views(data: Any, output_dir: Path) -> None:
    """Write clean and insertion-preserving FASTA views of one chain's MSA."""
    rows = [data.seq, *("".join(row) for row in data.raw)]
    clean_rows = ["".join(char for char in row if not char.islower()) for row in rows]
    insertion_rows = [analyzer_row for analyzer_row in rows]
    column_width = (
        max(len(chunk) for row in insertion_rows for chunk in _insertion_chunks(row))
        + 1
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_fasta(
        output_dir / "clean.fasta",
        clean_rows,
    )
    _write_fasta(
        output_dir / "aligned.fasta",
        [_aligned_row(row, column_width) for row in insertion_rows],
    )
    LOGGER.info(
        "Wrote MSA FASTA views to %s (column width: %d)",
        output_dir,
        column_width,
    )


def _insertion_chunks(row: str) -> list[str]:
    chunks: list[str] = []
    current: list[str] = []
    for char in row:
        if char.islower():
            current.append(char)
        else:
            chunks.append("".join(current))
            current = []
    chunks.append("".join(current))
    return chunks


def _aligned_row(row: str, column_width: int) -> str:
    """Keep only insertions, each padded to a fixed-width cell per column."""
    chunks = _insertion_chunks(row)
    parts = [chunk.ljust(column_width) for chunk in chunks[:-1]]
    parts.append(chunks[-1])
    return "".join(parts)


def _write_fasta(path: Path, rows: list[str]) -> None:
    path.write_text("\n".join(rows) + "\n")


def run(context: ResultsContext, output_dir: Path) -> None:
    """Generate covariance and indel plots for the configured SHS selections."""
    modules: dict[str, Any] = {}
    selected_ids = _selected_config_ids(context)
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_count = sum(
        _plot_config(config_id, context, output_dir, modules)
        for config_id in selected_ids
    )
    (output_dir / "summary.json").write_text(
        json.dumps({"config_ids": selected_ids, "plot_count": plot_count}, indent=2)
        + "\n"
    )
