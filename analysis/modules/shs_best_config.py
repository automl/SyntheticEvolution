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


def _load_module(name: str, path: Path) -> Any:
    if str(path.parent) not in sys.path:
        sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load SHS module from {path}")
    module = importlib.util.module_from_spec(spec)
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


def _plot_config(config_id: int, context: ResultsContext, output_dir: Path, analyzer: Any) -> int:
    config_dir = context.results_dir / "configs" / f"config_{config_id}" / "shs"
    task_dirs = sorted(path for path in config_dir.glob("rna_*") if path.is_dir())
    if not task_dirs:
        raise FileNotFoundError(f"No exported SHS tasks found for Config ID {config_id}")

    config_output = output_dir / f"config_{config_id}"
    plot_count = 0
    for task_dir in task_dirs:
        af3_input_path = task_dir / "af3_input.json"
        request_path = task_dir / "shs_generation_request.json"
        af3_input = json.loads(af3_input_path.read_text())
        request = json.loads(request_path.read_text())
        pair_map_module = _load_module(
            "analysis_pair_map",
            Path(__file__).parents[2] / "SHS-Generator-RNA" / "pair_map.py",
        )
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
        analyzer.plot_covariance(
            data.seq, features.covariance, show=False,
            out=str(task_output / "covariance.png"),
        )
        analyzer.plot_covariance_classification(
            data.seq, features.covariance, pairs=data.pairs, show=False,
            out=str(task_output / "covariance_classification.png"),
        )
        analyzer.plot_deletions_and_insertion(
            features.col_deletion_rate,
            features.col_insertion_mean,
            pairs=data.pairs,
            show=False,
            out=str(task_output / "deletion_insertion.png"),
        )
        plt.close("all")
        plot_count += 3
    return plot_count


def _write_msa_views(data: Any, output_dir: Path) -> None:
    """Write clean and insertion-preserving FASTA views of one chain's MSA."""
    rows = [data.seq, *("".join(row) for row in data.raw)]
    clean_rows = ["".join(char for char in row if not char.islower()) for row in rows]
    insertion_rows = [analyzer_row for analyzer_row in rows]
    insertion_widths = [
        max(len(chunks[index]) for chunks in map(_insertion_chunks, insertion_rows))
        + 1
        for index in range(len(data.seq) + 1)
    ]

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_fasta(
        output_dir / "clean.fasta",
        clean_rows,
    )
    _write_fasta(
        output_dir / "aligned.fasta",
        [_aligned_row(row, insertion_widths) for row in insertion_rows],
    )
    LOGGER.info(
        "Wrote MSA FASTA views to %s (maximum column width: %d)",
        output_dir,
        max(insertion_widths),
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


def _aligned_row(row: str, insertion_widths: list[int]) -> str:
    chunks = _insertion_chunks(row)
    aligned = [char for char in row if not char.islower()]
    parts = [
        chunks[index]
        + char
        + " " * (insertion_widths[index] - len(chunks[index]) - 1)
        for index, char in enumerate(aligned)
    ]
    parts.append(chunks[-1])
    return "".join(parts)


def _write_fasta(path: Path, rows: list[str]) -> None:
    path.write_text("\n".join(rows) + "\n")


def run(context: ResultsContext, output_dir: Path) -> None:
    """Generate covariance and indel plots for the configured SHS selections."""
    generator_dir = Path(__file__).parents[2] / "SHS-Generator-RNA"
    if str(generator_dir) not in sys.path:
        sys.path.insert(0, str(generator_dir))
    analyzer_path = Path(__file__).parents[1] / "SHS-Analyzer" / "shs_analyzer.py"
    analyzer = _load_module("analysis_shs_analyzer", analyzer_path)
    selected_ids = _selected_config_ids(context)
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_count = sum(
        _plot_config(config_id, context, output_dir, analyzer)
        for config_id in selected_ids
    )
    (output_dir / "summary.json").write_text(
        json.dumps({"config_ids": selected_ids, "plot_count": plot_count}, indent=2)
        + "\n"
    )
