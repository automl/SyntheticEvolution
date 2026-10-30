"""Load and validate an exported HPO results directory."""

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any


_CONFIG_ID_PATTERN = re.compile(r"^Config ID:\s*(\d+)\s*$", re.MULTILINE)
_SUMMARY_FILES = (
    "best_config_trajectory.txt",
    "best_config.txt",
    "short.csv",
)


@dataclass(frozen=True)
class ResultsContext:
    """Shared read-only input for analysis modules."""

    results_dir: Path
    output_dir: Path
    run_metadata: dict[str, Any]
    manifest: dict[str, Any]
    config_ids: list[int]
    trial_results: dict[int, dict[str, Any]]
    best_config_text: str
    trajectory_text: str
    short_csv_text: str


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read JSON object from {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def _required_file(results_dir: Path, relative_path: str) -> Path:
    path = results_dir / relative_path
    if not path.is_file():
        raise FileNotFoundError(f"Required analysis input is missing: {path}")
    return path


def _parse_config_ids(trajectory_text: str, trajectory_path: Path) -> list[int]:
    config_ids = [int(value) for value in _CONFIG_ID_PATTERN.findall(trajectory_text)]
    if not config_ids:
        raise ValueError(f"No Config IDs found in {trajectory_path}")
    if len(config_ids) != len(set(config_ids)):
        raise ValueError(f"Duplicate Config ID in {trajectory_path}")
    return config_ids


def load_results(results_dir: Path, output_dir: Path) -> ResultsContext:
    """Load the exported result contract without accessing the HPO workspace."""
    results_dir = results_dir.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    if not results_dir.is_dir():
        raise NotADirectoryError(f"Results directory does not exist: {results_dir}")

    run_metadata_path = _required_file(results_dir, "run_metadata.json")
    manifest_path = _required_file(results_dir, "manifest.json")
    trajectory_path = _required_file(
        results_dir, "summary/best_config_trajectory.txt"
    )
    best_config_path = _required_file(results_dir, "summary/best_config.txt")
    short_csv_path = _required_file(results_dir, "summary/short.csv")

    trajectory_text = trajectory_path.read_text()
    config_ids = _parse_config_ids(trajectory_text, trajectory_path)
    manifest = _read_json(manifest_path)
    manifest_config_ids = manifest.get("config_ids")
    if manifest_config_ids is not None and manifest_config_ids != config_ids:
        raise ValueError(
            "Manifest configuration IDs do not match the best-config trajectory"
        )

    trial_results: dict[int, dict[str, Any]] = {}
    for config_id in config_ids:
        relative_path = f"configs/config_{config_id}/trial_result.json"
        trial_path = _required_file(results_dir, relative_path)
        trial_results[config_id] = _read_json(trial_path)

    return ResultsContext(
        results_dir=results_dir,
        output_dir=output_dir,
        run_metadata=_read_json(run_metadata_path),
        manifest=manifest,
        config_ids=config_ids,
        trial_results=trial_results,
        best_config_text=best_config_path.read_text(),
        trajectory_text=trajectory_text,
        short_csv_text=short_csv_path.read_text(),
    )
