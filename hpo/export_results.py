#!/usr/bin/env python3
"""Copy the first portable snapshot of a completed HPO run.

The export directory is always ``<config directory>/results``. Keep this
script intentionally small while the exported result contract is evolving.
"""
import argparse
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

from hpo.pipeline.run import get_workspace_path, load_config

EXPORT_DIR_NAME = "results"


def _confirm_recreate(destination: Path) -> None:
    answer = input(
        f"Export results already exists at {destination}. "
        "Recreate it? [Y/n] "
    ).strip().lower()
    if answer not in {"", "y", "yes"}:
        raise SystemExit("Export cancelled.")
    shutil.rmtree(destination)


def export_results(config_path: Path) -> Path:
    """Copy the selected run summaries and trajectory trial results."""
    config_path = config_path.expanduser().resolve()
    config = load_config(config_path, mode="export")
    workspace = get_workspace_path(config)
    destination = config_path.parent / EXPORT_DIR_NAME
    if destination.exists():
        if not destination.is_dir():
            raise NotADirectoryError(
                f"Export output is not a directory: {destination}"
            )
        _confirm_recreate(destination)
    destination.mkdir(parents=True, exist_ok=True)

    neps_dir = workspace / "neps"
    summary_dir = neps_dir / "summary"
    summary_files = (
        "best_config_trajectory.txt",
        "best_config.txt",
        "short.csv",
    )
    copied: list[str] = []

    for filename in summary_files:
        source = summary_dir / filename
        if not source.is_file():
            raise FileNotFoundError(f"Expected summary file does not exist: {source}")
        target = destination / "summary" / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        copied.append(target.relative_to(destination).as_posix())

    run_metadata = workspace / "run_metadata.json"
    if not run_metadata.is_file():
        raise FileNotFoundError(f"Expected run metadata does not exist: {run_metadata}")
    shutil.copy2(run_metadata, destination / run_metadata.name)
    copied.append(run_metadata.name)

    config_ids = [
        int(match)
        for match in re.findall(
            r"^Config ID:\s*(\d+)\s*$",
            (summary_dir / "best_config_trajectory.txt").read_text(),
            re.MULTILINE,
        )
    ]
    if not config_ids:
        raise ValueError(
            f"No Config IDs found in {summary_dir / 'best_config_trajectory.txt'}"
        )

    for config_id in config_ids:
        source_dir = neps_dir / "configs" / f"config_{config_id}" / "artifacts"
        source = source_dir / "trial_result.json"
        if not source.is_file():
            raise FileNotFoundError(f"Expected trial result does not exist: {source}")
        target = destination / "configs" / f"config_{config_id}" / "trial_result.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        copied.append(target.relative_to(destination).as_posix())

        for task_dir in sorted(source_dir.glob("rna_*")):
            if not task_dir.is_dir():
                continue
            for filename in ("af3_input.json", "shs_generation_request.json"):
                artifact = task_dir / filename
                if not artifact.is_file():
                    raise FileNotFoundError(f"Expected SHS artifact does not exist: {artifact}")
                artifact_target = (
                    destination
                    / "configs"
                    / f"config_{config_id}"
                    / "shs"
                    / task_dir.name
                    / filename
                )
                artifact_target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(artifact, artifact_target)
                copied.append(artifact_target.relative_to(destination).as_posix())

    manifest = {
        "config": str(config_path),
        "workspace": str(workspace),
        "exported_at": datetime.now(timezone.utc).isoformat(),
        "config_ids": config_ids,
        "files": copied,
    }
    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Copy a compact result snapshot next to the selected config."
    )
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="YAML configuration file identifying the completed HPO run",
    )
    args = parser.parse_args()

    destination = export_results(args.config)
    print(destination)


if __name__ == "__main__":
    main()
