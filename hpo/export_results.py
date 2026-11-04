#!/usr/bin/env python3
"""Copy the first portable snapshot of a completed HPO run.

The export directory is always ``<config directory>/results``. Keep this
script intentionally small while the exported result contract is evolving.
"""
import argparse
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from hpo.pipeline import failures
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


def _warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def _export_errors(neps_dir: Path, destination: Path) -> dict[str, Any]:
    """Copy all error files and return the failure block for the manifest.

    The block is written also without failures, so it is clear where errors
    would be listed.
    """
    errors_dir = destination / "errors"
    files: list[str] = []
    reports: dict[str, dict[str, Any]] = {}
    stopped: set[str] = set()

    for path in sorted((neps_dir / "configs").glob(
        f"config_*/artifacts/{failures.FAILURE_REPORT}"
    )):
        config_name = path.parent.parent.name
        target = errors_dir / "reports" / f"{config_name}.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        files.append(target.relative_to(destination).as_posix())
        reports[config_name[len("config_"):]] = json.loads(path.read_text())

    # NePS logs every crash; it also covers trials without a failure report
    # (runs from before the reports existed, or crashes outside the pipeline).
    shared_errors = neps_dir / "shared_errors.jsonl"
    if shared_errors.is_file():
        target = errors_dir / shared_errors.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(shared_errors, target)
        files.append(target.relative_to(destination).as_posix())
        records, skipped = failures.read_jsonl(shared_errors)
        if skipped:
            _warn(f"Skipped {skipped} non-JSON line(s) in {shared_errors}")
        for record in records:
            trial_id = str(record.get("trial_id"))
            if trial_id in reports:
                continue
            # Trials stopped before they ran are not failures of their config.
            if failures.is_stop(str(record.get("err", ""))):
                stopped.add(trial_id)
                continue
            signal = re.search(r"Signals\.(SIG\w+)", str(record.get("err", "")))
            reports[trial_id] = dict(
                error_type=record.get("err_type"),
                signal=signal.group(1) if signal else None,
            )

    finished = sum(
        1 for path in (neps_dir / "configs").glob("config_*/report.yaml")
        if path.parent.name[len("config_"):] not in stopped
    )
    block = failures.summarize(reports, finished, files)
    if block["count"]:
        _warn(
            failures.summary_line(block)
            + f". See {errors_dir.relative_to(destination.parent)}/"
        )
    return block


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
        "failures": _export_errors(neps_dir, destination),
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
