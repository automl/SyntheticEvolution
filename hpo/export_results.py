#!/usr/bin/env python3
"""Copy the first portable snapshot of a completed HPO run or check run.

The export directory is always ``<config directory>/results``. Keep this
script intentionally small while the exported result contract is evolving.
"""
import argparse
import csv
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from hpo.pipeline import archive, failures, retention
from hpo.pipeline.run import get_workspace_path, load_config, shs_generation_request, storage_settings

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


def _export_shs(
    source_dir: Path,
    target_dir: Path,
    result: dict[str, Any],
    run_metadata: dict[str, Any],
    include_archives: bool,
) -> list[Path]:
    """Write ``shs/rna_<index>/{af3_input.json, shs_generation_request.json}``.

    Both are derived instead of stored: AF3's ``data.json`` in the RNA's
    archive is the AF3 input with the synthetic MSA, and the request is
    rebuilt from the dataset in ``run_metadata.json`` and the trial's seeds
    (results from before the seeds were recorded: the run's seeds).
    """
    archives = sorted(source_dir.glob(f"rna_*{archive.SUFFIX}"))
    if not archives:
        return _copy_legacy_shs(source_dir, target_dir)
    rows = run_metadata["dataset"]
    if len(archives) != len(rows):
        raise FileNotFoundError(
            f"Expected {len(rows)} RNA archives in {source_dir}, found {len(archives)}"
        )
    written = []
    for path in archives:
        task_name = path.name[: -len(archive.SUFFIX)]
        index = int(task_name[len("rna_"):])
        task_target = target_dir / "shs" / task_name
        task_target.mkdir(parents=True, exist_ok=True)
        (task_target / "af3_input.json").write_bytes(archive.read_data_json(path))
        seeds = [result.get(key, run_metadata["config"][key]) for key in ("shs_seed", "af3_seed")]
        request = shs_generation_request(rows[index], index, result["parameters"], *seeds)
        (task_target / "shs_generation_request.json").write_text(json.dumps(request, indent=2) + "\n")
        written += [task_target / "af3_input.json", task_target / "shs_generation_request.json"]
        if include_archives:
            copy = target_dir / "archives" / path.name
            copy.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, copy)
            written.append(copy)
    return written


def _copy_legacy_shs(source_dir: Path, target_dir: Path) -> list[Path]:
    """Runs from before the per-RNA archives kept the files in rna_<index>/."""
    written = []
    for task_dir in sorted(source_dir.glob("rna_*")):
        if not task_dir.is_dir():
            continue
        for filename in ("af3_input.json", "shs_generation_request.json"):
            artifact = task_dir / filename
            if not artifact.is_file():
                raise FileNotFoundError(f"Expected SHS artifact does not exist: {artifact}")
            artifact_target = target_dir / "shs" / task_dir.name / filename
            artifact_target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(artifact, artifact_target)
            written.append(artifact_target)
    return written


def _export_trials(
    config: dict[str, Any],
    workspace: Path,
    destination: Path,
    run_metadata: dict[str, Any],
    include_archives: bool,
) -> tuple[list[str], dict[str, Any]]:
    """Standalone run with ``trials``: each trial's result and SHS files.

    Returns the copied files and the failure block. The trial list comes from
    the given config (trials can be added after the run started); each result
    must match its trial's parameters and seeds.
    """
    copied: list[str] = []
    reports: dict[str, dict[str, Any]] = {}
    report_files: list[str] = []
    finished = 0
    for trial in config["trials"]:
        name = trial["name"]
        source_dir = workspace / "trials" / name
        source = source_dir / "trial_result.json"
        report = source_dir / failures.FAILURE_REPORT
        if source.is_file():
            result = json.loads(source.read_text())
            expected = dict(
                parameters={**config.get("fixed", {}), **trial.get("parameters", {})},
                **{key: trial.get(key, config[key]) for key in ("shs_seed", "af3_seed")},
            )
            if {key: result.get(key) for key in expected} != expected:
                raise ValueError(f"Trial {name} in {source_dir} does not match its entry in the config")
            target = destination / "trials" / name / "trial_result.json"
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            written = [target] + _export_shs(source_dir, target.parent, result, run_metadata, include_archives)
            copied += [path.relative_to(destination).as_posix() for path in written]
        elif report.is_file():
            target = destination / "errors" / "reports" / f"{name}.json"
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(report, target)
            report_files.append(target.relative_to(destination).as_posix())
            reports[name] = json.loads(report.read_text())
        else:
            _warn(f"Trial {name} has not finished; it is not exported")
            continue
        finished += 1

    block = failures.summarize(reports, finished, report_files)
    if block["count"]:
        _warn(failures.summary_line(block) + f". See {(destination / 'errors').relative_to(destination.parent)}/")
    return copied, block


def _trajectory_ids(summary_dir: Path) -> list[str]:
    return re.findall(
        r"^Config ID:\s*(\d+)\s*$",
        (summary_dir / "best_config_trajectory.txt").read_text(),
        re.MULTILINE,
    )


def prune_run(config_path: Path, dry_run: bool = False) -> dict[str, int]:
    """Keep archives only for the trajectory and the top k of a finished run.

    Applies the ``trajectory`` retention exactly, also to a run that used
    ``all``. Uses NePS's summary (``full.csv``, ``best_config_trajectory.txt``),
    so run it after the run has finished. Failed trials are never pruned.
    """
    config = load_config(config_path.expanduser().resolve(), mode="export")
    if "trials" in config:
        raise SystemExit("--prune applies to NePS runs only")
    top_k = storage_settings(config)["keep_top_k"]
    neps_dir = get_workspace_path(config) / "neps"
    with (neps_dir / "summary" / "full.csv").open(newline="") as stream:
        trials = {
            row["id"]: (float(row["time_sampled"]), float(row["objective_to_minimize"]))
            for row in csv.DictReader(stream)
            if row["state"] == "success"
        }
    kept = retention.kept_trials(trials, top_k)
    # NePS's own trajectory is never pruned, even if the rule above disagreed.
    trajectory = set(_trajectory_ids(neps_dir / "summary"))
    if trajectory - kept.keys():
        _warn(f"NePS trajectory configs kept although the retention rule differs: {sorted(trajectory - kept.keys())}")

    totals = dict(configs=0, files=0, bytes=0)
    for trial_id in sorted(trials.keys() - kept.keys() - trajectory):
        artifacts = neps_dir / "configs" / f"config_{trial_id}" / "artifacts"
        removed = retention.prune(artifacts, dry_run=dry_run)
        if removed["files"]:
            totals["configs"] += 1
            totals["files"] += removed["files"]
            totals["bytes"] += removed["bytes"]
            if not dry_run:
                retention.record(artifacts, None, pruned_by="export_results --prune")
    action = "Would remove" if dry_run else "Removed"
    print(
        f"{action} {totals['files']} archives ({totals['bytes'] / 1e6:.1f} MB) of "
        f"{totals['configs']} configs; kept {len(kept.keys() | trajectory)} of {len(trials)} "
        f"successful configs (trajectory and top {top_k})."
    )
    return totals


def export_results(config_path: Path, include_archives: bool = False) -> Path:
    """Copy the selected run summaries and trajectory trial results.

    ``include_archives`` also copies the complete per-RNA archives (AF3
    samples, confidences, DSSR output) of the exported configs.
    """
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

    run_metadata = workspace / "run_metadata.json"
    if not run_metadata.is_file():
        raise FileNotFoundError(f"Expected run metadata does not exist: {run_metadata}")
    shutil.copy2(run_metadata, destination / run_metadata.name)
    copied: list[str] = [run_metadata.name]
    run_metadata_data = json.loads(run_metadata.read_text())
    manifest = {
        "config": str(config_path),
        "workspace": str(workspace),
        "exported_at": datetime.now(timezone.utc).isoformat(),
    }

    # Check run (standalone with `trials`): results under trials/<name>/.
    if "trials" in config:
        files, manifest["failures"] = _export_trials(
            config, workspace, destination, run_metadata_data, include_archives,
        )
        manifest.update(trials=config["trials"], files=copied + files)
        (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        return destination

    neps_dir = workspace / "neps"
    summary_dir = neps_dir / "summary"
    summary_files = (
        "best_config_trajectory.txt",
        "best_config.txt",
        "short.csv",
        "full.csv",
    )

    for filename in summary_files:
        source = summary_dir / filename
        if not source.is_file():
            raise FileNotFoundError(f"Expected summary file does not exist: {source}")
        target = destination / "summary" / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        copied.append(target.relative_to(destination).as_posix())

    config_ids = [int(match) for match in _trajectory_ids(summary_dir)]
    if not config_ids:
        # No trial succeeded; the errors are then the most important export.
        _warn(
            f"No Config IDs found in {summary_dir / 'best_config_trajectory.txt'}: "
            "no trial succeeded. Exporting summaries and errors only."
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

        written = _export_shs(
            source_dir, target.parent, json.loads(source.read_text()),
            run_metadata_data, include_archives,
        )
        copied += [path.relative_to(destination).as_posix() for path in written]

    manifest.update(
        config_ids=config_ids,
        failures=_export_errors(neps_dir, destination),
        files=copied,
    )
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
        help="YAML configuration file identifying the completed HPO run or check run",
    )
    parser.add_argument(
        "--include-archives",
        action="store_true",
        help="Also copy the per-RNA archives (AF3 samples, confidences, DSSR output)",
    )
    parser.add_argument(
        "--prune",
        action="store_true",
        help="Instead of exporting: delete the archives of configs that are neither on "
             "the trajectory nor in the top storage.keep_top_k (after the run has finished)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="With --prune: only report what would be removed",
    )
    args = parser.parse_args()
    if args.dry_run and not args.prune:
        parser.error("--dry-run requires --prune")
    if args.prune:
        prune_run(args.config, dry_run=args.dry_run)
        return

    destination = export_results(args.config, include_archives=args.include_archives)
    print(destination)


if __name__ == "__main__":
    main()
