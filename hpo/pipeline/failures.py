"""Failure reports for trials and failure summaries for runs."""

import json
import os
import signal
import socket
import subprocess
import traceback
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Optional

from hpo.pipeline.run import write_json
from hpo.pipeline.signals import ControllerStop

FAILURE_REPORT = 'failure_report.json'
LOG_NAMES = ('shs_generation.log', 'af3.log', 'dssr.log')
LOG_TAIL_BYTES = 64 * 1024


def _relative(value: str, trial_directory: Path, scratch: Optional[Path] = None) -> str:
    """Shorten trial paths; absolute paths belong in run_metadata.json.

    Paths in the trial directory become relative, paths in the scratch
    directory (deleted after the trial) start with ``<scratch>/``.
    """
    replacements = {str(trial_directory) + os.sep: ''}
    if scratch is not None and scratch != trial_directory:
        replacements[str(scratch) + os.sep] = '<scratch>' + os.sep
    # Longest prefix first, in case one directory lies inside the other.
    for prefix in sorted(replacements, key=len, reverse=True):
        value = value.replace(prefix, replacements[prefix])
    return value


def _log_tail(path: Path) -> str:
    """Return at most the last LOG_TAIL_BYTES of a log; errors are at the end."""
    size = path.stat().st_size
    with path.open('rb') as stream:
        stream.seek(max(0, size - LOG_TAIL_BYTES))
        text = stream.read().decode('utf-8', errors='replace')
    if size > LOG_TAIL_BYTES:
        text = f'[... first {size - LOG_TAIL_BYTES} bytes truncated ...]\n' + text
    return text


def write_failure_report(
    trial_directory: Path,
    rows: list[dict[str, Any]],
    progress: dict[str, Any],
    error: BaseException,
    log_directory: Optional[Path] = None,
) -> dict[str, Any]:
    """Write everything needed to debug a failed trial into one file.

    Logs are read from ``log_directory/rna_<index>/`` (the trial's scratch
    directory; default: the trial directory) and embedded in the report.
    """
    if log_directory is None:
        log_directory = trial_directory
    rna_index = progress.get('rna_index')
    report: dict[str, Any] = dict(
        stage=progress.get('stage'),
        rna_index=rna_index,
        rna_id=rows[rna_index]['id'] if rna_index is not None else None,
        error_type=type(error).__name__,
        message=_relative(str(error), trial_directory, log_directory),
        traceback=''.join(
            traceback.format_exception(type(error), error, error.__traceback__)
        ),
        # A tool's own timeout, or Slurm's time limit / cancellation (ControllerStop).
        timed_out=isinstance(error, (subprocess.TimeoutExpired, ControllerStop)),
        command=None,
        returncode=None,
        signal=None,
        host=socket.gethostname(),
        slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        slurm_array_task_id=os.environ.get('SLURM_ARRAY_TASK_ID'),
        time=datetime.now(timezone.utc).isoformat(),
        logs={},
    )

    # CalledProcessError and TimeoutExpired carry the failed command.
    command = getattr(error, 'cmd', None)
    if command is not None:
        if isinstance(command, (list, tuple)):
            report['command'] = [_relative(str(arg), trial_directory, log_directory) for arg in command]
        else:
            report['command'] = _relative(str(command), trial_directory, log_directory)
    if isinstance(error, subprocess.CalledProcessError):
        report['returncode'] = error.returncode
        if error.returncode < 0:
            try:
                report['signal'] = signal.Signals(-error.returncode).name
            except ValueError:
                pass

    if rna_index is not None:
        task_dir = log_directory / f'rna_{rna_index:05d}'
        for name in LOG_NAMES:
            if (task_dir / name).is_file():
                report['logs'][name] = _log_tail(task_dir / name)

    write_json(trial_directory / FAILURE_REPORT, report)
    return report


def describe(report: dict[str, Any]) -> str:
    """One-line description of a failure, used in logs and export warnings."""
    cause = report.get('signal') or report.get('error_type') or 'unknown error'
    where = report.get('stage') or 'unknown stage'
    rna = report.get('rna_id')
    return f'{cause} in {where}' + (f' on {rna}' if rna is not None else '')


def read_jsonl(path: Path) -> tuple[list[dict[str, Any]], int]:
    """Read JSON lines, skipping lines that are not JSON objects.

    Returns the records and the number of skipped lines.
    """
    records, skipped = [], 0
    for line in path.read_text(errors='replace').splitlines():
        if not line.strip():
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            skipped += 1
            continue
        if isinstance(record, dict):
            records.append(record)
        else:
            skipped += 1
    return records, skipped


def summarize(
    reports: dict[str, dict[str, Any]],
    finished: Optional[int],
    files: Iterable[str] = (),
) -> dict[str, Any]:
    """Build the failure block for summaries; always present, also without failures.

    ``reports`` maps config IDs to failure reports (or to records with the same
    keys). ``finished`` is the number of finished trials, if known.
    """
    groups = Counter(
        (
            report.get('error_type'),
            report.get('signal'),
            report.get('stage'),
            report.get('rna_id'),
        )
        for report in reports.values()
    )
    count = len(reports)
    return dict(
        status='found' if count else 'none',
        count=count,
        finished=finished,
        rate=round(count / finished, 4) if finished else (None if count else 0.0),
        groups=[
            dict(
                error_type=error_type,
                signal=signal_name,
                stage=stage,
                rna_id=rna_id,
                count=group_count,
                description=describe(
                    dict(error_type=error_type, signal=signal_name, stage=stage, rna_id=rna_id)
                ),
            )
            for (error_type, signal_name, stage, rna_id), group_count in groups.most_common()
        ],
        files=sorted(files),
    )


def summary_line(block: dict[str, Any]) -> str:
    """Human-readable form of a failure block."""
    if block['status'] == 'none':
        return 'No failed trials.'
    rate = f" ({block['rate']:.1%})" if block['rate'] is not None else ''
    of = f" of {block['finished']}" if block['finished'] is not None else ''
    groups = '; '.join(f"{g['count']}x {g['description']}" for g in block['groups'])
    return f"{block['count']}{of} trials failed{rate}: {groups}"


# ---- Failure limit (NePS runs) ---------------------------------------------

STOPPED = 'STOPPED.json'
STOPPED_MARKER = '[stopped]'
FAILURE_LIMIT_EXIT_CODE = 3


def stop_time(run_root: Path) -> Optional[float]:
    """Time of the newest renamed stop file (STOPPED-<n>.json), if any.

    Renaming STOPPED.json restarts the run; the failure limit then counts only
    trials that finished after that stop.
    """
    times = []
    for path in run_root.glob('STOPPED-*.json'):
        try:
            times.append(float(json.loads(path.read_text())['time']))
        except (OSError, KeyError, TypeError, ValueError):
            raise ValueError(f'Cannot read the stop time from {path}') from None
    return max(times, default=None)


def check_not_stopped(run_root: Path) -> None:
    """Refuse to start a run that the failure limit has stopped."""
    path = run_root / STOPPED
    if not path.exists():
        return
    number = 1
    while (run_root / f'STOPPED-{number}.json').exists():
        number += 1
    try:
        reason = json.loads(path.read_text()).get('reason', 'unknown reason')
    except (OSError, ValueError):
        reason = 'unreadable stop file'
    raise RuntimeError(
        f'The run was stopped by the failure limit ({reason}). Inspect the failures, '
        f'then rename the stop file to restart:\n  mv {path} {run_root / f"STOPPED-{number}.json"}\n'
        'After the restart, the failure limit counts only trials that finish later.'
    )


def write_stop(run_root: Path, info: dict[str, Any]) -> bool:
    """Create STOPPED.json; return False if another controller created it first."""
    path = run_root / STOPPED
    temporary = run_root / f'.{STOPPED}.{os.getpid()}.{socket.gethostname()}.tmp'
    temporary.write_text(json.dumps(info, indent=2, sort_keys=True))
    try:
        # A hard link fails if the target exists: exclusive and never partial.
        os.link(temporary, path)
        return True
    except FileExistsError:
        return False
    finally:
        temporary.unlink(missing_ok=True)


def is_stop(message: Optional[str]) -> bool:
    """Whether a NePS error message belongs to a trial stopped before it ran."""
    return bool(message) and message.startswith(STOPPED_MARKER)


def failure_limit_exceeded(
    outcomes: list[tuple[float, bool]],
    settings: dict[str, Any],
) -> Optional[dict[str, Any]]:
    """Check the failure limit.

    ``outcomes`` holds ``(time_end, failed)`` of the finished trials that count
    (without stopped trials, and only those after the newest restart).
    Returns the stop information, or None while the limit holds.
    """
    if len(outcomes) < settings['min_trials']:
        return None
    ordered = [failed for _, failed in sorted(outcomes)]
    window = ordered[-settings['window']:] if settings['window'] else ordered
    rate = sum(window) / len(window)
    if rate <= settings['max_rate']:
        return None
    scope = f"last {len(window)}" if settings['window'] else 'all'
    return dict(
        reason=f"{rate:.1%} of the {scope} finished trials failed (max_rate {settings['max_rate']:.0%})",
        rate=round(rate, 4),
        failed=sum(window),
        counted=len(window),
        finished=len(ordered),
        settings=settings,
    )
