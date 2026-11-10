"""Run the NePS hyperparameter search for one shared RNA dataset."""
import argparse
import json
import os
import re
import socket
import subprocess
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any, Optional

from hpo.pipeline import failures, retention, signals
from hpo.pipeline.environment import NEPS_PACKAGE, NEPS_VERSION
from hpo.pipeline.run import load_config, prepare_run, storage_settings
from hpo.pipeline.trial import run_pipeline

import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)



def parse_slurm_duration(text: str) -> Optional[int]:
    """Seconds of a Slurm duration ([D-]HH:MM:SS or MM:SS); None for UNLIMITED etc."""
    match = re.fullmatch(r'(?:(\d+)-)?(?:(\d+):)?(\d+):(\d+)', text.strip())
    if match is None:
        return None
    days, hours, minutes, seconds = (int(group) if group else 0 for group in match.groups())
    return ((days * 24 + hours) * 60 + minutes) * 60 + seconds


def job_time_left() -> Optional[int]:
    """Seconds until this Slurm job's time limit; None outside Slurm or without a limit."""
    job = os.environ.get('SLURM_JOB_ID')
    if not job:
        return None
    try:
        result = subprocess.run(['squeue', '-h', '-j', job, '-o', '%L'],
                                capture_output=True, text=True, timeout=60, check=True)
    except (OSError, subprocess.SubprocessError):
        logger.warning('Could not read the time left of job %s; trials are not checked against it', job)
        return None
    return parse_slurm_duration(result.stdout)


def read_trials(neps_root: Path) -> dict:
    """All trials from NePS's own state (one cache file, read under NePS's lock)."""
    from neps.state import NePSState

    return NePSState.create_or_load(path=neps_root, load_only=True).lock_and_read_trials()


def error_message(trial) -> Optional[str]:
    report = trial.report
    return None if report is None or report.err is None else str(report.err)


def current_trial():
    """The NePS trial this controller is evaluating (id, metadata.time_sampled)."""
    from neps.runtime import get_in_progress_trial

    return get_in_progress_trial()


def apply_retention(neps_root: Path, trial_id: str, time_sampled: float, loss: float, top_k: int) -> None:
    """Prune the archives of configs that are no longer kept (``trajectory`` mode).

    Only configs that drop out because of this trial can change, so only
    those are touched: the trial itself, or configs it pushes off the
    trajectory or out of the top k.
    """
    from neps.state import State

    trials = {
        other_id: (trial.metadata.time_sampled, float(trial.report.objective_to_minimize))
        for other_id, trial in read_trials(neps_root).items()
        if trial.metadata.state == State.SUCCESS and trial.report is not None
        and trial.report.objective_to_minimize is not None
    }
    # The current trial is still `evaluating` in NePS's state.
    trials.pop(trial_id, None)
    before = retention.kept_trials(trials, top_k)
    trials[trial_id] = (time_sampled, loss)
    after = retention.kept_trials(trials, top_k)

    def artifacts(other_id: str) -> Path:
        return neps_root / 'configs' / f'config_{other_id}' / 'artifacts'

    for other_id in sorted((before.keys() | {trial_id}) - after.keys()):
        removed = retention.prune(artifacts(other_id))
        retention.record(artifacts(other_id), None, pruned_by=None if other_id == trial_id else trial_id)
        logger.info(
            'Pruned the archives of config_%s (%d files, %d bytes): not on the trajectory or in the top %d',
            other_id, removed['files'], removed['bytes'], top_k,
        )
    if trial_id in after:
        retention.record(artifacts(trial_id), after[trial_id])


def check_failure_limit(run_root: Path, neps_root: Path, settings: dict[str, Any]) -> Optional[dict]:
    """Return the stop information if the failure limit is exceeded.

    Trials stopped by the controller do not count, and after a restart
    (renamed stop file) only trials that finished later count.
    """
    from neps.state import State

    since = failures.stop_time(run_root)
    outcomes = []
    for trial in read_trials(neps_root).values():
        state = trial.metadata.state
        if state not in (State.SUCCESS, State.FAILED, State.CRASHED):
            continue
        if failures.is_stop(error_message(trial)):
            continue
        end = trial.metadata.time_end or 0.0
        if since is not None and end <= since:
            continue
        outcomes.append((end, state != State.SUCCESS))
    return failures.failure_limit_exceeded(outcomes, settings)


def stop_trial(trial_id: str, neps_root: Path, message: str, exit_code: int) -> None:
    """Report the current trial as crashed, then end the controller.

    NePS never recovers a trial left `evaluating`, so the trial is reported
    before exiting. SystemExit passes NePS's `except Exception` and ends
    `neps.run()`.
    """
    import neps

    neps.save_pipeline_results(
        # Without an explicit None, NePS stores 0, the best possible loss.
        {'exception': RuntimeError(message), 'objective_to_minimize': None},
        trial_id,
        neps_root,
    )
    log = logger.info if exit_code == 0 else logger.error
    log('Stopping the controller: %s (trial %s reported as crashed)', message, trial_id)
    raise SystemExit(exit_code)


def log_run_summary(neps_root: Path) -> None:
    """Log how many trials finished and failed, so failures cannot go unnoticed.

    Reads NePS's own trial state; stage and RNA come from the failure reports.
    """
    try:
        from neps.state import State

        trials = read_trials(neps_root)
        finished = {
            trial_id: trial for trial_id, trial in trials.items()
            if trial.metadata.state in (State.SUCCESS, State.FAILED, State.CRASHED)
        }
        # Trials stopped before they ran are not failures of their config.
        stopped = [
            trial_id for trial_id, trial in finished.items()
            if failures.is_stop(error_message(trial))
        ]
        for trial_id in stopped:
            del finished[trial_id]
        reports = {}
        for trial_id, trial in finished.items():
            if trial.metadata.state == State.SUCCESS:
                continue
            path = neps_root / 'configs' / f'config_{trial_id}' / 'artifacts' / failures.FAILURE_REPORT
            # NePS does not keep the exception type after reloading its state.
            reports[trial_id] = (
                json.loads(path.read_text()) if path.is_file()
                else dict(error_type='unknown (no failure report)')
            )
        block = failures.summarize(reports, len(finished))
    except Exception:
        logger.warning('Could not summarize the run in %s', neps_root, exc_info=True)
        return
    log = logger.error if block['count'] else logger.info
    line = failures.summary_line(block).rstrip('.') + '.'
    if stopped:
        line += f' {len(stopped)} trial(s) stopped by the controller.'
    log('Run summary: %s', line)


def main():
    """Load the search configuration and run the NePS controller."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--config',
        default='hpo/examples/helix.example.yaml',
        help='YAML configuration file for the dataset and NePS search',
    )
    args = parser.parse_args()

    # The search-space and neps.run API are tied to this tested package version.
    installed_neps_version = version(NEPS_PACKAGE)
    if installed_neps_version != NEPS_VERSION:
        raise RuntimeError(
            f'This implementation requires {NEPS_PACKAGE} {NEPS_VERSION}; '
            f'found {installed_neps_version}'
        )
    import neps

    config = load_config(args.config, mode='neps')

    run_root, dataset, fingerprint = prepare_run(config)
    try:
        failures.check_not_stopped(run_root)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    neps_root = run_root / "neps"
    stop_path = run_root / failures.STOPPED
    storage = storage_settings(config)
    signals.install()
    # Do not start a trial that cannot finish before the stop signal.
    time_left = job_time_left()
    deadline = None if time_left is None else time.monotonic() + time_left - signals.STOP_MARGIN_SECONDS
    trial_durations = []

    # Convert YAML search definitions into NePS parameter objects.
    space = {}
    for name, raw in config['generators'][config['optimize']]['search'].items():
        spec = dict(raw)
        kind = spec.pop('type')
        parameter_class = {
            'integer': neps.Integer,
            'float': neps.Float,
            'categorical': neps.Categorical,
        }[kind]
        space[name] = parameter_class(**spec)

    def evaluate_pipeline(pipeline_directory, **parameters):
        trial = current_trial()
        trial_id = trial.id

        # Stop before running: a signal arrived between trials, another
        # controller hit the failure limit, or this trial would exceed it.
        if signals.requested():
            stop_trial(
                trial_id, neps_root,
                f'{failures.STOPPED_MARKER} {signals.requested()} received before the trial started',
                signals.EXIT_CODE,
            )
        if not stop_path.exists():
            exceeded = check_failure_limit(run_root, neps_root, config['failures'])
            if exceeded is not None:
                stop = dict(
                    exceeded,
                    time=time.time(),
                    time_iso=datetime.now(timezone.utc).isoformat(),
                    trial_id=trial_id,
                    host=socket.gethostname(),
                    slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                )
                if failures.write_stop(run_root, stop):
                    logger.error('FAILURE LIMIT EXCEEDED: %s. Wrote %s', exceeded['reason'], stop_path)
        if stop_path.exists():
            reason = json.loads(stop_path.read_text()).get('reason', 'see ' + str(stop_path))
            stop_trial(
                trial_id, neps_root,
                f'{failures.STOPPED_MARKER} failure limit: {reason}',
                failures.FAILURE_LIMIT_EXIT_CODE,
            )
        if deadline is not None and trial_durations:
            remaining = deadline - time.monotonic()
            if remaining < max(trial_durations):
                stop_trial(
                    trial_id, neps_root,
                    f'{failures.STOPPED_MARKER} not enough time left for another trial '
                    f'({remaining:.0f} s before the stop signal, longest trial {max(trial_durations):.0f} s)',
                    0,
                )

        started = time.monotonic()
        try:
            with signals.trial_running():
                loss = run_pipeline(
                    config,
                    dataset,
                    pipeline_directory / 'artifacts',
                    fingerprint,
                    parameters,
                )
        except BaseException as error:
            # A time limit or cancellation counts as a failure of this trial;
            # the failure report was written by run_pipeline.
            if isinstance(error, signals.ControllerStop) or signals.requested():
                stop_trial(trial_id, neps_root, f'[timeout] {error}', signals.EXIT_CODE)
            # A failed trial also tells how long trials take.
            trial_durations.append(time.monotonic() - started)
            raise
        trial_durations.append(time.monotonic() - started)

        if storage['keep_artifacts'] == 'trajectory':
            try:
                apply_retention(neps_root, trial_id, trial.metadata.time_sampled, loss, storage['keep_top_k'])
            except Exception:
                # Never fail a finished trial over retention; --prune cleans up later.
                logger.exception('Retention failed after trial %s; archives are kept', trial_id)
        return loss

    if (neps_root / "pipeline_space.pkl").is_file():
        space_kwargs = {}
    else:
        space_kwargs = {"pipeline_space": space}

    optimizer_kwargs = {}
    if config.get('optimizer') is not None:
        optimizer_kwargs['optimizer'] = config['optimizer']

    try:
        neps.run(
            evaluate_pipeline= evaluate_pipeline,
            root_directory= neps_root,
            evaluations_to_spend= config.get('evaluations'),  # None: until the time is up
            ignore_errors= config['ignore_errors'],
            # A failed trial gets the worst loss (every objective is 1 - score),
            # so the optimizer moves away from it; without this, NePS's GP fills
            # in the median loss. The trial stays `crashed` and off the trajectory.
            objective_value_on_error=1.0,
            **space_kwargs,
            **optimizer_kwargs,
        )
    finally:
        log_run_summary(neps_root)


if __name__ == '__main__':
    main()
