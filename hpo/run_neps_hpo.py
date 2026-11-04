"""Run the NePS hyperparameter search for one shared RNA dataset."""
import argparse
import fcntl
import json
import os
import socket
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any, Optional

from hpo.pipeline import failures, signals
from hpo.pipeline.run import load_config, prepare_run
from hpo.pipeline.trial import run_pipeline

import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

NEPS_PACKAGE = 'neural-pipeline-search'
NEPS_VERSION = '0.16.0'


def read_trials(neps_root: Path) -> dict:
    """All trials from NePS's own state (one cache file, read under NePS's lock)."""
    from neps.state import NePSState

    return NePSState.create_or_load(path=neps_root, load_only=True).lock_and_read_trials()


def error_message(trial) -> Optional[str]:
    report = trial.report
    return None if report is None or report.err is None else str(report.err)


def current_trial_id() -> str:
    from neps.runtime import get_in_progress_trial

    return get_in_progress_trial().id


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
    logger.error('Stopping the controller: %s (trial %s reported as crashed)', message, trial_id)
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
    extra = f'. {len(stopped)} trial(s) stopped by the controller.' if stopped else ''
    log('Run summary: %s%s', failures.summary_line(block), extra)


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
    signals.install()

    # Convert YAML search definitions into NePS parameter objects.
    space = {}
    for name, raw in config['search'].items():
        spec = dict(raw)
        kind = spec.pop('type')
        parameter_class = {
            'integer': neps.Integer,
            'float': neps.Float,
            'categorical': neps.Categorical,
        }[kind]
        space[name] = parameter_class(**spec)

    def evaluate_pipeline(pipeline_directory, **parameters):
        trial_id = current_trial_id()

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

        try:
            with signals.trial_running():
                return run_pipeline(
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
            raise

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
            evaluations_to_spend= config['evaluations'],
            ignore_errors= config['ignore_errors'],
            **space_kwargs,
            **optimizer_kwargs,
        )
    finally:
        log_run_summary(neps_root)


if __name__ == '__main__':
    main()
