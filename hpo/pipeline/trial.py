"""Functions and helpers needed to execute a single trial."""

import json
import os
import re
import select
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Optional, Union

import logging
logger = logging.getLogger(__name__)

from hpo.pipeline.run import (
    Dataset, GENERATOR_DIRS, normalize_pairs, scratch_root, shs_generation_request, typed_parameters,
    write_json,
)
import hpo.pipeline.archive as archive
import hpo.pipeline.extraction as extraction
import hpo.pipeline.failures as failures
import hpo.pipeline.scoring as scoring


@contextmanager
def trial_scratch(root: Path) -> Iterator[Path]:
    """Private scratch directory for one trial, removed when the trial ends.

    Logs and temporary tool files live here, so they never reach the shared
    workspace. The Slurm job ID in the name identifies leftovers of killed jobs.
    """
    root.mkdir(parents=True, exist_ok=True)
    job = os.environ.get('SLURM_JOB_ID', 'local')
    path = Path(tempfile.mkdtemp(prefix=f'shs_hpo_{job}_', dir=root))
    try:
        yield path
    finally:
        shutil.rmtree(path, ignore_errors=True)


def keep_logs(scratch: Path, destination: Path) -> None:
    """Copy all logs from scratch (standalone trials keep their logs).

    Per-RNA logs keep their ``rna_<index>/`` folder; the trial's AF3 log
    (one for all RNAs) lands directly in ``destination``.
    """
    for log in sorted([*scratch.glob('*.log'), *scratch.glob('rna_*/*.log')]):
        target = destination / log.relative_to(scratch)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(log, target)


# This controller's AF3 worker (af3_worker.py) and the command that started it.
AF3_WORKER = Path(__file__).resolve().parent / 'af3_worker.py'
_af3_worker: Optional[tuple[list[str], subprocess.Popen]] = None


def stop_af3_worker() -> None:
    """Kill the AF3 worker, e.g. after a timeout or while it predicts."""
    global _af3_worker
    if _af3_worker is not None:
        _af3_worker[1].kill()
        _af3_worker[1].wait()
        _af3_worker = None


def run_af3(af3_inputs: dict[str, Path], output: Path, log_path: Path, config: dict[str, Any]) -> None:
    """Predict all inputs with this controller's persistent AF3 worker.

    A new AF3 process loads and compiles its model first (~37 s on an A100,
    more than a 36-nt prediction takes), so one worker serves every trial of
    the controller and compiles once per bucket. It starts with the first
    trial and again after a failure, timeout or stop. ``af3_inputs`` maps task
    names to inputs; AF3 writes ``output/<name>/`` and its log to
    ``log_path``. ``af3_timeout_seconds`` applies per input.

    The controller job loads the AF3 module (``controller.slurm``), which
    provides AF3_PYTHON and ALPHAFOLD_BIN_DIR.
    """
    global _af3_worker
    names = ('AF3_PYTHON', 'ALPHAFOLD_BIN_DIR')
    missing = [name for name in names if not os.environ.get(name)]
    if missing:
        raise RuntimeError(
            f"{', '.join(missing)} not set. Submit the pipeline through the "
            "controller Slurm script, which loads the AF3 module."
        )
    input_directory = output.parent / 'af3_inputs'
    input_directory.mkdir(parents=True, exist_ok=True)
    for name, path in af3_inputs.items():
        shutil.copy2(path, input_directory / f'{name}.json')
    command = [
        os.environ['AF3_PYTHON'], '-u', str(AF3_WORKER),
        f"--model_dir={Path(config['model_dir']).expanduser().resolve()}",
    ]
    # AF3 draws the diffusion noise for the padded size, so the buckets change
    # the predictions like a different seed would.
    if config.get('af3_buckets'):
        command.append(f"--buckets={','.join(map(str, config['af3_buckets']))}")
    if _af3_worker is not None and (_af3_worker[0] != command or _af3_worker[1].poll() is not None):
        stop_af3_worker()
    if _af3_worker is None:
        # Start-up output and errors go to the log of the trial that starts it.
        with log_path.open('a') as log:
            _af3_worker = command, subprocess.Popen(
                command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log, text=True,
            )
    process = _af3_worker[1]
    request = dict(input_dir=str(input_directory), output_dir=str(output), log=str(log_path))
    timeout = config['af3_timeout_seconds'] * len(af3_inputs)
    try:
        process.stdin.write(json.dumps(request) + '\n')
        process.stdin.flush()
        ready, _, _ = select.select([process.stdout], [], [], timeout)
        reply = process.stdout.readline() if ready else None
    except BaseException:
        # A stop signal or a dead worker: the model state is unknown.
        stop_af3_worker()
        raise
    if reply is None:
        stop_af3_worker()
        raise subprocess.TimeoutExpired(command, timeout)
    if reply != 'done\n':
        returncode = process.wait()
        _af3_worker = None
        raise subprocess.CalledProcessError(returncode, command)


def af3_predictions(output: Path, log_path: Path) -> dict[str, Path]:
    """Models of the inputs AF3 finished, also after AF3 failed or was stopped.

    An input counts as finished once AF3 logged that it is done and wrote
    exactly one nonempty top-level model (samples live one level deeper).
    """
    log = log_path.read_text(errors='replace') if log_path.is_file() else ''
    predictions = {}
    for name in re.findall(r'^Done processing fold input (\S+)\.$', log, re.MULTILINE):
        models = list((output / name).glob(f'*{archive.MODEL_SUFFIX}'))
        if len(models) == 1 and models[0].stat().st_size > 0:
            predictions[name] = models[0]
    return predictions


def run_pipeline(
    config: dict[str, Any],
    rows: Dataset,
    trial_directory: Union[str, Path],
    identity: str,
    parameters: Optional[dict[str, Any]] = None,
    keep_logs_in_trial: bool = False,
) -> float:
    """Run one parameter configuration across every RNA in the dataset.

    The workspace receives one archive per RNA and one trial record:
    ``trial_result.json`` on success; ``trial_metadata.json`` while running and
    after a failure, together with ``failure_report.json``. Everything else
    lives in a scratch directory. ``keep_logs_in_trial`` (standalone mode)
    copies the tool logs to ``<trial>/logs/``.
    """

    trial_directory = Path(trial_directory).resolve()
    trial_directory.mkdir(parents=True, exist_ok=True)
    logger.info('Preparing trial in %s', trial_directory)
    archive.remove_stale_temporaries(trial_directory)

    trial_metadata_path = trial_directory / 'trial_metadata.json'
    trial_result_path = trial_directory / 'trial_result.json'

    trial_parameters = typed_parameters(config, rows, parameters)
    # Seeds are part of the signature: standalone trials may set their own.
    trial_signature = dict(
        fingerprint=identity, parameters=trial_parameters,
        shs_seed=config['shs_seed'], af3_seed=config['af3_seed'],
    )

    # A finished trial keeps only its result, which also holds the signature.
    if trial_result_path.exists():
        result = json.loads(trial_result_path.read_text())
        if {key: result.get(key) for key in trial_signature} != trial_signature:
            raise ValueError('Trial parameters changed; choose a fresh trial directory')
        logger.info('Using existing result from %s', trial_result_path)
        return result['loss']

    if trial_metadata_path.exists():
        logger.info('Existing trial verified: %s', trial_directory)
        if json.loads(trial_metadata_path.read_text()) != trial_signature:
            raise ValueError('Trial parameters changed; choose a fresh trial directory')
    else:
        logger.info('New trial initialized')
        write_json(trial_metadata_path, trial_signature)

    # Stages record where they are, so a failure report can name the stage and RNA.
    progress: dict[str, Any] = {}
    with trial_scratch(scratch_root(config)) as scratch:
        try:
            loss = _run_stages(
                config, rows, trial_directory, scratch, trial_signature, progress
            )
        except BaseException as error:
            # BaseException: timeouts and cancellations must be reported as well.
            try:
                report = failures.write_failure_report(
                    trial_directory, rows, progress, error,
                    log_directory=scratch, storage=storage_usage(trial_directory),
                )
            except Exception:
                logger.exception('Could not write the failure report for %s', trial_directory)
            else:
                logger.error(
                    'Trial FAILED: %s. Report: %s',
                    failures.describe(report),
                    trial_directory / failures.FAILURE_REPORT,
                )
            raise
        finally:
            if keep_logs_in_trial:
                try:
                    keep_logs(scratch, trial_directory / 'logs')
                except Exception:
                    logger.exception('Could not copy the logs to %s', trial_directory / 'logs')

    # The result replaces the running record; an earlier failure report no
    # longer applies.
    trial_metadata_path.unlink(missing_ok=True)
    (trial_directory / failures.FAILURE_REPORT).unlink(missing_ok=True)
    return loss


def storage_usage(trial_directory: Path, exclude: tuple[str, ...] = ()) -> dict[str, int]:
    """Files and bytes the trial has published (trial directory, top level)."""
    files = [
        path for path in trial_directory.iterdir()
        if path.is_file() and not path.name.startswith('.') and path.name not in exclude
    ]
    return dict(files=len(files), bytes=sum(path.stat().st_size for path in files))


def _generate(
    config: dict[str, Any],
    row: dict[str, Any],
    index: int,
    trial_parameters: dict[str, dict[str, Any]],
    task_scratch: Path,
) -> Path:
    """Generate the AF3 input with the synthetic MSA in scratch.

    The generator and its parameters are chosen by the chain's type.
    """
    chain_type = row['type']
    request_path = task_scratch / 'shs_generation_request.json'
    af3_input_path = task_scratch / 'af3_input.json'
    write_json(request_path, shs_generation_request(
        row, index, trial_parameters, config['shs_seed'], config['af3_seed']
    ))

    # Run the generator directly in its SHS environment. It writes
    # shs_effective_parameters.json next to the request.
    with (task_scratch / 'shs_generation.log').open('w') as log:
        subprocess.run(
            [
                config['generators'][chain_type]['python'],
                str((GENERATOR_DIRS[chain_type] / 'shs_generator.py').resolve()),
                '--request-json', str(request_path),
                '--output-json', str(af3_input_path),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=config['generator_timeout_seconds'],
        )
    if not af3_input_path.is_file():
        raise RuntimeError(f"SHS generation completed without creating af3_input.json for {row['id']}")
    return af3_input_path


def _archive_predictions(
    pending: list[str],
    scratch: Path,
    trial_directory: Path,
    progress: dict[str, Any],
) -> list[str]:
    """Publish the archive of every prediction AF3 finished; return their names.

    Flattens AF3's ``<name>/`` folder (its data.json holds the synthetic MSA)
    and adds the generator's effective parameters.
    """
    predictions = af3_predictions(scratch / 'af3', scratch / 'af3.log')
    progress.update(stage='archive')
    archived = []
    for name in pending:
        if name not in predictions:
            continue
        content = scratch / name / 'archive'
        for path in predictions[name].parent.iterdir():
            shutil.move(str(path), str(content / path.name))
        effective = scratch / name / 'shs_effective_parameters.json'
        if effective.is_file():
            shutil.copy2(effective, content / effective.name)
        archive.publish(content, archive.path_for(trial_directory, name))
        archived.append(name)
    return archived


def _run_stages(
    config: dict[str, Any],
    rows: Dataset,
    trial_directory: Path,
    scratch: Path,
    trial_signature: dict[str, Any],
    progress: dict[str, Any],
) -> float:
    """Predict and score every RNA, then publish the trial result.

    1. Generate the AF3 inputs of all RNAs without a readable archive.
    2. Predict them with the AF3 worker and publish each archive; also after
       an AF3 failure or stop, for the inputs AF3 finished.
    3. Per chain, extract the pairs (DSSR for RNA, contacts for protein), add
       ``dssr.json`` or ``contacts.json`` to the archive and score.
    """
    trial_result_path = trial_directory / 'trial_result.json'
    trial_parameters = trial_signature['parameters']
    task_names = [f'rna_{index:05d}' for index in range(len(rows))]
    logger.info('Trial contains %d RNA rows', len(rows))

    ############################ 1. - generate SHS #############################
    af3_inputs = {}
    for index, (row, task_name) in enumerate(zip(rows, task_names)):
        progress.update(stage='generation', rna_index=index)
        content = scratch / task_name / 'archive'
        content.mkdir(parents=True)
        rna_archive = archive.path_for(trial_directory, task_name)
        if archive.model_name(rna_archive) is not None:
            logger.info('Reusing AF3 prediction %d/%d for %s', index + 1, len(rows), row['id'])
            archive.extract(rna_archive, content)
            continue
        if rna_archive.exists():
            logger.warning('Replacing unreadable archive %s', rna_archive)
            rna_archive.unlink()
        logger.info('Generating AF3 input %d/%d for %s', index + 1, len(rows), row['id'])
        af3_inputs[task_name] = _generate(config, row, index, trial_parameters, scratch / task_name)

    ################################# 2. - AF3 #################################
    if af3_inputs:
        pending = list(af3_inputs)
        progress.update(stage='af3', rna_index=None)
        logger.info('Running AF3 for %d RNA(s)', len(pending))
        try:
            run_af3(af3_inputs, scratch / 'af3', scratch / 'af3.log', config)
        except BaseException:
            # Keep what AF3 finished, then name the RNA it was working on.
            try:
                archived = _archive_predictions(pending, scratch, trial_directory, progress)
            except Exception:
                logger.exception('Could not archive the finished AF3 predictions')
                archived = []
            unfinished = [name for name in pending if name not in archived]
            progress.update(stage='af3', rna_index=task_names.index(unfinished[0]) if unfinished else None)
            raise
        archived = _archive_predictions(pending, scratch, trial_directory, progress)
        unfinished = [name for name in pending if name not in archived]
        if unfinished:
            progress.update(stage='af3', rna_index=task_names.index(unfinished[0]))
            raise RuntimeError(f'AF3 finished without a valid model for {unfinished[0]}; see af3.log')

    task_evaluations = []
    for index, (row, task_name) in enumerate(zip(rows, task_names)):
        task_scratch = scratch / task_name
        content = task_scratch / 'archive'
        rna_archive = archive.path_for(trial_directory, task_name)
        model_name = archive.model_name(rna_archive)

        ############################## 3. - EXTRACTION #############################
        progress.update(stage='extraction', rna_index=index)
        for stale in extraction.RESULT_FILES:
            (content / stale).unlink(missing_ok=True)
        try:
            extracted_pairs = extraction.extract_pairs(
                row['type'],
                content / model_name,
                content,
                row['sequence'],
                config['extraction'],
                scratch=task_scratch,
            )
        finally:
            # dssr.json is written before its pairs are validated; keep it for
            # analysis, also when the validation fails.
            if any((content / name).is_file() for name in extraction.RESULT_FILES):
                archive.publish(content, rna_archive)
        normalized_pairs = normalize_pairs(extracted_pairs, len(row['sequence']))

        ################################ 4. - SCORE #################################
        progress.update(stage='scoring')
        metrics = scoring.score_row(
            row,
            normalized_pairs,
            objective=config.get('scoring_objective', 'f1'),
            wl_iterations=config.get('wl_iterations', 5),
        )

        task_evaluation = dict(
            id=row['id'],
            **metrics,
            target_pairs=row['pairs'],
            **{key: row[key] for key in ('interactions', 'mutation_rates') if key in row},
            predicted_pairs=sorted(normalized_pairs),
            # The model inside the RNA's archive, relative to the trial.
            model=f'{rna_archive.name}/{model_name}',
        )
        task_evaluations.append(task_evaluation)

        logger.info(
            'Scored %s: loss=%.6f, objective=%s',
            row['id'],
            metrics['loss'],
            config.get('scoring_objective', 'f1'),
        )
    progress.update(stage='aggregation', rna_index=None)
    loss = scoring.aggregate(task_evaluations)
    write_json(
        trial_result_path,
        dict(
            **trial_signature,
            loss=loss,
            scores=task_evaluations,
            # Without the records that are removed once the result exists.
            storage=storage_usage(
                trial_directory, exclude=('trial_metadata.json', failures.FAILURE_REPORT)
            ),
        ),
    )
    # e.g. config_12/artifacts (NePS), <run_name>/standalone_trial or trials/<name>
    trial_name = f'{trial_directory.parent.name}/{trial_directory.name}'
    logger.info('Completed trial %s with aggregate loss %.6f', trial_name, loss)
    return loss
