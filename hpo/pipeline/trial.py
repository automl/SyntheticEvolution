"""Functions and helpers needed to execute a single trial."""

import json
import os
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Optional, Union

import logging
logger = logging.getLogger(__name__)

from hpo.pipeline.run import (
    Dataset, GENERATOR_DIR, normalize_pairs, scratch_root, shs_generation_request, write_json,
)
import hpo.pipeline.archive as archive
import hpo.pipeline.dssr as dssr
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
    """Copy all logs from scratch (standalone trials keep their logs)."""
    for log in sorted(scratch.glob('rna_*/*.log')):
        target = destination / log.parent.name / log.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(log, target)


def run_af3(af3_input: Path, output: Path, log_path: Path, config: dict[str, Any]) -> Path:
    """Run one AF3 prediction and return its top-level model.

    The controller job loads the AF3 module (``controller.slurm``), which
    provides AF3_PYTHON, ALPHAFOLD_BIN_DIR and ALPHAFOLD_DATABASES.
    """
    names = ('AF3_PYTHON', 'ALPHAFOLD_BIN_DIR', 'ALPHAFOLD_DATABASES')
    missing = [name for name in names if not os.environ.get(name)]
    if missing:
        raise RuntimeError(
            f"{', '.join(missing)} not set. Submit the pipeline through the "
            "controller Slurm script, which loads the AF3 module."
        )
    with log_path.open('w') as log:
        subprocess.run(
            [
                os.environ['AF3_PYTHON'],
                str(Path(os.environ['ALPHAFOLD_BIN_DIR']) / 'run_alphafold.py'),
                f'--json_path={af3_input}',
                f'--output_dir={output}',
                f"--model_dir={Path(config['model_dir']).expanduser().resolve()}",
                f"--db_dir={os.environ['ALPHAFOLD_DATABASES']}",
                '--norun_data_pipeline',
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=config['af3_timeout_seconds'],
        )
    # AF3 writes <output>/<name>/<name>_model.cif; samples live one level deeper.
    models = list(output.glob(f'*/*{archive.MODEL_SUFFIX}'))
    if len(models) != 1 or models[0].stat().st_size == 0:
        raise RuntimeError(f'Expected exactly one nonempty AF3 model in {output}, found {models}')
    return models[0]


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

    trial_parameters = {**config.get('fixed', {}), **(parameters or {})}
    trial_signature = dict(fingerprint=identity, parameters=trial_parameters)

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
    trial_parameters: dict[str, Any],
    task_scratch: Path,
) -> Path:
    """Generate the AF3 input with the synthetic MSA in scratch."""
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
                config['shs_python'],
                str((GENERATOR_DIR / 'shs_generator.py').resolve()),
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


def _run_stages(
    config: dict[str, Any],
    rows: Dataset,
    trial_directory: Path,
    scratch: Path,
    trial_signature: dict[str, Any],
    progress: dict[str, Any],
) -> float:
    """Predict and score every RNA in turn, then publish the trial result.

    Per RNA: generate and predict in scratch, publish the archive right after
    AF3 (so a later failure cannot lose the prediction), then run DSSR and add
    ``dssr.json`` to the archive. An existing readable archive skips
    generation and AF3.
    """
    trial_result_path = trial_directory / 'trial_result.json'
    trial_parameters = trial_signature['parameters']

    logger.info('Trial contains %d RNA rows', len(rows))
    task_evaluations = []
    for index, row in enumerate(rows):
        task_name = f'rna_{index:05d}'
        task_scratch = scratch / task_name
        content = task_scratch / 'archive'
        content.mkdir(parents=True)
        rna_archive = archive.path_for(trial_directory, task_name)

        ######################## 1.-2. - generate SHS, AF3 ########################
        progress.update(stage='generation', rna_index=index)
        model_name = archive.model_name(rna_archive)
        if model_name is not None:
            logger.info('Reusing AF3 prediction %d/%d for %s', index + 1, len(rows), row['id'])
            archive.extract(rna_archive, content)
        else:
            if rna_archive.exists():
                logger.warning('Replacing unreadable archive %s', rna_archive)
                rna_archive.unlink()
            logger.info('Generating AF3 input %d/%d for %s', index + 1, len(rows), row['id'])
            af3_input_path = _generate(config, row, index, trial_parameters, task_scratch)

            progress.update(stage='af3')
            logger.info('Running AF3 prediction %d/%d for %s', index + 1, len(rows), row['id'])
            model = run_af3(af3_input_path, task_scratch / 'af3', task_scratch / 'af3.log', config)

            # Flatten AF3's <name>/ folder; its data.json holds the synthetic MSA.
            for path in model.parent.iterdir():
                shutil.move(str(path), str(content / path.name))
            effective = task_scratch / 'shs_effective_parameters.json'
            if effective.is_file():
                shutil.copy2(effective, content / effective.name)

            progress.update(stage='archive')
            archive.publish(content, rna_archive)
            model_name = model.name

        ################################ 3. - DSSR #################################
        progress.update(stage='dssr')
        (content / 'dssr.json').unlink(missing_ok=True)
        try:
            dssr_predicted_pairs = dssr.run_dssr(
                content / model_name,
                content,
                row['sequence'],
                config['dssr'],
                config['dssr_timeout_seconds'],
                scratch=task_scratch,
            )
        finally:
            # dssr.json is written before its pairs are validated; keep it for
            # analysis, also when the validation fails.
            if (content / 'dssr.json').is_file():
                archive.publish(content, rna_archive)
        normalized_pairs = normalize_pairs(dssr_predicted_pairs, len(row['sequence']))

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
    logger.info('Completed trial %s with aggregate loss %.6f', trial_directory.name, loss)
    return loss
