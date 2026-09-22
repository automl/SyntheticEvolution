"""Functions and helpers needed to execute a single trial."""

import json
import subprocess
from pathlib import Path
from typing import Any, Optional, Union
import os
import shutil

import logging
logger = logging.getLogger(__name__)

from hpo.pipeline.run import Dataset, PIPELINE_DIR, GENERATOR_DIR, normalize_pairs, write_json
import hpo.pipeline.dssr as dssr
import hpo.pipeline.scoring as scoring


def completed_af3_model(output_directory: Union[str, Path]) -> Optional[Path]:
    """Return the validated AF3 model recorded by a success marker."""

    output_directory = Path(output_directory).resolve()
    marker = output_directory / "af3_success.json"

    if not marker.is_file():
        return None

    try:
        marker_data = json.loads(marker.read_text())
        model = Path(marker_data["model"]).resolve()
    except (
        OSError,
        KeyError,
        TypeError,
        json.JSONDecodeError,
    ):
        return None

    try:
        model.relative_to(output_directory)
    except ValueError:
        return None

    if not model.is_file() or model.stat().st_size == 0:
        return None

    return model


def run_af3_tasks(tasks: list[dict[str, str]], config: dict[str, Any]) -> None:
    """Run all incomplete AF3 tasks inside the controller allocation."""

    af3_python = os.environ.get("AF3_PYTHON")
    if not af3_python:
        raise RuntimeError("AF3_PYTHON is not set. Submit the pipeline through the controller SLURM script that loads the AF3 module.")

    af3_python_path = Path(af3_python)
    if not af3_python_path.is_file():
        raise FileNotFoundError(f"AF3 Python does not exist: {af3_python_path}")

    if not tasks:
        raise ValueError("No AF3 tasks were generated")

    af3_task_script = PIPELINE_DIR / "af3_task.py"
    if not af3_task_script.is_file():
        raise FileNotFoundError(af3_task_script)

    completed = 0
    reused = 0

    for index, task in enumerate(tasks, start=1):
        input_path = Path(task["input"]).resolve()
        output_directory = Path(task["output"]).resolve()
        model_directory = Path(task["model_dir"]).expanduser().resolve()

        if completed_af3_model(output_directory) is not None:
            reused += 1
            logger.info("Reusing AF3 prediction %d/%d", index, len(tasks))
            continue

        # An output without a valid success marker is incomplete or stale.
        # Remove only this task's AF3 output, never the complete task directory.
        if output_directory.exists():
            shutil.rmtree(output_directory)

        output_directory.parent.mkdir(parents=True, exist_ok=True)
        log_path = output_directory.parent / "af3.log"

        logger.info("Running AF3 prediction %d/%d for %s", index, len(tasks), input_path)

        with log_path.open("w") as log:
            subprocess.run(
                [
                    str(af3_python_path),
                    str(af3_task_script),
                    "--input",
                    str(input_path),
                    "--output",
                    str(output_directory),
                    "--model-dir",
                    str(model_directory),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=config["af3_timeout_seconds"],
            )

        model = completed_af3_model(output_directory)
        if model is None:
            raise RuntimeError(f"AF3 exited successfully but did not publish a valid success marker for {input_path}. Inspect {log_path}.")
        completed += 1
        
    logger.info(
        "All %d AF3 tasks are ready: %d executed, %d reused",
        len(tasks),
        completed,
        reused,
    )


def run_pipeline(
    config: dict[str, Any],
    rows: Dataset,
    trial_directory: Union[str, Path],
    identity: str,
    parameters: Optional[dict[str, Any]] = None
) -> float:
    """Run one parameter configuration across every RNA in the dataset."""

    trial_directory = Path(trial_directory).resolve()
    trial_directory.mkdir(parents=True, exist_ok=True)
    logger.info('Preparing trial in %s', trial_directory)

    # Stable files and directories used throughout this trial.
    trial_metadata_path = trial_directory / 'trial_metadata.json'
    trial_result_path = trial_directory / 'trial_result.json'
    af3_tasks_path = trial_directory / 'af3_tasks.json'

    trial_parameters = {**config.get('fixed', {}), **(parameters or {})}
    trial_signature = dict(fingerprint=identity, parameters=trial_parameters)

    if trial_metadata_path.exists():
        logger.info('Existing trial verified: %s', trial_directory)
        if json.loads(trial_metadata_path.read_text()) != trial_signature:
            raise ValueError('Trial parameters changed; choose a fresh trial directory')
    else:
        logger.info('New trial initialized')
        write_json(trial_metadata_path, trial_signature)

    if trial_result_path.exists():
        logger.info('Using existing result from %s', trial_result_path)
        return json.loads(trial_result_path.read_text())['loss']

    logger.info('Trial contains %d RNA rows', len(rows))
    ############################# 1. - generate SHS ###########################
    logger.info('Generating AF3 inputs and request for %d RNA rows', len(rows))
    af3_tasks = []
    for index, row in enumerate(rows):
        task_dir = trial_directory / f"rna_{index:05d}"
        task_dir.mkdir(exist_ok=True)
        request_path = task_dir / "shs_generation_request.json"
        af3_input_path = task_dir / "af3_input.json"
        af3_output_directory = task_dir / "af3"

        shs_generation_request = dict(
            **row,   # contains (id=, sequence=, pairs=)
            task_name=f'rna_{index:05d}', 
            parameters=trial_parameters,
            shs_seed=config['shs_seed'],
            af3_seed=config['af3_seed']
        )
        write_json(request_path, shs_generation_request)

        # Run the generator directly in its SHS environment.
        if not af3_input_path.is_file():
            logger.info("Generating AF3 input %d/%d for %s", index + 1, len(rows), row["id"])

            with (task_dir / 'shs_generation.log').open('w') as log:
                subprocess.run(
                    [
                        config['shs_python'],
                        str((GENERATOR_DIR / 'shs_generator.py').resolve()),
                        '--request-json', str(request_path.resolve()),
                        '--output-json', str(af3_input_path.resolve()),
                    ],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                    timeout=config['generator_timeout_seconds'],
                )
    ################################# 2. - AF3 #################################
        if not af3_input_path.is_file():
            raise RuntimeError(f"SHS generation completed without creating {af3_input_path}")
        af3_tasks.append(
            {   "input": str(af3_input_path.resolve()),
                "output": str(af3_output_directory.resolve()),
                "model_dir": str(
                    Path(config["model_dir"]).expanduser().resolve()),
            })

    write_json(af3_tasks_path, af3_tasks)
    run_af3_tasks(af3_tasks, config)

    ################################ 3. - DSSR #################################
    logger.info('Scoring %d predicted structures', len(rows))
    task_evaluations = []
    for index, row in enumerate(rows):
        task_dir = trial_directory / f'rna_{index:05d}'
        predicted_model = completed_af3_model(task_dir / "af3")
        if predicted_model is None:
            raise RuntimeError(f"No valid AF3 prediction found for {row['id']}")
        
        dssr_predicted_pairs = dssr.run_dssr(
            predicted_model,
            task_dir / 'dssr',
            row['sequence'],
            config['dssr'],
            config['dssr_timeout_seconds']
        )
        normalized_pairs = normalize_pairs(dssr_predicted_pairs, len(row['sequence']))
    
        ################################ 4. - SCORE #################################
        metrics = scoring.score_row(row, normalized_pairs)

        task_evaluation = dict(
            id=row['id'], 
            **metrics, 
            target_pairs=row['pairs'],
            **{key: row[key] for key in ('interactions', 'mutation_rates') if key in row}, 
            predicted_pairs=sorted(normalized_pairs), 
            model=str(predicted_model)
        )
        write_json(task_dir / 'base_pair_evaluation.json', task_evaluation)

        task_evaluations.append(task_evaluation)
        logger.info('Scored %s: loss=%.6f, f1=%.6f', row['id'], metrics['loss'], metrics['f1'])
    loss = scoring.aggregate(task_evaluations)
    write_json(
        trial_result_path, 
        dict(
            loss=loss, 
            scores=task_evaluations, 
            parameters=trial_parameters
            )
        )
    logger.info('Completed trial %s with aggregate loss %.6f', trial_directory.name, loss)
    return loss
