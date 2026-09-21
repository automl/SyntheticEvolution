"""One shared configuration across every dataset row; synchronous Slurm controller.

Only one controller may use a run directory. Job IDs and submission intent are
persisted before waiting. On ambiguous submission, stop rather than duplicate work.
"""
import csv
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any, Optional, Union
import yaml
import logging
import os
import shutil
logger = logging.getLogger(__name__)

import hpo.pipeline.dssr as dssr
import hpo.pipeline.scoring as scoring

PIPELINE_DIR = Path(__file__).resolve().parent  # SyntheticEvolution/hpo/pipeline
HPO_DIR = PIPELINE_DIR.parent                   # SyntheticEvolution/hpo
ROOT = HPO_DIR.parent                           # SyntheticEvolution
GENERATOR_DIR = ROOT / 'SHS-Generator'

Dataset = list[dict[str, Union[str, list[tuple[int, int]]]]]


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True))
    temporary.replace(path)


def repo_path(value):
    """Convert path to a global path. 
    If it is not absolute use repo root as starting point."""
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    else:
        return (ROOT / path).resolve() 


def command(args):
    try:
        result = subprocess.run(
            args,
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        )
    except subprocess.CalledProcessError as exc:
        logger.error(
            'Command failed: %s\n%s',
            args,
            (exc.stderr or exc.stdout or 'No error output').strip(),
        )
        raise
    return result.stdout.strip()


def load_config(path, *, mode, submitting_controller=False):
    """Load YAML and validate settings required by the selected entry point."""
    if mode not in {'standalone', 'neps', 'export'}:
        raise ValueError("mode must be 'standalone', 'neps', or 'export'")

    def require_text(settings, key, prefix=''):
        value = settings.get(key)
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f'{prefix}{key} must be a nonempty string')
        return value

    config = yaml.safe_load(repo_path(path).read_text())
    if not isinstance(config, dict):
        raise ValueError('Config must be a YAML mapping')


    def require_positive(settings, key, prefix='', *, integer=False):
        value = settings.get(key)
        expected_type = int if integer else (int, float)
        if (
            isinstance(value, bool)
            or not isinstance(value, expected_type)
            or not math.isfinite(value)
            or value <= 0
        ):
            kind = 'integer' if integer else 'number'
            raise ValueError(f'{prefix}{key} must be a positive {kind}')
        return value

    def require_mapping(key):
        value = config.get(key)
        if not isinstance(value, dict):
            raise ValueError(f'{key} must be a YAML mapping')
        return value

    # Settings shared by standalone and NePS trials.
    run_name = require_text(config, 'run_name')
    if not re.fullmatch(r'[A-Za-z0-9_-]+', run_name):
        raise ValueError(
            'run_name must contain only letters, digits, underscore or hyphen'
        )

    require_text(config, 'workspace_name')

    # Export reads saved results; execution settings are not needed.
    if mode == 'export':
        return config

    for key in (
        'input', 'shs_python', 'module', 'model_dir', 'dssr',
    ):
        require_text(config, key)

    for key in (
        'generator_timeout_seconds', 'dssr_timeout_seconds', 'af3_timeout_seconds'
    ):
        require_positive(config, key)

    for key in ('shs_seed', 'af3_seed'):
        value = config.get(key)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f'{key} must be a nonnegative integer')

    # Optional parameter sections must have the expected structure.
    for key in ('fixed', 'search'):
        if key in config and not isinstance(config[key], dict):
            raise ValueError(
                f'{key} must be a YAML mapping; use {{}} for an empty section'
            )

    fixed = config.get('fixed', {})
    search = config.get('search', {})
    overlap = fixed.keys() & search.keys()
    if overlap:
        raise ValueError(
            'Parameters cannot appear in both fixed and search: '
            + ', '.join(sorted(map(str, overlap)))
        )

    if mode == 'standalone':
        if 'search' in config:
            logger.warning(
                "Standalone mode: ignoring 'search'. Using 'fixed' values "
                "and generator defaults for all remaining parameters."
            )
    elif mode == 'neps':
        if not search:
            raise ValueError('NePS mode requires a nonempty search section')
        require_positive(config, 'evaluations', integer=True)
        if config.get('optimizer') is not None:
            require_text(config, 'optimizer')
        if not isinstance(config.get('ignore_errors'), bool):
            raise ValueError('ignore_errors must be a YAML boolean: true or false')
        require_positive

    # These settings are used to launch the controller through Slurm.
    if submitting_controller:
        if mode == 'export':
            raise ValueError('Cannot submit an export as a controller job')

        controller = require_mapping('controller')

        if mode == 'neps':
            require_text(config, 'neps_python')
            require_positive(controller, "count", "controller.", integer=True)

        if mode == 'standalone' and "count" in controller:
            logger.warning("Standalone mode: ignoring 'controller.count'")

        for key in ("partition", "memory", "gres", "time"):
            require_text(controller, key, "controller.")
        require_positive(controller, "cpus", "controller.", integer=True)
    return config


def normalize_pairs(pairs: list[list[int]], size: int) -> set[tuple[int, int]]:
    """Validate RNA base pairs and return unique pairs in canonical order."""
    result = set()
    for pair in pairs:
        if len(pair) != 2 or any(type(x) is not int for x in pair):
            raise ValueError('Pairs must contain exactly two integer positions')
        a, b = pair
        if not (0 <= a < size and 0 <= b < size) or a == b:
            raise ValueError(f'Invalid pair {pair} for length {size}')
        result.add(tuple(sorted((a, b))))
    return result

def create_dataset(path) -> Dataset:
    """Create a dataset from a path.
    ### Input:
    @path : single file or dir of files.
    ### returns
    dataset = list of dicts for each line with

        id : identifying string
        sequence : RNA sequence str containing only ACGU
        pairs : sorted list of unique base pairs in canonical order
    ### Raises
    ValueError for empty dataset or duplicate ids or invalid sequences
    """
    path = repo_path(path)
    files = sorted(path.glob('*.csv')) if path.is_dir() else [path]
    rows, seen = [], set()
    for file in files:
        with file.open(newline='') as stream:
            for row in csv.DictReader(stream):
                name, seq = row['id'].strip(), row['sequence'].strip().upper()
                if not name or name in seen or not seq or set(seq) - set('ACGU'):
                    raise ValueError(f'Invalid/duplicate ID or RNA sequence in {file}: {name}')
                seen.add(name)
                pairs = sorted(normalize_pairs(json.loads(row['pairs']), len(seq)))
                rows.append(dict(id=name, sequence=seq, pairs=pairs))
    if not rows:
        raise ValueError('Dataset is empty')
    return rows


def get_workspace_path(config) -> Path: 
    """Get absolut workspace path based on config."""
    raw = config.get('workspace_path') or command(['ws_find', config['workspace_name']])
    path = Path(raw).expanduser()
    if not path.is_absolute() or not path.is_dir():
        raise ValueError(f'Workspace must exist and have an absolute path: {raw}')
    return path.resolve() / 'shs_hpo' / config['run_name']


def fingerprint(config, rows):
    """Create an identifying fingerprint of the run. This prevents accidentally 
    reusing old results for a meaningfully different experiment. It includes:

    - Dataset rows
    - Scientific configuration (everything but 'evaluations', 'controller', 'poll_seconds', 'wait_timeout_seconds')
    - Relevant source code ('PIPELINE_DIR/*.py', 'generator/*.py')
    """
    code = {}
    for file in sorted(PIPELINE_DIR.glob('*.py')):
        code[str(file.relative_to(ROOT))] = hashlib.sha256(file.read_bytes()).hexdigest()
    for file in sorted(repo_path(GENERATOR_DIR).glob('*.py')):
        code[str(file)] = hashlib.sha256(file.read_bytes()).hexdigest()
    scientific = {k:v for k,v in config.items() if k not in ('evaluations', 'controller', 'poll_seconds', 'wait_timeout_seconds')}
    return hashlib.sha256(json.dumps([scientific, rows, code], sort_keys=True).encode()).hexdigest()


def prepare_run(config) -> tuple[Path, Dataset, str]:
    """Validate inputs, create the run workspace, create persist run metadata.

    Returns:
        A tuple containing the workspace path, validated dataset, and run
        fingerprint.

    Raises:
        ValueError: If the run name, dataset, or existing run metadata is invalid.
        FileNotFoundError: If a required generator dependency is missing.
    """
    logger.info('Preparing run <%s>.', config['run_name'])
    rows = create_dataset(config['input'])

    for name in ('shs_generator.py', 'pair_map.py', 'json_generator.py'):
        if not (repo_path(GENERATOR_DIR) / name).is_file():
            raise FileNotFoundError(GENERATOR_DIR / name)

    ws_root = get_workspace_path(config)
    ws_root.mkdir(parents=True, exist_ok=True)

    identity = fingerprint(config, rows)

    guard_path = ws_root / 'run_metadata.json'
    if guard_path.exists():
        logger.info('Existing run verified: %s', ws_root)
        if json.loads(guard_path.read_text())['fingerprint'] != identity:
            raise ValueError('Run inputs/code/config changed. Choose a new run_name.')
    else:
        logger.info('New run initialized: %s', ws_root)
        write_json(guard_path, dict(fingerprint=identity, config=config, dataset=rows))
    return ws_root, rows, identity


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
        target_pairs = map(tuple,row['pairs'])
        metrics = scoring.score_pairs(target_pairs, normalized_pairs)

        task_evaluation = dict(
            id=row['id'], 
            **metrics, 
            target_pairs=row['pairs'], 
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
