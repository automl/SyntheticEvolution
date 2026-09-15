"""One shared configuration across every dataset row; synchronous Slurm controller.

Only one controller may use a run directory. Job IDs and submission intent are
persisted before waiting. On ambiguous submission, stop rather than duplicate work.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import time
from pathlib import Path
from typing import Any, Optional, Union
import yaml
import logging
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
    return subprocess.run(args, check=True, capture_output=True, text=True, timeout=60).stdout.strip()


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
        'generator_timeout_seconds', 'poll_seconds',
        'wait_timeout_seconds', 'dssr_timeout_seconds',
    ):
        require_positive(config, key)

    for key in ('shs_seed', 'af3_seed'):
        value = config.get(key)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f'{key} must be a nonnegative integer')

    gpu = require_mapping('gpu')
    for key in ('partition', 'memory', 'gres', 'time'):
        require_text(gpu, key, 'gpu.')
    for key in ('cpus', 'concurrency'):
        require_positive(gpu, key, 'gpu.', integer=True)

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
    else:
        if not search:
            raise ValueError('NePS mode requires a nonempty search section')
        require_positive(config, 'evaluations', integer=True)
        require_text(config, 'optimizer')

    # These settings are used to launch the controller through Slurm.
    if submitting_controller:
        if mode == 'export':
            raise ValueError('Cannot submit an export as a controller job')

        if mode == 'neps':
            require_text(config, 'neps_python')

        controller = require_mapping('controller')
        for key in ('partition', 'memory', 'time'):
            require_text(controller, key, 'controller.')
        require_positive(controller, 'cpus', 'controller.', integer=True)
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
    - Relevant source code ('PIPELINE_DIR/*.py', 'HPO_DIR/slurm/*.slurm', 'generator/*.py')
    """
    code = {}
    for file in sorted(PIPELINE_DIR.glob('*.py')) + sorted((HPO_DIR / 'slurm').glob('*.slurm')):
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


def parse_accounting(accounting: str) -> dict[str, tuple[str, str]]:
    """Parse ``sacct`` output into job IDs mapped to state and exit code."""
    records = {}
    for line in accounting.splitlines():
        fields = line.split('|')
        if len(fields) < 3 or not fields[1].split():
            continue
        records[fields[0].strip()] = (
            fields[1].split()[0].rstrip('+'),
            fields[2].strip(),
        )
    return records


def wait_job(job: int, count: int, config):
    """Wait until every task in a Slurm array completes successfully.

    Polls ``sacct`` for all expected array elements and returns when each task
    has state ``COMPLETED`` with exit code ``0:0``. Raises ``RuntimeError`` if
    a task fails or completes with a non-zero exit code, and raises
    ``TimeoutError`` when the configured wait limit is exceeded.
    """
    started = time.monotonic()
    while True:
        # sacct must list every array element with successful state AND exit code.
        accounting = command(['sacct', '-j', job, '--noheader', '--parsable2', '--format=JobID%64,State%40,ExitCode'])
        records = parse_accounting(accounting)
        expected = [records.get(f'{job}_{i}') for i in range(count)]
        failed = [x for x in expected if x and x[0] in ('FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE')]
        if any(x and x[0] == 'COMPLETED' and x[1] != '0:0' for x in expected):
            raise RuntimeError(f'Array {job} has a nonzero completed-task exit code')
        if failed:
            raise RuntimeError(f'Array {job} failed: {failed}; inspect logs. No partial score returned.')
        if all(x == ('COMPLETED','0:0') for x in expected):
            logger.info('AF3 array %s finished: %d/%d tasks completed successfully.', job, count, count)
            return
        if time.monotonic() - started > config['wait_timeout_seconds']:
            raise TimeoutError(f'Waiting for array {job} timed out; job may still run. Saved ID permits reattachment.')
        completed_count = sum(x == ('COMPLETED', '0:0') for x in expected)
        logger.info('Waiting for AF3 array %s: %d/%d tasks completed.', job, completed_count, count)
        time.sleep(config['poll_seconds'])


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
    af3_job_id_path = trial_directory / 'af3_job_id.json'
    af3_tasks_path = trial_directory / 'af3_tasks.json'
    af3_submission_intent_path = trial_directory / 'af3_submission_intent.json'
    logs_directory = trial_directory / 'logs'

    trial_parameters = {**config.get('fixed', {}), **(parameters or {})}
    trial_signature = dict(fingerprint=identity, parameters=trial_parameters)

    if trial_metadata_path.exists():
        logger.info('Existing trial verified: %s', trial_directory)
        if json.loads(trial_metadata_path.read_text()) != trial_signature:
            raise ValueError('Trial parameters changed; choose a fresh trial directory')
    else:
        logger.info('New run initialized: %s', trial_directory)
        write_json(trial_metadata_path, trial_signature)

    if trial_result_path.exists():
        logger.info('Using existing result from %s', trial_result_path)
        return json.loads(trial_result_path.read_text())['loss']

    logger.info('Trial contains %d RNA rows', len(rows))
    if af3_job_id_path.exists():
        job = json.loads(af3_job_id_path.read_text())['job_id']
        logger.info('Found saved AF3 array job ID %s. Skipping shs generation and AF3 job submission.', job)
    elif af3_submission_intent_path.exists():
        raise RuntimeError(
            f'Found AF3 submission-intent file at "{af3_submission_intent_path}", '
            f'but the AF3 job ID file was not found at "{af3_job_id_path}". '
            f'Check the sbatch command in "{af3_submission_intent_path}", '
            f'AF3 stdout/stderr logs in "{logs_directory}", and the job status in Slurm.'
        )
    else:
        ###################### 1. - generate SHS (and AF3 request) ######################
        logger.info('Generating AF3 inputs and request for %d RNA rows', len(rows))
        af3_tasks = []
        for index, row in enumerate(rows):
            logger.debug('Generating input for %s (%d/%d)', row['id'], index + 1, len(rows))

            # Each single RNA chain with its structure is a task with its own files.
            task_dir = trial_directory / f'rna_{index:05d}'
            task_dir.mkdir(exist_ok=True)
            shs_generation_request_path = task_dir / 'shs_generation_request.json'
            af3_input_path = task_dir / 'af3_input.json'

            shs_generation_request = dict(**row,   # contains (id=, sequence=, pairs=)
                           task_name=f'rna_{index:05d}', 
                           parameters=trial_parameters,
                           shs_seed=config['shs_seed'],
                           af3_seed=config['af3_seed']
            )
            write_json(shs_generation_request_path, shs_generation_request)
            # Run the generator directly in its SHS environment.
            with (task_dir / 'shs_generation.log').open('w') as log:
                subprocess.run(
                    [
                        config['shs_python'],
                        str((GENERATOR_DIR / 'shs_generator.py').resolve()),
                        '--request-json', str(shs_generation_request_path.resolve()),
                        '--output-json', str(af3_input_path.resolve()),
                    ],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                    timeout=config['generator_timeout_seconds'],
                )
            af3_tasks.append(dict(input=str(af3_input_path), 
                                  output=str(task_dir / 'af3'),
                                  model_dir=str(Path(config['model_dir']).expanduser())
                                  ))
        write_json(af3_tasks_path, af3_tasks)
        logger.info('Wrote Slurm task list with %d tasks to %s', len(af3_tasks), af3_tasks_path)

        ################################# 2. - AF3 #################################
        gpu_node_settings = config['gpu']
        logs_directory.mkdir(exist_ok=True)
        job_name = 'shs-' + trial_directory.name
        array_range = f'0-{len(rows) - 1}%{gpu_node_settings["concurrency"]}'
        output_log = logs_directory / '%A_%a.out'
        error_log = logs_directory / '%A_%a.err'
        slurm_script = HPO_DIR / 'slurm/af3-inference-array.slurm'
        af3_task_script = HPO_DIR / 'af3_task.py'

        args = [
            'sbatch',
            '--parsable',
            '--job-name=' + job_name,
            '--array=' + array_range,
            '--partition=' + gpu_node_settings['partition'],
            '--cpus-per-task=' + str(gpu_node_settings['cpus']),
            '--mem=' + gpu_node_settings['memory'],
            '--gres=' + gpu_node_settings['gres'],
            '--time=' + gpu_node_settings['time'],
            '--output=' + str(output_log),
            '--error=' + str(error_log),
            str(slurm_script),
            str(af3_task_script),
            str(af3_tasks_path),
            config['module'],
        ]
        write_json(af3_submission_intent_path, dict(command=args))

        logger.info('Submitting AF3 array for %d tasks', len(af3_tasks))
        job = command(args).split(';')[0]
        if not job.isdigit():
            raise RuntimeError('Unrecognized sbatch job ID; inspect submission intent')
        write_json(af3_job_id_path, dict(job_id=job))
        logger.info('Submitted AF3 array with job ID %s', job)

    job = json.loads(af3_job_id_path.read_text())['job_id']

    logger.info('Waiting for AF3 array job with ID %s to complete', job)
    wait_job(job, len(rows), config)

    ################################ 3. - DSSR #################################
    logger.info('Scoring %d predicted structures', len(rows))
    task_evaluations = []
    for index, row in enumerate(rows):
        task_dir = trial_directory / f'rna_{index:05d}'
        predicted_model = json.loads((task_dir / 'af3/af3_success.json').read_text())['model']
        if not Path(predicted_model).is_file():
            raise FileNotFoundError(predicted_model)
        
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
            id=row['id'], **metrics, 
            target_pairs=row['pairs'], 
            predicted_pairs=sorted(normalized_pairs), 
            model=predicted_model
        )
        write_json(task_dir / 'base_pair_evaluation.json', task_evaluation)

        task_evaluations.append(task_evaluation)
        logger.info('Scored %s: loss=%.6f, f1=%.6f', row['id'], metrics['loss'], metrics['f1'])
    loss = scoring.aggregate(task_evaluations)
    write_json(trial_result_path, dict(loss=loss, scores=task_evaluations, parameters=trial_parameters, job_id=job))
    logger.info('Completed trial %s with aggregate loss %.6f', trial_directory.name, loss)
    return loss
