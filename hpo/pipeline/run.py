"""Shared helpers for a whole run and run validation."""

import csv
import hashlib
import json
import subprocess
import math
import re
import yaml
from typing import Union
from pathlib import Path

import logging
logger = logging.getLogger(__name__)

PIPELINE_DIR = Path(__file__).resolve().parent  # SyntheticEvolution/hpo/pipeline
HPO_DIR = PIPELINE_DIR.parent                   # SyntheticEvolution/hpo
ROOT = HPO_DIR.parent                           # SyntheticEvolution
GENERATOR_DIR = ROOT / 'SHS-Generator'

Dataset = list[dict[str, Union[str, list[tuple[int, int]]]]]

def repo_path(value):
    """Convert path to a global path. 
    If it is not absolute use repo root as starting point."""
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    else:
        return (ROOT / path).resolve()


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


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True))
    temporary.replace(path)


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
    - Scientific configuration
    - Relevant source code ('PIPELINE_DIR/*.py', 'generator/*.py')
    """
    non_scientific_configurations = (
        'controller', 'evaluations', 'ignore_errors',
        'dssr_timeout_seconds', 'af3_timeout_seconds', 'generator_timeout_seconds',
        'controller_python', 'neps_python', 'shs_python', 
    )

    code = {}
    for file in sorted(PIPELINE_DIR.glob('*.py')):
        code[str(file.relative_to(ROOT))] = hashlib.sha256(file.read_bytes()).hexdigest()
    for file in sorted(repo_path(GENERATOR_DIR).glob('*.py')):
        code[str(file)] = hashlib.sha256(file.read_bytes()).hexdigest()
    scientific = {k:v for k,v in config.items() if k not in non_scientific_configurations}
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

    if check_existing_run(ws_root, identity):
        logger.info('Existing run verified: %s', ws_root)
    else:
        write_json(
            guard_path,
            dict(fingerprint=identity, config=config, dataset=rows),
        )
        logger.info('New run initialized: %s', ws_root)

    return ws_root, rows, identity


def check_existing_run(ws_root: Path, identity: str) -> bool:
    """Validate an existing run's fingerprint without modifying files.

    Returns True if run metadata exists and matches, otherwise False.
    Raises ValueError if metadata is invalid or the fingerprint differs.
    """
    guard_path = ws_root / 'run_metadata.json'
    if not guard_path.exists():
        return False

    try:
        metadata = json.loads(guard_path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ValueError(f'Invalid run metadata: {guard_path}') from exc

    if not isinstance(metadata, dict) or 'fingerprint' not in metadata:
        raise ValueError(f'Missing fingerprint in run metadata: {guard_path}')

    if metadata['fingerprint'] != identity:
        raise ValueError(
            f'Run inputs/code/config changed for {ws_root}. '
            'Choose a new run_name.'
        )

    return True