"""Shared helpers for a whole run and run validation."""

import csv
import hashlib
import json
import subprocess
import math
import re
import yaml
from typing import Any
from pathlib import Path
import tempfile

from hpo.pipeline.signals import STOP_MARGIN_SECONDS

import logging
logger = logging.getLogger(__name__)

PIPELINE_DIR = Path(__file__).resolve().parent  # SyntheticEvolution/hpo/pipeline
HPO_DIR = PIPELINE_DIR.parent                   # SyntheticEvolution/hpo
ROOT = HPO_DIR.parent                           # SyntheticEvolution
GENERATOR_DIR = ROOT / 'SHS-Generator-RNA'

Dataset = list[dict[str, Any]]

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

    buckets = config.get('af3_buckets')
    if buckets is not None and (
        not isinstance(buckets, list) or not buckets
        or any(isinstance(b, bool) or not isinstance(b, int) or b <= 0 for b in buckets)
        or buckets != sorted(set(buckets))
    ):
        raise ValueError('af3_buckets must be a nonempty, strictly increasing list of positive integers')

    storage_settings(config)

    scoring_objective = config.get('scoring_objective', 'f1')
    if scoring_objective not in {'f1', 'f1_soft', 'wl'}:
        raise ValueError("scoring_objective must be 'f1', 'f1_soft' or 'wl'")

    if 'wl_iterations' in config:
        require_positive(config, 'wl_iterations', integer=True)

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
        # Without evaluations, each controller samples until its time is up.
        if config.get('evaluations') is not None:
            require_positive(config, 'evaluations', integer=True)
        if config.get('optimizer') is not None:
            require_text(config, 'optimizer')
        if not isinstance(config.get('ignore_errors'), bool):
            raise ValueError('ignore_errors must be a YAML boolean: true or false')
        config['failures'] = failure_limit_settings(config.get('failures'))

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

        # Slurm asks the controller to stop STOP_MARGIN_SECONDS before the
        # limit; a shorter job would be stopped before it does any work.
        seconds = slurm_seconds(controller['time'])
        if seconds is None:
            raise ValueError(
                "controller.time must be a Slurm time limit such as '30', "
                "'00:30:00' or '1-00:00:00'"
            )
        if seconds < STOP_MARGIN_SECONDS + MIN_WORK_SECONDS:
            raise ValueError(
                f'controller.time must be at least '
                f'{(STOP_MARGIN_SECONDS + MIN_WORK_SECONDS) // 60} minutes: the last '
                f'{STOP_MARGIN_SECONDS} s are reserved for a clean stop'
            )
    return config


# Minimum wall time left for trials after the stop margin.
MIN_WORK_SECONDS = 300


def slurm_seconds(value):
    """Seconds in a Slurm time limit, math.inf for no limit, None if invalid.

    Formats as in ``sbatch --time``: minutes, minutes:seconds,
    hours:minutes:seconds, days-hours, days-hours:minutes and
    days-hours:minutes:seconds.
    """
    value = value.strip()
    if value.lower() in {'infinite', 'unlimited'}:
        return math.inf
    match = re.fullmatch(r'(?:(\d+)-)?(\d+)(?::(\d+))?(?::(\d+))?', value)
    if match is None:
        return None
    days, *parts = match.groups()
    numbers = [int(part) for part in parts if part is not None]
    if days is None:
        units = {1: (60,), 2: (60, 1), 3: (3600, 60, 1)}[len(numbers)]
    else:
        units = (3600, 60, 1)
    return int(days or 0) * 86400 + sum(n * unit for n, unit in zip(numbers, units))

STORAGE_DEFAULTS = dict(scratch=None, keep_artifacts='all', keep_top_k=5)


def storage_settings(config):
    """Validate the optional storage section and return it with defaults.

    ``keep_artifacts``: ``all`` keeps the archives of every config (the
    default, so a forgotten setting never loses data); ``trajectory`` keeps
    them only for configs on NePS's best-config trajectory or among the
    ``keep_top_k`` lowest losses. Failed trials always keep their archives.
    """
    storage = config.get('storage')
    if storage is None:
        storage = {}
    if not isinstance(storage, dict):
        raise ValueError('storage must be a YAML mapping')
    unknown = storage.keys() - STORAGE_DEFAULTS.keys()
    if unknown:
        raise ValueError('Unknown storage settings: ' + ', '.join(sorted(map(str, unknown))))
    settings = {**STORAGE_DEFAULTS, **{k: v for k, v in storage.items() if v is not None}}
    scratch = settings['scratch']
    if scratch is not None and (not isinstance(scratch, str) or not Path(scratch).is_absolute()):
        raise ValueError('storage.scratch must be an absolute path')
    if settings['keep_artifacts'] not in ('all', 'trajectory'):
        raise ValueError("storage.keep_artifacts must be 'all' or 'trajectory'")
    top_k = settings['keep_top_k']
    if isinstance(top_k, bool) or not isinstance(top_k, int) or top_k < 0:
        raise ValueError('storage.keep_top_k must be a nonnegative integer')
    return settings


def scratch_root(config) -> Path:
    """Directory under which each trial gets its private scratch directory.

    Defaults to /dev/shm (RAM): on Helix, $TMPDIR is on GPFS and would put the
    same metadata load on the shared filesystem as the workspace.
    """
    scratch = storage_settings(config).get('scratch')
    if scratch is not None:
        return Path(scratch)
    shm = Path('/dev/shm')
    return shm if shm.is_dir() else Path(tempfile.gettempdir())


FAILURE_LIMIT_DEFAULTS = dict(min_trials=5, max_rate=0.20, window=50)


def failure_limit_settings(raw):
    """Validate the NePS failure limit and fill in defaults.

    ``max_rate: 1.0`` disables the limit; ``window: null`` uses all trials.
    """
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ValueError('failures must be a YAML mapping')
    unknown = raw.keys() - FAILURE_LIMIT_DEFAULTS.keys()
    if unknown:
        raise ValueError('Unknown failures settings: ' + ', '.join(sorted(map(str, unknown))))
    settings = {**FAILURE_LIMIT_DEFAULTS, **raw}

    def positive_integer(value):
        return isinstance(value, int) and not isinstance(value, bool) and value > 0

    if not positive_integer(settings['min_trials']):
        raise ValueError('failures.min_trials must be a positive integer')
    if settings['window'] is not None and not positive_integer(settings['window']):
        raise ValueError('failures.window must be a positive integer or null')
    rate = settings['max_rate']
    if isinstance(rate, bool) or not isinstance(rate, (int, float)) or not 0 < rate <= 1:
        raise ValueError('failures.max_rate must be a number in (0, 1]; 1.0 disables the limit')
    settings['max_rate'] = float(rate)
    return settings


def write_json(path, value):
    """Atomically publish JSON using a unique temporary file per writer."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode='w',
            encoding='utf-8',
            dir=path.parent,
            prefix=f'.{path.name}.',
            suffix='.tmp',
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            json.dump(value, stream, indent=2, sort_keys=True)

        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


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
        if not isinstance(pair, (list, tuple)) or len(pair) != 2 or any(type(x) is not int for x in pair):
            raise ValueError('Pairs must contain exactly two integer positions')
        a, b = pair
        if not (0 <= a < size and 0 <= b < size) or a == b:
            raise ValueError(f'Invalid pair {pair} for length {size}')
        result.add(tuple(sorted((a, b))))
    return result


def normalize_interactions(raw, size):
    """Validate and retain strengths, including explicit zero-strength entries."""
    if not isinstance(raw, list):
        raise ValueError('interactions must be a JSON list')
    result = {}
    for entry in raw:
        if not isinstance(entry, list) or len(entry) not in (2, 3):
            raise ValueError('Interactions require two indices and optional strength')
        key = next(iter(normalize_pairs([entry[:2]], size)))
        strength = entry[2] if len(entry) == 3 else 1.0
        validate_probability(strength, 'interaction strength')
        if key in result and result[key] != strength:
            raise ValueError(f'Conflicting interaction strengths for pair {key}')
        result[key] = float(strength)
    return [[a, b, strength] for (a, b), strength in sorted(result.items())]


def validate_probability(value, name):
    # Keep submission validation independent of the SHS environment (NumPy).
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f'{name} must be a finite number between 0 and 1')


def shs_generation_request(row, index, parameters, shs_seed, af3_seed):
    """The generator's input for one RNA; the exporter rebuilds it the same way."""
    return dict(
        **row,   # contains (id=, sequence=, pairs=)
        task_name=f'rna_{index:05d}',
        parameters=parameters,
        shs_seed=shs_seed,
        af3_seed=af3_seed,
    )


def create_dataset(path) -> Dataset:
    """Create a dataset from a path.
    ### Input:
    @path : single file or dir of files.
    ### returns
    dataset = list of dicts for each line with

        id : identifying string
        sequence : RNA sequence str containing only ACGU
        pairs : sorted scoring targets; derived from positive interactions if absent
        interactions, mutation_rates : optional validated generator inputs
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
                item = dict(id=name, sequence=seq)
                interaction_text = (row.get('interactions') or '').strip()
                pair_text = (row.get('pairs') or '').strip()
                if interaction_text:
                    item['interactions'] = normalize_interactions(json.loads(interaction_text), len(seq))
                if pair_text:
                    raw_pairs = json.loads(pair_text)
                    if not isinstance(raw_pairs, list):
                        raise ValueError('pairs must be a JSON list')
                    pairs = sorted(normalize_pairs(raw_pairs, len(seq)))
                elif interaction_text:
                    pairs = [(a, b) for a, b, strength in item['interactions'] if strength > 0]
                else:
                    raise ValueError(f'Provide pairs or interactions for {name}; use [] for no pairs')
                item['pairs'] = pairs
                rates_text = (row.get('mutation_rates') or '').strip()
                if rates_text:
                    rates = json.loads(rates_text)
                    if not isinstance(rates, list) or len(rates) != len(seq):
                        raise ValueError('mutation_rates must have exactly one rate per RNA position')
                    for rate in rates:
                        validate_probability(rate, 'mutation rate')
                    item['mutation_rates'] = rates
                rows.append(item)
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
        'storage', 'failures',
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