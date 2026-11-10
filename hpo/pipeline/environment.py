"""Check before submitting that everything the controller job uses exists.

A wrong interpreter path or module name otherwise only fails after the job
has waited in the GPU queue. Every problem is collected and reported at once.
"""

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

NEPS_PACKAGE = 'neural-pipeline-search'
NEPS_VERSION = '0.16.0'

# Runs in the checked interpreter: which modules are importable (without
# importing them) and which package versions are installed.
_PROBE = '''
import importlib.metadata, importlib.util, json, sys
def version(name):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None
modules, packages = sys.argv[1].split(','), sys.argv[2].split(',')
print(json.dumps(dict(
    missing=[name for name in modules if name and importlib.util.find_spec(name) is None],
    versions={name: version(name) for name in packages if name},
)))
'''


def _check_executable(key: str, value: Any) -> list[str]:
    path = Path(str(value)).expanduser()
    if not path.is_file():
        return [f'{key}: {value} does not exist']
    if not os.access(path, os.X_OK):
        return [f'{key}: {value} is not executable']
    return []


def _check_interpreter(key: str, python: str, modules: list[str], packages: dict[str, str]) -> list[str]:
    """Start the interpreter once and check its modules and package versions."""
    try:
        result = subprocess.run(
            [str(Path(python).expanduser()), '-c', _PROBE, ','.join(modules), ','.join(packages)],
            capture_output=True, text=True, timeout=120, check=True,
        )
        found = json.loads(result.stdout)
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        return [f'{key}: {python} does not start: {error}']
    problems = [f'{key}: module {name} is not installed in {python}' for name in found['missing']]
    for name, expected in packages.items():
        installed = found['versions'].get(name)
        if installed != expected:
            problems.append(f'{key}: {name} {expected} required in {python}, found {installed or "none"}')
    return problems


def _check_module(module: str) -> list[str]:
    """Find the environment module's file (Lmod: <MODULEPATH entry>/<name>[.lua])."""
    roots = [Path(entry) for entry in os.environ.get('MODULEPATH', '').split(':') if entry]
    if not roots:
        return [f'module: MODULEPATH is not set, cannot find {module}; submit from a Helix login node']
    if any((root / module).exists() or (root / f'{module}.lua').exists() for root in roots):
        return []
    return [f'module: {module} not found in MODULEPATH']


def check_environment(config: dict[str, Any], mode: str) -> list[str]:
    """Problems that would make the controller job fail; empty if none.

    Every generator with a ``python`` configured is checked, and DSSR if
    ``extraction.dssr`` is set. Protein extraction needs biopython and numpy in
    the controller environment.
    """
    controller_key = 'neps_python' if mode == 'neps' else 'controller_python'
    generators = {
        chain_type: settings['python']
        for chain_type, settings in config['generators'].items()
        if settings.get('python')
    }
    problems = _check_executable(controller_key, config[controller_key])
    for chain_type, python in generators.items():
        problems += _check_executable(f'generators.{chain_type}.python', python)

    dssr = config['extraction'].get('dssr')
    if dssr and shutil.which(str(Path(dssr).expanduser())) is None:
        problems.append(f'extraction.dssr: {dssr} is neither an executable path nor on PATH')

    model_dir = Path(config['model_dir']).expanduser()
    if not model_dir.is_dir() or not any(model_dir.iterdir()):
        problems.append(f'model_dir: {config["model_dir"]} is missing or empty')

    problems += _check_module(config['module'])

    # Import checks only for interpreters that exist.
    if not any(problem.startswith(f'{controller_key}:') for problem in problems):
        modules = ['yaml']
        if config.get('scoring_objective', 'f1') == 'wl':
            modules += ['numpy', 'grakel']
        if 'protein' in generators:
            modules += ['numpy', 'Bio']
        packages = {NEPS_PACKAGE: NEPS_VERSION} if mode == 'neps' else {}
        problems += _check_interpreter(
            controller_key, config[controller_key], list(dict.fromkeys(modules)), packages
        )
    for chain_type, python in generators.items():
        key = f'generators.{chain_type}.python'
        if not any(problem.startswith(f'{key}:') for problem in problems):
            problems += _check_interpreter(key, python, ['numpy'], {})
    return problems


def require_environment(config: dict[str, Any], mode: str) -> None:
    """Stop the submission with a list of all problems found."""
    problems = check_environment(config, mode)
    if problems:
        raise SystemExit(
            f'Submission cancelled: {len(problems)} problem(s) with the environment:\n'
            + '\n'.join(f'  - {problem}' for problem in problems)
        )
