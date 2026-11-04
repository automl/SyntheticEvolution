"""Environment check before submitting: interpreters, packages, DSSR, models, AF3 module."""
import sys

import pytest

from hpo.pipeline import environment


def executable(path, text='#!/bin/sh\n'):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(0o755)
    return str(path)


@pytest.fixture
def setup(tmp_path, monkeypatch):
    """A working environment; the interpreters are the Python running the tests."""
    (tmp_path / 'models' / 'af3.bin').parent.mkdir()
    (tmp_path / 'models' / 'af3.bin').write_text('weights')
    modules = tmp_path / 'modulefiles'
    (modules / 'bio' / 'alphafold').mkdir(parents=True)
    (modules / 'bio' / 'alphafold' / '3.0.1.lua').write_text('')
    monkeypatch.setenv('MODULEPATH', f'/nonexistent:{modules}')
    executable(tmp_path / 'bin' / 'x3dna-dssr')
    monkeypatch.setenv('PATH', f"{tmp_path / 'bin'}:/usr/bin:/bin")
    # The test interpreter has no NePS; the NePS mode test fakes the probe.
    return dict(controller_python=sys.executable, neps_python=sys.executable, shs_python=sys.executable,
                dssr='x3dna-dssr', model_dir=str(tmp_path / 'models'), module='bio/alphafold/3.0.1')


def test_working_environment_has_no_problems(setup):
    assert environment.check_environment(setup, 'standalone') == []


def test_every_problem_is_reported_at_once(setup, tmp_path, monkeypatch):
    # The error from the measured run: a mixed-up interpreter path.
    setup['controller_python'] = '/home/user/.conda/controller/synEvo/bin/python'
    setup['shs_python'] = executable(tmp_path / 'not-executable', '')
    (tmp_path / 'not-executable').chmod(0o644)
    setup['dssr'] = 'missing-dssr'
    setup['model_dir'] = str(tmp_path / 'empty')
    (tmp_path / 'empty').mkdir()
    setup['module'] = 'bio/alphafold/9.9.9'
    problems = environment.check_environment(setup, 'standalone')
    assert problems == [
        'controller_python: /home/user/.conda/controller/synEvo/bin/python does not exist',
        f'shs_python: {tmp_path / "not-executable"} is not executable',
        'dssr: missing-dssr is neither an executable path nor on PATH',
        f'model_dir: {tmp_path / "empty"} is missing or empty',
        'module: bio/alphafold/9.9.9 not found in MODULEPATH',
    ]
    with pytest.raises(SystemExit, match='5 problem'):
        environment.require_environment(setup, 'standalone')


def test_interpreter_without_the_needed_modules_is_reported(setup, tmp_path):
    # A Python that exists and starts, but lacks yaml/numpy (the -I flag hides site-packages).
    setup['shs_python'] = executable(tmp_path / 'isolated-python', f'#!/bin/sh\nexec {sys.executable} -I -S "$@"\n')
    assert environment.check_environment(setup, 'standalone') == [
        f'shs_python: module numpy is not installed in {setup["shs_python"]}',
    ]


def test_neps_mode_requires_the_pinned_neps_version(setup, monkeypatch):
    found = {}
    def probe(key, python, modules, packages):
        found[key] = (modules, packages)
        return []
    monkeypatch.setattr(environment, '_check_interpreter', probe)
    setup['scoring_objective'] = 'wl'
    assert environment.check_environment(setup, 'neps') == []
    assert found['neps_python'] == (['yaml', 'numpy', 'grakel'], {'neural-pipeline-search': '0.16.0'})
    assert found['shs_python'] == (['numpy'], {})


def test_wrong_neps_version_is_reported(setup):
    problems = environment.check_environment(setup, 'neps')
    assert problems == [f'neps_python: neural-pipeline-search 0.16.0 required in {sys.executable}, found none']


def test_module_cannot_be_checked_outside_the_cluster(setup, monkeypatch):
    monkeypatch.delenv('MODULEPATH')
    assert environment.check_environment(setup, 'standalone') == [
        'module: MODULEPATH is not set, cannot find bio/alphafold/3.0.1; submit from a Helix login node',
    ]
