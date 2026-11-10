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
    return dict(controller_python=sys.executable, neps_python=sys.executable,
                generators=dict(rna=dict(python=sys.executable)),
                extraction=dict(dssr='x3dna-dssr'),
                model_dir=str(tmp_path / 'models'), module='bio/alphafold/3.0.1')


def test_working_environment_has_no_problems(setup):
    assert environment.check_environment(setup, 'standalone') == []


def test_every_problem_is_reported_at_once(setup, tmp_path, monkeypatch):
    # The error from the measured run: a mixed-up interpreter path.
    setup['controller_python'] = '/home/user/.conda/controller/synEvo/bin/python'
    setup['generators']['rna']['python'] = executable(tmp_path / 'not-executable', '')
    (tmp_path / 'not-executable').chmod(0o644)
    setup['extraction']['dssr'] = 'missing-dssr'
    setup['model_dir'] = str(tmp_path / 'empty')
    (tmp_path / 'empty').mkdir()
    setup['module'] = 'bio/alphafold/9.9.9'
    problems = environment.check_environment(setup, 'standalone')
    assert problems == [
        'controller_python: /home/user/.conda/controller/synEvo/bin/python does not exist',
        f'generators.rna.python: {tmp_path / "not-executable"} is not executable',
        'extraction.dssr: missing-dssr is neither an executable path nor on PATH',
        f'model_dir: {tmp_path / "empty"} is missing or empty',
        'module: bio/alphafold/9.9.9 not found in MODULEPATH',
    ]
    with pytest.raises(SystemExit, match='5 problem'):
        environment.require_environment(setup, 'standalone')


def test_interpreter_without_the_needed_modules_is_reported(setup, tmp_path):
    # A Python that exists and starts, but lacks yaml/numpy (the -I flag hides site-packages).
    python = executable(tmp_path / 'isolated-python', f'#!/bin/sh\nexec {sys.executable} -I -S "$@"\n')
    setup['generators']['rna']['python'] = python
    assert environment.check_environment(setup, 'standalone') == [
        f'generators.rna.python: module numpy is not installed in {python}',
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
    assert found['generators.rna.python'] == (['numpy'], {})


def test_wrong_neps_version_is_reported(setup):
    problems = environment.check_environment(setup, 'neps')
    assert problems == [f'neps_python: neural-pipeline-search 0.16.0 required in {sys.executable}, found none']


def test_module_cannot_be_checked_outside_the_cluster(setup, monkeypatch):
    monkeypatch.delenv('MODULEPATH')
    assert environment.check_environment(setup, 'standalone') == [
        'module: MODULEPATH is not set, cannot find bio/alphafold/3.0.1; submit from a Helix login node',
    ]


def test_protein_generator_adds_its_interpreter_and_the_biopython_requirement(setup, monkeypatch):
    found = {}
    def probe(key, python, modules, packages):
        found[key] = modules
        return []
    monkeypatch.setattr(environment, '_check_interpreter', probe)
    setup['generators']['protein'] = dict(python=sys.executable)
    assert environment.check_environment(setup, 'standalone') == []
    assert found['generators.protein.python'] == ['numpy']
    assert found['controller_python'] == ['yaml', 'numpy', 'Bio']


def test_generator_without_python_and_missing_dssr_setting_are_skipped(setup):
    setup['generators']['protein'] = dict(python=None)
    del setup['extraction']['dssr']
    assert environment.check_environment(setup, 'standalone') == []


def test_missing_protein_interpreter_is_reported(setup, tmp_path):
    setup['generators']['protein'] = dict(python=str(tmp_path / 'nope'))
    assert environment.check_environment(setup, 'standalone') == [
        f'generators.protein.python: {tmp_path / "nope"} does not exist',
    ]
