"""Inspect sbatch commands without invoking Slurm; include shell-sensitive paths."""
import importlib
import shlex
import sys
from unittest.mock import Mock

import pytest


@pytest.mark.parametrize('mode,python_key,job_name', [
    ('neps', 'neps_python', 'shs-hpo'), ('standalone', 'shs_python', 'shs-standalone')])
def test_controller_submission_resources_and_wrapped_command(tmp_path, monkeypatch, capsys, mode, python_key, job_name):
    suffix = 'neps_hpo' if mode == 'neps' else 'standalone_trial'
    module = importlib.import_module('hpo.submit_' + suffix)
    config_path = tmp_path / "config space 'quote' $literal.yaml"
    root = tmp_path / 'repo with spaces'
    workspace = tmp_path / 'workspace'
    config = dict(neps_python='/env/neps space/bin/python', shs_python='/env/shs space/bin/python',
                  controller=dict(partition='cpu-single', cpus=2, memory='5G', time='03:00:00'))
    load = Mock(return_value=config)
    command = Mock(return_value='123')
    monkeypatch.setattr(module, 'load_config', load)
    monkeypatch.setattr(module, 'get_workspace_path', lambda cfg: workspace)
    monkeypatch.setattr(module, 'repo_path', lambda path: config_path)
    monkeypatch.setattr(module, 'ROOT', root)
    monkeypatch.setattr(module, 'command', command)
    monkeypatch.setattr(sys, 'argv', ['submit', '--config', str(config_path)])
    module.main()
    load.assert_called_once_with(str(config_path), mode=mode, submitting_controller=True)
    command.assert_called_once()
    args = command.call_args.args[0]
    assert args[:2] == ['sbatch', '--parsable']
    for flag in ['--job-name=' + job_name, '--chdir=' + str(root), '--nodes=1', '--ntasks=1',
                 '--partition=cpu-single', '--cpus-per-task=2', '--mem=5G', '--time=03:00:00',
                 '--output=' + str(workspace / 'controller_logs/%j.out'),
                 '--error=' + str(workspace / 'controller_logs/%j.err')]:
        assert flag in args
    wrapped = next(arg[len('--wrap='):] for arg in args if arg.startswith('--wrap='))
    assert shlex.split(wrapped) == ['exec', config[python_key], '-u', '-m', 'hpo.run_' + suffix,
                                    '--config', str(config_path)]
    assert (workspace / 'controller_logs').is_dir()
    assert capsys.readouterr().out.strip() == '123'
