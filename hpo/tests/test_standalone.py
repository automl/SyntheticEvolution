"""CLI behavior beyond --help; no external process is started."""
import json
import sys
from unittest.mock import Mock

import pytest

from hpo import run_standalone_trial as entry


def test_validate_only_does_not_prepare_workspace_or_run_trial(monkeypatch, capsys):
    rows = [dict(id='rna', sequence='ACGU', pairs=[])]
    monkeypatch.setattr(sys, 'argv', ['standalone', '--config', 'input.yaml', '--validate-only'])
    load = Mock(return_value={'input': 'data.csv'})
    dataset = Mock(return_value=rows)
    prepare = Mock(side_effect=AssertionError('validation must not prepare workspace'))
    run = Mock(side_effect=AssertionError('validation must not run trial'))
    monkeypatch.setattr(entry, 'load_config', load)
    monkeypatch.setattr(entry, 'create_dataset', dataset)
    monkeypatch.setattr(entry, 'prepare_run', prepare)
    monkeypatch.setattr(entry, 'run_pipeline', run)
    entry.main()
    load.assert_called_once_with('input.yaml', mode='standalone')
    dataset.assert_called_once_with('data.csv')
    assert json.loads(capsys.readouterr().out) == rows
    prepare.assert_not_called()
    run.assert_not_called()


def test_standalone_passes_prepared_dataset_and_prints_loss(tmp_path, monkeypatch, capsys):
    config = {'optimize': 'rna', 'generators': {'rna': {'fixed': {'N': 20}}}}
    rows = [dict(id='rna', sequence='ACGU', pairs=[])]
    monkeypatch.setattr(sys, 'argv', ['standalone', '--config', 'input.yaml'])
    monkeypatch.setattr(entry, 'load_config', Mock(return_value=config))
    prepare = Mock(return_value=(tmp_path, rows, 'identity'))
    run = Mock(return_value=0.25)
    monkeypatch.setattr(entry, 'prepare_run', prepare)
    monkeypatch.setattr(entry, 'run_pipeline', run)
    monkeypatch.setattr(entry.signals, 'install', lambda: None)
    entry.main()
    prepare.assert_called_once_with(config)
    # Standalone trials keep their logs; file count does not matter for one trial.
    run.assert_called_once_with(config, rows, tmp_path / 'standalone_trial', 'identity', None,
                                keep_logs_in_trial=True)
    assert capsys.readouterr().out.strip() == 'standalone_trial 0.25'


def run_trials(tmp_path, monkeypatch, run):
    config = {'fixed': {'N': 20}, 'shs_seed': 42, 'af3_seed': 42, 'trials': [
        {'name': 'best', 'base': 7, 'parameters': {'N': 8}},
        {'name': 'best_s1', 'shs_seed': 1, 'af3_seed': 2, 'parameters': {'N': 8}},
        {'name': 'nomsa', 'af3_seed': 3, 'parameters': {'N': 1}},
    ]}
    rows = [dict(id='rna', sequence='ACGU', pairs=[])]
    monkeypatch.setattr(sys, 'argv', ['standalone', '--config', 'input.yaml'])
    monkeypatch.setattr(entry, 'load_config', Mock(return_value=config))
    monkeypatch.setattr(entry, 'prepare_run', Mock(return_value=(tmp_path, rows, 'identity')))
    monkeypatch.setattr(entry, 'run_pipeline', run)
    monkeypatch.setattr(entry.signals, 'install', lambda: None)
    entry.main()


def test_each_trial_runs_in_its_directory_with_its_parameters_and_seeds(tmp_path, monkeypatch, capsys):
    run = Mock(side_effect=[0.1, 0.2, 0.3])
    run_trials(tmp_path, monkeypatch, run)
    calls = [(c.args[0]['shs_seed'], c.args[0]['af3_seed'], c.args[2], c.args[4], c.kwargs)
             for c in run.call_args_list]
    # Seeds default to the run's; logs of many trials are not kept.
    assert calls == [
        (42, 42, tmp_path / 'trials' / 'best', {'N': 8}, {'keep_logs_in_trial': False}),
        (1, 2, tmp_path / 'trials' / 'best_s1', {'N': 8}, {'keep_logs_in_trial': False}),
        (42, 3, tmp_path / 'trials' / 'nomsa', {'N': 1}, {'keep_logs_in_trial': False}),
    ]
    assert capsys.readouterr().out.split('\n')[:3] == ['best 0.1', 'best_s1 0.2', 'nomsa 0.3']


def test_a_failed_trial_does_not_stop_the_others(tmp_path, monkeypatch):
    run = Mock(side_effect=[0.1, RuntimeError('dssr crashed'), 0.3])
    with pytest.raises(SystemExit, match='Failed trials: best_s1'):
        run_trials(tmp_path, monkeypatch, run)
    assert run.call_count == 3


def test_a_stop_signal_ends_the_run(tmp_path, monkeypatch):
    run = Mock(side_effect=[0.1, entry.signals.ControllerStop('SIGUSR1'), 0.3])
    with pytest.raises(SystemExit) as stop:
        run_trials(tmp_path, monkeypatch, run)
    assert (stop.value.code, run.call_count) == (entry.signals.EXIT_CODE, 2)


def test_a_signal_between_trials_stops_before_the_next(tmp_path, monkeypatch):
    run = Mock(return_value=0.1)
    monkeypatch.setattr(entry.signals, 'requested', Mock(side_effect=[None, 'SIGUSR1', 'SIGUSR1']))
    with pytest.raises(SystemExit) as stop:
        run_trials(tmp_path, monkeypatch, run)
    assert (stop.value.code, run.call_count) == (entry.signals.EXIT_CODE, 1)
