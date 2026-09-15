"""CLI behavior beyond --help; no external process is started."""
import json
import sys
from unittest.mock import Mock

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
    config = {'fixed': {'N': 20}}
    rows = [dict(id='rna', sequence='ACGU', pairs=[])]
    monkeypatch.setattr(sys, 'argv', ['standalone', '--config', 'input.yaml'])
    monkeypatch.setattr(entry, 'load_config', Mock(return_value=config))
    prepare = Mock(return_value=(tmp_path, rows, 'identity'))
    run = Mock(return_value=0.25)
    monkeypatch.setattr(entry, 'prepare_run', prepare)
    monkeypatch.setattr(entry, 'run_pipeline', run)
    entry.main()
    prepare.assert_called_once_with(config)
    run.assert_called_once_with(config, rows, tmp_path / 'standalone_trial', 'identity')
    assert capsys.readouterr().out.strip() == '0.25'
