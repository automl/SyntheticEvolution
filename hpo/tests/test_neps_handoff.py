"""Test our NePS adapter with a fake NePS module, not optimizer internals."""
import sys
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from hpo import run_neps_hpo as entry


@pytest.fixture
def neps_case(tmp_path, monkeypatch):
    config = dict(search={'N': dict(type='integer', lower=1, upper=30),
                          'rate': dict(type='float', lower=0.0, upper=1.0),
                          'approach': dict(type='categorical', choices=['a', 'b'])},
                  evaluations=2, ignore_errors=False, optimizer='random_search')
    rows = [dict(id='rna', sequence='ACGU', pairs=[[0, 3]])]
    fake = SimpleNamespace(Integer=Mock(return_value='integer-space'),
                           Float=Mock(return_value='float-space'),
                           Categorical=Mock(return_value='category-space'), run=Mock())
    monkeypatch.setitem(sys.modules, 'neps', fake)
    monkeypatch.setattr(sys, 'argv', ['run_neps_hpo', '--config', 'config.yaml'])
    monkeypatch.setattr(entry, 'version', lambda name: '0.16.0')
    load = Mock(return_value=config)
    monkeypatch.setattr(entry, 'load_config', load)
    monkeypatch.setattr(entry, 'prepare_run', lambda cfg: (tmp_path, rows, 'fingerprint'))
    # An explicit signature catches misplaced positional arguments as well as values.
    calls = []
    def pipeline(config_arg, rows_arg, trial_directory, identity, parameters=None):
        assert config_arg is config
        assert rows_arg is rows
        assert identity == 'fingerprint'
        calls.append((trial_directory, parameters))
        return 0.25
    monkeypatch.setattr(entry, 'run_pipeline', pipeline)
    return config, fake, calls, tmp_path, load


def test_neps_space_and_each_evaluation_reach_trial(neps_case):
    config, fake, calls, root, load = neps_case
    def search(**kwargs):
        assert kwargs['root_directory'] == root / 'neps'
        assert kwargs['evaluations_to_spend'] == 2
        assert kwargs['optimizer'] == 'random_search'
        assert kwargs['ignore_errors'] is False
        assert kwargs['pipeline_space'] == dict(N='integer-space', rate='float-space', approach='category-space')
        evaluate = kwargs['evaluate_pipeline']
        assert evaluate(root / 'evaluation_1', N=10, rate=0.1, approach='a') == 0.25
        assert evaluate(root / 'evaluation_2', N=20, rate=0.2, approach='b') == 0.25
    fake.run.side_effect = search
    entry.main()
    load.assert_called_once_with('config.yaml', mode='neps')
    fake.Integer.assert_called_once_with(lower=1, upper=30)
    fake.Float.assert_called_once_with(lower=0.0, upper=1.0)
    fake.Categorical.assert_called_once_with(choices=['a', 'b'])
    assert calls == [(root / 'evaluation_1/artifacts', dict(N=10, rate=0.1, approach='a')),
                     (root / 'evaluation_2/artifacts', dict(N=20, rate=0.2, approach='b'))]
    assert config['search']['N']['type'] == 'integer'  # Building space must not mutate YAML data.


def test_neps_failure_is_not_replaced_with_a_numeric_loss(neps_case, monkeypatch):
    _, fake, _, root, _ = neps_case
    monkeypatch.setattr(entry, 'run_pipeline', Mock(side_effect=RuntimeError('trial failed')))
    fake.run.side_effect = lambda **kw: kw['evaluate_pipeline'](root / 'evaluation', N=2)
    with pytest.raises(RuntimeError, match='trial failed'):
        entry.main()


def test_wrong_neps_version_stops_before_loading_run(neps_case, monkeypatch):
    _, fake, _, _, load = neps_case
    monkeypatch.setattr(entry, 'version', lambda name: '0.15.0')
    with pytest.raises(RuntimeError, match='0.16.0'):
        entry.main()
    load.assert_not_called()
    fake.run.assert_not_called()


def test_second_controller_cannot_start_same_search(neps_case):
    # A real local advisory lock verifies the duplicate-controller guard.
    _, fake, _, root, _ = neps_case
    with (root / 'controller.lock').open('w') as lock:
        entry.fcntl.flock(lock, entry.fcntl.LOCK_EX | entry.fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError):
            entry.main()
    fake.run.assert_not_called()
