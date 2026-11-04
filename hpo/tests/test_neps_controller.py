"""The NePS controller with a fake NePS module: handoff to the trial and run summary.

The failure limit and the stop signals are tested in ``test_failures.py``
and ``test_signals.py``.
"""
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
                  evaluations=2, ignore_errors=False, optimizer='random_search',
                  failures=dict(min_trials=5, max_rate=0.2, window=50))
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
    # Signal handlers and NePS's trial state are covered by their own tests.
    monkeypatch.setattr(entry.signals, 'install', lambda: None)
    monkeypatch.setattr(entry, 'current_trial', lambda: SimpleNamespace(id='1', metadata=SimpleNamespace(time_sampled=1.0)))
    monkeypatch.setattr(entry, 'check_failure_limit', lambda *args: None)
    monkeypatch.setattr(entry, 'job_time_left', lambda: None)  # not in a Slurm job
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


# ---- Handoff to the trial ---------------------------------------------------

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


def test_neps_failure_reaches_neps_which_scores_it_with_the_worst_loss(neps_case, monkeypatch):
    _, fake, _, root, _ = neps_case
    monkeypatch.setattr(entry, 'run_pipeline', Mock(side_effect=RuntimeError('trial failed')))
    def search(**kwargs):
        assert kwargs['objective_value_on_error'] == 1.0
        kwargs['evaluate_pipeline'](root / 'evaluation', N=2)
    fake.run.side_effect = search
    with pytest.raises(RuntimeError, match='trial failed'):
        entry.main()


def test_wrong_neps_version_stops_before_loading_run(neps_case, monkeypatch):
    _, fake, _, _, load = neps_case
    monkeypatch.setattr(entry, 'version', lambda name: '0.15.0')
    with pytest.raises(RuntimeError, match='0.16.0'):
        entry.main()
    load.assert_not_called()
    fake.run.assert_not_called()


# ---- Run summary ------------------------------------------------------------

def fake_neps_state(monkeypatch, states):
    """Fake NePS state; ``states`` maps trial IDs to a state or (state, error message)."""
    from enum import Enum
    State = Enum('State', {name.upper(): name for name in ('success', 'failed', 'crashed', 'evaluating')})
    trials = {}
    for trial_id, state in states.items():
        state, err = state if isinstance(state, tuple) else (state, None)
        report = SimpleNamespace(err=None if err is None else Exception(err))
        trials[trial_id] = SimpleNamespace(
            metadata=SimpleNamespace(state=State[state.upper()]), report=report)
    loaded = SimpleNamespace(lock_and_read_trials=lambda: trials)
    NePSState = SimpleNamespace(create_or_load=Mock(return_value=loaded))
    monkeypatch.setitem(sys.modules, 'neps.state', SimpleNamespace(NePSState=NePSState, State=State))
    return NePSState


def test_run_summary_reports_failures_from_reports(tmp_path, monkeypatch, caplog):
    import json
    NePSState = fake_neps_state(monkeypatch, {
        '1': 'success', '2': 'crashed', '3': 'crashed', '4': 'evaluating',
        '5': ('crashed', '[stopped] failure limit: ...'),
    })
    report = tmp_path / 'configs' / 'config_2' / 'artifacts' / 'failure_report.json'
    report.parent.mkdir(parents=True)
    report.write_text(json.dumps(dict(error_type='CalledProcessError', signal='SIGSEGV',
                                      stage='dssr', rna_id='rna_b')))
    entry.log_run_summary(tmp_path)
    NePSState.create_or_load.assert_called_once_with(path=tmp_path, load_only=True)
    [record] = [r for r in caplog.records if 'Run summary' in r.getMessage()]
    assert record.levelname == 'ERROR'
    assert record.getMessage() == ('Run summary: 2 of 3 trials failed (66.7%): '
                                   '1x SIGSEGV in dssr on rna_b; 1x unknown (no failure report) in unknown stage'
                                   '. 1 trial(s) stopped by the controller.')


def test_run_summary_without_failures_is_info(tmp_path, monkeypatch, caplog):
    fake_neps_state(monkeypatch, {'1': 'success'})
    caplog.set_level('INFO')
    entry.log_run_summary(tmp_path)
    [record] = [r for r in caplog.records if 'Run summary' in r.getMessage()]
    assert (record.levelname, record.getMessage()) == ('INFO', 'Run summary: No failed trials.')


def test_run_summary_never_breaks_the_controller(tmp_path, monkeypatch, caplog):
    monkeypatch.setitem(sys.modules, 'neps.state', None)  # import fails
    entry.log_run_summary(tmp_path)
    assert 'Could not summarize' in caplog.text


def test_run_summary_is_logged_when_neps_fails(neps_case, monkeypatch):
    _, fake, _, root, _ = neps_case
    summary = Mock()
    monkeypatch.setattr(entry, 'log_run_summary', summary)
    fake.run.side_effect = RuntimeError('worker failed')
    with pytest.raises(RuntimeError):
        entry.main()
    summary.assert_called_once_with(root / 'neps')


# ---- Retention (trajectory mode) --------------------------------------------

def finished_run(tmp_path, monkeypatch, losses):
    """NePS state with successful trials sampled at time = id; all have archives."""
    import json
    from enum import Enum
    State = Enum('State', {'SUCCESS': 'success', 'CRASHED': 'crashed', 'EVALUATING': 'evaluating'})
    trials = {}
    for trial_id, loss in losses.items():
        state = State.SUCCESS if loss is not None else State.CRASHED
        trials[trial_id] = SimpleNamespace(
            metadata=SimpleNamespace(state=state, time_sampled=float(trial_id)),
            report=SimpleNamespace(objective_to_minimize=loss, err=None))
        artifacts = tmp_path / 'configs' / f'config_{trial_id}' / 'artifacts'
        artifacts.mkdir(parents=True)
        (artifacts / 'rna_00000.tar.gz').write_bytes(b'archive')
        if loss is not None:
            (artifacts / 'trial_result.json').write_text(json.dumps(dict(loss=loss)))
    loaded = SimpleNamespace(lock_and_read_trials=lambda: trials)
    monkeypatch.setitem(sys.modules, 'neps.state', SimpleNamespace(
        NePSState=SimpleNamespace(create_or_load=Mock(return_value=loaded)), State=State))


def archives(tmp_path):
    return sorted(path.parent.parent.name for path in tmp_path.glob('configs/*/artifacts/*.tar.gz'))


def retention_record(tmp_path, trial_id):
    import json
    return json.loads((tmp_path / 'configs' / f'config_{trial_id}' / 'artifacts' / 'trial_result.json').read_text())['retention']


def test_retention_prunes_a_trial_that_is_not_kept(tmp_path, monkeypatch):
    # Trials 1-3 are on the trajectory; trial 4 (crashed) is never pruned.
    finished_run(tmp_path, monkeypatch, {'1': .5, '2': .4, '3': .3, '4': None, '5': .9})
    entry.apply_retention(tmp_path, '5', 5.0, .9, top_k=0)
    assert archives(tmp_path) == ['config_1', 'config_2', 'config_3', 'config_4']
    assert retention_record(tmp_path, '5') == dict(kept=False, reasons=[])


def test_retention_prunes_configs_pushed_off_by_an_earlier_sampled_trial(tmp_path, monkeypatch):
    # Trial 2 was sampled before 3 but finishes later, with a better loss.
    finished_run(tmp_path, monkeypatch, {'1': .5, '2': .2, '3': .3})
    entry.apply_retention(tmp_path, '2', 2.0, .2, top_k=0)
    assert archives(tmp_path) == ['config_1', 'config_2']
    assert retention_record(tmp_path, '3') == dict(kept=False, reasons=[], pruned_after='2')
    assert retention_record(tmp_path, '2') == dict(kept=True, reasons=['trajectory'])


@pytest.mark.parametrize('mode', ['all', 'trajectory'])
def test_retention_runs_only_in_trajectory_mode_and_never_fails_a_trial(neps_case, monkeypatch, mode):
    config, fake, _, root, _ = neps_case
    config['storage'] = {'keep_artifacts': mode}
    apply = Mock(side_effect=RuntimeError('broken state'))
    monkeypatch.setattr(entry, 'apply_retention', apply)
    fake.run.side_effect = lambda **kw: kw['evaluate_pipeline'](root / 'config_1', N=2)
    entry.main()
    if mode == 'trajectory':
        apply.assert_called_once_with(root / 'neps', '1', 1.0, 0.25, 5)
    else:
        apply.assert_not_called()


# ---- Time limit -------------------------------------------------------------

@pytest.mark.parametrize('text,seconds', [
    ('1-02:03:04', 93784), ('02:03:04', 7384), ('23:59', 1439), ('UNLIMITED', None), ('', None),
])
def test_slurm_time_left_is_parsed(text, seconds):
    assert entry.parse_slurm_duration(text) == seconds


def test_no_trial_starts_that_cannot_finish_before_the_stop_signal(neps_case, monkeypatch):
    config, fake, calls, root, _ = neps_case
    del config['evaluations']  # no budget: sample until the time is up
    clock = [0.0]
    monkeypatch.setattr(entry.time, 'monotonic', lambda: clock[0])
    # 1000 s before the stop signal: room for two 500 s trials.
    monkeypatch.setattr(entry, 'job_time_left', lambda: 1000 + entry.signals.STOP_MARGIN_SECONDS)
    def trial(*args, **kwargs):
        clock[0] += 500
        return 0.25
    monkeypatch.setattr(entry, 'run_pipeline', trial)
    fake.save_pipeline_results = Mock()
    def search(**kwargs):
        assert kwargs['evaluations_to_spend'] is None
        for index in range(1, 10):
            kwargs['evaluate_pipeline'](root / f'config_{index}', N=index)
    fake.run.side_effect = search
    with pytest.raises(SystemExit) as stop:
        entry.main()
    assert stop.value.code == 0
    # Two trials ran; the third was sampled but stopped before it started.
    assert clock[0] == 1000
    [call] = fake.save_pipeline_results.call_args_list
    result = call.args[0]
    assert str(result['exception']).startswith('[stopped] not enough time left')
    assert result['objective_to_minimize'] is None
