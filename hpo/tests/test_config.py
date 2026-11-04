"""Configuration loading and the validation each entry point needs."""
import copy
from pathlib import Path

import pytest
import yaml

import hpo.pipeline.run as run
from hpo.pipeline.run import load_config

FIXTURE = Path(__file__).parent / 'fixtures' / 'valid_standalone_config.yaml'
CONTROLLER = dict(count=1, partition='gpu-single', cpus=1, memory='20G', gres='gpu:A100:1', time='04:00:00')


@pytest.fixture
def load(tmp_path, monkeypatch):
    """Write a config as YAML and load it like an entry point does."""
    # Resolve test configs locally, independently of the repository layout.
    monkeypatch.setattr(run, 'repo_path', lambda value: Path(value))
    path = tmp_path / 'config.yaml'

    def load(config, mode='standalone', **kwargs):
        path.write_text(yaml.safe_dump(config))
        return load_config(path, mode=mode, **kwargs)
    return load


@pytest.fixture
def standalone():
    return yaml.safe_load(FIXTURE.read_text())


@pytest.fixture
def neps(standalone):
    return dict(
        standalone,
        evaluations=2,
        optimizer='bayesian_optimization',
        search={'N': {'type': 'integer', 'lower': 2, 'upper': 128}},
        # Explicit, so loading returns the config unchanged.
        failures={'min_trials': 5, 'max_rate': 0.2, 'window': 50},
    )


# ---- Shared settings --------------------------------------------------------

@pytest.mark.parametrize('value', [None, ['not', 'a', 'mapping']])
def test_config_must_be_a_mapping(load, value):
    with pytest.raises(ValueError, match='YAML mapping'):
        load(value)


@pytest.mark.parametrize('value', ['missing', '', None])
def test_required_value_cannot_be_missing_empty_or_null(load, standalone, value):
    if value == 'missing':
        del standalone['dssr']
    else:
        standalone['dssr'] = value
    with pytest.raises(ValueError, match='dssr'):
        load(standalone)


@pytest.mark.parametrize('name', ['', '../other-run', 'run/name', 'run name', None])
def test_run_name_must_be_a_safe_folder_name(load, standalone, name):
    with pytest.raises(ValueError, match='run_name'):
        load(dict(standalone, run_name=name))


@pytest.mark.parametrize('mode', ['standalone', 'neps'])
def test_fixed_is_optional(load, standalone, neps, mode):
    config = neps if mode == 'neps' else standalone
    config.pop('fixed', None)
    assert load(config, mode=mode) == config


def test_scratch_is_optional_and_must_be_absolute(load, standalone):
    loaded = load(dict(standalone, storage={'scratch': '/dev/shm/custom'}))
    assert run.scratch_root(loaded) == Path('/dev/shm/custom')
    for storage in ({'scratch': 'relative/path'}, {'scratch': 3}, {'typo': '/x'}, [],
                    {'keep_artifacts': 'none'}, {'keep_top_k': -1}, {'keep_top_k': True}):
        with pytest.raises(ValueError, match='storage'):
            load(dict(standalone, storage=storage))


def test_storage_defaults_keep_every_config(standalone):
    # A forgotten setting must never cost a big run its archives.
    assert run.storage_settings(standalone) == dict(scratch=None, keep_artifacts='all', keep_top_k=5)
    assert run.storage_settings({'storage': {'scratch': None, 'keep_artifacts': 'trajectory'}}) == dict(
        scratch=None, keep_artifacts='trajectory', keep_top_k=5)


@pytest.mark.parametrize('shm_exists', [True, False])
def test_scratch_defaults_to_shared_memory(monkeypatch, shm_exists):
    monkeypatch.setattr(run.Path, 'is_dir', lambda self: shm_exists)
    expected = Path('/dev/shm') if shm_exists else Path(run.tempfile.gettempdir())
    assert run.scratch_root({}) == expected


# ---- Standalone -------------------------------------------------------------

def test_standalone_needs_no_neps_settings(load, standalone):
    assert load(standalone) == standalone


def test_standalone_warns_about_search_and_keeps_fixed(load, neps, caplog):
    neps['fixed'] = {'mutation_rate_paired': 0.2}
    loaded = load(neps)
    assert "ignoring 'search'" in caplog.text
    assert loaded['fixed'] == neps['fixed']


# ---- NePS -------------------------------------------------------------------

def test_neps_needs_no_submission_settings(load, neps):
    assert load(neps, mode='neps') == neps


@pytest.mark.parametrize('search', ['missing', {}])
def test_neps_requires_a_nonempty_search(load, neps, search):
    if search == 'missing':
        del neps['search']
    else:
        neps['search'] = search
    with pytest.raises(ValueError, match='search'):
        load(neps, mode='neps')


def test_af3_buckets_are_optional_and_strictly_increasing(load, standalone):
    assert 'af3_buckets' not in load(standalone)
    assert load(dict(standalone, af3_buckets=[64, 128]))['af3_buckets'] == [64, 128]
    for value in ([], [128, 64], [64, 64], [0, 64], [64.0], [True], '64,128'):
        with pytest.raises(ValueError, match='af3_buckets'):
            load(dict(standalone, af3_buckets=value))


def test_evaluations_are_optional_but_must_be_positive(load, neps):
    del neps['evaluations']
    assert 'evaluations' not in load(neps, mode='neps')
    for value in (0, -1, 1.5, True):
        with pytest.raises(ValueError, match='evaluations'):
            load(dict(neps, evaluations=value), mode='neps')


def test_parameter_cannot_be_both_fixed_and_searched(load, neps):
    neps['fixed'] = {'N': 20}
    with pytest.raises(ValueError, match='both fixed and search'):
        load(neps, mode='neps')


def test_failure_limit_defaults_and_partial_settings(load, neps):
    del neps['failures']
    assert load(neps, mode='neps')['failures'] == {'min_trials': 5, 'max_rate': 0.2, 'window': 50}
    neps['failures'] = {'max_rate': 1, 'window': None}
    assert load(neps, mode='neps')['failures'] == {'min_trials': 5, 'max_rate': 1.0, 'window': None}


@pytest.mark.parametrize('failures', [
    [], {'max_rate': 0}, {'max_rate': 1.5}, {'max_rate': True},
    {'window': 0}, {'min_trials': 2.5}, {'typo': 1},
])
def test_invalid_failure_limit_is_rejected(load, neps, failures):
    neps['failures'] = failures
    with pytest.raises(ValueError, match='failures'):
        load(neps, mode='neps')


# ---- Controller submission --------------------------------------------------

def test_submission_requires_the_controller_settings(load, neps):
    neps.pop('controller')
    with pytest.raises(ValueError, match='controller'):
        load(neps, mode='neps', submitting_controller=True)
    neps['controller'] = CONTROLLER
    with pytest.raises(ValueError, match='neps_python'):
        load(neps, mode='neps', submitting_controller=True)
    neps['neps_python'] = '/example/neps/bin/python'
    assert load(neps, mode='neps', submitting_controller=True) == neps


@pytest.mark.parametrize('key,value', [
    ('gres', 'missing'), ('gres', ''), ('gres', None), ('gres', 1),
    ('count', 0), ('count', -1), ('count', 1.5), ('count', True), ('count', '2'),
])
def test_invalid_controller_resources_are_rejected(load, neps, key, value):
    neps['neps_python'] = '/example/neps/bin/python'
    neps['controller'] = copy.deepcopy(CONTROLLER)
    if value == 'missing':
        del neps['controller'][key]
    else:
        neps['controller'][key] = value
    with pytest.raises(ValueError, match=rf'controller\.{key}'):
        load(neps, mode='neps', submitting_controller=True)


@pytest.mark.parametrize('time,seconds', [
    ('10', 600), ('10:30', 630), ('01:02:03', 3723), ('2-3', 183600),
    ('1-02:03', 93780), ('1-02:03:04', 93784), ('UNLIMITED', float('inf')),
    ('', None), ('1:2:3:4', None), ('10m', None), ('-5', None),
])
def test_slurm_time_limits_are_parsed(time, seconds):
    assert run.slurm_seconds(time) == seconds


@pytest.mark.parametrize('time,accepted', [
    ('00:07:00', True), ('7', True), ('00:06:59', False), ('00:02:00', False), ('ten', False),
])
def test_controller_time_must_leave_room_after_the_stop_margin(load, neps, time, accepted):
    neps['neps_python'] = '/example/neps/bin/python'
    neps['controller'] = dict(CONTROLLER, time=time)
    if accepted:
        load(neps, mode='neps', submitting_controller=True)
    else:
        with pytest.raises(ValueError, match=r'controller\.time'):
            load(neps, mode='neps', submitting_controller=True)


# ---- Export -----------------------------------------------------------------

def test_export_only_requires_the_run_location(load):
    config = {'run_name': 'test-run', 'workspace_name': 'alphafold'}
    assert load(config, mode='export') == config
