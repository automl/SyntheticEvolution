"""One trial: per-RNA flow, archives, records, logs, failures and resuming.

The generator, AF3 and DSSR are replaced by fakes that write what the real
tools write; no tool installation is needed.
"""
import json
import subprocess
import tarfile
from pathlib import Path
from unittest.mock import Mock

import pytest

import hpo.pipeline.run as run
import hpo.pipeline.trial as trial


@pytest.fixture
def trial_case(tmp_path, monkeypatch):
    """A two-RNA trial whose fake tools write what the real tools write."""
    config = {
        "fixed": {"N": 20},
        "shs_python": "/test/shs/bin/python",
        "shs_seed": 11,
        "af3_seed": 2,
        "model_dir": str(tmp_path / "models"),
        "generator_timeout_seconds": 30,
        "af3_timeout_seconds": 60,
        "dssr": "test-dssr",
        "dssr_timeout_seconds": 10,
        "storage": {"scratch": str(tmp_path / "scratch")},
        "module": "bio/alphafold/3.0.1",
    }

    rows = [
        {"id": "paired", "sequence": "ACGU", "pairs": [[0, 3]],
         "interactions": [[0, 3, 0.3]], "mutation_rates": [0.2] * 4},
        {"id": "unpaired", "sequence": "AAAA", "pairs": []},
    ]

    directory = tmp_path / "trial with spaces"
    requests = []

    def generate(args, **kwargs):
        assert args[:2] == [
            config["shs_python"],
            str((run.GENERATOR_DIR / "shs_generator.py").resolve()),
        ]
        assert args[2] == "--request-json"
        assert args[4] == "--output-json"
        assert kwargs["check"] is True
        assert kwargs["timeout"] == 30
        assert kwargs["stderr"] == subprocess.STDOUT
        assert not kwargs["stdout"].closed

        request_path = Path(args[3])
        output_path = Path(args[5])
        # Generator inputs and outputs live in scratch, not in the workspace.
        assert Path(config["storage"]["scratch"]) in request_path.parents
        request = json.loads(request_path.read_text())
        requests.append(request)
        kwargs["stdout"].write(f"generated {request['task_name']}\n")
        output_path.write_text(json.dumps({"name": request["task_name"], "msa": "synthetic"}))
        request_path.with_name("shs_effective_parameters.json").write_text('{"N": 20}')

    def run_af3(af3_input_path, output, log_path, cfg):
        assert cfg is config
        af3_input = json.loads(af3_input_path.read_text())
        name = af3_input["name"]
        prediction = output / name
        (prediction / "seed-42_sample-0").mkdir(parents=True)
        model = prediction / f"{name}_model.cif"
        model.write_text(f"model of {name}")
        (prediction / f"{name}_data.json").write_text(json.dumps(af3_input))
        (prediction / "seed-42_sample-0" / "model.cif").write_text("sample")
        log_path.write_text(f"af3 ran {name}\n")
        return model

    def fake_dssr(model, directory, sequence, exe, timeout, scratch):
        assert (exe, timeout) == ("test-dssr", 10)
        (scratch / "dssr.log").write_text(f"dssr ran {model.name}\n")
        (directory / "dssr.json").write_text(json.dumps({"model": model.name}))
        # RNA 1 is perfect; RNA 2 has one spurious pair: mean loss 0.5.
        return [(0, 3)]

    generation = Mock(side_effect=generate)
    af3 = Mock(side_effect=run_af3)
    dssr = Mock(side_effect=fake_dssr)

    monkeypatch.setattr(trial.subprocess, "run", generation)
    monkeypatch.setattr(trial, "run_af3", af3)
    monkeypatch.setattr(trial.dssr, "run_dssr", dssr)

    return config, rows, directory, requests, generation, af3, dssr


def run_pipeline(case, **kwargs):
    config, rows, directory, *_ = case
    return trial.run_pipeline(config, rows, directory, 'identity', {'mutation_rate_paired': 0.3}, **kwargs)


def members(archive):
    with tarfile.open(archive) as tar:
        return sorted(tar.getnames())


def member_text(archive, name):
    with tarfile.open(archive) as tar:
        return tar.extractfile(name).read().decode()


def workspace_files(directory):
    return sorted(path.relative_to(directory).as_posix() for path in directory.rglob('*'))


def scratch_root(case):
    return Path(case[0]['storage']['scratch'])


def read_failure_report(directory):
    return json.loads((directory / 'failure_report.json').read_text())


def fail_dssr_on_second(case, failure):
    dssr = case[6]
    score = dssr.side_effect
    def fail_on_second(model, directory, sequence, exe, timeout, scratch):
        if model.name == 'rna_00001_model.cif':
            (scratch / 'dssr.log').write_text('segfault soon\n')
            raise failure
        return score(model, directory, sequence, exe, timeout, scratch)
    dssr.side_effect = fail_on_second


# ---- Successful trial -------------------------------------------------------

def test_fresh_trial_publishes_one_archive_per_rna(trial_case):
    config, rows, directory, requests, generation, af3, dssr = trial_case
    assert run_pipeline(trial_case) == pytest.approx(0.5)

    # The workspace holds only the archives and the result.
    assert workspace_files(directory) == ['rna_00000.tar.gz', 'rna_00001.tar.gz', 'trial_result.json']
    assert members(directory / 'rna_00000.tar.gz') == [
        'dssr.json', 'rna_00000_data.json', 'rna_00000_model.cif',
        'seed-42_sample-0/model.cif', 'shs_effective_parameters.json',
    ]
    assert json.loads(member_text(directory / 'rna_00001.tar.gz', 'dssr.json')) == {'model': 'rna_00001_model.cif'}
    assert json.loads(member_text(directory / 'rna_00000.tar.gz', 'rna_00000_data.json'))['msa'] == 'synthetic'

    for index, request in enumerate(requests):
        assert request == dict(rows[index], task_name=f'rna_{index:05d}',
                               parameters={'N': 20, 'mutation_rate_paired': 0.3},
                               shs_seed=11, af3_seed=2)
    assert (generation.call_count, af3.call_count, dssr.call_count) == (2, 2, 2)
    # DSSR reads the model from scratch.
    for index, call in enumerate(dssr.call_args_list):
        model = call.args[0]
        assert model.name == f'rna_{index:05d}_model.cif'
        assert scratch_root(trial_case) in model.parents
        assert call.args[2] == rows[index]['sequence']

    result = json.loads((directory / 'trial_result.json').read_text())
    # The result replaces trial_metadata.json and carries the trial signature.
    assert (result['fingerprint'], result['parameters']) == ('identity', requests[0]['parameters'])
    assert [s['id'] for s in result['scores']] == ['paired', 'unpaired']
    assert [s['loss'] for s in result['scores']] == [0, 1]
    # Optional dataset columns are recorded with the scores.
    assert (result['scores'][0]['interactions'], result['scores'][0]['mutation_rates']) == (
        rows[0]['interactions'], rows[0]['mutation_rates'])
    assert 'interactions' not in result['scores'][1]
    assert [s['model'] for s in result['scores']] == [
        'rna_00000.tar.gz/rna_00000_model.cif', 'rna_00001.tar.gz/rna_00001_model.cif',
    ]
    archive_bytes = sum((directory / f'rna_{i:05d}.tar.gz').stat().st_size for i in range(2))
    assert result['storage'] == {'files': 2, 'bytes': archive_bytes}
    assert result['loss'] == 0.5
    assert list(scratch_root(trial_case).iterdir()) == []

    # Cached completion performs no further external work.
    assert run_pipeline(trial_case) == 0.5
    assert (generation.call_count, af3.call_count, dssr.call_count) == (2, 2, 2)


def test_changed_identity_rejects_cached_result(trial_case):
    run_pipeline(trial_case)
    with pytest.raises(ValueError):
        trial.run_pipeline(trial_case[0], trial_case[1], trial_case[2], 'different',
                           {'mutation_rate_paired': 0.3})
    with pytest.raises(ValueError):
        trial.run_pipeline(trial_case[0], trial_case[1], trial_case[2], 'identity',
                           {'mutation_rate_paired': 0.4})


def test_without_fixed_or_sampled_parameters_the_generator_defaults_apply(trial_case):
    config, rows, directory, requests, *_ = trial_case
    del config['fixed']
    assert trial.run_pipeline(config, rows, directory, 'identity') == pytest.approx(0.5)
    assert [request['parameters'] for request in requests] == [{}, {}]


@pytest.mark.parametrize('fails', [False, True], ids=['success', 'failure'])
def test_standalone_trial_keeps_logs(trial_case, fails):
    directory = trial_case[2]
    if fails:
        fail_dssr_on_second(trial_case, RuntimeError('DSSR failed'))
    try:
        run_pipeline(trial_case, keep_logs_in_trial=True)
    except RuntimeError:
        assert fails
    logs = [name for name in workspace_files(directory) if name.endswith('.log')]
    assert logs == [f'logs/rna_{i:05d}/{tool}.log' for i in range(2) for tool in ('af3', 'dssr', 'shs_generation')]
    assert (directory / 'logs/rna_00000/shs_generation.log').read_text() == 'generated rna_00000\n'
    assert list(scratch_root(trial_case).iterdir()) == []


# ---- Failures ---------------------------------------------------------------

@pytest.mark.parametrize('failure', [subprocess.CalledProcessError(1, 'generator'),
                                     subprocess.TimeoutExpired('generator', 30)])
def test_generation_failure_never_submits(trial_case, failure):
    trial_case[4].side_effect = failure
    with pytest.raises(type(failure)):
        run_pipeline(trial_case)
    trial_case[5].assert_not_called()
    assert workspace_files(trial_case[2]) == ['failure_report.json', 'trial_metadata.json']
    report = read_failure_report(trial_case[2])
    assert (report['stage'], report['rna_index'], report['rna_id']) == ('generation', 0, 'paired')
    assert report['timed_out'] is isinstance(failure, subprocess.TimeoutExpired)
    assert report['storage'] == {'files': 1, 'bytes': (trial_case[2] / 'trial_metadata.json').stat().st_size}


def test_dssr_crash_keeps_published_archives_and_reports(trial_case):
    directory = trial_case[2]
    dssr = trial_case[6]
    def crash_on_second(model, directory, sequence, exe, timeout, scratch):
        if model.name == 'rna_00000_model.cif':
            (directory / 'dssr.json').write_text('{}')
            return [(0, 3)]
        (scratch / 'dssr.log').write_text('dssr output\n')
        raise subprocess.CalledProcessError(-11, ['/opt/x3dna-dssr', f'-i={model}', '--json'])
    dssr.side_effect = crash_on_second
    with pytest.raises(subprocess.CalledProcessError):
        run_pipeline(trial_case)

    # AF3 output of both RNAs survives: the archive is published before DSSR.
    assert workspace_files(directory) == [
        'failure_report.json', 'rna_00000.tar.gz', 'rna_00001.tar.gz', 'trial_metadata.json',
    ]
    assert 'dssr.json' in members(directory / 'rna_00000.tar.gz')
    assert 'dssr.json' not in members(directory / 'rna_00001.tar.gz')
    assert 'rna_00001_model.cif' in members(directory / 'rna_00001.tar.gz')

    report = read_failure_report(directory)
    assert (report['stage'], report['rna_index'], report['rna_id']) == ('dssr', 1, 'unpaired')
    assert (report['error_type'], report['returncode'], report['signal']) == ('CalledProcessError', -11, 'SIGSEGV')
    assert report['timed_out'] is False
    # Scratch paths are marked; only run_metadata.json holds absolute paths.
    assert report['command'] == ['/opt/x3dna-dssr', '-i=<scratch>/rna_00001/archive/rna_00001_model.cif', '--json']
    assert str(directory) not in report['message']
    assert report['logs'] == {
        'shs_generation.log': 'generated rna_00001\n', 'af3.log': 'af3 ran rna_00001\n',
        'dssr.log': 'dssr output\n',
    }
    assert 'Traceback' in report['traceback']
    assert {'host', 'slurm_job_id', 'slurm_array_task_id', 'time'} <= report.keys()
    assert report['storage']['files'] == 3
    assert list(scratch_root(trial_case).iterdir()) == []


def test_rejected_dssr_output_is_kept_in_the_archive(trial_case):
    directory = trial_case[2]
    dssr = trial_case[6]
    score = dssr.side_effect
    dssr.side_effect = lambda *args, **kwargs: (score(*args, **kwargs), [(0, 99)])[1]
    with pytest.raises(ValueError):
        run_pipeline(trial_case)
    assert json.loads(member_text(directory / 'rna_00000.tar.gz', 'dssr.json')) == {'model': 'rna_00000_model.cif'}
    assert read_failure_report(directory)['stage'] == 'dssr'


def test_af3_without_valid_prediction_publishes_nothing(trial_case):
    trial_case[5].side_effect = RuntimeError('Expected exactly one nonempty AF3 model')
    with pytest.raises(RuntimeError):
        run_pipeline(trial_case)
    assert not list(trial_case[2].glob('*.tar.gz'))
    assert read_failure_report(trial_case[2])['stage'] == 'af3'


def test_af3_failure_names_the_rna_and_keeps_the_log_tail(trial_case):
    af3 = trial_case[5]
    finish_af3 = af3.side_effect
    def fail_second(af3_input_path, output, log_path, cfg):
        if 'rna_00001' in str(af3_input_path):
            log_path.write_text('x' * 70000 + 'CUDA error')
            raise subprocess.CalledProcessError(1, ['af3'])
        return finish_af3(af3_input_path, output, log_path, cfg)
    af3.side_effect = fail_second
    with pytest.raises(subprocess.CalledProcessError):
        run_pipeline(trial_case)
    report = read_failure_report(trial_case[2])
    assert (report['stage'], report['rna_id'], report['signal']) == ('af3', 'unpaired', None)
    log = report['logs']['af3.log']
    # Only the end of a long log is kept, where the error is.
    assert log.startswith('[... first 4474 bytes truncated ...]') and log.endswith('CUDA error')
    # The trial stops at the first failing RNA; the first archive stays.
    assert [path.name for path in trial_case[2].glob('*.tar.gz')] == ['rna_00000.tar.gz']


def test_interruption_is_reported_like_a_failure(trial_case):
    # BaseExceptions (Ctrl-C, Slurm stops) also get a report and a clean scratch.
    fail_dssr_on_second(trial_case, KeyboardInterrupt())
    with pytest.raises(KeyboardInterrupt):
        run_pipeline(trial_case)
    report = read_failure_report(trial_case[2])
    assert (report['error_type'], report['logs']['dssr.log']) == ('KeyboardInterrupt', 'segfault soon\n')
    assert list(scratch_root(trial_case).iterdir()) == []


# ---- Resuming ---------------------------------------------------------------

def test_resume_reuses_archives_and_replaces_the_failure_record(trial_case):
    config, rows, directory, requests, generation, af3, dssr = trial_case
    score = dssr.side_effect
    fail_dssr_on_second(trial_case, RuntimeError('DSSR failed'))
    with pytest.raises(RuntimeError):
        run_pipeline(trial_case)
    assert (generation.call_count, af3.call_count) == (2, 2)

    dssr.side_effect = score
    assert run_pipeline(trial_case) == pytest.approx(0.5)
    # Both AF3 predictions came from the archives.
    assert (generation.call_count, af3.call_count) == (2, 2)
    assert workspace_files(directory) == ['rna_00000.tar.gz', 'rna_00001.tar.gz', 'trial_result.json']
    assert 'dssr.json' in members(directory / 'rna_00001.tar.gz')


def test_unreadable_archive_and_stale_temporaries_are_replaced(trial_case):
    directory = trial_case[2]
    directory.mkdir()
    (directory / 'rna_00000.tar.gz').write_bytes(b'truncated')
    (directory / '.rna_00001.tar.gz.node.123.tmp').write_bytes(b'partial')
    assert run_pipeline(trial_case) == pytest.approx(0.5)
    assert trial_case[4].call_count == 2
    assert 'rna_00000_model.cif' in members(directory / 'rna_00000.tar.gz')
    assert workspace_files(directory) == ['rna_00000.tar.gz', 'rna_00001.tar.gz', 'trial_result.json']


# ---- Scratch and AF3 --------------------------------------------------------

def test_trial_scratch_is_private_and_removed(tmp_path, monkeypatch):
    monkeypatch.setenv('SLURM_JOB_ID', '4242')
    with trial.trial_scratch(tmp_path / 'root') as first, trial.trial_scratch(tmp_path / 'root') as second:
        assert first != second
        assert first.parent == tmp_path / 'root'
        assert first.name.startswith('shs_hpo_4242_')
        (first / 'rna_00000').mkdir()
        (first / 'rna_00000' / 'af3.log').write_text('log')
    assert list((tmp_path / 'root').iterdir()) == []


@pytest.fixture
def af3_environment(tmp_path, monkeypatch):
    for name, value in (('AF3_PYTHON', '/af3/bin/python'), ('ALPHAFOLD_BIN_DIR', '/af3/bin'),
                        ('ALPHAFOLD_DATABASES', '/af3/db')):
        monkeypatch.setenv(name, value)
    return dict(model_dir=str(tmp_path / 'models'), af3_timeout_seconds=10)


def test_af3_runs_with_the_module_environment_and_returns_the_model(tmp_path, monkeypatch, af3_environment):
    output, log = tmp_path / 'af3', tmp_path / 'af3.log'
    def predict(args, **kwargs):
        assert args == ['/af3/bin/python', '/af3/bin/run_alphafold.py', f'--json_path={tmp_path / "in.json"}',
                        f'--output_dir={output}', f'--model_dir={tmp_path / "models"}',
                        '--db_dir=/af3/db', '--norun_data_pipeline']
        assert (kwargs['check'], kwargs['timeout'], kwargs['stdout'].name) == (True, 10, str(log))
        kwargs['stdout'].write('af3 output\n')
        for path in ('rna_00000/rna_00000_model.cif', 'rna_00000/seed-42_sample-0/model.cif'):
            (output / path).parent.mkdir(parents=True, exist_ok=True)
            (output / path).write_text('model')
    monkeypatch.setattr(trial.subprocess, 'run', predict)
    model = trial.run_af3(tmp_path / 'in.json', output, log, af3_environment)
    assert model == output / 'rna_00000' / 'rna_00000_model.cif'
    assert log.read_text() == 'af3 output\n'


@pytest.mark.parametrize('models', [[], ['a/a_model.cif', 'b/b_model.cif'], ['a/a_model.cif:empty']])
def test_af3_needs_exactly_one_nonempty_model(tmp_path, monkeypatch, af3_environment, models):
    output = tmp_path / 'af3'
    def predict(args, **kwargs):
        for entry in models:
            path, _, empty = entry.partition(':')
            (output / path).parent.mkdir(parents=True, exist_ok=True)
            (output / path).write_text('' if empty else 'model')
    monkeypatch.setattr(trial.subprocess, 'run', predict)
    with pytest.raises(RuntimeError, match='exactly one nonempty AF3 model'):
        trial.run_af3(tmp_path / 'in.json', output, tmp_path / 'af3.log', af3_environment)


def test_af3_requires_the_module_environment(tmp_path, monkeypatch, af3_environment):
    monkeypatch.delenv('ALPHAFOLD_DATABASES')
    with pytest.raises(RuntimeError, match='ALPHAFOLD_DATABASES not set'):
        trial.run_af3(tmp_path / 'in.json', tmp_path / 'af3', tmp_path / 'af3.log', af3_environment)
