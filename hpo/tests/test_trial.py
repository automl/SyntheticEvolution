"""One trial: per-RNA flow, archives, records, logs, failures and resuming.

The generator, AF3 and DSSR are replaced by fakes that write what the real
tools write; no tool installation is needed.
"""
import json
import os
import re
import signal
import subprocess
import sys
import tarfile
from pathlib import Path
from unittest.mock import Mock

import pytest

import hpo.pipeline.run as run
import hpo.pipeline.trial as trial
from hpo.pipeline import signals


def fake_prediction(name, input_path, output, log):
    """What AF3 writes for one finished input, including its log lines."""
    log.write(f"Processing fold input {name}\n")
    af3_input = json.loads(input_path.read_text())
    prediction = output / name
    (prediction / "seed-42_sample-0").mkdir(parents=True)
    (prediction / f"{name}_model.cif").write_text(f"model of {name}")
    (prediction / f"{name}_data.json").write_text(json.dumps(af3_input))
    (prediction / "seed-42_sample-0" / "model.cif").write_text("sample")
    log.write(f"Done processing fold input {name}.\n")


@pytest.fixture
def trial_case(tmp_path, monkeypatch):
    """A two-RNA trial whose fake tools write what the real tools write."""
    config = {
        "optimize": "rna",
        "generators": {"rna": {"python": "/test/shs/bin/python", "fixed": {"N": 20}}},
        "extraction": {"dssr": "test-dssr", "dssr_timeout_seconds": 10},
        "shs_seed": 11,
        "af3_seed": 2,
        "model_dir": str(tmp_path / "models"),
        "generator_timeout_seconds": 30,
        "af3_timeout_seconds": 60,
        "storage": {"scratch": str(tmp_path / "scratch")},
        "module": "bio/alphafold/3.0.1",
    }

    rows = [
        {"id": "paired", "type": "rna", "sequence": "ACGU", "pairs": [[0, 3]],
         "interactions": [[0, 3, 0.3]], "mutation_rates": [0.2] * 4},
        {"id": "unpaired", "type": "rna", "sequence": "AAAA", "pairs": []},
    ]

    directory = tmp_path / "trial with spaces"
    requests = []

    def generate(args, **kwargs):
        assert args[:2] == [
            config["generators"]["rna"]["python"],
            str((run.GENERATOR_DIRS["rna"] / "shs_generator.py").resolve()),
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

    def run_af3(af3_inputs, output, log_path, cfg):
        assert cfg is config
        with log_path.open("w") as log:
            for name, input_path in af3_inputs.items():
                fake_prediction(name, input_path, output, log)

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
    # The real extractor dispatch runs; only the DSSR binary wrapper is replaced.
    monkeypatch.setattr(trial.extraction, "run_dssr", dssr)

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
        # The generator gets neither the chain type nor the other row bookkeeping.
        assert request == dict({k: v for k, v in rows[index].items() if k != 'type'},
                               task_name=f'rna_{index:05d}',
                               parameters={'N': 20, 'mutation_rate_paired': 0.3},
                               shs_seed=11, af3_seed=2)
    # One AF3 process predicts all RNAs, so the model is compiled once.
    assert (generation.call_count, af3.call_count, dssr.call_count) == (2, 1, 2)
    assert list(af3.call_args.args[0]) == ['rna_00000', 'rna_00001']
    # DSSR reads the model from scratch.
    for index, call in enumerate(dssr.call_args_list):
        model = call.args[0]
        assert model.name == f'rna_{index:05d}_model.cif'
        assert scratch_root(trial_case) in model.parents
        assert call.args[2] == rows[index]['sequence']

    result = json.loads((directory / 'trial_result.json').read_text())
    # The result replaces trial_metadata.json and carries the trial signature.
    assert (result['fingerprint'], result['parameters'], result['shs_seed'], result['af3_seed']) == (
        'identity', {'rna': requests[0]['parameters']}, 11, 2)
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
    assert (generation.call_count, af3.call_count, dssr.call_count) == (2, 1, 2)


def test_changed_identity_rejects_cached_result(trial_case):
    run_pipeline(trial_case)
    with pytest.raises(ValueError):
        trial.run_pipeline(trial_case[0], trial_case[1], trial_case[2], 'different',
                           {'mutation_rate_paired': 0.3})
    with pytest.raises(ValueError):
        trial.run_pipeline(trial_case[0], trial_case[1], trial_case[2], 'identity',
                           {'mutation_rate_paired': 0.4})
    # Standalone trials may set their own seeds.
    with pytest.raises(ValueError):
        trial.run_pipeline(dict(trial_case[0], shs_seed=12), trial_case[1], trial_case[2], 'identity',
                           {'mutation_rate_paired': 0.3})


def test_without_fixed_or_sampled_parameters_the_generator_defaults_apply(trial_case):
    config, rows, directory, requests, *_ = trial_case
    del config['generators']['rna']['fixed']
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
    assert logs == ['logs/af3.log'] + [
        f'logs/rna_{i:05d}/{tool}.log' for i in range(2) for tool in ('dssr', 'shs_generation')]
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
    assert (report['stage'], report['rna_index'], report['rna_id']) == ('extraction', 1, 'unpaired')
    assert (report['error_type'], report['returncode'], report['signal']) == ('CalledProcessError', -11, 'SIGSEGV')
    assert report['timed_out'] is False
    # Scratch paths are marked; only run_metadata.json holds absolute paths.
    assert report['command'] == ['/opt/x3dna-dssr', '-i=<scratch>/rna_00001/archive/rna_00001_model.cif', '--json']
    assert str(directory) not in report['message']
    # The trial's AF3 log covers all RNAs.
    assert report['logs']['shs_generation.log'] == 'generated rna_00001\n'
    assert report['logs']['dssr.log'] == 'dssr output\n'
    assert report['logs']['af3.log'].endswith('Done processing fold input rna_00001.\n')
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
    assert read_failure_report(directory)['stage'] == 'extraction'


def fail_af3_on_second(case, failure, log_text=''):
    """AF3 finishes the first input, then fails while predicting the second."""
    def run_af3(af3_inputs, output, log_path, cfg):
        with log_path.open('w') as log:
            first, second = af3_inputs.items()
            fake_prediction(*first, output, log)
            log.write(f'Processing fold input {second[0]}\n{log_text}')
        raise failure
    case[5].side_effect = run_af3


@pytest.mark.parametrize('failure', [subprocess.CalledProcessError(1, ['af3']), KeyboardInterrupt()],
                         ids=['af3-error', 'interrupted'])
def test_af3_failure_keeps_finished_predictions_and_names_the_rna(trial_case, failure):
    fail_af3_on_second(trial_case, failure, log_text='x' * 70000 + 'CUDA error')
    with pytest.raises(type(failure)):
        run_pipeline(trial_case)
    # The first prediction is archived although the AF3 process failed.
    assert [path.name for path in trial_case[2].glob('*.tar.gz')] == ['rna_00000.tar.gz']
    report = read_failure_report(trial_case[2])
    assert (report['stage'], report['rna_id'], report['error_type']) == ('af3', 'unpaired', type(failure).__name__)
    log = report['logs']['af3.log']
    # Only the end of a long log is kept, where the error is.
    assert log.startswith('[... first ') and log.endswith('CUDA error')
    assert len(log) < 64 * 1024 + 100
    assert list(scratch_root(trial_case).iterdir()) == []


def test_af3_without_a_model_fails_that_rna(trial_case):
    def skip_second(af3_inputs, output, log_path, cfg):
        with log_path.open('w') as log:
            fake_prediction('rna_00000', af3_inputs['rna_00000'], output, log)
    trial_case[5].side_effect = skip_second
    with pytest.raises(RuntimeError, match='rna_00001'):
        run_pipeline(trial_case)
    assert [path.name for path in trial_case[2].glob('*.tar.gz')] == ['rna_00000.tar.gz']
    report = read_failure_report(trial_case[2])
    assert (report['stage'], report['rna_id']) == ('af3', 'unpaired')


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
    assert (generation.call_count, af3.call_count) == (2, 1)

    dssr.side_effect = score
    assert run_pipeline(trial_case) == pytest.approx(0.5)
    # Both AF3 predictions came from the archives.
    assert (generation.call_count, af3.call_count) == (2, 1)
    assert workspace_files(directory) == ['rna_00000.tar.gz', 'rna_00001.tar.gz', 'trial_result.json']
    assert 'dssr.json' in members(directory / 'rna_00001.tar.gz')


def test_resume_after_af3_failure_predicts_only_the_missing_rnas(trial_case):
    config, rows, directory, requests, generation, af3, dssr = trial_case
    run_af3 = af3.side_effect
    fail_af3_on_second(trial_case, subprocess.CalledProcessError(1, ['af3']))
    with pytest.raises(subprocess.CalledProcessError):
        run_pipeline(trial_case)
    af3.side_effect = run_af3
    assert run_pipeline(trial_case) == pytest.approx(0.5)
    assert list(af3.call_args.args[0]) == ['rna_00001']
    assert generation.call_count == 3  # rna_00001 is generated again, rna_00000 is not


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


# Stands in for af3_worker.py: same protocol, no AF3. Inputs named "*crash"
# make it exit, "*hang" makes it sleep.
FAKE_WORKER = r"""
import json, os, sys, time
from pathlib import Path
reply = os.fdopen(os.dup(1), 'w', buffering=1)
os.dup2(2, 1)
print('starting', sys.argv[1:])
for line in sys.stdin:
    request = json.loads(line)
    with open(request['log'], 'a') as log:
        os.dup2(log.fileno(), 1)
        os.dup2(log.fileno(), 2)
    print('worker', os.getpid())
    for path in sorted(Path(request['input_dir']).glob('*.json')):
        name = json.loads(path.read_text())['name']
        if name.endswith('crash'):
            sys.exit(3)
        if name.endswith('hang'):
            time.sleep(60)
        (Path(request['output_dir']) / name).mkdir(parents=True)
        (Path(request['output_dir']) / name / f'{name}_model.cif').write_text('model')
        print(f'Done processing fold input {name}.')
    reply.write('done\n')
"""


@pytest.fixture
def af3_environment(tmp_path, monkeypatch):
    worker = tmp_path / 'fake_worker.py'
    worker.write_text(FAKE_WORKER)
    monkeypatch.setattr(trial, 'AF3_WORKER', worker)
    monkeypatch.setenv('AF3_PYTHON', sys.executable)
    monkeypatch.setenv('ALPHAFOLD_BIN_DIR', '/af3/bin')
    yield dict(model_dir=str(tmp_path / 'models'), af3_timeout_seconds=10)
    trial.stop_af3_worker()


def predict(tmp_path, trial_name, names, config):
    """One trial's AF3 call in its own scratch directory; returns its log."""
    scratch = tmp_path / trial_name
    inputs = {}
    for name in names:
        inputs[name] = scratch / name / 'af3_input.json'
        inputs[name].parent.mkdir(parents=True)
        inputs[name].write_text(json.dumps({'name': name}))
    trial.run_af3(inputs, scratch / 'af3', scratch / 'af3.log', config)
    return (scratch / 'af3.log').read_text()


def worker_pid(log):
    return re.search(r'^worker (\d+)$', log, re.MULTILINE).group(1)


def test_one_af3_worker_serves_every_trial(tmp_path, af3_environment):
    first = predict(tmp_path, 'trial_1', ['rna_00000', 'rna_00001'], af3_environment)
    second = predict(tmp_path, 'trial_2', ['rna_00000'], af3_environment)
    # The trial that starts the worker also gets its start-up output.
    assert f"starting ['--model_dir={tmp_path / 'models'}']" in first
    assert 'starting' not in second
    assert worker_pid(first) == worker_pid(second)
    assert sorted(trial.af3_predictions(tmp_path / 'trial_1' / 'af3', tmp_path / 'trial_1' / 'af3.log')) \
        == ['rna_00000', 'rna_00001']


def test_changed_af3_settings_start_a_new_worker(tmp_path, af3_environment):
    first = predict(tmp_path, 'trial_1', ['rna_00000'], af3_environment)
    second = predict(tmp_path, 'trial_2', ['rna_00000'], dict(af3_environment, af3_buckets=[64, 128]))
    assert worker_pid(first) != worker_pid(second)
    assert "'--buckets=64,128'" in second


def test_crashed_af3_worker_fails_the_trial_and_is_replaced(tmp_path, af3_environment):
    with pytest.raises(subprocess.CalledProcessError) as error:
        predict(tmp_path, 'trial_1', ['rna_00000', 'rna_00001_crash'], af3_environment)
    assert error.value.returncode == 3
    # What AF3 finished before the crash is still in the log.
    assert 'Done processing fold input rna_00000.' in (tmp_path / 'trial_1' / 'af3.log').read_text()
    assert 'starting' in predict(tmp_path, 'trial_2', ['rna_00000'], af3_environment)


def test_af3_timeout_applies_per_input_and_kills_the_worker(tmp_path, af3_environment):
    config = dict(af3_environment, af3_timeout_seconds=0.25)
    with pytest.raises(subprocess.TimeoutExpired) as error:
        predict(tmp_path, 'trial_1', ['rna_00000', 'rna_00001_hang'], config)
    assert error.value.timeout == 0.5
    assert trial._af3_worker is None


def test_stop_signal_kills_the_af3_worker(tmp_path, af3_environment, monkeypatch):
    monkeypatch.setattr(signals, '_requested', None)
    previous = signal.signal(signal.SIGUSR1, signals._handle)
    signal.signal(signal.SIGALRM, lambda *args: os.kill(os.getpid(), signal.SIGUSR1))
    try:
        signal.setitimer(signal.ITIMER_REAL, 0.5)
        with pytest.raises(signals.ControllerStop), signals.trial_running():
            predict(tmp_path, 'trial_1', ['rna_00000', 'rna_00001_hang'], af3_environment)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, signal.SIG_DFL)
        signal.signal(signal.SIGUSR1, previous)
    log = (tmp_path / 'trial_1' / 'af3.log').read_text()
    assert 'Done processing fold input rna_00000.' in log
    assert trial._af3_worker is None
    with pytest.raises(ProcessLookupError):
        os.kill(int(worker_pid(log)), 0)


@pytest.mark.parametrize('models,done,found', [
    (['a/a_model.cif'], ['a'], ['a']),
    (['a/a_model.cif', 'b/b_model.cif'], ['a'], ['a']),             # b still running
    (['a/a_model.cif', 'a/seed-42_sample-0/model.cif'], ['a'], ['a']),
    (['a/a_model.cif:empty'], ['a'], []),
    (['a/a_model.cif', 'a/x_model.cif'], ['a'], []),
    ([], ['a'], []),
])
def test_af3_predictions_need_the_done_line_and_one_nonempty_model(tmp_path, models, done, found):
    output = tmp_path / 'af3'
    for entry in models:
        path, _, empty = entry.partition(':')
        (output / path).parent.mkdir(parents=True, exist_ok=True)
        (output / path).write_text('' if empty else 'model')
    log = tmp_path / 'af3.log'
    log.write_text(''.join(f'Processing fold input {n}\nDone processing fold input {n}.\n' for n in done))
    assert sorted(trial.af3_predictions(output, log)) == found


def test_af3_requires_the_module_environment(tmp_path, monkeypatch, af3_environment):
    monkeypatch.delenv('ALPHAFOLD_BIN_DIR')
    with pytest.raises(RuntimeError, match='ALPHAFOLD_BIN_DIR not set'):
        trial.run_af3({}, tmp_path / 'af3', tmp_path / 'af3.log', af3_environment)


# ---- Mixed chain types ------------------------------------------------------

def test_mixed_dataset_uses_each_chain_types_generator_parameters_and_extractor(trial_case, monkeypatch):
    """Only the optimized generator is sampled; the protein chain gets its own tools and archive file."""
    config, rows, directory, requests, generation, af3, dssr = trial_case
    config['generators']['protein'] = {'python': '/test/protein/bin/python', 'fixed': {'N': 5}}
    rows[1] = {'id': 'protein', 'type': 'protein', 'sequence': 'MKTA', 'pairs': [[0, 3]]}
    commands = []

    def generate(args, **kwargs):
        commands.append(args[:2])
        request = json.loads(Path(args[3]).read_text())
        requests.append(request)
        Path(args[5]).write_text(json.dumps({'name': request['task_name']}))
    generation.side_effect = generate
    contacts = Mock(return_value=[(0, 3)])
    monkeypatch.setattr(trial.extraction, 'contact_pairs', contacts)

    assert run_pipeline(trial_case) == 0

    assert commands == [
        ['/test/shs/bin/python', str(run.GENERATOR_DIRS['rna'] / 'shs_generator.py')],
        ['/test/protein/bin/python', str(run.GENERATOR_DIRS['protein'] / 'shs_generator.py')],
    ]
    assert requests[0]['parameters'] == {'N': 20, 'mutation_rate_paired': 0.3}
    assert requests[1]['parameters'] == {'N': 5}   # the protein generator is not optimized
    assert 'type' not in requests[1]
    # DSSR only for the RNA, contacts only for the protein; each file lands in its own archive.
    assert dssr.call_count == 1
    assert contacts.call_args.args[1] == 'MKTA'
    assert contacts.call_args.kwargs == {'atom': 'CA', 'cutoff': 8.0, 'min_separation': 6}
    assert 'dssr.json' in members(directory / 'rna_00000.tar.gz')
    assert 'contacts.json' not in members(directory / 'rna_00000.tar.gz')
    assert 'contacts.json' in members(directory / 'rna_00001.tar.gz')
    assert 'dssr.json' not in members(directory / 'rna_00001.tar.gz')
    result = json.loads((directory / 'trial_result.json').read_text())
    assert result['parameters'] == {'rna': requests[0]['parameters'], 'protein': {'N': 5}}


def test_parameters_are_stored_only_for_chain_types_in_the_dataset(trial_case):
    config, rows, directory, *_ = trial_case
    config['generators']['protein'] = {'python': '/test/protein/bin/python', 'fixed': {'N': 5}}
    run_pipeline(trial_case)
    result = json.loads((directory / 'trial_result.json').read_text())
    assert list(result['parameters']) == ['rna']
