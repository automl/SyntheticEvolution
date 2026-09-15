"""AF3 worker output validation and scheduler failure handling, without a GPU."""
import json
import subprocess
import sys
from unittest.mock import Mock

import pytest
from hpo.pipeline import af3_task, trial


@pytest.mark.parametrize('outcome', ['success', 'missing', 'empty', 'multiple', 'failed'])
def test_worker_selects_array_task_and_publishes_only_valid_output(tmp_path, monkeypatch, outcome):
    tasks = [dict(input=str(tmp_path / f'{i}.json'), output=str(tmp_path / f'out{i}'),
                  model_dir='/test/models') for i in range(2)]
    manifest = tmp_path / 'tasks.json'
    manifest.write_text(json.dumps(tasks))
    monkeypatch.setattr(sys, 'argv', ['af3_task.py', str(manifest)])
    monkeypatch.setenv('SLURM_ARRAY_TASK_ID', '1')
    monkeypatch.setenv('ALPHAFOLD_BIN_DIR', '/test/af3')
    monkeypatch.setenv('ALPHAFOLD_DATABASES', '/test/databases')
    output = tmp_path / 'out1'
    def execute(args, **kwargs):
        assert args == [sys.executable, '/test/af3/run_alphafold.py',
                        '--json_path=' + tasks[1]['input'], '--output_dir=' + tasks[1]['output'],
                        '--model_dir=/test/models', '--db_dir=/test/databases', '--norun_data_pipeline']
        assert kwargs == {'check': True}
        if outcome == 'failed':
            raise subprocess.CalledProcessError(1, args)
        if outcome != 'missing':
            folder = output / 'rna'
            folder.mkdir()
            (folder / 'rna_model.cif').write_text('' if outcome == 'empty' else 'model')
            # Per-sample models must not be mistaken for the selected top model.
            sample = folder / 'seed-1_sample-0'
            sample.mkdir()
            (sample / 'model.cif').write_text('sample')
            if outcome == 'multiple':
                (folder / 'other_model.cif').write_text('another')
    fake = Mock(side_effect=execute)
    monkeypatch.setattr(af3_task.subprocess, 'run', fake)
    if outcome == 'success':
        af3_task.main()
        marker = json.loads((output / 'af3_success.json').read_text())
        assert marker == {'model': str((output / 'rna/rna_model.cif').resolve())}
        assert not (output / 'af3_success.tmp').exists()
    else:
        with pytest.raises((RuntimeError, subprocess.CalledProcessError)):
            af3_task.main()
        assert not (output / 'af3_success.json').exists()
    fake.assert_called_once()
    assert not (tmp_path / 'out0').exists()


@pytest.mark.parametrize('state,exit_code', [('FAILED', '1:0'), ('CANCELLED', '0:15'),
    ('TIMEOUT', '0:0'), ('OUT_OF_MEMORY', '0:0'), ('NODE_FAIL', '0:0'),
    ('PREEMPTED', '0:0'), ('BOOT_FAIL', '0:0'), ('DEADLINE', '0:0'), ('COMPLETED', '1:0')])
def test_wait_rejects_failed_array_element(monkeypatch, state, exit_code):
    monkeypatch.setattr(trial, 'command', Mock(return_value=f'42_0|COMPLETED|0:0\n42_1|{state}|{exit_code}\n'))
    sleep = Mock(side_effect=AssertionError('must fail immediately'))
    monkeypatch.setattr(trial.time, 'sleep', sleep)
    with pytest.raises(RuntimeError):
        trial.wait_job('42', 2, dict(wait_timeout_seconds=10, poll_seconds=1))


def test_wait_requires_every_array_element_not_parent_or_batch(monkeypatch):
    command = Mock(side_effect=[
        '42|COMPLETED|0:0\n42_0|COMPLETED|0:0\n42_1.batch|COMPLETED|0:0',
        '42_0|COMPLETED|0:0\n42_1|COMPLETED|0:0'])
    monkeypatch.setattr(trial, 'command', command)
    monkeypatch.setattr(trial.time, 'monotonic', lambda: 0)
    sleep = Mock()
    monkeypatch.setattr(trial.time, 'sleep', sleep)
    trial.wait_job('42', 2, dict(wait_timeout_seconds=10, poll_seconds=1))
    assert command.call_count == 2
    sleep.assert_called_once_with(1)


def test_wait_timeout_does_not_sleep_in_real_time(monkeypatch):
    monkeypatch.setattr(trial, 'command', Mock(return_value=''))
    monkeypatch.setattr(trial.time, 'monotonic', Mock(side_effect=[0, 11]))
    monkeypatch.setattr(trial.time, 'sleep', Mock(side_effect=AssertionError('unexpected sleep')))
    with pytest.raises(TimeoutError):
        trial.wait_job('42', 2, dict(wait_timeout_seconds=10, poll_seconds=1))
