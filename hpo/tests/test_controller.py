"""Trial orchestration tests: no generator, Slurm, AF3 or DSSR installation needed."""
import json
import subprocess
from pathlib import Path
import tempfile
import unittest
import pytest
from unittest.mock import Mock, patch

from hpo.pipeline import trial
from hpo.pipeline.dssr import parse_pairs


FIXTURE = json.loads(
	(Path(__file__).parent / 'fixtures' / 'dssr_example.json').read_text()
)
SEQ = 'AGUAGUAGUAGUAGA'


@pytest.fixture
def trial_case(tmp_path, monkeypatch):
    config = dict(fixed={'N': 20}, shs_python='/test/shs/bin/python',
                  shs_seed=11, af3_seed=2, model_dir=str(tmp_path / 'models'),
                  generator_timeout_seconds=30, dssr='test-dssr',
                  dssr_timeout_seconds=10, module='bio/alphafold/3.0.1',
                  gpu=dict(partition='gpu-single', cpus=8, memory='20G',
                           gres='gpu:A100:1', time='00:15:00', concurrency=1))
    rows = [dict(id='paired', sequence='ACGU', pairs=[[0, 3]]),
            dict(id='unpaired', sequence='AAAA', pairs=[])]
    directory = tmp_path / 'trial with spaces'
    requests = []

    def generate(args, **kwargs):
        assert args[:2] == [config['shs_python'], str((trial.GENERATOR_DIR / 'shs_generator.py').resolve())]
        assert args[2] == '--request-json' and args[4] == '--output-json'
        assert kwargs['check'] is True
        assert kwargs['timeout'] == 30
        assert kwargs['stderr'] == subprocess.STDOUT
        assert not kwargs['stdout'].closed
        request = json.loads(Path(args[3]).read_text())
        requests.append(request)
        Path(args[5]).write_text(json.dumps({'name': request['task_name']}))

    def finish(job, count, cfg):
        assert (job, count, cfg) == ('42', 2, config)
        for task in json.loads((directory / 'af3_tasks.json').read_text()):
            model = Path(task['output']) / 'prediction' / 'rna_model.cif'
            model.parent.mkdir(parents=True, exist_ok=True)
            model.write_text('test model')
            trial.write_json(Path(task['output']) / 'af3_success.json', {'model': str(model)})

    generation = Mock(side_effect=generate)
    submit = Mock(return_value='42;helix')
    wait = Mock(side_effect=finish)
    # RNA 1 is perfect; RNA 2 has one spurious pair. Mean loss must be 0.5.
    dssr = Mock(side_effect=lambda model, directory, sequence, exe, timeout: [(0, 3)])
    monkeypatch.setattr(trial.subprocess, 'run', generation)
    monkeypatch.setattr(trial, 'command', submit)
    monkeypatch.setattr(trial, 'wait_job', wait)
    monkeypatch.setattr(trial.dssr, 'run_dssr', dssr)
    return config, rows, directory, requests, generation, submit, wait, dssr


def run(case):
    config, rows, directory, *_ = case
    return trial.run_pipeline(config, rows, directory, 'identity', {'mutation_rate_paired': 0.3})


def test_fresh_trial_wires_two_rnas_and_averages_scores(trial_case):
    config, rows, directory, requests, generation, submit, wait, dssr = trial_case
    assert run(trial_case) == pytest.approx(0.5)
    assert generation.call_count == 2
    for index, request in enumerate(requests):
        assert request == dict(rows[index], task_name=f'rna_{index:05d}',
                               parameters={'N': 20, 'mutation_rate_paired': 0.3},
                               shs_seed=11, af3_seed=2)
    tasks = json.loads((directory / 'af3_tasks.json').read_text())
    assert len(tasks) == 2
    for index, task in enumerate(tasks):
        assert task == dict(input=str(directory / f'rna_{index:05d}' / 'af3_input.json'),
                            output=str(directory / f'rna_{index:05d}' / 'af3'),
                            model_dir=config['model_dir'])
    args = submit.call_args.args[0]
    assert args[0] == 'sbatch'
    for flag in ['--parsable', '--array=0-1%1', '--partition=gpu-single',
                 '--cpus-per-task=8', '--mem=20G', '--gres=gpu:A100:1', '--time=00:15:00']:
        assert flag in args
    assert args[-4:] == [str(trial.HPO_DIR / 'slurm/af3-inference-array.slurm'),
                         str(trial.HPO_DIR / 'af3_task.py'), str(directory / 'af3_tasks.json'), config['module']]
    assert dssr.call_count == 2
    for index, call in enumerate(dssr.call_args_list):
        assert call.args[1:] == (directory / f'rna_{index:05d}' / 'dssr', rows[index]['sequence'], 'test-dssr', 10)
    result = json.loads((directory / 'trial_result.json').read_text())
    assert result['job_id'] == '42'
    assert result['parameters'] == requests[0]['parameters']
    assert [s['id'] for s in result['scores']] == ['paired', 'unpaired']
    assert [s['loss'] for s in result['scores']] == [0, 1]
    assert result['loss'] == 0.5
    # Cached completion performs no further external work.
    assert run(trial_case) == 0.5
    assert (generation.call_count, submit.call_count, wait.call_count, dssr.call_count) == (2, 1, 1, 2)


@pytest.mark.parametrize('failure', [subprocess.CalledProcessError(1, 'generator'),
                                     subprocess.TimeoutExpired('generator', 30)])
def test_generation_failure_never_submits(trial_case, failure):
    trial_case[4].side_effect = failure
    with pytest.raises(type(failure)):
        run(trial_case)
    trial_case[5].assert_not_called()
    trial_case[6].assert_not_called()
    assert not (trial_case[2] / 'trial_result.json').exists()


def test_timeout_resume_reattaches_without_regeneration(trial_case):
    wait = trial_case[6]
    finish = wait.side_effect
    wait.side_effect = TimeoutError('still pending')
    with pytest.raises(TimeoutError):
        run(trial_case)
    assert not (trial_case[2] / 'trial_result.json').exists()
    wait.side_effect = finish
    assert run(trial_case) == 0.5
    assert trial_case[4].call_count == 2
    trial_case[5].assert_called_once()


@pytest.mark.parametrize('reply', ['not a job id', ''])
def test_uncertain_submission_blocks_automatic_resubmission(trial_case, reply):
    trial_case[5].return_value = reply
    with pytest.raises(RuntimeError):
        run(trial_case)
    assert (trial_case[2] / 'af3_submission_intent.json').exists()
    assert not (trial_case[2] / 'af3_job_id.json').exists()
    with pytest.raises(RuntimeError):
        run(trial_case)
    trial_case[5].assert_called_once()
    assert trial_case[4].call_count == 2


@pytest.mark.parametrize('failure', [RuntimeError('DSSR failed'),
                                     subprocess.TimeoutExpired('dssr', 10)])
def test_second_rna_failure_does_not_publish_partial_loss(trial_case, failure):
    trial_case[7].side_effect = [[(0, 3)], failure]
    with pytest.raises(type(failure)):
        run(trial_case)
    assert (trial_case[2] / 'rna_00000/base_pair_evaluation.json').exists()
    assert not (trial_case[2] / 'trial_result.json').exists()


@pytest.mark.parametrize('broken', ['missing_marker', 'missing_model', 'invalid_pairs'])
def test_invalid_prediction_cannot_be_scored_as_success(trial_case, broken):
    wait = trial_case[6]
    finish = wait.side_effect
    def damage(*args):
        finish(*args)
        marker = trial_case[2] / 'rna_00000/af3/af3_success.json'
        if broken == 'missing_marker':
            marker.unlink()
        elif broken == 'missing_model':
            Path(json.loads(marker.read_text())['model']).unlink()
    wait.side_effect = damage
    if broken == 'invalid_pairs':
        trial_case[7].side_effect = None
        trial_case[7].return_value = [(0, 99)]
    with pytest.raises((FileNotFoundError, ValueError)):
        run(trial_case)
    assert not (trial_case[2] / 'trial_result.json').exists()


def test_changed_identity_rejects_cached_result(trial_case):
    run(trial_case)
    with pytest.raises(ValueError):
        trial.run_pipeline(trial_case[0], trial_case[1], trial_case[2], 'different',
                           {'mutation_rate_paired': 0.3})


def test_submission_command_failure_preserves_intent_for_manual_recovery(trial_case):
    trial_case[5].side_effect = subprocess.TimeoutExpired('sbatch', 60)
    with pytest.raises(subprocess.TimeoutExpired):
        run(trial_case)
    assert (trial_case[2] / 'af3_submission_intent.json').exists()
    assert not (trial_case[2] / 'af3_job_id.json').exists()
    with pytest.raises(RuntimeError):
        run(trial_case)
    trial_case[5].assert_called_once()


class DSSRTests(unittest.TestCase):
	def test_dssr_json_parser(self):
		predicted = parse_pairs(FIXTURE, SEQ)
		self.assertEqual(
			predicted,
			[(0, 14), (1, 13), (2, 12), (3, 11), (4, 10), (5, 9)],
		)

	def test_mapping_rejects_missing_and_mismatch(self):
		with self.assertRaises(ValueError):
			parse_pairs(FIXTURE, 'A' * len(SEQ))
		broken = dict(FIXTURE, nts=FIXTURE['nts'][:-1])
		with self.assertRaises(ValueError):
			parse_pairs(broken, SEQ)

	def test_resume_and_failure_no_partial_loss(self):
		with tempfile.TemporaryDirectory() as tmp:
			directory = Path(tmp)
			config = dict(fixed={}, dssr='fake', dssr_timeout_seconds=10)
			rows = [dict(
				id='one',
				sequence=SEQ,
				pairs=[[0, 14], [1, 13], [7, 12]],
			)]
			trial.write_json(
				directory / 'trial_metadata.json',
				dict(fingerprint='abc', parameters={}),
			)
			trial.write_json(
				directory / 'af3_job_id.json',
				dict(job_id='7'),
			)
			model = directory / 'rna_00000/af3/model.cif'
			model.parent.mkdir(parents=True)
			model.write_text('mock model')
			trial.write_json(
				model.parent / 'af3_success.json',
				dict(model=str(model)),
			)

			with (
				patch.object(trial, 'wait_job') as wait,
				patch.object(
					trial.dssr,
					'run_dssr',
					return_value=parse_pairs(FIXTURE, SEQ),
				),
			):
				loss = trial.run_pipeline(config, rows, directory, 'abc')
				self.assertAlmostEqual(loss, 5 / 9)
				self.assertAlmostEqual(
					trial.run_pipeline(config, rows, directory, 'abc'),
					loss,
				)
				self.assertEqual(wait.call_count, 1)

			with self.assertRaises(ValueError):
				trial.run_pipeline(config, rows, directory, 'abc', {'N': 3})

			(directory / 'trial_result.json').unlink()
			with patch.object(trial, 'wait_job', side_effect=RuntimeError('failed')):
				with self.assertRaises(RuntimeError):
					trial.run_pipeline(config, rows, directory, 'abc')
			self.assertFalse((directory / 'trial_result.json').exists())


if __name__ == '__main__':
	unittest.main()
