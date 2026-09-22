"""Trial orchestration tests: no generator, Slurm, AF3 or DSSR installation needed."""
import json
import subprocess
from pathlib import Path
import tempfile
import unittest
import pytest
from unittest.mock import Mock, patch

import hpo.pipeline.trial as trial
import hpo.pipeline.run as run
from hpo.pipeline.dssr import parse_pairs


FIXTURE = json.loads(
    (Path(__file__).parent / 'fixtures' / 'dssr_example.json').read_text()
)
SEQ = 'AGUAGUAGUAGUAGA'


@pytest.fixture
def trial_case(tmp_path, monkeypatch):
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
        "module": "bio/alphafold/3.0.1",
        "controller": {
            "partition": "gpu-single",
            "cpus": 1,
            "memory": "20G",
            "gres": "gpu:A100:1",
            "time": "24:00:00",
            "count": 1,
        },
    }

    rows = [
        {
            "id": "paired",
            "sequence": "ACGU",
            "pairs": [[0, 3]],
        },
        {
            "id": "unpaired",
            "sequence": "AAAA",
            "pairs": [],
        },
    ]

    directory = tmp_path / "trial with spaces"
    requests = []

    def generate(args, **kwargs):
        assert args[:2] == [
            config["shs_python"],
            str(
                (
                    run.GENERATOR_DIR / "shs_generator.py"
                ).resolve()
            ),
        ]
        assert args[2] == "--request-json"
        assert args[4] == "--output-json"

        assert kwargs["check"] is True
        assert kwargs["timeout"] == 30
        assert kwargs["stderr"] == subprocess.STDOUT
        assert not kwargs["stdout"].closed

        request_path = Path(args[3])
        output_path = Path(args[5])

        request = json.loads(request_path.read_text())
        requests.append(request)

        output_path.write_text(
            json.dumps({"name": request["task_name"]})
        )

    def run_af3(tasks, cfg):
        assert cfg is config
        assert len(tasks) == len(rows)

        # run_pipeline should persist exactly the task collection passed to
        # run_af3_tasks().
        saved_tasks = json.loads(
            (directory / "af3_tasks.json").read_text()
        )
        assert saved_tasks == tasks

        for task in tasks:
            input_path = Path(task["input"])
            output_directory = Path(task["output"])

            assert input_path.is_file()
            assert task["model_dir"] == str(
                Path(config["model_dir"]).resolve()
            )

            model = (
                output_directory
                / "prediction"
                / "rna_model.cif"
            )
            model.parent.mkdir(parents=True, exist_ok=True)
            model.write_text("test model")

            run.write_json(
                output_directory / "af3_success.json",
                {"model": str(model.resolve())},
            )

    generation = Mock(side_effect=generate)
    af3 = Mock(side_effect=run_af3)

    # RNA 1 is perfect; RNA 2 has one spurious pair.
    # The aggregate mean loss must therefore be 0.5.
    dssr = Mock(
        side_effect=lambda model, directory, sequence, exe, timeout: [(0, 3)]
    )

    monkeypatch.setattr(trial.subprocess, "run", generation)
    monkeypatch.setattr(trial, "run_af3_tasks", af3)
    monkeypatch.setattr(trial.dssr, "run_dssr", dssr)

    return (
        config,
        rows,
        directory,
        requests,
        generation,
        af3,
        dssr,
    )


def run_pipeline(case):
    config, rows, directory, *_ = case
    return trial.run_pipeline(config, rows, directory, 'identity', {'mutation_rate_paired': 0.3})


def test_fresh_trial_wires_two_rnas_and_averages_scores(trial_case):
    config, rows, directory, requests, generation, af3, dssr = trial_case
    assert run_pipeline(trial_case) == pytest.approx(0.5)
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
    af3.assert_called_once_with(tasks, config)
    assert [Path(task['input']) for task in af3.call_args.args[0]] == [
        directory / 'rna_00000' / 'af3_input.json',
        directory / 'rna_00001' / 'af3_input.json',
    ]
    assert not (directory / 'af3_job_id.json').exists()
    assert not (directory / 'af3_submission_intent.json').exists()
    assert dssr.call_count == 2
    for index, call in enumerate(dssr.call_args_list):
        assert call.args[1:] == (directory / f'rna_{index:05d}' / 'dssr', rows[index]['sequence'], 'test-dssr', 10)
    result = json.loads((directory / 'trial_result.json').read_text())
    assert result['parameters'] == requests[0]['parameters']
    assert [s['id'] for s in result['scores']] == ['paired', 'unpaired']
    assert [s['loss'] for s in result['scores']] == [0, 1]
    assert result['loss'] == 0.5
    # Cached completion performs no further external work.
    assert run_pipeline(trial_case) == 0.5
    assert (generation.call_count, af3.call_count, dssr.call_count) == (2, 1, 2)


@pytest.mark.parametrize('failure', [subprocess.CalledProcessError(1, 'generator'),
                                     subprocess.TimeoutExpired('generator', 30)])
def test_generation_failure_never_submits(trial_case, failure):
    trial_case[4].side_effect = failure
    with pytest.raises(type(failure)):
        run_pipeline(trial_case)
    trial_case[5].assert_not_called()
    assert not (trial_case[2] / 'trial_result.json').exists()


@pytest.mark.parametrize('failure', [RuntimeError('DSSR failed'),
                                     subprocess.TimeoutExpired('dssr', 10)])
def test_second_rna_failure_does_not_publish_partial_loss(trial_case, failure):
    (config, rows, directory, requests, generation, af3, dssr) = trial_case
    dssr.side_effect = [[(0, 3)], failure]
    with pytest.raises(type(failure)):
        run_pipeline(trial_case)
    assert (trial_case[2] / 'rna_00000/base_pair_evaluation.json').exists()
    assert not (trial_case[2] / 'trial_result.json').exists()


@pytest.mark.parametrize(
    "broken, expected_exception",
    [
        ("missing_marker", RuntimeError),
        ("missing_model", RuntimeError),
        ("invalid_pairs", ValueError),
    ],
)
def test_invalid_prediction_cannot_be_scored_as_success(trial_case, broken, expected_exception):
    (config, rows, directory, requests, generation, af3, dssr) = trial_case
    finish_af3 = af3.side_effect
    def damage(*args, **kwargs):
        finish_af3(*args, **kwargs)
        marker = directory / 'rna_00000/af3/af3_success.json'
        if broken == 'missing_marker':
            marker.unlink()
        elif broken == 'missing_model':
            marker_data = json.loads(marker.read_text())
            Path(marker_data["model"]).unlink()
    af3.side_effect = damage
    if broken == 'invalid_pairs':
        dssr.side_effect = None
        dssr.return_value = [(0, 99)]
    with pytest.raises(expected_exception):
        run_pipeline(trial_case)
    assert not (directory / 'trial_result.json').exists()


def test_changed_identity_rejects_cached_result(trial_case):
    run_pipeline(trial_case)
    with pytest.raises(ValueError):
        trial.run_pipeline(trial_case[0], trial_case[1], trial_case[2], 'different',
                           {'mutation_rate_paired': 0.3})


@pytest.mark.parametrize(
    'parameters',
    [{}, {'mutation_rate_paired': 0.3}],
    ids=['generator-defaults', 'sampled-parameters'],
)
def test_pipeline_works_without_fixed(trial_case, parameters):
    config, rows, directory, requests, generation, af3, dssr = trial_case
    del config['fixed']

    if parameters:
        loss = trial.run_pipeline(
            config, rows, directory, 'identity', parameters
        )
    else:
        loss = trial.run_pipeline(
            config, rows, directory, 'identity'
        )

    assert loss == pytest.approx(0.5)
    assert generation.call_count == len(rows)
    assert len(requests) == len(rows)
    assert all(
        request['parameters'] == parameters
        for request in requests
    )
    af3.assert_called_once()
    assert dssr.call_count == len(rows)

    result = json.loads(
        (directory / 'trial_result.json').read_text()
    )
    assert result['parameters'] == parameters
    assert result['loss'] == pytest.approx(0.5)


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

            config = {
                "fixed": {},
                "shs_python": "/test/shs/bin/python",
                "shs_seed": 11,
                "af3_seed": 2,
                "model_dir": "/test/models",
                "generator_timeout_seconds": 30,
                "af3_timeout_seconds": 60,
                "dssr": "fake",
                "dssr_timeout_seconds": 10,
            }

            rows = [
                {
                    "id": "one",
                    "sequence": SEQ,
                    "pairs": [
                        [0, 14],
                        [1, 13],
                        [7, 12],
                    ],
                }
            ]

            run.write_json(
                directory / "trial_metadata.json",
                {
                    "fingerprint": "abc",
                    "parameters": {},
                },
            )

            task_directory = directory / "rna_00000"
            task_directory.mkdir(parents=True)

            # The SHS stage has already completed, so resumption must not try
            # to invoke the generator.
            af3_input = task_directory / "af3_input.json"
            af3_input.write_text('{"name": "rna_00000"}')

            # The AF3 stage has already completed and published a valid marker.
            af3_output = task_directory / "af3"
            model = af3_output / "model.cif"
            model.parent.mkdir(parents=True)
            model.write_text("mock model")

            run.write_json(
                af3_output / "af3_success.json",
                {"model": str(model.resolve())},
            )

            with (
                patch.object(
                    trial,
                    "run_af3_tasks",
                ) as af3,
                patch.object(
                    trial.dssr,
                    "run_dssr",
                    return_value=parse_pairs(FIXTURE, SEQ),
                ) as dssr,
            ):
                loss = trial.run_pipeline(
                    config,
                    rows,
                    directory,
                    "abc",
                )
                self.assertAlmostEqual(loss, 5 / 9)

                # The second call reuses the completed trial result.
                self.assertAlmostEqual(
                    trial.run_pipeline(
                        config,
                        rows,
                        directory,
                        "abc",
                    ),
                    loss,
                )

                # AF3 orchestration and DSSR run only on the first call.
                self.assertEqual(af3.call_count, 1)
                self.assertEqual(dssr.call_count, 1)

            # run_pipeline should recreate the task manifest.
            tasks = json.loads(
                (directory / "af3_tasks.json").read_text()
            )
            self.assertEqual(
                tasks,
                [
                    {
                        "input": str(af3_input.resolve()),
                        "output": str(af3_output.resolve()),
                        "model_dir": str(
                            Path(config["model_dir"])
                            .expanduser()
                            .resolve()
                        ),
                    }
                ],
            )

            # The existing trial directory cannot be reused with different
            # parameters.
            with self.assertRaises(ValueError):
                trial.run_pipeline(
                    config,
                    rows,
                    directory,
                    "abc",
                    {"N": 3},
                )

            # Remove only the completed trial result. The next invocation must
            # enter the pipeline again.
            (directory / "trial_result.json").unlink()

            with patch.object(
                trial,
                "run_af3_tasks",
                side_effect=RuntimeError("failed"),
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "failed",
                ):
                    trial.run_pipeline(
                        config,
                        rows,
                        directory,
                        "abc",
                    )

            # A failed pipeline invocation must not publish a partial result.
            self.assertFalse(
                (directory / "trial_result.json").exists()
            )


if __name__ == '__main__':
    unittest.main()
