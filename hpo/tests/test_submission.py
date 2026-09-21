"""Inspect controller sbatch commands without invoking Slurm."""

import importlib
import sys
from unittest.mock import Mock

import pytest


@pytest.mark.parametrize(
    "mode,python_key,job_name",
    [
        ("neps", "neps_python", "shs-hpo"),
        (
            "standalone",
            "controller_python",
            "shs-standalone",
        ),
    ],
)
def test_controller_submission_resources_and_command(
    tmp_path,
    monkeypatch,
    capsys,
    mode,
    python_key,
    job_name,
):
    suffix = (
        "neps_hpo"
        if mode == "neps"
        else "standalone_trial"
    )
    submission_module = importlib.import_module(
        f"hpo.submit_{suffix}"
    )

    config_path = (
        tmp_path / "config space 'quote' $literal.yaml"
    )
    root = tmp_path / "repo with spaces"
    workspace = tmp_path / "workspace"

    config = {
        "module": "bio/alphafold/3.0.1",
        "neps_python": "/env/neps space/bin/python",
        "shs_python": "/env/shs space/bin/python",
        "controller_python": (
            "/env/controller space/bin/python"
        ),
        "controller": {
            "partition": "gpu-single",
            "cpus": 2,
            "memory": "20G",
            "gres": "gpu:A100:1",
            "time": "03:00:00",
        },
    }

    load_config = Mock(return_value=config)
    command = Mock(return_value="123")

    monkeypatch.setattr(
        submission_module,
        "load_config",
        load_config,
    )
    monkeypatch.setattr(
        submission_module,
        "get_workspace_path",
        lambda cfg: workspace,
    )
    monkeypatch.setattr(
        submission_module,
        "repo_path",
        lambda path: config_path,
    )
    monkeypatch.setattr(submission_module, "ROOT", root)
    monkeypatch.setattr(
        submission_module,
        "command",
        command,
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["submit", "--config", str(config_path)],
    )

    submission_module.main()

    load_config.assert_called_once_with(
        str(config_path),
        mode=mode,
        submitting_controller=True,
    )
    command.assert_called_once()

    args = command.call_args.args[0]

    assert args[:2] == ["sbatch", "--parsable"]
    assert args.count("sbatch") == 1

    # The controller is the only submitted job.
    assert not any(
        "af3-inference-array.slurm" in argument
        for argument in args
    )
    assert not any(
        argument.startswith("--array=")
        for argument in args
    )

    expected_flags = [
        f"--job-name={job_name}",
        f"--chdir={root}",
        "--nodes=1",
        "--ntasks=1",
        "--partition=gpu-single",
        "--cpus-per-task=2",
        "--mem=20G",
        "--gres=gpu:A100:1",
        "--time=03:00:00",
        (
            "--output="
            + str(workspace / "controller_logs/%j.out")
        ),
        (
            "--error="
            + str(workspace / "controller_logs/%j.err")
        ),
    ]

    for flag in expected_flags:
        assert flag in args

    # Shell-sensitive paths must remain separate list elements.
    # No shell command should be constructed with --wrap.
    assert not any(
        argument.startswith("--wrap=")
        for argument in args
    )

    controller_script = (
        root / "hpo" / "slurm" / "controller.slurm"
    )
    controller_module = f"hpo.run_{suffix}"

    assert args[-5:] == [
        str(controller_script),
        config["module"],
        config[python_key],
        controller_module,
        str(config_path),
    ]

    assert (
        workspace / "controller_logs"
    ).is_dir()

    assert capsys.readouterr().out.strip() == "123"