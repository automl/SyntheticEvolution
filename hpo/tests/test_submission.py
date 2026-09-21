"""Inspect controller sbatch commands without invoking Slurm."""

import importlib
import sys
from unittest.mock import Mock

import pytest


import os
from pathlib import Path
import subprocess

from hpo.pipeline.trial import ROOT


def test_controller_script_starts_requested_module(tmp_path):
    script = ROOT / "hpo/slurm/controller.slurm"
    assert script.is_file(), f"Missing controller script: {script}"

    # Bash reads BASH_ENV before executing the script.
    # Provide a fake module command without loading cluster software.
    bash_env = tmp_path / "bash_env.sh"
    bash_env.write_text(
        'module() {\n'
        '    [[ "$#" -eq 2 && "$1" == "load" '
        '&& "$2" == "bio/alphafold/3.0.1" ]] || return 1\n'
        '    export TEST_MODULE_LOADED=1\n'
        '}\n'
    )

    controller_python = tmp_path / "controller python"
    controller_python.write_text(
        "#!/bin/bash\n"
        "set -euo pipefail\n"
        '[[ "${TEST_MODULE_LOADED:-}" == "1" ]]\n'
        '[[ -n "${AF3_PYTHON:-}" ]]\n'
        'printf "%s\\n" "$@"\n'
    )
    controller_python.chmod(0o755)

    config_path = tmp_path / "config space 'quote' $literal.yaml"
    config_path.write_text("{}\n")

    env = os.environ.copy()
    env["BASH_ENV"] = str(bash_env)

    # Do not let an inherited variable mask a missing shell assignment.
    env.pop("controller_module", None)
    env.pop("AF3_PYTHON", None)
    env.pop("TEST_MODULE_LOADED", None)

    result = subprocess.run(
        [
            "bash",
            str(script),
            "bio/alphafold/3.0.1",
            str(controller_python),
            "hpo.run_standalone_trial",
            str(config_path),
        ],
        env=env,
        capture_output=True,
        text=True,
        timeout=10,
    )

    assert result.returncode == 0, (
        f"Controller script failed:\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert result.stdout.splitlines() == [
        "-u",
        "-m",
        "hpo.run_standalone_trial",
        "--config",
        str(config_path),
    ]

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