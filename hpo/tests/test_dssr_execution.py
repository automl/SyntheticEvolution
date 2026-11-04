"""DSSR subprocess contract and error propagation; the binary is replaced."""
import json
import subprocess
from pathlib import Path
from unittest.mock import Mock

import pytest
from hpo.pipeline import dssr


@pytest.mark.parametrize('outcome', ['success', 'exit_failure', 'timeout', 'invalid_json', 'missing_output'])
def test_dssr_execution_and_output_errors(tmp_path, monkeypatch, outcome):
    model = tmp_path / 'input model.cif'
    model.write_text('model')
    directory = tmp_path / 'rna dir'
    scratch = tmp_path / 'scratch' / 'rna_00000'
    log_path = scratch / 'dssr.log'
    document = {'test_document': True}
    parse = Mock(return_value=[(0, 3)])
    monkeypatch.setattr(dssr, 'parse_pairs', parse)
    work_directories = []
    def execute(args, **kwargs):
        assert args[0] == '/test/x3dna-dssr'
        assert '-i=' + str(model.resolve()) in args
        assert '--json' in args
        output = Path(next(arg[3:] for arg in args if arg.startswith('-o=')))
        # DSSR runs and writes its JSON in a private directory inside scratch.
        work_directory = Path(kwargs['cwd'])
        work_directories.append(work_directory)
        assert output.parent == work_directory
        assert work_directory.parent == scratch.resolve()
        assert kwargs['check'] is True
        assert kwargs['timeout'] == 7
        assert kwargs['stderr'] == subprocess.STDOUT
        assert kwargs['stdout'].name == str(log_path.resolve())
        assert not kwargs['stdout'].closed
        (work_directory / 'dssr-pairs.pdb').write_text('side file')
        if outcome == 'exit_failure':
            raise subprocess.CalledProcessError(1, args)
        if outcome == 'timeout':
            raise subprocess.TimeoutExpired(args, 7)
        if outcome == 'invalid_json':
            output.write_text('{broken json')
        elif outcome == 'success':
            output.write_text(json.dumps(document))
    command = Mock(side_effect=execute)
    monkeypatch.setattr(dssr.subprocess, 'run', command)
    if outcome == 'success':
        assert dssr.run_dssr(model, directory, 'ACGU', '/test/x3dna-dssr', 7, scratch=scratch) == [(0, 3)]
        parse.assert_called_once_with(document, 'ACGU')
        assert json.loads((directory / 'dssr.json').read_text()) == document
    else:
        expected = {'exit_failure': subprocess.CalledProcessError,
                    'timeout': subprocess.TimeoutExpired, 'invalid_json': ValueError,
                    'missing_output': FileNotFoundError}[outcome]
        with pytest.raises(expected):
            dssr.run_dssr(model, directory, 'ACGU', '/test/x3dna-dssr', 7, scratch=scratch)
        parse.assert_not_called()
    command.assert_called_once()
    # The workspace gets only the JSON, and only when DSSR produced valid JSON.
    # The log stays in scratch; side files are removed in every case.
    kept = ['dssr.json'] if outcome == 'success' else []
    assert sorted(path.name for path in directory.iterdir()) == kept
    assert sorted(path.name for path in scratch.iterdir()) == ['dssr.log']
    assert not work_directories[0].exists()


def test_dssr_json_is_kept_when_mapping_is_rejected(tmp_path, monkeypatch):
    model = tmp_path / 'model.cif'
    model.write_text('model')
    def execute(args, **kwargs):
        output = Path(next(arg[3:] for arg in args if arg.startswith('-o=')))
        output.write_text(json.dumps({'nts': []}))
    monkeypatch.setattr(dssr.subprocess, 'run', execute)
    with pytest.raises(ValueError):
        dssr.run_dssr(model, tmp_path / 'rna', 'ACGU', '/test/x3dna-dssr', 7)
    assert json.loads((tmp_path / 'rna' / 'dssr.json').read_text()) == {'nts': []}
