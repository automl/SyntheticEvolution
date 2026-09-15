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
    directory = tmp_path / 'dssr results'
    document = {'test_document': True}
    parse = Mock(return_value=[(0, 3)])
    monkeypatch.setattr(dssr, 'parse_pairs', parse)
    def execute(args, **kwargs):
        assert args[0] == '/test/x3dna-dssr'
        assert '-i=' + str(model.resolve()) in args
        assert '--json' in args
        output = Path(next(arg[3:] for arg in args if arg.startswith('-o=')))
        assert output.parent == directory.resolve()
        assert Path(kwargs['cwd']).resolve() == directory.resolve()
        assert kwargs['check'] is True
        assert kwargs['timeout'] == 7
        assert kwargs['stderr'] == subprocess.STDOUT
        assert not kwargs['stdout'].closed
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
        assert dssr.run_dssr(model, directory, 'ACGU', '/test/x3dna-dssr', 7) == [(0, 3)]
        parse.assert_called_once_with(document, 'ACGU')
    else:
        expected = {'exit_failure': subprocess.CalledProcessError,
                    'timeout': subprocess.TimeoutExpired, 'invalid_json': ValueError,
                    'missing_output': FileNotFoundError}[outcome]
        with pytest.raises(expected):
            dssr.run_dssr(model, directory, 'ACGU', '/test/x3dna-dssr', 7)
        parse.assert_not_called()
    command.assert_called_once()
