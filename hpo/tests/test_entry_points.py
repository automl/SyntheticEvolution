"""Every entry point starts and shows its help without touching the cluster."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('module', [
    'hpo.run_standalone_trial',
    'hpo.run_neps_hpo',
    'hpo.submit_standalone_trial',
    'hpo.submit_neps_hpo',
    'hpo.export_results',
    'hpo.make_checks',
])
def test_entry_point_shows_help(module):
    result = subprocess.run([sys.executable, '-m', module, '--help'],
                            cwd=ROOT, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert 'usage:' in result.stdout
