"""Run-level helpers: fingerprint, existing-run check and atomic JSON writes."""
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier

import pytest

from hpo.pipeline.run import check_existing_run, fingerprint, write_json

ROWS = [dict(id='rna', sequence='ACGU', pairs=[(0, 3)])]


# ---- Fingerprint ------------------------------------------------------------

def test_fingerprint_covers_dataset_inputs():
    variants = [
        [dict(ROWS[0], interactions=[[0, 3, .2]], mutation_rates=[.1] * 4)],
        [dict(ROWS[0], interactions=[[0, 3, .3]], mutation_rates=[.1] * 4)],
        [dict(ROWS[0], interactions=[[0, 3, .2]], mutation_rates=[.2] * 4)],
    ]
    assert len({fingerprint({}, rows) for rows in variants}) == 3


def test_fingerprint_covers_scientific_settings_only():
    base = dict(shs_seed=1, fixed={'N': 8})
    identity = fingerprint(base, ROWS)
    assert fingerprint(dict(base, shs_seed=2), ROWS) != identity
    # Resources, timeouts, retention and failure handling may change within a run.
    for key, value in [('controller', {'count': 4}), ('evaluations', 99), ('ignore_errors', True),
                       ('af3_timeout_seconds', 1), ('storage', {'scratch': '/x'}),
                       ('failures', {'max_rate': 1.0})]:
        assert fingerprint(dict(base, **{key: value}), ROWS) == identity, key


# ---- Existing runs ----------------------------------------------------------

def test_new_run_has_no_metadata(tmp_path):
    assert check_existing_run(tmp_path, 'expected') is False
    assert list(tmp_path.iterdir()) == []


def test_matching_run_is_verified_without_changes(tmp_path):
    metadata = tmp_path / 'run_metadata.json'
    metadata.write_text(json.dumps({'fingerprint': 'expected'}))
    before = metadata.read_bytes()
    assert check_existing_run(tmp_path, 'expected') is True
    assert metadata.read_bytes() == before


@pytest.mark.parametrize('contents', ['{"fingerprint": "different"}', '{}', '[]', 'invalid json'])
def test_changed_or_invalid_run_is_rejected_without_changes(tmp_path, contents):
    metadata = tmp_path / 'run_metadata.json'
    metadata.write_text(contents)
    with pytest.raises(ValueError):
        check_existing_run(tmp_path, 'expected')
    assert metadata.read_text() == contents


# ---- Atomic JSON ------------------------------------------------------------

def test_write_json_handles_simultaneous_writers(tmp_path, monkeypatch):
    target = tmp_path / 'run_metadata.json'
    metadata = {'fingerprint': 'same-run'}
    ready_to_replace = Barrier(2)
    original_replace = Path.replace

    def synchronized_replace(source, destination):
        if Path(destination) == target:
            # Both temporary files are written before either rename.
            ready_to_replace.wait(timeout=5)
        return original_replace(source, destination)

    monkeypatch.setattr(Path, 'replace', synchronized_replace)
    with ThreadPoolExecutor(max_workers=2) as workers:
        for future in [workers.submit(write_json, target, metadata) for _ in range(2)]:
            future.result(timeout=10)

    assert json.loads(target.read_text()) == metadata
    assert list(tmp_path.iterdir()) == [target]  # no temporary files left
