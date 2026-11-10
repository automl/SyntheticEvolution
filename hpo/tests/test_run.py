"""Run-level helpers: fingerprint, existing-run check and atomic JSON writes."""
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier

import pytest

from hpo.pipeline.run import (
    check_existing_run, fingerprint, shs_generation_request, typed_parameters, validate_config_for_dataset,
    write_json,
)

ROWS = [dict(id='rna', chains=[dict(id='A', type='rna', sequence='ACGU')],
             type='rna', sequence='ACGU', pairs=[(0, 3)])]
CONFIG = dict(
    shs_seed=1, optimize='rna',
    generators=dict(rna=dict(python='/env/rna/python', fixed={'N': 8})),
    extraction=dict(dssr='/env/dssr', dssr_timeout_seconds=300),
)


# ---- Fingerprint ------------------------------------------------------------

def test_fingerprint_covers_dataset_inputs():
    variants = [
        [dict(ROWS[0], interactions=[[0, 3, .2]], mutation_rates=[.1] * 4)],
        [dict(ROWS[0], interactions=[[0, 3, .3]], mutation_rates=[.1] * 4)],
        [dict(ROWS[0], interactions=[[0, 3, .2]], mutation_rates=[.2] * 4)],
    ]
    assert len({fingerprint(CONFIG, rows) for rows in variants}) == 3


def test_fingerprint_covers_scientific_settings_only():
    base = CONFIG
    identity = fingerprint(base, ROWS)
    assert fingerprint(dict(base, shs_seed=2), ROWS) != identity
    # Generator parameters and the protein contact definition change the results ...
    fixed = dict(base, generators=dict(rna=dict(python='/env/rna/python', fixed={'N': 9})))
    assert fingerprint(fixed, ROWS) != identity
    contacts = dict(base, extraction=dict(base['extraction'], protein={'cutoff': 9.0}))
    assert fingerprint(contacts, ROWS) != identity
    # ... executable paths and DSSR's timeout do not.
    moved = dict(
        base, generators=dict(rna=dict(python='/other/python', fixed={'N': 8})),
        extraction=dict(dssr='/other/dssr', dssr_timeout_seconds=1),
    )
    assert fingerprint(moved, ROWS) == identity
    # Resources, timeouts, retention and failure handling may change within a run.
    for key, value in [('controller', {'count': 4}), ('evaluations', 99), ('ignore_errors', True),
                       ('af3_timeout_seconds', 1), ('storage', {'scratch': '/x'}),
                       ('failures', {'max_rate': 1.0}), ('trials', [{'name': 'a'}])]:
        assert fingerprint(dict(base, **{key: value}), ROWS) == identity, key


# ---- Generation request and per-type settings -------------------------------

def test_generation_request_keeps_only_generator_fields_and_the_rows_type_parameters():
    row = dict(ROWS[0], interactions=[[0, 3, .2]], mutation_rates=[.1] * 4)
    request = shs_generation_request(row, 3, {'rna': {'N': 8}, 'protein': {'N': 2}}, 11, 22)
    assert request == dict(
        id='rna', sequence='ACGU', pairs=[(0, 3)], interactions=[[0, 3, .2]], mutation_rates=[.1] * 4,
        task_name='rna_00003', parameters={'N': 8}, shs_seed=11, af3_seed=22,
    )
    protein = dict(row, type='protein')
    assert shs_generation_request(protein, 0, {'rna': {'N': 8}, 'protein': {'N': 2}}, 1, 2)['parameters'] == {'N': 2}


def test_typed_parameters_override_only_the_optimized_generators_fixed_values():
    config = dict(CONFIG, generators=dict(
        rna=dict(python='/env/rna/python', fixed={'N': 8, 'a': 1}),
        protein=dict(python='/env/protein/python', fixed={'N': 2}),
        ))
    protein = dict(ROWS[0], type='protein')
    assert typed_parameters(config, ROWS, {'N': 64}) == {'rna': {'N': 64, 'a': 1}}
    assert typed_parameters(config, ROWS + [protein], {'N': 64}) == {'rna': {'N': 64, 'a': 1}, 'protein': {'N': 2}}
    assert typed_parameters(config, [protein]) == {'protein': {'N': 2}}   # nothing sampled for the other type


def test_dataset_chain_types_need_their_generator_and_extractor():
    protein = dict(ROWS[0], type='protein')
    config = dict(CONFIG, extraction={})
    with pytest.raises(ValueError, match='generators.protein'):
        validate_config_for_dataset(CONFIG, [protein])
    with pytest.raises(ValueError, match='extraction.dssr'):
        validate_config_for_dataset(config, ROWS)
    # Protein chains need no DSSR.
    protein_config = dict(CONFIG, generators=dict(protein=dict(python='/env/protein/python')), extraction={})
    validate_config_for_dataset(protein_config, [protein])


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
