"""Dataset input: CSV rows, pairs, interactions and per-position mutation rates."""
import csv
import json

import pytest

from hpo.pipeline.run import create_dataset, normalize_pairs


def write_csv(path, rows, columns=('id', 'sequence', 'pairs')):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)
    return path


def dataset(tmp_path, **values):
    """One RNA ``example``/``ACGU`` with the given optional columns (JSON-encoded)."""
    row = dict(id='example', sequence='ACGU', **values)
    encoded = [json.dumps(v) if isinstance(v, list) else v for v in row.values()]
    return create_dataset(write_csv(tmp_path / 'input.csv', [encoded], columns=list(row)))


# ---- Pairs ------------------------------------------------------------------

def test_pairs_are_normalized_and_deduplicated():
    assert normalize_pairs([[0, 1], [1, 0]], 2) == {(0, 1)}


@pytest.mark.parametrize('pairs', [[[0, 0]], [[-1, 2]], [[0, 4]], [[0.0, 1]], [[True, 1]], [[0, 1, 2]]])
def test_invalid_pairs_are_rejected(pairs):
    with pytest.raises(ValueError):
        normalize_pairs(pairs, 3)


# ---- CSV files --------------------------------------------------------------

def test_csv_normalizes_case_whitespace_and_reversed_pairs(tmp_path):
    path = write_csv(tmp_path / 'input.csv', [(' one ', ' acgu ', '[[3,0],[0,3]]')])
    assert create_dataset(path) == [dict(id='one', sequence='ACGU', pairs=[(0, 3)])]


def test_folder_input_reads_every_csv(tmp_path):
    for index in range(2):
        write_csv(tmp_path / f'{index}.csv', [(f'rna{index}', 'ACGU', '[[0,3]]')])
    assert [row['id'] for row in create_dataset(tmp_path)] == ['rna0', 'rna1']


def test_duplicate_id_across_files_is_rejected(tmp_path):
    for name in ('a.csv', 'b.csv'):
        write_csv(tmp_path / name, [('same', 'ACGU', '[]')])
    with pytest.raises(ValueError):
        create_dataset(tmp_path)


@pytest.mark.parametrize('row', [
    ('', 'ACGU', '[]'), ('one', '', '[]'), ('one', 'ACNT', '[]'),
    ('one', 'ACGU', 'not-json'), ('one', 'ACGU', '[[0,4]]'), ('one', 'ACGU', '[[true,3]]'),
])
def test_invalid_row_is_rejected(tmp_path, row):
    with pytest.raises(ValueError):
        create_dataset(write_csv(tmp_path / 'input.csv', [row]))


def test_empty_csv_is_rejected(tmp_path):
    with pytest.raises(ValueError):
        create_dataset(write_csv(tmp_path / 'empty.csv', []))


# ---- Optional columns: interactions and mutation rates ----------------------

def test_interactions_define_pairs_and_keep_zero_strengths(tmp_path):
    [row] = dataset(tmp_path, interactions=[[3, 0, .2], [1, 2, 0]], mutation_rates=[0, .3, .5, 1])
    assert row['pairs'] == [(0, 3)]  # only positive strengths are scoring targets
    assert row['interactions'] == [[0, 3, .2], [1, 2, 0]]
    assert row['mutation_rates'] == [0, .3, .5, 1]


def test_explicit_pairs_take_precedence_over_interactions(tmp_path):
    [row] = dataset(tmp_path, pairs=[[0, 3]], interactions=[])
    assert (row['pairs'], row['interactions']) == ([(0, 3)], [])


def test_blank_optional_columns_are_ignored(tmp_path):
    [row] = dataset(tmp_path, pairs=[[0, 3]], interactions='', mutation_rates='')
    assert set(row) == {'id', 'sequence', 'pairs'}


@pytest.mark.parametrize('values', [
    dict(interactions=[[0, 4, 1]]), dict(interactions=[[0, 0, 1]]),
    dict(interactions=[[0, 3, -1]]), dict(interactions=[[0, 3, float('nan')]]),
    dict(interactions=[[0, 3, True]]), dict(interactions=[[False, 3, 1]]),
    dict(interactions=[[0, 3, 1.1]]), dict(interactions=[[0, 3, float('inf')]]),
    dict(interactions=[[0, 3, .1], [3, 0, .2]]),
    dict(pairs=[], mutation_rates=[]), dict(pairs=[], mutation_rates=[0, 0, 0, 2]),
    dict(pairs=[], mutation_rates=[False, 0, 0, 0]), dict(pairs='null'),
])
def test_invalid_optional_columns_are_rejected(tmp_path, values):
    with pytest.raises(ValueError):
        dataset(tmp_path, **values)
