"""CSV boundary cases beyond the existing pair-normalization tests."""
import csv
import pytest
from hpo.pipeline.trial import create_dataset


def write_csv(path, rows, columns=('id', 'sequence', 'pairs')):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


@pytest.mark.parametrize('row', [('', 'ACGU', '[]'), ('one', '', '[]'),
    ('one', 'ACNT', '[]'), ('one', 'ACGU', 'not-json'), ('one', 'ACGU', '[[0,4]]'),
    ('one', 'ACGU', '[[true,3]]')])
def test_bad_csv_row_rejected(tmp_path, row):
    path = tmp_path / 'input.csv'
    write_csv(path, [row])
    with pytest.raises(ValueError):
        create_dataset(path)


def test_empty_csv_rejected(tmp_path):
    path = tmp_path / 'empty.csv'
    write_csv(path, [])
    with pytest.raises(ValueError):
        create_dataset(path)


def test_duplicate_id_across_files_rejected(tmp_path):
    for name in ['a.csv', 'b.csv']:
        write_csv(tmp_path / name, [('same', 'ACGU', '[]')])
    with pytest.raises(ValueError):
        create_dataset(tmp_path)


def test_csv_normalizes_case_whitespace_and_reversed_pairs(tmp_path):
    path = tmp_path / 'input.csv'
    write_csv(path, [(' one ', ' acgu ', '[[3,0],[0,3]]')])
    assert create_dataset(path) == [dict(id='one', sequence='ACGU', pairs=[(0, 3)])]
