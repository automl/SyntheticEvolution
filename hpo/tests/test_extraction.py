"""Extraction: DSSR (RNA) parsing and execution, protein contacts and dispatch."""
import sys
import json
import subprocess
from pathlib import Path
from unittest.mock import Mock

import pytest

from hpo.pipeline import extraction
from hpo.pipeline.extraction import parse_pairs

FIXTURE = json.loads((Path(__file__).parent / 'fixtures' / 'dssr_example.json').read_text())
SEQUENCE = 'AGUAGUAGUAGUAGA'


# ---- Parsing ----------------------------------------------------------------

def test_pairs_are_mapped_through_validated_nucleotide_ids():
    assert parse_pairs(FIXTURE, SEQUENCE) == [(0, 14), (1, 13), (2, 12), (3, 11), (4, 10), (5, 9)]


@pytest.mark.parametrize('data,sequence', [
    (FIXTURE, 'A' * len(SEQUENCE)),                       # base mismatch
    (dict(FIXTURE, nts=FIXTURE['nts'][:-1]), SEQUENCE),   # missing residue
], ids=['mismatched-sequence', 'missing-residue'])
def test_unexpected_residues_are_rejected(data, sequence):
    with pytest.raises(ValueError):
        parse_pairs(data, sequence)


# ---- Running DSSR -----------------------------------------------------------


@pytest.mark.parametrize('outcome', ['success', 'exit_failure', 'timeout', 'invalid_json', 'missing_output'])
def test_dssr_runs_in_scratch_and_propagates_errors(tmp_path, monkeypatch, outcome):
    model = tmp_path / 'input model.cif'
    model.write_text('model')
    directory = tmp_path / 'rna dir'
    scratch = tmp_path / 'scratch' / 'rna_00000'
    log_path = scratch / 'dssr.log'
    document = {'test_document': True}
    parse = Mock(return_value=[(0, 3)])
    monkeypatch.setattr(extraction, 'parse_pairs', parse)
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
    monkeypatch.setattr(extraction.subprocess, 'run', command)
    if outcome == 'success':
        assert extraction.run_dssr(model, directory, 'ACGU', '/test/x3dna-dssr', 7, scratch=scratch) == [(0, 3)]
        parse.assert_called_once_with(document, 'ACGU')
        assert json.loads((directory / 'dssr.json').read_text()) == document
    else:
        expected = {'exit_failure': subprocess.CalledProcessError,
                    'timeout': subprocess.TimeoutExpired, 'invalid_json': ValueError,
                    'missing_output': FileNotFoundError}[outcome]
        with pytest.raises(expected):
            extraction.run_dssr(model, directory, 'ACGU', '/test/x3dna-dssr', 7, scratch=scratch)
        parse.assert_not_called()
    command.assert_called_once()
    # The workspace gets only the JSON, and only when DSSR produced valid JSON.
    # The log stays in scratch; side files are removed in every case.
    kept = ['dssr.json'] if outcome == 'success' else []
    assert sorted(path.name for path in directory.iterdir()) == kept
    assert sorted(path.name for path in scratch.iterdir()) == ['dssr.log']
    assert not work_directories[0].exists()


def test_dssr_json_is_kept_when_its_pairs_are_rejected(tmp_path, monkeypatch):
    model = tmp_path / 'model.cif'
    model.write_text('model')
    def execute(args, **kwargs):
        output = Path(next(arg[3:] for arg in args if arg.startswith('-o=')))
        output.write_text(json.dumps({'nts': []}))
    monkeypatch.setattr(extraction.subprocess, 'run', execute)
    with pytest.raises(ValueError):
        extraction.run_dssr(model, tmp_path / 'rna', 'ACGU', '/test/x3dna-dssr', 7)
    assert json.loads((tmp_path / 'rna' / 'dssr.json').read_text()) == {'nts': []}


# ---- extraction dispatch and protein contacts -------------------------------

def write_cif(path, residues):
    """Write a minimal mmCIF with one CA per (resname, (x, y, z)) in chain A."""
    columns = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id',
               'label_comp_id', 'label_asym_id', 'label_seq_id', 'pdbx_PDB_ins_code',
               'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv',
               'auth_seq_id', 'auth_asym_id', 'pdbx_PDB_model_num']
    lines = ['data_test', 'loop_'] + [f'_atom_site.{c}' for c in columns]
    for i, (name, (x, y, z)) in enumerate(residues, start=1):
        lines.append(f'ATOM {i} C CA . {name} A {i} ? {x:.3f} {y:.3f} {z:.3f} 1.0 0.0 {i} A 1')
    path.write_text('\n'.join(lines) + '\n')
    return path


@pytest.fixture
def helix_like(tmp_path):
    # Residues on a line, 4 A apart: |i-j| = 1 -> 4 A, 2 -> 8 A (boundary), 3 -> 12 A.
    names = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY']
    return write_cif(tmp_path / 'model.cif',
                     [(n, (4.0 * i, 0.0, 0.0)) for i, n in enumerate(names)]), 'ACDEFG'


def test_contact_pairs_cutoff_is_inclusive_and_respects_min_separation(helix_like):
    pytest.importorskip('Bio')
    model, sequence = helix_like
    assert extraction.contact_pairs(model, sequence, cutoff=8.0, min_separation=2) == [
        (0, 2), (1, 3), (2, 4), (3, 5)]
    assert extraction.contact_pairs(model, sequence, cutoff=11.9, min_separation=2) == [
        (0, 2), (1, 3), (2, 4), (3, 5)]
    assert extraction.contact_pairs(model, sequence, cutoff=12.0, min_separation=3) == [
        (0, 3), (1, 4), (2, 5)]
    assert extraction.contact_pairs(model, sequence, cutoff=7.9, min_separation=1) == [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]


def test_contact_pairs_rejects_sequence_mismatch_and_wrong_length(helix_like):
    pytest.importorskip('Bio')
    model, sequence = helix_like
    with pytest.raises(ValueError):
        extraction.contact_pairs(model, 'AAAAAA')
    with pytest.raises(ValueError):
        extraction.contact_pairs(model, sequence[:-1])
    with pytest.raises(ValueError):
        extraction.contact_pairs(model, sequence, chain='B')
    with pytest.raises(ValueError):
        extraction.contact_pairs(model, sequence, atom='CB')


def test_extract_pairs_dispatches_by_chain_type(helix_like, tmp_path, monkeypatch):
    pytest.importorskip('Bio')
    model, sequence = helix_like
    settings = {'dssr': '/test/dssr', 'dssr_timeout_seconds': 5,
                'protein': {'cutoff': 8.0, 'min_separation': 2}}
    directory = tmp_path / 'extraction'
    assert extraction.extract_pairs('protein', model, directory, sequence, settings) == [
        (0, 2), (1, 3), (2, 4), (3, 5)]
    saved = json.loads((directory / 'contacts.json').read_text())
    assert saved['options'] == {'atom': 'CA', 'cutoff': 8.0, 'min_separation': 2}

    run = Mock(return_value=[(0, 3)])
    monkeypatch.setattr(extraction, 'run_dssr', run)
    assert extraction.extract_pairs('rna', model, directory, 'ACGU', settings) == [(0, 3)]
    run.assert_called_once_with(model, directory, 'ACGU', '/test/dssr', 5, scratch=None)
    run.reset_mock()
    scratch = tmp_path / 'scratch'
    extraction.extract_pairs('rna', model, directory, 'ACGU', settings, scratch=scratch)
    run.assert_called_once_with(model, directory, 'ACGU', '/test/dssr', 5, scratch=scratch)

    with pytest.raises(ValueError):
        extraction.extract_pairs('dna', model, directory, 'ACGT', settings)


def test_protein_extraction_without_biopython_names_the_environment(helix_like, monkeypatch):
    model, sequence = helix_like
    monkeypatch.setitem(sys.modules, 'Bio.PDB.MMCIFParser', None)
    with pytest.raises(ImportError, match='controller environment'):
        extraction.contact_pairs(model, sequence)
