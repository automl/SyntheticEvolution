"""Scoring objectives (F1, soft F1, Weisfeiler-Lehman, TM-score) and the dataset loss."""
import os

import numpy as np
import pytest

from hpo.pipeline.scoring import USALIGN, aggregate, score_pairs, score_row, score_soft_pairs, score_tm, score_wl


# ---- F1 ---------------------------------------------------------------------

def test_f1_counts_pairs_and_loss_is_one_minus_f1():
    predicted = [(0, 14), (1, 13), (2, 12), (3, 11), (4, 10), (5, 9)]
    score = score_pairs({(0, 14), (1, 13), (7, 12)}, predicted)
    assert (score['tp'], score['fp'], score['fn']) == (2, 4, 1)
    assert score['loss'] == pytest.approx(5 / 9)


@pytest.mark.parametrize('target,predicted,loss', [
    ([], [], 0),          # nothing to find, nothing predicted
    ([(0, 1)], [], 1),    # missed pair
    ([], [(0, 1)], 1),    # spurious pair
])
def test_f1_empty_cases(target, predicted, loss):
    assert score_pairs(target, predicted)['loss'] == loss


# ---- Soft F1 ----------------------------------------------------------------

HELIX = {(0, 20), (1, 19), (2, 18)}


@pytest.mark.parametrize('slip,f1', [(0, 1), (1, .75), (2, .5), (3, .25), (4, 0), (6, 0)])
def test_soft_f1_drops_gradually_with_slippage(slip, f1):
    slipped = {(i, j - slip) for i, j in HELIX}
    assert score_soft_pairs(HELIX, slipped)['f1'] == pytest.approx(f1)


@pytest.mark.parametrize('predicted', [
    set(), {(0, 20)}, {(0, 20), (1, 19), (2, 18), (5, 9)}, {(0, 12), (5, 9)}, {(3, 15)},
])
def test_soft_f1_equals_f1_without_slipped_pairs(predicted):
    assert score_soft_pairs(HELIX, predicted)['f1'] == pytest.approx(score_pairs(HELIX, predicted)['f1'])


def test_soft_f1_matches_each_pair_at_most_once():
    # One predicted pair next to two target pairs earns credit only once.
    assert score_soft_pairs({(0, 20), (1, 20)}, {(0, 19)})['credit'] == .75


def test_soft_f1_objective():
    row = dict(sequence='A' * 21, pairs=[list(p) for p in HELIX])
    assert score_row(row, {(i, j - 1) for i, j in HELIX}, objective='f1_soft')['loss'] == pytest.approx(.25)
    assert score_soft_pairs([], [])['f1'] == 1


def test_scoring_uses_the_explicit_pairs_not_the_generator_interactions():
    row = dict(sequence='ACGU', pairs=[[0, 3], [1, 2]], interactions=[[0, 1, .2]])
    assert score_row(row, {(0, 3), (0, 1)}) == dict(tp=1, fp=1, fn=1, f1=.5, loss=.5)


def test_unknown_objective_is_rejected():
    with pytest.raises(ValueError):
        score_row(dict(sequence='ACGU', pairs=[]), set(), objective='not-an-objective')


# ---- Weisfeiler-Lehman ------------------------------------------------------

@pytest.mark.parametrize('pairs', [{(0, 3)}, set()], ids=['paired', 'unpaired'])
def test_wl_identical_structures_have_zero_loss(pairs):
    metrics = score_row(dict(sequence='ACGU', pairs=list(pairs)), pairs, objective='wl')
    assert (metrics['wl'], metrics['loss']) == (pytest.approx(1.0), pytest.approx(0.0))


def test_wl_is_symmetric():
    forward = score_wl({(0, 3)}, {(0, 2)}, length=4)
    reverse = score_wl({(0, 2)}, {(0, 3)}, length=4)
    assert forward == pytest.approx(reverse)


@pytest.mark.parametrize('n_iter', [1, 5])
def test_wl_matches_the_reference_grakel_implementation(n_iter):
    """Compare against GraKeL directly, including the iteration count."""
    from grakel import Graph
    from grakel.kernels import VertexHistogram, WeisfeilerLehman

    target, predicted, length = {(0, 4)}, {(0, 3)}, 5

    def graph(pairs):
        matrix = np.zeros((length, length), dtype=int)
        for a, b in pairs:
            matrix[a, b] = matrix[b, a] = 1
        return Graph(initialization_object=matrix, node_labels={i: str(i) for i in range(length)})

    kernel = WeisfeilerLehman(n_iter=n_iter, normalize=True, base_graph_kernel=VertexHistogram)
    kernel.fit_transform([graph(target)])
    expected = float(kernel.transform([graph(predicted)])[0][0])

    actual = score_wl(target, predicted, length=length, n_iter=n_iter)
    assert actual['wl'] == pytest.approx(expected)
    assert actual['loss'] == pytest.approx(1.0 - expected)
    assert actual['wl'] < 1.0


# ---- TM-score ---------------------------------------------------------------

needs_usalign = pytest.mark.skipif(
    not os.access(USALIGN, os.X_OK) or USALIGN.read_bytes()[:4] != b'\x7fELF',
    reason=f'no Linux USalign at {USALIGN} (see hpo/README.md, Scoring)',
)
SEQUENCE = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ'
THREE = dict(zip('ACDEFGHIKLMNPQRSTVWY', 'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU '
                                         'MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split()))
# A random C-alpha trace (3.8 A steps): unlike a helix, a shifted copy does not fit it.
_steps = np.random.default_rng(0).normal(size=(len(SEQUENCE), 3))
TRACE = np.cumsum(3.8 * _steps / np.linalg.norm(_steps, axis=1, keepdims=True), axis=0)


def ca_pdb(path, sequence=SEQUENCE, coordinates=TRACE, numbers=None):
    """A C-alpha-only PDB file of chain A."""
    numbers = numbers or range(1, len(sequence) + 1)
    path.write_text(''.join(
        f'ATOM  {serial:5d}  CA  {THREE[code]} A{number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n'
        for serial, (code, number, (x, y, z)) in enumerate(zip(sequence, numbers, coordinates), 1)
    ) + 'END\n')
    return path


@needs_usalign
def test_tm_of_the_reference_itself_is_one(tmp_path):
    reference = ca_pdb(tmp_path / 'reference.pdb')
    metrics = score_row(dict(reference=str(reference)), set(), objective='tm',
                        model=ca_pdb(tmp_path / 'model.pdb'))
    assert (metrics['tm'], metrics['loss']) == (pytest.approx(1.0), pytest.approx(0.0))


@needs_usalign
def test_tm_matches_residues_by_sequence_not_by_number(tmp_path):
    # PDB numbering starting at 3, low-pLDDT residues 10-15 dropped: still a perfect model.
    kept = [i for i in range(len(SEQUENCE)) if not 10 <= i < 16]
    reference = ca_pdb(tmp_path / 'reference.pdb', ''.join(SEQUENCE[i] for i in kept),
                       TRACE[kept], [i + 3 for i in kept])
    assert score_tm(ca_pdb(tmp_path / 'model.pdb'), reference)['tm'] == pytest.approx(1.0)


@needs_usalign
def test_tm_punishes_a_register_shift(tmp_path):
    # Every residue sits where its fourth neighbour belongs. A structural alignment
    # would match the shifted residues and score this close to 1.
    shifted = np.concatenate([TRACE[4:], TRACE[-1] + np.arange(1, 5)[:, None] * 3.8 * np.eye(3)[0]])
    metrics = score_tm(ca_pdb(tmp_path / 'model.pdb', coordinates=shifted), ca_pdb(tmp_path / 'reference.pdb'))
    assert metrics['tm'] < 0.5


@needs_usalign
def test_tm_rejects_a_reference_of_another_sequence(tmp_path):
    other = ca_pdb(tmp_path / 'reference.pdb', 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVW')
    with pytest.raises(ValueError, match='does not match'):
        score_tm(ca_pdb(tmp_path / 'model.pdb'), other)


# ---- Dataset loss -----------------------------------------------------------

def test_dataset_loss_weights_each_rna_equally():
    large = score_pairs([(i, i + 10) for i in range(10)], [(i, i + 10) for i in range(10)])
    small = score_pairs([(0, 1)], [])
    assert aggregate([large, small]) == pytest.approx(0.5)


@pytest.mark.parametrize('scores', [
    [],
    [{'loss': 0.2}, {'loss': float('nan')}],
    [{'loss': 0.2}, {'loss': float('inf')}],
    [{'loss': 0.2}, {'loss': -float('inf')}],
], ids=['empty', 'nan', 'inf', '-inf'])
def test_dataset_loss_rejects_empty_or_nonfinite_scores(scores):
    with pytest.raises(ValueError):
        aggregate(scores)
