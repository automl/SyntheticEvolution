"""Scoring objectives (F1, Weisfeiler-Lehman) and the dataset loss."""
import numpy as np
import pytest

from hpo.pipeline.scoring import aggregate, score_pairs, score_row, score_wl


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
