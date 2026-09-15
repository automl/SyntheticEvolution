"""Dataset weighting and non-finite losses, supplementing existing F1 tests."""
import pytest
from hpo.pipeline.scoring import aggregate, score_pairs


def test_dataset_loss_weights_each_rna_equally():
    large = score_pairs([(i, i + 10) for i in range(10)], [(i, i + 10) for i in range(10)])
    small = score_pairs([(0, 1)], [])
    assert aggregate([large, small]) == pytest.approx(0.5)


@pytest.mark.parametrize('bad', [float('nan'), float('inf'), -float('inf')])
def test_nonfinite_loss_rejected(bad):
    with pytest.raises(ValueError):
        aggregate([{'loss': 0.2}, {'loss': bad}])


def test_empty_target_with_spurious_prediction_has_full_loss():
    assert score_pairs([], [(0, 1)]) == dict(tp=0, fp=1, fn=0, f1=0.0, loss=1.0)
