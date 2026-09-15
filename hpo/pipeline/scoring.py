"""All-pair F1, ignoring interaction types. Positions are zero-based.

Pair order and duplicates do not matter. Both sets empty => F1=1.
An execution failure must raise upstream; it is never an empty prediction.
Replace score_pairs for a new per-RNA objective, aggregate for dataset weighting.
"""
import math

def score_pairs(target, predicted):
    """Compare predicted base pairs with the target structure using F1.

    Pair order, duplicate pairs, and interaction types are ignored. Both
    inputs are converted to sets of zero-based ``(i, j)`` pairs. The counts
    are defined as:

    * ``tp``: pairs present in both target and prediction
    * ``fp``: predicted pairs absent from the target
    * ``fn``: target pairs missing from the prediction

    The F1 score is calculated as ``2 * tp / (2 * tp + fp + fn)``. If both
    sets are empty, the prediction is considered perfect and ``f1`` is ``1``.
    The returned loss is ``1 - f1``, so lower loss is better.
    """
    target, predicted = set(target), set(predicted)

    # Set operations make ordering and duplicate pairs irrelevant.
    tp = len(target & predicted)
    fp, fn = len(predicted - target), len(target - predicted)

    # An empty target and prediction contain no disagreement, so score them
    # as a perfect match instead of dividing by zero.
    f1 = 2 * tp / (2 * tp + fp + fn) if target or predicted else 1.0
    return dict(tp=tp, fp=fp, fn=fn, f1=f1, loss=1.0-f1)


def aggregate(scores):
    """Return the arithmetic mean of finite per-RNA losses.

    Every score dictionary must contain a finite ``loss`` value. The dataset
    objective is therefore ``sum(losses) / number_of_RNAs`` and gives every
    RNA equal weight. Empty or non-finite input is rejected.
    """
    values = [s['loss'] for s in scores]
    if not values or not all(math.isfinite(x) for x in values):
        raise ValueError('Require a finite score for every RNA')

    # Equal weighting prevents RNAs with more base pairs from dominating the
    # dataset-level objective.
    return sum(values) / len(values)
