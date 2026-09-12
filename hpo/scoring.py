"""All-pair F1, ignoring interaction types. Positions are zero-based.

Pair order and duplicates do not matter. Both sets empty => F1=1.
An execution failure must raise upstream; it is never an empty prediction.
Replace score_pairs for a new per-RNA objective, aggregate for dataset weighting.
"""
import math


def normalize_pairs(pairs, size):
    result = set()
    for pair in pairs:
        if len(pair) != 2 or any(type(x) is not int for x in pair):
            raise ValueError('Pairs must contain exactly two integer positions')
        a, b = pair
        if not (0 <= a < size and 0 <= b < size) or a == b:
            raise ValueError(f'Invalid pair {pair} for length {size}')
        result.add(tuple(sorted((a, b))))
    return result


def score_pairs(target, predicted):
    target, predicted = set(target), set(predicted)
    tp = len(target & predicted)
    fp, fn = len(predicted - target), len(target - predicted)
    f1 = 2 * tp / (2 * tp + fp + fn) if target or predicted else 1.0
    return dict(tp=tp, fp=fp, fn=fn, f1=f1, loss=1.0-f1)


def aggregate(scores):
    values = [s['loss'] for s in scores]
    if not values or not all(math.isfinite(x) for x in values):
        raise ValueError('Require a finite score for every RNA')
    return sum(values) / len(values)
