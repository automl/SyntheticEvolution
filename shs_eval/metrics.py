"""Base-pair comparison metrics, matching the definitions already used in
evaluate_secondary_structure_data.py so numbers stay comparable with the published evaluation."""

from typing import Dict, Iterable, Sequence, Set, Tuple

Pair = Tuple[int, int]


def _normalise(pairs: Iterable[Pair]) -> Set[Pair]:
    return {(min(int(a), int(b)), max(int(a), int(b))) for a, b in pairs if int(a) != int(b)}


def pair_metrics(predicted: Iterable[Pair], reference: Iterable[Pair],
                 length: int) -> Dict[str, float]:
    """Precision/recall/F1/MCC over the set of base pairs.

    True negatives are counted over the upper triangle of the contact matrix, which is what
    pairs2mat-based scoring in the published evaluation effectively does.
    """
    pred, ref = _normalise(predicted), _normalise(reference)
    tp = len(pred & ref)
    fp = len(pred - ref)
    fn = len(ref - pred)
    possible = length * (length - 1) // 2
    tn = max(possible - tp - fp - fn, 0)

    eps = 1e-8
    precision = tp / (tp + fp + eps)
    recall = tp / (tp + fn + eps)
    f1 = 2 * tp / (2 * tp + fp + fn + eps)
    denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    mcc = (tp * tn - fp * fn) / (denom + eps)
    return {
        "tp": float(tp), "fp": float(fp), "fn": float(fn), "tn": float(tn),
        "precision": float(precision), "recall": float(recall),
        "f1": float(f1), "mcc": float(mcc),
        "n_pred": float(len(pred)), "n_ref": float(len(ref)),
    }


def canonical_fraction(sequence: str, pairs: Iterable[Pair]) -> float:
    """Share of the given pairs that are Watson-Crick or wobble in this sequence."""
    canonical = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
    pairs = list(_normalise(pairs))
    if not pairs:
        return float("nan")
    hits = 0
    for i, j in pairs:
        if i < len(sequence) and j < len(sequence):
            hits += (sequence[i].upper(), sequence[j].upper()) in canonical
    return hits / len(pairs)


def summarise(rows: Sequence[Dict[str, float]], keys: Sequence[str]) -> Dict[str, float]:
    """Mean of each key over rows, skipping missing values."""
    out = {}
    for key in keys:
        values = [r[key] for r in rows if key in r and r[key] == r[key]]
        out[key] = sum(values) / len(values) if values else float("nan")
    return out
