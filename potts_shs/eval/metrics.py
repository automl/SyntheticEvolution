"""Base-pair comparison metrics, matching the definitions already used in
evaluate_secondary_structure_data.py so numbers stay comparable with the published evaluation."""

from typing import Dict, Iterable, List, Sequence, Set, Tuple

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

    if not ref:
        # 72 of the 919 benchmark rows carry no pairs at all. Against an empty reference these
        # ratios are undefined rather than zero, and scoring them as zero drags every mean down.
        # The counts stay, so the rows remain visible in scores.csv.
        nan = float("nan")
        precision = recall = f1 = mcc = nan

    return {
        "tp": float(tp), "fp": float(fp), "fn": float(fn), "tn": float(tn),
        "precision": float(precision), "recall": float(recall),
        "f1": float(f1), "mcc": float(mcc),
        "n_pred": float(len(pred)), "n_ref": float(len(ref)),
    }


CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}


def _is_canonical(sequence: str, i: int, j: int) -> bool:
    if i >= len(sequence) or j >= len(sequence):
        return False
    return (sequence[i].upper(), sequence[j].upper()) in CANONICAL


def canonical_fraction(sequence: str, pairs: Iterable[Pair]) -> float:
    """Share of the given pairs that are Watson-Crick or wobble in this sequence."""
    pairs = list(_normalise(pairs))
    if not pairs:
        return float("nan")
    return sum(_is_canonical(sequence, i, j) for i, j in pairs) / len(pairs)


def filter_canonical(sequence: str, pairs: Iterable[Pair]) -> List[Pair]:
    """Keep only the Watson-Crick and wobble pairs, by sequence rather than by DSSR's label.

    Used to bring a reference pair list into the same universe as a canonical-only DSSR readout.
    """
    return sorted(p for p in _normalise(pairs) if _is_canonical(sequence, *p))


def summarise(rows: Sequence[Dict[str, float]], keys: Sequence[str]) -> Dict[str, float]:
    """Mean of each key over rows, skipping missing values."""
    out = {}
    for key in keys:
        values = [r[key] for r in rows if key in r and r[key] == r[key]]
        out[key] = sum(values) / len(values) if values else float("nan")
    return out
