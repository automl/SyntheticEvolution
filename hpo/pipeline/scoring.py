"""Scoring objectives for RNA structure prediction.

The default objective is the existing pair-level F1 loss.  The optional
``wl`` objective is an implementation of Weisfeiler-Lehman (WL) graph-kernel score
based on evaluate_secondary_structure_data.py in this repo:

* each RNA position is a graph vertex;
* an edge represents a predicted/target base pair;
* vertices use their position as the default label;
* the normalized GraKeL ``WeisfeilerLehman`` kernel is used with a
  ``VertexHistogram`` base kernel and five WL iterations;
* the kernel similarity is converted to a minimization loss as ``1 - score``.

GraKeL is imported lazily so that the existing F1 objective remains usable in
environments where the optional WL dependency is not installed.
"""

import math
from typing import Iterable


SCORING_OBJECTIVES = ("f1", "wl")
DEFAULT_WL_ITERATIONS = 5


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

    tp = len(target & predicted)
    fp, fn = len(predicted - target), len(target - predicted)
    f1 = 2 * tp / (2 * tp + fp + fn) if target or predicted else 1.0
    return dict(tp=tp, fp=fp, fn=fn, f1=f1, loss=1.0 - f1)


def _pairs_to_graph_matrix(pairs: Iterable[tuple[int, int]], length: int):
    """Build the symmetric base-pair adjacency matrix used by the reference."""
    import numpy as np

    matrix = np.zeros((length, length), dtype=int)
    for a, b in pairs:
        matrix[a, b] = 1
        matrix[b, a] = 1
    return matrix


def _mat_to_graph(matrix, node_labels=None):
    """Convert an adjacency matrix to the GraKeL graph used by the reference."""
    try:
        from grakel import Graph
    except ImportError as exc:
        raise RuntimeError(
            "The 'wl' scoring objective requires GraKeL. Install the project "
            "dependency before running trials with scoring_objective: wl."
        ) from exc

    if node_labels is None:
        node_labels = {i: str(i) for i in range(matrix.shape[0])}

    return Graph(
        initialization_object=matrix.astype(int),
        node_labels=node_labels,
    )


def _get_wl_kernel(*, n_iter=DEFAULT_WL_ITERATIONS, normalize=True):
    """Create the same GraKeL WL kernel as the reference implementation."""
    try:
        from grakel.kernels import WeisfeilerLehman, VertexHistogram
    except ImportError as exc:
        raise RuntimeError(
            "The 'wl' scoring objective requires GraKeL. Install the project "
            "dependency before running trials with scoring_objective: wl."
        ) from exc

    return WeisfeilerLehman(
        n_iter=n_iter,
        normalize=normalize,
        base_graph_kernel=VertexHistogram,
    )


def score_wl(target, predicted, length, *, n_iter=DEFAULT_WL_ITERATIONS):
    """Score two RNA structures with the normalized WL graph kernel.

    This mirrors the reference's ``graph_distance_score_from_matrices``:
    the target graph is fitted first and the predicted graph is then passed to
    ``transform``.  The returned ``wl`` value is therefore a kernel similarity
    (higher is better), while ``loss`` converts it to the minimization form
    expected by the HPO pipeline.
    """
    target = set(target)
    predicted = set(predicted)

    target_graph = _mat_to_graph(_pairs_to_graph_matrix(target, length))
    predicted_graph = _mat_to_graph(_pairs_to_graph_matrix(predicted, length))

    kernel = _get_wl_kernel(n_iter=n_iter, normalize=True)
    kernel.fit_transform([target_graph])
    similarity = float(kernel.transform([predicted_graph])[0][0])

    return {
        "wl": similarity,
        "loss": 1.0 - similarity,
    }


def aggregate(scores):
    """Return the arithmetic mean of finite per-RNA losses.

    Every score dictionary must contain a finite ``loss`` value. The dataset
    objective is ``sum(losses) / number_of_RNAs`` and gives every RNA equal
    weight. Empty or non-finite input is rejected.
    """
    values = [s["loss"] for s in scores]
    if not values or not all(math.isfinite(x) for x in values):
        raise ValueError("Require a finite score for every RNA")
    return sum(values) / len(values)


def score_row(row, predicted, *, objective="f1", wl_iterations=DEFAULT_WL_ITERATIONS):
    """Score one dataset row using the selected objective.

    ``objective='f1'`` preserves the existing behavior. ``objective='wl'``
    uses the reference GraKeL Weisfeiler-Lehman graph kernel.
    """
    if objective == "f1":
        return score_pairs(map(tuple, row["pairs"]), predicted)

    if objective == "wl":
        return score_wl(
            map(tuple, row["pairs"]),
            predicted,
            len(row["sequence"]),
            n_iter=wl_iterations,
        )

    raise ValueError(
        f"Unknown scoring objective {objective!r}; expected one of "
        f"{', '.join(SCORING_OBJECTIVES)}"
    )
