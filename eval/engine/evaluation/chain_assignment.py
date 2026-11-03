"""
Shared one-to-one chain assignment for the per-chain "*_AllChains" RNA metrics.

Both extractEvalMetrics.py (RNABench env, TM) and ost_calcRNA_AllChains.py (ost env, LDDT) build
an all-vs-all ref x pred score matrix and must reduce it to one value per predicted chain. The old
code took the greedy per-column max, which lets several predicted chains claim the SAME reference
chain (a one-to-one violation) and systematically over-scores multichain / homomeric complexes.

This module computes a proper one-to-one maximum-weight assignment (Hungarian) so every predicted
chain is matched to a distinct reference chain, and returns the matched pairs plus the assigned
(raw) value per predicted chain. To keep TM / LDDT / ICS referring to the SAME reference partner,
the TM step (which runs first) persists the assignment as a chain map that the OST step reuses;
if that map is absent, each step falls back to assigning on its own matrix.

Kept dependency-light on purpose: uses scipy.optimize.linear_sum_assignment when available (the
same routine the plotting pipeline used), otherwise an exact small-n / greedy-descending fallback,
so it works in both conda envs without adding a hard dependency.
"""

from itertools import permutations


def _assign_scipy(matrix):
    """Optimal max-weight one-to-one assignment via scipy. Returns list of (row, col) or None."""
    try:
        import numpy as np
        from scipy.optimize import linear_sum_assignment
    except Exception:
        return None
    arr = np.array(matrix, dtype=float)
    if arr.size == 0:
        return []
    # linear_sum_assignment MINIMISES cost; negate to maximise the score. NaNs -> very low score.
    cost = -np.nan_to_num(arr, nan=-1e9)
    rows, cols = linear_sum_assignment(cost)
    return list(zip(rows.tolist(), cols.tolist()))


def _assign_fallback(matrix):
    """Dependency-free assignment: exact (permutations) for small n, greedy-descending otherwise."""
    n_ref = len(matrix)
    n_pred = len(matrix[0]) if n_ref else 0
    if n_ref == 0 or n_pred == 0:
        return []

    def _val(i, j):
        v = matrix[i][j]
        try:
            v = float(v)
        except (TypeError, ValueError):
            return 0.0
        return v if v == v else 0.0  # NaN -> 0.0

    k = min(n_ref, n_pred)
    # Exact search is cheap while the smaller dimension is tiny (RNA chain counts almost always are).
    if k <= 7:
        if n_ref <= n_pred:
            rows = list(range(n_ref))
            best, best_pairs = None, None
            for cols in permutations(range(n_pred), n_ref):
                total = sum(_val(i, cols[idx]) for idx, i in enumerate(rows))
                if best is None or total > best:
                    best, best_pairs = total, list(zip(rows, cols))
            return best_pairs or []
        else:
            cols = list(range(n_pred))
            best, best_pairs = None, None
            for rows in permutations(range(n_ref), n_pred):
                total = sum(_val(rows[idx], j) for idx, j in enumerate(cols))
                if best is None or total > best:
                    best, best_pairs = total, list(zip(rows, cols))
            return best_pairs or []

    # Large n: greedy-descending. Not guaranteed optimal but still strictly one-to-one.
    cells = sorted(((_val(i, j), i, j) for i in range(n_ref) for j in range(n_pred)),
                   reverse=True)
    used_r, used_c, pairs = set(), set(), []
    for _v, i, j in cells:
        if i in used_r or j in used_c:
            continue
        used_r.add(i); used_c.add(j); pairs.append((i, j))
        if len(pairs) == k:
            break
    return pairs


def assign_pairs(matrix):
    """
    Maximum-weight one-to-one assignment of a ref x pred score matrix.

    matrix: list-of-lists, matrix[ref_i][pred_j] = score (0.0 / NaN allowed).
    Returns a list of (ref_i, pred_j) index pairs, one per matched chain (min(n_ref, n_pred) pairs).
    """
    if not matrix or not matrix[0]:
        return []
    pairs = _assign_scipy(matrix)
    if pairs is None:
        pairs = _assign_fallback(matrix)
    return pairs


def assigned_by_pred_chain(matrix, ref_chains, pred_chains, positive_only=True):
    """
    Reduce a ref x pred matrix to the assigned (raw) value per PREDICTED chain, plus the chain map.

    ref_chains / pred_chains: sequences whose element [0] is the chain id (e.g. (chain_id, length)),
    matching the matrix row/column order.

    Returns (values, chain_map):
      values   -> list of single-key dicts keyed by predicted chain id, e.g. [{"B": 0.47}, {"C": 0.56}],
                  carrying the RAW assigned score (no rounding). Predicted chains left unmatched, or
                  matched to a non-positive score when positive_only, are omitted.
      chain_map -> {predicted_chain_id: reference_chain_id} for the matched pairs, so a later step can
                  report other metrics against the SAME reference partner.
    """
    values, chain_map = [], {}
    if not ref_chains or not pred_chains:
        return values, chain_map
    ref_ids = [c[0] if isinstance(c, (list, tuple)) else c for c in ref_chains]
    pred_ids = [c[0] if isinstance(c, (list, tuple)) else c for c in pred_chains]
    by_pred = {}
    for i, j in assign_pairs(matrix):
        if i >= len(ref_ids) or j >= len(pred_ids):
            continue
        by_pred[j] = i
    for j, pred_id in enumerate(pred_ids):
        if j not in by_pred:
            continue
        i = by_pred[j]
        v = matrix[i][j]
        try:
            v = float(v)
        except (TypeError, ValueError):
            continue
        if v != v:  # NaN
            continue
        if positive_only and not (v > 0):
            continue
        values.append({pred_id: v})
        chain_map[pred_id] = ref_ids[i]
    return values, chain_map


def assigned_from_map(matrix, ref_chains, pred_chains, chain_map, positive_only=True):
    """
    Reduce a ref x pred matrix using a PRE-COMPUTED {pred_chain: ref_chain} assignment (e.g. from the
    TM step), so this metric refers to the same reference partner. Falls back to assign_by_pred_chain
    for any predicted chain not covered by the map. Returns the list of single-key dicts (raw values).
    """
    if not chain_map:
        vals, _ = assigned_by_pred_chain(matrix, ref_chains, pred_chains, positive_only)
        return vals
    ref_ids = [c[0] if isinstance(c, (list, tuple)) else c for c in ref_chains]
    pred_ids = [c[0] if isinstance(c, (list, tuple)) else c for c in pred_chains]
    ref_index = {rid: i for i, rid in enumerate(ref_ids)}
    out = []
    for j, pred_id in enumerate(pred_ids):
        ref_id = chain_map.get(pred_id)
        if ref_id is None or ref_id not in ref_index:
            continue
        i = ref_index[ref_id]
        if i >= len(matrix) or j >= len(matrix[i]):
            continue
        v = matrix[i][j]
        try:
            v = float(v)
        except (TypeError, ValueError):
            continue
        if v != v:
            continue
        if positive_only and not (v > 0):
            continue
        out.append({pred_id: v})
    return out
