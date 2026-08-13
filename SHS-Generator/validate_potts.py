"""Offline validation of the substitution models: checks that realised mutation rates match the
requested ones and measures how cleanly a covariation readout recovers the seeded base pairs."""

import argparse
import itertools
import logging
import random
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

HERE_DIR = Path(__file__).resolve().parent
if str(HERE_DIR) not in sys.path:
    sys.path.insert(0, str(HERE_DIR))

import shs_generator
from pair_map import PairMap

BASES = ("A", "C", "G", "U")
CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE = {("G", "U"), ("U", "G")}

APPROACHES = ["potts", "watson_crick_cov", "watson_crick", "covariance", "original", "none"]

# A hairpin, a three-hairpin junction, and a pseudoknot, so the comparison is not one structure
# deep. Every stem is fully complementary, so a non-canonical pair in the output is the model's
# doing rather than the query's.
CASES: List[Tuple[str, str, str]] = [
    ("hairpin",
     "GGCGCGAUAACGCGCC",
     "((((((....))))))"),
    ("junction",
     "GGCUUCGGCCGGCUUCGGCCGGCAUUCGAUGCCAAAA",
     "(((....)))(((....)))((((.....))))...."),
    ("pseudoknot",
     "GGCGAAAACGAACGCCAAACGU",
     "((((...[[[..))))...]]]"),
]


def make_args(approach: str, seq: str, structure: str, n: int, seed: int,
              indels: bool, coupling: float, wobble: float,
              rate_matching: bool) -> argparse.Namespace:
    argv = [
        "shs_generator.py",
        "--rna-seq", seq,
        "--protein-seq", "M",
        "--structure", structure,
        "--pair-mutation-approach", approach,
        "-N", str(n),
        "--seed", str(seed),
        "--potts-coupling", str(coupling),
        "--potts-wobble", str(wobble),
    ]
    if not rate_matching:
        argv.append("--potts-no-rate-matching")
    if not indels:
        for flag in ("--stem_single_insertion_prob", "--stem_long_insertion_prob",
                     "--stem_single_deletion_prob", "--stem_pair_deletion_prob",
                     "--loop_single_insertion_prob", "--loop_single_deletion_prob",
                     "--loop_long_insertion_prob", "--loop_long_deletion_prob"):
            argv += [flag, "0"]
    saved, sys.argv = sys.argv, argv
    try:
        return shs_generator.parse_args()
    finally:
        sys.argv = saved


def generate(approach: str, seq: str, structure: str, n: int, seed: int, indels: bool,
             coupling: float, wobble: float, rate_matching: bool) -> Tuple[List[str], PairMap]:
    args = make_args(approach, seq, structure, n, seed, indels, coupling, wobble, rate_matching)
    generator = shs_generator.MsaGenerator(args)
    generator.pair_map = generator.get_structure(seq)
    return generator.generate_msa(seq), generator.pair_map


# -- measurements ----------------------------------------------------------

def strip_insertions(row: str) -> str:
    """Drop lowercase insertion columns, leaving the alignment columns of the query."""
    return "".join(c for c in row if not c.islower())


def mutation_rates(msa: Sequence[str], query: str) -> np.ndarray:
    """Per-column fraction of non-query, non-gap rows among the non-gap rows."""
    rows = np.array([list(r) for r in msa[1:]])
    query_arr = np.array(list(query))
    kept = rows != "-"
    differs = (rows != query_arr) & kept
    denom = kept.sum(axis=0)
    return np.divide(differs.sum(axis=0), denom, out=np.zeros(len(query)), where=denom > 0)


def pair_composition(msa: Sequence[str], pairs: Sequence[Tuple[int, int]]) -> Dict[str, float]:
    canonical = wobble = other = 0
    for row in msa[1:]:
        for i, j in pairs:
            combo = (row[i], row[j])
            if "-" in combo:
                continue
            if combo in CANONICAL:
                canonical += 1
            elif combo in WOBBLE:
                wobble += 1
            else:
                other += 1
    total = canonical + wobble + other
    if total == 0:
        return {"canonical": float("nan"), "wobble": float("nan"), "non_pair": float("nan")}
    return {"canonical": canonical / total, "wobble": wobble / total, "non_pair": other / total}


def mutual_information(msa: Sequence[str]) -> np.ndarray:
    """Column-pair mutual information over the five-letter alphabet, in bits."""
    rows = np.array([list(r) for r in msa[1:]])
    n, length = rows.shape
    alphabet = BASES + ("-",)
    lookup = {b: k for k, b in enumerate(alphabet)}
    codes = np.vectorize(lookup.get)(rows)
    mi = np.zeros((length, length))
    for i in range(length):
        for j in range(i + 1, length):
            joint = np.zeros((len(alphabet), len(alphabet)))
            np.add.at(joint, (codes[:, i], codes[:, j]), 1.0)
            joint /= n
            px = joint.sum(axis=1)
            py = joint.sum(axis=0)
            nz = joint > 0
            expected = np.outer(px, py)
            mi[i, j] = mi[j, i] = float(np.sum(joint[nz] * np.log2(joint[nz] / expected[nz])))
    return mi


def auroc(scores: Sequence[float], labels: Sequence[bool]) -> float:
    """Mann-Whitney AUROC with ties handled by mid-ranks."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    positives, negatives = labels.sum(), (~labels).sum()
    if positives == 0 or negatives == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1, dtype=float)
    sorted_scores = scores[order]
    start = 0
    for end in range(1, len(scores) + 1):
        if end == len(scores) or sorted_scores[end] != sorted_scores[start]:
            ranks[order[start:end]] = ranks[order[start:end]].mean()
            start = end
    return float((ranks[labels].sum() - positives * (positives + 1) / 2) / (positives * negatives))


def pair_auroc(mi: np.ndarray, pairs: Sequence[Tuple[int, int]], length: int) -> float:
    truth = {(min(i, j), max(i, j)) for i, j in pairs}
    scores, labels = [], []
    for i, j in itertools.combinations(range(length), 2):
        scores.append(mi[i, j])
        labels.append((i, j) in truth)
    return auroc(scores, labels)


# -- driver ----------------------------------------------------------------

def evaluate(approach: str, name: str, seq: str, structure: str, depth: int, args) -> Dict[str, float]:
    msa, pair_map = generate(approach, seq, structure, depth, args.seed, args.indels,
                             args.coupling, args.wobble, not args.no_rate_matching)
    pairs = sorted(pair_map.unique_pairs)
    paired = sorted({p for pair in pairs for p in pair})
    unpaired = [i for i in range(len(seq)) if i not in set(paired)]

    aligned = [strip_insertions(row) for row in msa]
    rates = mutation_rates(aligned, seq)
    composition = pair_composition(aligned, pairs)
    mi = mutual_information(aligned)

    query_canonical = (sum((seq[i], seq[j]) in CANONICAL or (seq[i], seq[j]) in WOBBLE
                           for i, j in pairs) / len(pairs)) if pairs else float("nan")

    return {
        "case": name,
        "depth": depth,
        "approach": approach,
        "query_canonical": query_canonical,
        "rate_paired": float(np.mean(rates[paired])) if paired else float("nan"),
        "rate_unpaired": float(np.mean(rates[unpaired])) if unpaired else float("nan"),
        "canonical": composition["canonical"],
        "wobble": composition["wobble"],
        "non_pair": composition["non_pair"],
        "mi_auroc": pair_auroc(mi, pairs, len(seq)),
        "mi_at_pairs": float(np.mean([mi[i, j] for i, j in pairs])) if pairs else float("nan"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depths", type=int, nargs="+", default=[20, 100, 800],
                        help="MSA depths to evaluate. 20 and 100 are the manuscript's operating "
                             "regime; the largest is used for the rate-matching check")
    parser.add_argument("--repeats", type=int, default=5,
                        help="Independent seeds per configuration; results are averaged")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--coupling", type=float, default=6.0)
    parser.add_argument("--wobble", type=float, default=0.85)
    parser.add_argument("--no-rate-matching", action="store_true")
    parser.add_argument("--indels", action="store_true",
                        help="Leave the indel layer at its defaults instead of switching it off")
    parser.add_argument("--approaches", nargs="+", default=APPROACHES)
    parser.add_argument("--tolerance", type=float, default=0.02,
                        help="Allowed deviation of the realised paired mutation rate from the requested one")
    args = parser.parse_args()

    logging.getLogger().setLevel(logging.WARNING)
    requested = 0.2  # --mutation-rate-paired / --mutation-rate-unpaired defaults

    metrics = ["rate_paired", "rate_unpaired", "canonical", "wobble", "non_pair",
               "mi_auroc", "mi_at_pairs", "query_canonical"]
    rows = []
    for depth in args.depths:
        for name, seq, structure in CASES:
            for approach in args.approaches:
                runs = []
                for repeat in range(args.repeats):
                    args.seed = args.seed + repeat
                    random.seed(args.seed)
                    runs.append(evaluate(approach, name, seq, structure, depth, args))
                    args.seed -= repeat
                averaged = {"case": name, "depth": depth, "approach": approach}
                for metric in metrics:
                    averaged[metric] = float(np.nanmean([r[metric] for r in runs]))
                rows.append(averaged)

    header = (f"{'depth':>6}  {'case':<12}{'approach':<18}{'mut/pair':>9}{'mut/loop':>9}"
              f"{'canon':>8}{'wobble':>8}{'non-pair':>9}{'MI AUROC':>10}{'MI@pairs':>10}")
    print(header)
    print("-" * len(header))
    last = None
    for r in rows:
        if last is not None and (r["case"], r["depth"]) != last:
            print()
        last = (r["case"], r["depth"])
        print(f"{r['depth']:>6}  {r['case']:<12}{r['approach']:<18}{r['rate_paired']:>9.3f}"
              f"{r['rate_unpaired']:>9.3f}{r['canonical']:>8.3f}{r['wobble']:>8.3f}"
              f"{r['non_pair']:>9.3f}{r['mi_auroc']:>10.3f}{r['mi_at_pairs']:>10.3f}")

    print(f"\nrequested mutation rate: {requested}   (mean of {args.repeats} seeds per row)")
    for name, _, _ in CASES:
        share = next(r["query_canonical"] for r in rows if r["case"] == name)
        if share < 1.0:
            print(f"NOTE {name}: only {share:.0%} of the seeded pairs are canonical or wobble in the "
                  f"query, so its 'canon' column is bounded well below 1")
    deepest = max(args.depths)
    failures = [r for r in rows
                if r["approach"] == "potts" and r["depth"] == deepest
                and abs(r["rate_paired"] - requested) > args.tolerance]
    if failures:
        for r in failures:
            print(f"FAIL rate matching: {r['case']} paired rate {r['rate_paired']:.3f} "
                  f"deviates from {requested} by more than {args.tolerance}")
        return 1
    if "potts" in args.approaches:
        print("OK rate matching: realised paired mutation rates within tolerance of the requested rate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
