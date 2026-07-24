"""Recover generator parameters from a synthetic MSA."""

import argparse
import logging
import json
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

VALID_BASES = list("AGCU-")


# ---------------------------------------------------------------------------
# Layer 1: parsing -> MsaData
# ---------------------------------------------------------------------------

def _split_row(row: str) -> Tuple[List[str], List[int]]:
    """Split one MSA row into (aligned_columns, insertions)."""
    aligned: List[str] = []
    insertions: List[int] = []
    sum = 0
    for c in row:
        if c.islower():
            sum += 1
        else:
            aligned.append(c)
            insertions.append(sum)
            sum = 0
    return aligned, insertions


def _extract_fasta(content: str) -> str:
    """Return the FASTA MSA text from either raw FASTA or an AF3 JSON string."""
    try:
        data = json.loads(content)
        if not isinstance(data, dict):
            raise json.decoder.JSONDecodeError("not an object", content, 0)
    except json.decoder.JSONDecodeError:
        return content
    all_fasta = [msa for c in data.get("sequences", []) if "rna" in c and (msa := c["rna"].get("unpairedMsa"))]
    if len(all_fasta) == 0:
        logging.error("No unpaired MSA found for any rna in the json input")
        raise ValueError("Invalid Input, aborting.")
    if len(all_fasta) > 1:
        logging.warning("Input json contains multiple rna chains with unpaired msa, using the first one")
    return all_fasta[0]


@dataclass(frozen=True, eq=False)
class MsaData:
    """Immutable parsed MSA."""
    seq: str
    aligned: np.ndarray
    insertions: Tuple[Tuple[Tuple[int, str], ...], ...]
    pairs: Optional[Dict[int, int]] = None

    @classmethod
    def from_path(cls, path: str, pairs: Optional[Dict[int, int]] = None) -> "MsaData":
        with open(path, "r") as f:
            content = f.read()
        return cls.from_text(content, pairs=pairs)

    @classmethod
    def from_text(cls, content: str, pairs: Optional[Dict[int, int]] = None) -> "MsaData":
        fasta = _extract_fasta(content)
        # Sequence lines only (drop headers); first is the seq, rest samples.
        lines = fasta.splitlines()[1::2]
        if len(lines) < 2:
            logging.error("Recovered MSA is empty (need a query and >=1 sample)")
            raise ValueError("Invalid Input, aborting.")
        seq_list, *sample_rows = lines
        seq = ''.join(seq_list)

        aligned_rows: List[List[str]] = []
        insertions: List[Tuple[Tuple[int, str], ...]] = []
        for row in sample_rows:
            cols, ins = _split_row(row)
            aligned_rows.append(cols)
            insertions.append(tuple(ins))

        if not all(len(r) == len(seq) for r in aligned_rows):
            logging.error("Sequences in MSA have different aligned lengths or the json was invalid")
            raise ValueError("Invalid Input, aborting.")
        aligned = np.array(aligned_rows)
        if not np.all(np.isin(aligned, VALID_BASES)):
            bad = np.unique(aligned[~np.isin(aligned, VALID_BASES)])
            logging.error("Invalid characters in MSA: %s", bad)
            raise ValueError("Invalid Input, aborting.")

        return cls(seq=seq, aligned=aligned,
                   insertions=tuple(insertions), pairs=pairs)


# ---------------------------------------------------------------------------
# Layer 2: derived features (computed at most once, shared by all estimators)
# ---------------------------------------------------------------------------

class Features:
    """Lazily-cached features of an MsaData."""

    def __init__(self, data: MsaData) -> None:
        self.data = data

    @cached_property
    def seq_arr(self) -> np.ndarray:
        return np.array(list(self.data.seq))

    @cached_property
    def diff_mask(self) -> np.ndarray:
        """(N, L) bool: sample column differs from the seq."""
        return self.data.aligned != self.seq_arr

    @cached_property
    def deletion_mask(self) -> np.ndarray:
        """(N, L) bool: column is a deletion ('-')."""
        return self.data.aligned == '-'

    @cached_property
    def mutation_mask(self) -> np.ndarray:
        """(N, L) bool: a real base->base substitution (differs AND not a gap)."""
        return self.diff_mask & ~self.deletion_mask

    @cached_property
    def col_deletion_rate(self) -> np.ndarray:
        """Per-column fraction of rows that got deleted"""
        return self.deletion_mask.mean(axis=0)

    @cached_property
    def col_mut_rate(self) -> np.ndarray:
        """Per-column fraction of rows that differ from the seq"""
        return self.diff_mask.mean(axis=0)

    @cached_property
    def covariance(self) -> np.ndarray:
        """Covariance of the mutation mask matrix, with the
        diagonal replaced by each column's mutation rate."""
        cov = np.cov(self.diff_mask.astype(int), rowvar=False)
        np.fill_diagonal(cov, self.col_mut_rate)
        return cov


# ---------------------------------------------------------------------------
# Layer 3: estimators (one per recovered parameter) + driver
# ---------------------------------------------------------------------------

@dataclass
class Estimate:
    """A single recovered parameter: its value plus whatever supporting detail
    the estimator wants to expose (sample sizes, confidence, per-column data)."""
    name: str
    value: float
    detail: Dict[str, Any] = field(default_factory=dict)


# An estimator reads the shared Features and the dict of already-recovered
# results (so pair-conditional estimators can depend on a recovered structure),
# and returns one Estimate. Register new parameters by adding to ESTIMATORS.
Estimator = Callable[[Features, Dict[str, Estimate]], Estimate]


def est_n(feat: Features, results: Dict[str, Estimate]) -> Estimate:
    name = "n_sequences"
    value = float(feat.data.aligned.shape[0] + 1) # N includes the query but the aligned data doesn't
    return Estimate(name, value)


def est_mutation_rate_unpaired(feat: Features, results: Dict[str, Estimate]) -> Estimate:
    """Mutation rate at unpaired positions -> recovers --mutation-rate-unpaired."""
    pairs = feat.data.pairs or {}
    unpaired_cols = [i for i in range(len(feat.data.seq)) if i not in pairs and bool(pairs)]
    mutations = feat.mutation_mask[:, unpaired_cols]
    no_mutation = ~feat.deletion_mask[:, unpaired_cols]
    denom = int(no_mutation.sum())
    rate = float(mutations.sum() / denom) if denom else float("nan")
    detail: Dict[str, Any] = {"unpaired_cols": len(unpaired_cols), "structure_known": bool(pairs)}
    return Estimate("mutation_rate_unpaired", rate, detail)


def est_mutation_rate_paired(feat: Features, results: Dict[str, Estimate]) -> Estimate:
    """Mutation rate at paired (stem) positions -> recovers --mutation-rate-paired.

    Note that mutation rates are underestimated when using the watson-crick pair mutation
    method as one of the bases can stay the same during mutation"""
    pairs = feat.data.pairs or {}
    paired_cols = [i for i in range(len(feat.data.seq)) if i in pairs]
    mutations = feat.mutation_mask[:, paired_cols]
    not_deleted = ~feat.deletion_mask[:, paired_cols]
    denom = int(not_deleted.sum())
    rate = float(mutations.sum() / denom) if denom else float("nan")
    detail: Dict[str, Any] = {"paired_cols": len(paired_cols), "structure_known": bool(pairs)}
    return Estimate("mutation_rate_paired", rate, detail)


# Registry, ordered so that structure-dependent estimators can read the results
# of earlier ones. Append new estimators here.
ESTIMATORS: List[Estimator] = [
    est_n,
    est_mutation_rate_unpaired,
    est_mutation_rate_paired,
]


def analyze_data(data: MsaData,
                 estimators: Optional[Sequence[Estimator]] = None) -> Dict[str, Estimate]:
    """analyze() for an already-parsed MsaData (e.g. straight from the notebook
    where the generator handed you its pairs)."""
    feat = Features(data)
    results: Dict[str, Estimate] = {}
    for est in (estimators if estimators is not None else ESTIMATORS):
        e = est(feat, results)
        results[e.name] = e
    return results


# ---------------------------------------------------------------------------
# Plotting / CLI
# ---------------------------------------------------------------------------

def plot_covariance(seq: str, cov: np.ndarray, title: str = "Recovered covariance and mutation rate",
                    show: bool = True, out: Optional[str] = None,
                    vmin: Optional[float] = None, vmax: Optional[float] = None,
                    show_values=False):
    """Heatmap of the recovered covariance matrix (see Features.covariance).
    
    Off-diagonal cell (i, j) is the covariance between positions i and j across the 
    MSA. The main diagonal is NOT a variance but each column's mutation rate"""

    fig, ax = plt.subplots(figsize=(cov.shape[0], cov.shape[1]))
    im = ax.imshow(cov, cmap="viridis", vmin=vmin, vmax=vmax)

    if show_values:
        for i in range(cov.shape[0]):
            for j in range(cov.shape[1]):
                ax.text(j, i, f"{round(cov[i, j], 2)}\n{seq[i]+seq[j]}",
                        ha="center", va="center", color="w")
    ax.set_title(title)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        logging.info("Figure written to: %s", out)
    if show:
        plt.show()
    return fig

def plot_covariance_classification(seq: str, cov: np.ndarray, pairs: Optional[Dict[int, int]] = None,
                                   title: str = "Spread of recovered covariance",
                                   show: bool = True, out: Optional[str] = None):
    """Scatter of every off-diagonal covariance value, colored by base-pair status.

    Each point is one possible combination of (i, j). If the pairs are known,
    paired combinations are colored green and unpaired combinations red.
    """
    L = cov.shape[0]
    i, j = np.triu_indices(L, k=1)
    x = np.array([seq[a] + seq[b] for a, b in zip(i, j)])
    y = cov[i, j]

    fig, ax = plt.subplots()
    if pairs:
        is_pair = np.array([pairs.get(i) == j for i, j in zip(i, j)])
        ax.scatter(x[~is_pair], y[~is_pair], c="red", label="Unpaired")
        ax.scatter(x[is_pair], y[is_pair], c="green", label="Paired")
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax.scatter(x, y, c="black")

    ax.set_xlabel("base pair")
    ax.set_ylabel("recovered covariance")
    ax.set_title(title)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        logging.info("Figure written to: %s", out)
    if show:
        plt.show()
    return fig

def plot_deletion_rate(seq: str, del_rate: np.ndarray, pairs: Optional[Dict[int, int]] = None,
                                   title: str = "Deletion rate for each position",
                                   show: bool = True, out: Optional[str] = None):
    """Bar plot of the deletion rates for each position, colored by base-pair status."""
    x = np.arange(len(del_rate))
    y = np.asarray(del_rate)
    print(x)
    print(y)

    fig, ax = plt.subplots()
    if pairs:
        paired_mask = np.isin(x, list(pairs.keys()))
        ax.bar(x[paired_mask], y[paired_mask], color="green", label="paired")
        ax.bar(x[~paired_mask], y[~paired_mask], color="red", label="unpaired")
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax.bar(x, y, color="black")

    ax.set_xlabel("position in the sequence")
    ax.set_ylabel("deletion rate")
    ax.set_title(title)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        logging.info("Figure written to: %s", out)
    if show:
        plt.show()
    return fig


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Recover generator parameters from a synthetic MSA.")
    p.add_argument("-i", "--input", help="Path to AF3 JSON or FASTA MSA file.", required=True)
    p.add_argument("-o", "--out", help="Output PNG path for the covariance heatmap.")
    p.add_argument("-s", "--show", action="store_true", help="Display the figure.")
    p.add_argument("-t", "--title", default="", help="Figure title.")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    data = MsaData.from_path(args.input)
    feat = Features(data)

    report = analyze_data(data)
    logging.info("Recovered parameters:")
    for name, est in report.items():
        logging.info("  %-24s = %-8s %s", name, round(est.value, 4), est.detail or "")

    title = args.title or "Recovered (covariation)"
    plot_covariance(data.seq, feat.covariance, title=title, show=args.show, out=args.out)


if __name__ == "__main__":
    main()
