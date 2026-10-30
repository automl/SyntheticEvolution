"""Recover generator parameters from a synthetic MSA."""

import argparse
import logging
import json
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt

from pair_map import PairMap

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
    insertions.append(sum)
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
    raw: np.ndarray
    aligned: np.ndarray
    insertions: np.ndarray
    pairs: Optional[PairMap] = None

    @classmethod
    def from_path(cls, path: str, pairs: Optional[PairMap] = None) -> "MsaData":
        with open(path, "r") as f:
            content = f.read()
        return cls.from_text(content, pairs=pairs)

    @classmethod
    def from_text(cls, content: str, pairs: Optional[PairMap] = None) -> "MsaData":
        fasta = _extract_fasta(content)
        # Sequence lines only (drop headers); first is the seq, rest samples.
        lines = fasta.splitlines()[1::2]
        if len(lines) < 2:
            logging.error("Recovered MSA is empty (need a query and >=1 sample)")
            raise ValueError("Invalid Input, aborting.")
        seq_list, *sample_rows = lines
        seq = ''.join(seq_list)
        raw = np.array(sample_rows)

        aligned_rows: List[List[str]] = []
        insertion_rows: List[List[int]] = []
        for row in sample_rows:
            cols, ins = _split_row(row)
            aligned_rows.append(cols)
            insertion_rows.append(ins)

        if not all(len(r) == len(seq) for r in aligned_rows):
            logging.error("Sequences in MSA have different aligned lengths or the json was invalid")
            raise ValueError("Invalid Input, aborting.")
        aligned = np.array(aligned_rows)
        insertions = np.array(insertion_rows)
        if not np.all(np.isin(aligned, VALID_BASES)):
            bad = np.unique(aligned[~np.isin(aligned, VALID_BASES)])
            logging.error("Invalid characters in MSA: %s", bad)
            raise ValueError("Invalid Input, aborting.")

        return cls(seq=seq, raw=raw, aligned=aligned, insertions=insertions, pairs=pairs)


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
    def col_insertion_mean(self) -> np.ndarray:
        """(N + 1) Per-column average number of insertions in the gap before it"""
        return self.data.insertions.mean(axis=0)

    @cached_property
    def col_mut_rate(self) -> np.ndarray:
        """Per-column fraction of rows that differ from the seq"""
        return self.diff_mask.mean(axis=0, where=~self.deletion_mask)

    @cached_property
    def covariance(self) -> np.ndarray:
        """Covariance of the mutation mask matrix, with the
        diagonal replaced by each column's mutation rate."""
        cov = np.cov(self.mutation_mask.astype(int), rowvar=False)
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
    unpaired_cols = [i for i in range(len(feat.data.seq)) if not pairs.is_paired(i) and bool(pairs)]
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
    paired_cols = [i for i in range(len(feat.data.seq)) if pairs.is_paired(i)]
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
    fig, ax = plt.subplots(figsize=(20, 20)) if show_values else plt.subplots()
    im = ax.imshow(cov, cmap="viridis", vmin=vmin, vmax=vmax)

    if show_values:
        for i in range(cov.shape[0]):
            for j in range(cov.shape[1]):
                ax.text(j, i, f"{round(cov[i, j], 2)}\n{seq[i]+seq[j]}",
                        ha="center", va="center", color="w")
    fig.colorbar(im)
    ax.set_title(title)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        logging.info("Figure written to: %s", out)
    if show:
        plt.show()
    return fig

def plot_covariance_classification(seq: str, cov: np.ndarray, pairs: Optional[PairMap] = None,
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
        interactions = np.array([pairs.interaction(i, j) for i, j in zip(i, j)])
        is_pair = interactions > 0
        ax.scatter(x[~is_pair], y[~is_pair], c="red", label="Unpaired", s=5)
        ax.scatter(x[is_pair], y[is_pair], c="green", label="Paired", s=interactions[is_pair] * 200 + 5)
        for xp, yp, strength in zip(x[is_pair], y[is_pair], interactions[is_pair]):
            ax.annotate(f"{strength:.2f}", (xp, yp), textcoords="offset points", xytext=(6, 3), fontsize=6)
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

def plot_deletions_and_insertion(deletions: np.ndarray, insertions: np.ndarray,
                                    pairs: Optional[PairMap] = None,
                                   title: str = "Deletion rate and average insertion count for each position",
                                   ylable: str = "deletion rate / average insertion count",
                                   show: bool = True, out: Optional[str] = None):
    """Bar plot of the deletion rates for each position, colored by base-pair status."""
    x_del = np.arange(len(deletions))
    y_del = np.asarray(deletions)
    x_ins = np.arange(-0.5, len(insertions)-0.5)
    y_ins = np.asarray(insertions)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(x_ins, y_ins, color="navy", width=0.3, label="insertions")
    if pairs:
        paired_mask = np.array([pairs.is_paired(i) for i, _ in enumerate(deletions)])
        ax.bar(x_del[paired_mask], y_del[paired_mask], color="tab:orange", width=0.7, label="deletions - paired")
        ax.bar(x_del[~paired_mask], y_del[~paired_mask], color="crimson", width=0.7, label="deletions - unpaired")
    else:
        ax.bar(x_del, y_del, color="tomato", width=0.7, label="deletions")

    plt.xticks(x_del)

    ax.legend()
    ax.set_xlabel("position in the sequence")
    ax.set_ylabel(ylable)
    ax.set_title(title)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        logging.info("Figure written to: %s", out)
    if show:
        plt.show()
    return fig


def plot_mutation_rate_comparison(input_rates: np.ndarray, output_rates: np.ndarray,
                                  title: str = "Input vs recovered mutation rates",
                                  show: bool = True, out: Optional[str] = None):
    """Compare input and recovered per-position mutation rates.

    Returns the figure together with summary accuracy metrics.
    """
    input_rates = np.asarray(input_rates, dtype=float)
    output_rates = np.asarray(output_rates, dtype=float)

    length = min(len(input_rates), len(output_rates))
    if length == 0:
        raise ValueError("Need at least one mutation rate to compare.")

    input_rates = input_rates[:length]
    output_rates = output_rates[:length]

    mae = float(np.mean(np.abs(input_rates - output_rates)))
    rmse = float(np.sqrt(np.mean((input_rates - output_rates) ** 2)))
    if length > 1 and np.std(input_rates) > 0 and np.std(output_rates) > 0:
        pearson_r = float(np.corrcoef(input_rates, output_rates)[0, 1])
    else:
        pearson_r = float("nan")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    positions = np.arange(length)

    axes[0].plot(positions, input_rates, marker="o", label="Input")
    axes[0].plot(positions, output_rates, marker="o", label="Recovered")
    axes[0].set_title(title)
    axes[0].set_xlabel("Position")
    axes[0].set_ylabel("Mutation rate")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].scatter(input_rates, output_rates, alpha=0.85)
    lo = float(min(input_rates.min(), output_rates.min()))
    hi = float(max(input_rates.max(), output_rates.max()))
    axes[1].plot([lo, hi], [lo, hi], "k--", linewidth=1)
    axes[1].set_title("Per-position recovery")
    axes[1].set_xlabel("Input mutation rate")
    axes[1].set_ylabel("Recovered mutation rate")
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        logging.info("Figure written to: %s", out)
    if show:
        plt.show()
    return fig, {"mae": mae, "rmse": rmse, "pearson_r": pearson_r}


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
