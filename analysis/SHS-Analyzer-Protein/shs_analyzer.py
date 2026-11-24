"""Recover generator parameters from a synthetic MSA."""

import argparse
import logging
import json
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import importlib
import sys
import warnings

SHS_DIR = Path(__file__).resolve().parents[2] / "SHS-Generator-Protein"
if str(SHS_DIR) not in sys.path:
    sys.path.insert(0, str(SHS_DIR))

from pair_map import PairMap, AMINO_ACIDS

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Real sequences can hold ambiguous or rare residues (B, J, O, U, X, Z); the
# features only compare residues with the query and '-', so they are kept.
VALID_BASES = list(AMINO_ACIDS + 'BJOUXZ') + ['-']


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
    # Native AF3 JSON stores MSAs per chain type; protein chains use "protein",
    # so accept both (protein first, since this is the protein analyzer).
    all_fasta = []
    for chain_type in ("protein", "rna"):
        all_fasta += [msa for c in data.get("sequences", [])
                      if chain_type in c and (msa := c[chain_type].get("unpairedMsa"))]
    if len(all_fasta) == 0:
        logging.error("No unpaired MSA found for any protein or rna chain in the json input")
        raise ValueError("Invalid Input, aborting.")
    if len(all_fasta) > 1:
        logging.warning("Input json contains multiple chains with unpaired msa, using the first one")
    return all_fasta[0]


@dataclass(frozen=True, eq=False)
class MsaData:
    """Immutable parsed Protein."""
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
        """Per-column fraction of rows that differ from the seq.
        Columns deleted in every row have no data and come out as NaN."""
        with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
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
    """Mutation rate at unpaired positions -> recovers --mutation-rate-unpaired.

    Without a known structure, all positions are treated as unpaired, which
    blends the paired and unpaired rates and is reported via detail."""
    pairs = feat.data.pairs
    if pairs is not None:
        unpaired_cols = [i for i in range(len(feat.data.seq)) if not pairs.is_paired(i)]
    else:
        unpaired_cols = list(range(len(feat.data.seq)))
    mutations = feat.mutation_mask[:, unpaired_cols]
    no_mutation = ~feat.deletion_mask[:, unpaired_cols]
    denom = int(no_mutation.sum())
    rate = float(mutations.sum() / denom) if denom else float("nan")
    detail: Dict[str, Any] = {"unpaired_cols": len(unpaired_cols), "structure_known": pairs is not None}
    return Estimate("mutation_rate_unpaired", rate, detail)


def est_mutation_rate_paired(feat: Features, results: Dict[str, Estimate]) -> Estimate:
    """Mutation rate at paired (contact) positions -> recovers --mutation-rate-paired.

    Without a known structure, no paired columns exist and the estimate is NaN.
    Note that mutation rates are underestimated for the covariance approach:
    a partner can mutate back to the original residue."""
    pairs = feat.data.pairs
    if pairs is not None:
        paired_cols = [i for i in range(len(feat.data.seq)) if pairs.is_paired(i)]
    else:
        paired_cols = []
    mutations = feat.mutation_mask[:, paired_cols]
    not_deleted = ~feat.deletion_mask[:, paired_cols]
    denom = int(not_deleted.sum())
    rate = float(mutations.sum() / denom) if denom else float("nan")
    detail: Dict[str, Any] = {"paired_cols": len(paired_cols), "structure_known": pairs is not None}
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


def _finish(fig, out: Optional[str], show: bool) -> None:
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        logging.info("Figure written to: %s", out)
    if show:
        plt.show()


def plot_recovered_pairs(cov: np.ndarray, pairs: PairMap,
                         title: str = "Input vs recovered pairs",
                         show: bool = True, out: Optional[str] = None):
    """Compare the input pairs with the pairs recovered from the MSA.

    The k position pairs with the highest covariance (Features.covariance) are
    the recovered pairs, k = number of input pairs. Left: the covariance (lower
    triangle) and which pairs were recovered (upper triangle). Right: covariance
    per input group, a connected set of input contacts (PairMap.partners);
    pairs that are not input pairs are "no group". Returns the figure and the
    fraction of the k recovered pairs that are input pairs.
    """
    length = cov.shape[0]
    i, j = np.triu_indices(length, k=1)
    values = cov[i, j]
    is_input = np.array([pairs.is_pair(a, b) for a, b in zip(i, j)])
    k = int(is_input.sum())
    if k == 0:
        raise ValueError("Need at least one input pair to compare with.")
    # Each group is named by its smallest position; "no group" comes last.
    component = np.array([min(a, *pairs.partners(a)) if p else -1 for a, p in zip(i, is_input)])
    anchors = np.unique(component[is_input])
    group = np.where(is_input, np.searchsorted(anchors, component), len(anchors))
    x = group + np.random.default_rng(0).uniform(-0.35, 0.35, len(group))
    recovered = np.zeros(len(values), dtype=bool)
    recovered[np.argsort(-values)[:k]] = True
    precision = float((is_input & recovered).sum() / k)
    groups = [
        (is_input & recovered, dict(c="tab:green", marker="s"), "input pair, recovered"),
        (is_input & ~recovered, dict(facecolors="none", edgecolors="tab:orange", marker="s"), "input pair, missed"),
        (~is_input & recovered, dict(c="tab:red", marker="x"), "recovered, not an input pair"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    lower = np.where(np.tri(length, k=-1, dtype=bool), cov, np.nan)
    im = axes[0].imshow(lower, cmap="viridis", vmin=0)
    fig.colorbar(im, ax=axes[0], label="covariance")
    for mask, style, label in groups:
        axes[0].scatter(j[mask], i[mask], s=1500 / length, label=label, **style)
    axes[0].set_title(f"{title}\nlower: covariance · upper: top-{k} recovered vs input")
    axes[0].set_xlabel("position")
    axes[0].set_ylabel("position")

    rest = ~is_input & ~recovered
    axes[1].scatter(x[rest], values[rest], c="lightgray", s=8, label="other position pairs")
    for mask, style, label in groups:
        axes[1].scatter(x[mask], values[mask], s=30, label=label, **style)
    axes[1].axhline(np.sort(values)[-k], color="k", linestyle="--", linewidth=1, label=f"top-{k} cutoff")
    axes[1].set_title(f"{precision:.0%} of the recovered pairs are input pairs")
    axes[1].set_xticks(range(len(anchors) + 1), [f"group {n + 1}" for n in range(len(anchors))] + ["no group"])
    axes[1].set_xlabel("input group (connected input contacts)")
    axes[1].set_ylabel("covariance")
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    axes[1].grid(True, axis="y", alpha=0.3)
    _finish(fig, out, show)
    return fig, {"k": k, "precision": precision}


def _codes(residues: np.ndarray) -> np.ndarray:
    """Index in AMINO_ACIDS per residue; -1 for gaps and non-canonical residues."""
    codes = np.full(residues.shape, -1)
    for k, a in enumerate(AMINO_ACIDS):
        codes[residues == a] = k
    return codes


def _background(seq: str, kind: str) -> np.ndarray:
    """The generator's substitution_background: 'uniform', or 'query' composition +1 per residue."""
    if kind not in ("uniform", "query"):
        raise ValueError(f"Unknown substitution background: {kind}")
    counts = np.ones(len(AMINO_ACIDS))
    if kind == "query":
        counts += [seq.count(a) for a in AMINO_ACIDS]
    return counts / counts.sum()


def _draws(background: np.ndarray, original: int) -> np.ndarray:
    """A substitution never keeps the original residue: background without it, renormalized."""
    q = background.copy()
    q[original] = 0
    return q / q.sum()


def plot_residue_class_preferences(data: MsaData, table: np.ndarray, classes: Dict[str, str],
                                   background: str = "uniform", min_expected: float = 5,
                                   title: str = "Residue pair preferences at the input contacts",
                                   show: bool = True, out: Optional[str] = None):
    """Input residue pair table (MutationParameters.residue_table()) against the MSA.

    For every input contact i < j and every row where j was substituted (the
    generator draws j after i), counts the residue pair (a_i, b_j), and the count
    expected without preferences: the background without j's query residue.
    The recovered energy is -ln(observed / expected) (+0.5 each), negative = favourable as in
    the table, and scaled by the contacts' interaction strengths. A constant
    added to a table row changes no draw, so both are compared row-centered.
    Cells expected fewer than min_expected times are left out. Residues are
    ordered by class. Returns the figure and the Pearson r of the comparison.
    """
    order = [AMINO_ACIDS.index(a) for members in classes.values() for a in members]
    bg = _background(data.seq, background)
    codes = _codes(data.aligned)
    query = _codes(np.array(list(data.seq)))
    observed = np.zeros((len(AMINO_ACIDS), len(AMINO_ACIDS)))
    expected = np.zeros_like(observed)
    for i, j in data.pairs.unique_pairs:
        rows = (codes[:, i] >= 0) & (codes[:, j] >= 0) & (codes[:, j] != query[j])
        if query[j] < 0 or not rows.any():
            continue
        np.add.at(observed, (codes[rows, i], codes[rows, j]), 1)
        expected += np.outer(np.bincount(codes[rows, i], minlength=len(AMINO_ACIDS)), _draws(bg, query[j]))
    recovered = -np.log((observed + 0.5) / (expected + 0.5))
    recovered[expected < min_expected] = np.nan
    with warnings.catch_warnings():  # rows without any kept cell stay NaN
        warnings.simplefilter("ignore", RuntimeWarning)
        recovered -= np.nanmean(recovered, axis=1, keepdims=True)
        centered = np.where(np.isnan(recovered), np.nan, table)
        centered -= np.nanmean(centered, axis=1, keepdims=True)
    recovered, centered = recovered[np.ix_(order, order)], centered[np.ix_(order, order)]
    ok = np.isfinite(recovered) & np.isfinite(centered)
    pearson_r = float(np.corrcoef(centered[ok], recovered[ok])[0, 1]) if ok.sum() > 2 and centered[ok].std() > 0 else float("nan")

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    labels = [AMINO_ACIDS[k] for k in order]
    names = list(classes)
    starts = np.cumsum([0] + [len(m) for m in classes.values()])
    table = table[np.ix_(order, order)]
    for ax, matrix, name in ((axes[0], table, "input table"),
                             (axes[1], recovered, "recovered from the MSA (row-centered)")):
        limit = np.nanmax(np.abs(matrix)) or 1
        im = ax.imshow(matrix, cmap="RdBu", vmin=-limit, vmax=limit)
        fig.colorbar(im, ax=ax, label="energy (negative = favourable)")
        ax.set_xticks(range(len(labels)), labels)
        ax.set_yticks(range(len(labels)), labels)
        for b in starts[1:-1] - 0.5:
            ax.axhline(b, color="k", linewidth=0.5)
            ax.axvline(b, color="k", linewidth=0.5)
        ax.set_title(name)
    # Name each class block after its r_<class>_<class> parameter.
    limit = np.abs(table).max() or 1
    for r, row_class in enumerate(names):
        for c, col_class in enumerate(names):
            block = table[starts[r]:starts[r + 1], starts[c]:starts[c + 1]]
            first, second = sorted((row_class, col_class), key=names.index)
            axes[0].text((starts[c] + starts[c + 1] - 1) / 2, (starts[r] + starts[r + 1] - 1) / 2,
                         f"{first}-\n{second}", ha="center", va="center", fontsize=7,
                         color="white" if abs(np.median(block)) > 0.5 * limit else "black")
    axes[1].set_ylabel("residue of the earlier contact member")
    axes[1].set_xlabel("residue drawn for the later member")
    axes[2].scatter(centered[ok], recovered[ok], s=10, alpha=0.6)
    axes[2].axhline(0, color="gray", linewidth=0.5)
    axes[2].axvline(0, color="gray", linewidth=0.5)
    axes[2].set_title(f"{title}\nPearson r = {pearson_r:.2f} over {int(ok.sum())} residue pairs")
    axes[2].set_xlabel("input energy (row-centered)")
    axes[2].set_ylabel("recovered energy (row-centered)")
    axes[2].grid(True, alpha=0.3)
    _finish(fig, out, show)
    return fig, {"pearson_r": pearson_r, "residue_pairs": int(ok.sum()), "draws": int(observed.sum())}


def plot_substitution_background(data: MsaData,
                                 title: str = "Composition of substituted residues",
                                 show: bool = True, out: Optional[str] = None):
    """Substituted residues against both substitution_background settings.

    Uses the positions without a lower-index contact partner, where every
    approach draws from the background alone. A substitution never keeps the
    original residue, so the expectation leaves it out per position. Returns the
    figure and the mean log-likelihood per substitution under each setting; the
    higher one fits the MSA better.
    """
    codes = _codes(data.aligned)
    query = _codes(np.array(list(data.seq)))
    backgrounds = {kind: _background(data.seq, kind) for kind in ("uniform", "query")}
    observed = np.zeros(len(AMINO_ACIDS))
    expected = {kind: np.zeros(len(AMINO_ACIDS)) for kind in backgrounds}
    loglik = dict.fromkeys(backgrounds, 0.0)
    for i in range(len(data.seq)):
        if query[i] < 0 or (data.pairs is not None and any(p < i for p in data.pairs.direct_partners(i))):
            continue
        drawn = codes[:, i][(codes[:, i] >= 0) & (codes[:, i] != query[i])]
        counts = np.bincount(drawn, minlength=len(AMINO_ACIDS))
        observed += counts
        for kind, bg in backgrounds.items():
            q = _draws(bg, query[i])
            expected[kind] += len(drawn) * q
            loglik[kind] += counts[q > 0] @ np.log(q[q > 0])
    n = int(observed.sum())
    if n == 0:
        raise ValueError("No substitutions at positions drawn from the background alone.")
    loglik = {kind: value / n for kind, value in loglik.items()}

    fig, ax = plt.subplots(figsize=(14, 4))
    x = np.arange(len(AMINO_ACIDS))
    ax.bar(x, observed / n, color="lightgray", label=f"observed ({n} substitutions)")
    for kind, marker in (("uniform", "o"), ("query", "x")):
        ax.scatter(x, expected[kind] / n, marker=marker, zorder=3,
                   label=f"expected for '{kind}' (log-likelihood {loglik[kind]:.4f})")
    ax.set_xticks(x, list(AMINO_ACIDS))
    ax.set_xlabel("substituted residue")
    ax.set_ylabel("fraction of substitutions")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    _finish(fig, out, show)
    # With a query that holds every residue equally often, both settings are the same.
    best = "either" if np.isclose(*loglik.values()) else max(loglik, key=loglik.get)
    return fig, {"substitutions": n, "loglik": loglik, "best": best}


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
