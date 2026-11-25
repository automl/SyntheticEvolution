"""Recover generator parameters from a synthetic MSA."""

import argparse
import logging
import json
from collections import Counter
from dataclasses import dataclass, field, replace
from functools import cached_property
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt

from pair_map import PairMap, build_pair_map

logger = logging.getLogger("analysis.shs_analyzer")

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
        logger.error("No unpaired MSA found for any rna in the json input")
        raise ValueError("Invalid Input, aborting.")
    if len(all_fasta) > 1:
        logger.warning("Input json contains multiple rna chains with unpaired msa, using the first one")
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
            logger.error("Recovered MSA is empty (need a query and >=1 sample)")
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
            logger.error("Sequences in MSA have different aligned lengths or the json was invalid")
            raise ValueError("Invalid Input, aborting.")
        aligned = np.array(aligned_rows)
        insertions = np.array(insertion_rows)
        if not np.all(np.isin(aligned, VALID_BASES)):
            bad = np.unique(aligned[~np.isin(aligned, VALID_BASES)])
            logger.error("Invalid characters in MSA: %s", bad)
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
    def paired(self) -> np.ndarray:
        """(L,) bool: position has a partner in the known structure."""
        return np.array([self.data.pairs.is_paired(i) for i in range(len(self.data.seq))])

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
    def col_insertion_rate(self) -> np.ndarray:
        """(L + 1) Per-gap fraction of rows with an insertion in the gap before the column"""
        return (self.data.insertions > 0).mean(axis=0)

    @cached_property
    def col_mut_rate(self) -> np.ndarray:
        """Per-column fraction of rows that differ from the seq"""
        return self.diff_mask.mean(axis=0, where=~self.deletion_mask)

    @cached_property
    def covariance(self) -> np.ndarray:
        """Covariance of the mutation mask of two columns over the rows where both
        are present (a deleted base is neither mutated nor kept). The diagonal is
        replaced by each column's mutation rate."""
        x = self.mutation_mask.astype(float)  # 0 where deleted
        present = (~self.deletion_mask).astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            n = present.T @ present
            mean = (x.T @ present) / n  # [i, j]: rate of i in the rows where i and j are present
            cov = (x.T @ x) / n - mean * mean.T
        np.fill_diagonal(cov, self.col_mut_rate)
        return cov


# ---------------------------------------------------------------------------
# Layer 3: estimators (each recovers one or more parameters) + driver
# ---------------------------------------------------------------------------

@dataclass
class Estimate:
    """A single recovered parameter: its value plus whatever supporting detail
    the estimator wants to expose (sample sizes, confidence, per-column data)."""
    name: str
    value: float
    detail: Dict[str, Any] = field(default_factory=dict)


# An estimator reads the shared Features and returns the Estimates of the
# parameters it recovers. Register new estimators by adding to ESTIMATORS.
Estimator = Callable[[Features], List[Estimate]]

WATSON_CRICK = {"A": "U", "C": "G", "G": "C", "U": "A"}
WOBBLE = {"G": "U", "U": "G"}


def _div(a, b) -> float:
    return float(a / b) if b else np.nan


def est_n(feat: Features) -> List[Estimate]:
    # N includes the query but the aligned data doesn't
    return [Estimate("n_sequences", float(len(feat.data.aligned) + 1))]


def est_mutation_rates(feat: Features) -> List[Estimate]:
    """Mutation rate of the kept unpaired and paired positions.

    'original' underestimates the paired rate, because a mutated pair can keep one
    of its bases (e.g. GC -> GU)."""
    out = []
    for name, cols in (("mutation_rate_unpaired", ~feat.paired), ("mutation_rate_paired", feat.paired)):
        kept = ~feat.deletion_mask[:, cols]
        rate = _div(feat.mutation_mask[:, cols].sum(), kept.sum())
        out.append(Estimate(name, rate, {"columns": int(cols.sum())}))
    return out


def est_insertions(feat: Features) -> List[Estimate]:
    """Insertion rates of stem gaps (both neighbours paired) and loop gaps (all
    others, including both ends).

    A gap gets a long insertion (>= 2 bases) with the long rate, otherwise a single
    one with the single rate. Long lengths are drawn up to max(int(L * fraction), 2),
    so the longest insertion / L is a lower bound of the max fraction."""
    paired, L = feat.paired, len(feat.paired)
    stem = np.zeros(L + 1, bool)
    stem[1:L] = paired[:-1] & paired[1:]
    out = []
    for kind, gaps in (("stem", stem), ("loop", ~stem)):
        ins = feat.data.insertions[:, gaps]
        n_long = int((ins >= 2).sum())
        longest = int(ins.max(initial=0))
        out += [
            Estimate(f"{kind}_long_insertion_prob", _div(n_long, ins.size), {"gaps": int(gaps.sum())}),
            Estimate(f"{kind}_single_insertion_prob", _div((ins == 1).sum(), ins.size - n_long)),
            Estimate(f"{kind}_max_insertion_fraction", longest / L if longest >= 2 else np.nan,
                     {"longest": longest}),
        ]
    return out


def est_loop_deletions(feat: Features) -> List[Estimate]:
    """Long and single deletion rates of unpaired positions.

    A long deletion run ends at the next paired position. So no run is active at
    position i if i is the first position or follows a paired or a kept position.
    For two unpaired positions i, i+1 after such a start,
    q = P(i deleted) = long + (1 - long) * single and
    r = P(i and i+1 deleted) = long + (1 - long) * single * q,
    which gives long = (r - q^2) / (1 - q).

    loop_max_deletion_fraction is not recovered: a new long deletion can start
    inside a running one, so runs can be longer than the maximum."""
    paired, dm = feat.paired, feat.deletion_mask
    no_run = np.ones_like(dm)
    no_run[:, 1:] = paired[:-1] | ~dm[:, :-1]
    cols = ~paired[:-1] & ~paired[1:]  # i and i+1 unpaired
    start = no_run[:, :-1][:, cols]
    first = start & dm[:, :-1][:, cols]
    both = first & dm[:, 1:][:, cols]
    q, r = _div(first.sum(), start.sum()), _div(both.sum(), start.sum())
    long = _div(r - q * q, 1 - q)
    return [Estimate("loop_long_deletion_prob", long, {"samples": int(start.sum())}),
            Estimate("loop_single_deletion_prob", _div(q - long, 1 - long))]


def est_stem_deletions(feat: Features) -> List[Estimate]:
    """Single and pair deletion rates from the basic pairs.

    A pair (a, b) with interaction I loses both bases with pair * I, otherwise each
    base independently with single. So P(only a) / P(neither) = single / (1 - single)
    and P(neither) = (1 - pair * I) * (1 - single)^2."""
    pairs, dm = feat.data.pairs, feat.deletion_mask
    neither = np.array([(~dm[:, a] & ~dm[:, b]).sum() for a, b in pairs.basic_pairs])
    one = sum((dm[:, a] ^ dm[:, b]).sum() for a, b in pairs.basic_pairs)
    x = _div(one / 2, neither.sum())
    single = x / (1 + x)
    pair_times_i = 1 - neither / len(dm) / (1 - single) ** 2
    pair = _div(pair_times_i.sum(), sum(pairs.interaction(a, b) for a, b in pairs.basic_pairs))
    return [Estimate("stem_single_deletion_prob", single, {"basic_pairs": len(neither)}),
            Estimate("stem_pair_deletion_prob", pair)]


def _pair_weights(partner: str, wc, wobble) -> Dict[str, Any]:
    """Weight of each replacement base given the partner base, as in the generator:
    the canonical bases share wc (G-U wobble gets wc * wobble), the mismatches share
    1 - wc, plus a floor of 1e-5."""
    canonical = {WATSON_CRICK[partner]: wc}
    if partner in WOBBLE:
        canonical = {WATSON_CRICK[partner]: wc * (1 - wobble), WOBBLE[partner]: wc * wobble}
    return {x: canonical.get(x, (1 - wc) / (4 - len(canonical))) + 1e-5 for x in "ACGU"}


def est_pair_preferences(feat: Features, steps: int = 101) -> List[Estimate]:
    """pair_wc_prob and pair_wobble_prob (watson_crick and watson_crick_cov only).

    Maximum likelihood on a grid. Uses every row where the second base b of a basic
    pair (a, b) mutated and a is kept: the new base of b is drawn by a's base,
    excluding b's original base. Assumes interaction 1 (weaker interactions only
    change the 1e-5 floor's share)."""
    seq, aligned = feat.data.seq, feat.data.aligned
    counts: Counter = Counter()
    for a, b in feat.data.pairs.basic_pairs:
        rows = feat.mutation_mask[:, b] & ~feat.deletion_mask[:, a]
        counts.update((p, seq[b], c) for p, c in zip(aligned[rows, a], aligned[rows, b]))
    events = sum(counts.values())
    if not events:
        return [Estimate("pair_wc_prob", np.nan, {"events": 0}), Estimate("pair_wobble_prob", np.nan)]
    wc, wobble = np.meshgrid(np.linspace(0, 1, steps), np.linspace(0, 1, steps), indexing="ij")
    log_likelihood = np.zeros_like(wc)
    for (partner, old, new), n in counts.items():
        w = _pair_weights(partner, wc, wobble)
        log_likelihood += n * np.log(w[new] / sum(w[x] for x in "ACGU" if x != old))
    best = np.unravel_index(np.argmax(log_likelihood), wc.shape)
    return [Estimate("pair_wc_prob", float(wc[best]), {"events": events}),
            Estimate("pair_wobble_prob", float(wobble[best]))]


def est_interaction_strength(feat: Features) -> List[Estimate]:
    """Mean interaction strength I of the basic pairs (covariance and
    watson_crick_cov only).

    The second base b of a pair (a, b) mutates with m_b + I m_b (1 - m_a) / (2 m_a)
    if a mutated and with m_b (1 - I / 2) otherwise, so
    cov(a, b) = I m_b (1 - m_a) / 2. Valid while the first probability is <= 1."""
    rate, cov = feat.col_mut_rate, feat.covariance
    values = [_div(2 * cov[a, b], rate[b] * (1 - rate[a])) for a, b in feat.data.pairs.basic_pairs]
    return [Estimate("interaction_strength", float(np.mean(values)) if values else np.nan,
                     {"basic_pairs": len(values)})]


# All estimators except est_n need the structure (MsaData.pairs).
ESTIMATORS: List[Estimator] = [
    est_n,
    est_mutation_rates,
    est_insertions,
    est_loop_deletions,
    est_stem_deletions,
    est_pair_preferences,
    est_interaction_strength,
]


def analyze_data(data: MsaData,
                 estimators: Optional[Sequence[Estimator]] = None) -> Dict[str, Estimate]:
    """Run the estimators on an already-parsed MsaData (e.g. straight from the
    notebook where the generator handed you its pairs)."""
    if estimators is None:
        estimators = ESTIMATORS
        if data.pairs is None:
            logger.warning("No structure given: only n_sequences is recovered")
            estimators = [est_n]
    feat = Features(data)
    return {e.name: e for est in estimators for e in est(feat)}


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
        logger.info("Figure written to: %s", out)
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
        logger.info("Figure written to: %s", out)
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
        logger.info("Figure written to: %s", out)
    if show:
        plt.show()
    return fig


def plot_indels(deletion_rate: np.ndarray, insertion_rate: np.ndarray,
                pairs: Optional[PairMap] = None,
                title: str = "Deletion rate per position and insertion rate per gap",
                show: bool = True, out: Optional[str] = None):
    """Deletion rate per position (top) and insertion rate per gap (bottom) on a shared x-axis.

    Gap k lies between positions k - 1 and k, so it is drawn at k - 0.5. With a
    structure, green marks paired positions and stem gaps (both neighbours paired),
    red marks unpaired positions and loop gaps."""
    L = len(deletion_rate)
    paired = np.array([bool(pairs) and pairs.is_paired(i) for i in range(L)])
    stem = np.zeros(L + 1, bool)
    stem[1:L] = paired[:-1] & paired[1:]

    fig, (ax_del, ax_ins) = plt.subplots(2, 1, sharex=True, figsize=(10, 5))
    for ax, x, y, is_stem, width, names in (
            (ax_del, np.arange(L), np.asarray(deletion_rate), paired, 0.8, ("paired", "unpaired")),
            (ax_ins, np.arange(L + 1) - 0.5, np.asarray(insertion_rate), stem, 0.5, ("stem gap", "loop gap"))):
        if pairs:
            ax.bar(x[is_stem], y[is_stem], width, color="green", label=names[0])
            ax.bar(x[~is_stem], y[~is_stem], width, color="red", label=names[1])
            ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1))
        else:
            ax.bar(x, y, width, color="gray")
        ax.grid(axis="y", alpha=0.3)
    ax_del.set_ylabel("deletion rate")
    ax_ins.set_ylabel("insertion rate")
    ax_ins.set_xlabel("position in the sequence (gaps between positions)")
    ax_ins.set_xticks(np.arange(0, L, 5))
    ax_ins.set_xticks(np.arange(L), minor=True)
    ax_del.set_title(title)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        logger.info("Figure written to: %s", out)
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
        logger.info("Figure written to: %s", out)
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
    p.add_argument("--structure", help="Dot-bracket or zero-based JSON pair list of the query; "
                                       "needed for every parameter except n_sequences.")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    data = MsaData.from_path(args.input)
    if args.structure:
        try:
            structure = json.loads(args.structure)
        except json.JSONDecodeError:
            structure = args.structure
        data = replace(data, pairs=build_pair_map(data.seq, structure))
    feat = Features(data)

    report = analyze_data(data)
    logger.info("Recovered parameters:")
    for name, est in report.items():
        logger.info("  %-32s = %-8s %s", name, round(est.value, 4), est.detail or "")

    title = args.title or "Recovered (covariation)"
    plot_covariance(data.seq, feat.covariance, title=title, show=args.show, out=args.out)


if __name__ == "__main__":
    main()