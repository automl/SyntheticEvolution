"""Generator defaults and validation of global mutation settings.

MutationParameters holds global algorithm settings. PairMap is the ONLY source
of per-position mutation rates; paired/unpaired rates are used when building it.
This module performs no file I/O and does not depend on argparse or the pipeline.
"""
from dataclasses import dataclass, fields
from numbers import Integral
# Re-export existing helpers for compatibility with older callers.
from pair_map import (
    build_pair_map, probability, validate_sequence,
    DEFAULT_MUTATION_RATE_PAIRED, DEFAULT_MUTATION_RATE_UNPAIRED,
)

APPROACHES = ('watson_crick', 'covariance', 'watson_crick_cov', 'original', 'none')


def validate_seed(seed, name='seed', allow_none=True):
    if seed is None and allow_none:
        return
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError(f'{name} must be a nonnegative integer')


@dataclass(frozen=True)
class MutationParameters:
    """Global generation settings; N includes the unchanged query row.

    Defaults match the existing generator. Per-position rates live in PairMap.
    The existing minimum long insertion/deletion length of 2 is retained.
    """
    N: int = 20
    pair_mutation_approach: str = 'watson_crick_cov'
    stem_single_insertion_prob: float = 0.05
    stem_long_insertion_prob: float = 0.01
    stem_single_deletion_prob: float = 0.01
    stem_pair_deletion_prob: float = 0.01
    loop_single_insertion_prob: float = 0.1
    loop_single_deletion_prob: float = 0.1
    loop_long_insertion_prob: float = 0.02
    loop_long_deletion_prob: float = 0.02
    wobble_prob: float = 0.1
    max_insertion_fraction: float = 0.1
    max_deletion_fraction: float = 0.1


    # Partner base -> candidate base weights; normalized separately for each row.
    pair_weight_A_A: float = 0.0
    pair_weight_A_C: float = 0.0
    pair_weight_A_G: float = 0.0
    pair_weight_A_U: float = 1.0
    pair_weight_C_A: float = 0.0
    pair_weight_C_C: float = 0.0
    pair_weight_C_G: float = 1.0
    pair_weight_C_U: float = 0.0
    pair_weight_G_A: float = 0.0
    pair_weight_G_C: float = 0.75
    pair_weight_G_G: float = 0.0
    pair_weight_G_U: float = 0.25
    pair_weight_U_A: float = 0.75
    pair_weight_U_C: float = 0.0
    pair_weight_U_G: float = 0.25
    pair_weight_U_U: float = 0.0

    def pair_mutation_probabilities(self):
        """Resolve partner-conditioned preferences; these are not final frequencies."""
        matrix = {}
        for base in "ACGU":
            weights = {other: getattr(self, f"pair_weight_{base}_{other}") for other in "ACGU"}
            total = sum(weights.values())
            if total <= 0:
                raise ValueError(f"Pair weight row {base} must have a positive sum")
            matrix[base] = {other: weight / total for other, weight in weights.items()}
        return matrix

    def __post_init__(self):
        if isinstance(self.N, bool) or not isinstance(self.N, Integral) or self.N < 1:
            raise ValueError('N must be an integer >= 1')
        if self.pair_mutation_approach not in APPROACHES:
            raise ValueError(f'Unknown pair_mutation_approach: {self.pair_mutation_approach}')
        for field in fields(self):
            if field.name not in ('N', 'pair_mutation_approach'):
                probability(getattr(self, field.name), field.name)
        self.pair_mutation_probabilities()


def split_parameters(values):
    """Split a request's parameter dict into algorithm settings and structure rates.

    Missing values receive defaults; unknown keys fail. This is shared by the CLI
    and notebook, so request JSON need not expose internal dataclass organization.
    """
    if not isinstance(values, dict):
        raise ValueError('parameters must be a JSON object')
    values = dict(values)
    paired = probability(values.pop('mutation_rate_paired', DEFAULT_MUTATION_RATE_PAIRED), 'mutation_rate_paired')
    unpaired = probability(values.pop('mutation_rate_unpaired', DEFAULT_MUTATION_RATE_UNPAIRED), 'mutation_rate_unpaired')
    unknown = set(values) - {f.name for f in fields(MutationParameters)}
    if unknown:
        raise ValueError(f'Unknown generator parameters: {sorted(unknown)}')
    return MutationParameters(**values), paired, unpaired
