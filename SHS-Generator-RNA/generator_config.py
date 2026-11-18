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
WATSON_CRICK = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
WOBBLE = {'G': 'U', 'U': 'G'}


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
    loop_max_insertion_fraction: float = 0.1
    loop_max_deletion_fraction: float = 0.1
    stem_max_insertion_fraction: float = 0.1

    # Replacement base for a mutated pair member, given its (mutated) partner
    # (watson_crick and watson_crick_cov only). With pair_wc_prob the new base
    # pairs canonically with the partner, otherwise it is a mismatch (uniform over
    # the non-canonical bases). For G and U partners, a canonical replacement is
    # the G-U wobble with pair_wobble_prob and Watson-Crick otherwise. The defaults
    # give the former default (pure WC for A/C partners, 3:1 WC:wobble for G/U).
    pair_wc_prob: float = 1.0
    pair_wobble_prob: float = 0.25

    def pair_mutation_probabilities(self):
        """Partner base -> probability of each replacement base.

        These are the sampling preferences, not final MSA frequencies: the current
        base is excluded when sampling, and multiplets mix several partner rows.
        """
        wc, wobble = self.pair_wc_prob, self.pair_wobble_prob
        matrix = {}
        for partner in 'ACGU':
            canonical = {WATSON_CRICK[partner]: wc}
            if partner in WOBBLE:
                canonical = {WATSON_CRICK[partner]: wc * (1 - wobble), WOBBLE[partner]: wc * wobble}
            mismatches = [base for base in 'ACGU' if base not in canonical]
            matrix[partner] = {base: canonical.get(base, (1 - wc) / len(mismatches)) for base in 'ACGU'}
        return matrix

    def __post_init__(self):
        if isinstance(self.N, bool) or not isinstance(self.N, Integral) or self.N < 1:
            raise ValueError('N must be an integer >= 1')
        if self.pair_mutation_approach not in APPROACHES:
            raise ValueError(f'Unknown pair_mutation_approach: {self.pair_mutation_approach}')
        for field in fields(self):
            if field.name not in ('N', 'pair_mutation_approach'):
                probability(getattr(self, field.name), field.name)


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
