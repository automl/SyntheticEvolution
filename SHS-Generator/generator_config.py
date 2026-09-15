"""Generator defaults, validation, and construction of a PairMap from raw input.

MutationParameters holds global algorithm settings. PairMap is the ONLY source
of per-position mutation rates; paired/unpaired rates are used when building it.
This module performs no file I/O and does not depend on argparse or the pipeline.
"""
from dataclasses import dataclass, fields
import math
from numbers import Integral, Real
import numpy as np
from pair_map import PairMap

APPROACHES = ('watson_crick', 'covariance', 'watson_crick_cov', 'original', 'none')
DEFAULT_MUTATION_RATE_PAIRED = 0.2
DEFAULT_MUTATION_RATE_UNPAIRED = 0.2


def probability(value, name):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f'{name} must be a finite number between 0 and 1')
    return float(value)


def validate_seed(seed, name='seed', allow_none=True):
    if seed is None and allow_none:
        return
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError(f'{name} must be a nonnegative integer')


def validate_sequence(sequence):
    if not isinstance(sequence, str) or not sequence:
        raise ValueError('RNA sequence must be a nonempty string')
    sequence = sequence.upper()
    if set(sequence) - set('ACGU'):
        raise ValueError('RNA sequence must contain only A, C, G and U')
    return sequence


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


def build_pair_map(sequence, structure, mutation_rate_paired=DEFAULT_MUTATION_RATE_PAIRED,
                   mutation_rate_unpaired=DEFAULT_MUTATION_RATE_UNPAIRED,
                   mutation_rates=None):
    """Validate raw structure and build PairMap without silently dropping bad pairs.

    Accept dot-bracket notation, two-element pairs, or three-element weighted
    interactions. Positions are zero-based. Empty lists mean no interactions.
    A supplied per-position rate list overrides the two scalar rates.
    """
    sequence = validate_sequence(sequence)
    size = len(sequence)
    paired = probability(mutation_rate_paired, 'mutation_rate_paired')
    unpaired = probability(mutation_rate_unpaired, 'mutation_rate_unpaired')
    if isinstance(structure, str):
        if len(structure) != size:
            raise ValueError('Dot-bracket length must equal RNA length')
        stacks = {c: [] for c in '([{<'}
        closing = dict(zip(')]}>', '([{<'))
        raw = []
        for i, char in enumerate(structure):
            if char in stacks:
                stacks[char].append(i)
            elif char in closing:
                stack = stacks[closing[char]]
                if not stack:
                    raise ValueError(f'Unmatched closing bracket at position {i}')
                raw.append((stack.pop(), i))
            elif char != '.':
                raise ValueError(f'Invalid dot-bracket character: {char}')
        if any(stacks.values()):
            raise ValueError('Unmatched opening bracket in structure')
    elif isinstance(structure, dict):
        raw = list(structure.items())
    elif isinstance(structure, (list, tuple, np.ndarray)):
        raw = structure
    else:
        raise ValueError('Structure must be dot-bracket text or a list of pairs')

    interactions = {}
    for pair in raw:
        if not isinstance(pair, (list, tuple, np.ndarray)) or len(pair) not in (2, 3):
            raise ValueError('Each pair must contain two indices and optional strength')
        a, b = pair[:2]
        if any(isinstance(x, bool) or not isinstance(x, Integral) for x in (a, b)):
            raise ValueError('Pair indices must be integers')
        if not (0 <= a < size and 0 <= b < size) or a == b:
            raise ValueError(f'Invalid pair {(a, b)} for sequence length {size}')
        strength = probability(pair[2], 'interaction strength') if len(pair) == 3 else 1.0
        key = tuple(sorted((int(a), int(b))))
        if key in interactions and interactions[key] != strength:
            raise ValueError(f'Conflicting interaction strengths for pair {key}')
        interactions[key] = strength

    # A zero strength explicitly means no interaction; do not let PairMap clamp it.
    triples = [(a, b, strength) for (a, b), strength in interactions.items() if strength > 0]
    if mutation_rates is None:
        paired_positions = {i for a, b, _ in triples for i in (a, b)}
        rates = [paired if i in paired_positions else unpaired for i in range(size)]
    else:
        if not isinstance(mutation_rates, (list, tuple, np.ndarray)) or len(mutation_rates) != size:
            raise ValueError('mutation_rates must have exactly one rate per RNA position')
        rates = [probability(rate, f'mutation_rates[{i}]') for i, rate in enumerate(mutation_rates)]
    return PairMap.from_interactions(size, triples, rates)
