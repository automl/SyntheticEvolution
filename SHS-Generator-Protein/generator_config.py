"""Generator defaults and validation of global mutation settings.

MutationParameters holds global algorithm settings. PairMap is the ONLY source
of per-position mutation rates; paired/unpaired rates are used when building it.
This module performs no file I/O and does not depend on argparse or the pipeline.
Run it as a script to print the HPO search entries of the residue-class table
(see search_space()).
"""
from dataclasses import dataclass, fields
import json
import math
from numbers import Integral, Real
import numpy as np
# Re-export existing helpers for compatibility with older callers.
from pair_map import (
    build_pair_map, probability, validate_sequence, AMINO_ACIDS,
    DEFAULT_MUTATION_RATE_PAIRED, DEFAULT_MUTATION_RATE_UNPAIRED,
)

APPROACHES = ('covariance', 'none', 'residue_class_pref', 'covariance_residue_class_pref')
BACKGROUNDS = ('uniform', 'query')

# Residue-class pair table. Classes in this order; every residue in exactly one.
RESIDUE_CLASSES = {'POS': 'KRH', 'NEG': 'DE', 'POL': 'STNQC', 'SMALL': 'GAVLIMP', 'BIG': 'FWY'}
AROMATIC = 'HFWY'  # o_ARO_ARO; includes His, which the POS class separates from F/W/Y
TABLE_BOUND = 3.0


def validate_seed(seed, name='seed', allow_none=True):
    if seed is None and allow_none:
        return
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError(f'{name} must be a nonnegative integer')


def is_table_parameter(name):
    return name.startswith(('r_', 'o_'))


@dataclass(frozen=True)
class MutationParameters:
    """Global generation settings; N includes the unchanged query row.

    Defaults match the existing generator. Per-position rates live in PairMap.
    The existing minimum long insertion/deletion length of 2 is retained.

    substitution_background sets the distribution every substitution (not
    insertion) is drawn from: 'uniform' over the other 19 residues (the original
    behaviour) or 'query', f_a = (count of a in the query + 1) / (L + 20).

    The r_*/o_* fields are residue-class pair energies used by the residue_class
    approaches; see residue_table(). Negative = favourable, bounded to [-3, 3],
    default 0. The approaches never change the meaning of the mutation rates.
    """
    N: int = 20
    pair_mutation_approach: str = 'covariance'
    substitution_background: str = 'uniform'
    stem_single_insertion_prob: float = 0.05
    stem_long_insertion_prob: float = 0.01
    stem_single_deletion_prob: float = 0.01
    stem_pair_deletion_prob: float = 0.01
    loop_single_insertion_prob: float = 0.1
    loop_single_deletion_prob: float = 0.1
    loop_long_insertion_prob: float = 0.02
    loop_long_deletion_prob: float = 0.02
    loop_max_insertion_fraction: float = 0.1
    loop_max_deletion_fraction: float = 0.1
    stem_max_insertion_fraction: float = 0.1

    # Class pair energies r_{C1}_{C2}, C1 not later than C2 in RESIDUE_CLASSES.
    r_POS_POS: float = 0.0
    r_POS_NEG: float = 0.0
    r_POS_POL: float = 0.0
    r_POS_SMALL: float = 0.0
    r_POS_BIG: float = 0.0
    r_NEG_NEG: float = 0.0
    r_NEG_POL: float = 0.0
    r_NEG_SMALL: float = 0.0
    r_NEG_BIG: float = 0.0
    r_POL_POL: float = 0.0
    r_POL_SMALL: float = 0.0
    r_POL_BIG: float = 0.0
    r_SMALL_SMALL: float = 0.0
    r_SMALL_BIG: float = 0.0
    r_BIG_BIG: float = 0.0
    # Overrides added on top of the class value.
    o_CYS_CYS: float = 0.0  # Cys-Cys cell only
    o_ARO_ARO: float = 0.0  # all pairs within AROMATIC
    o_GLY_GLY: float = 0.0  # Gly-Gly cell only

    def residue_table(self):
        """Symmetric 20x20 pair energies in AMINO_ACIDS order; negative = favourable.

        Block-constant by design: within two classes every pair has the same
        class value (Leu-Leu equals Leu-Ile), except for the three overrides.
        Gaps are not part of the table.
        """
        names = list(RESIDUE_CLASSES)
        residue_class = {aa: name for name, members in RESIDUE_CLASSES.items() for aa in members}
        table = np.empty((len(AMINO_ACIDS), len(AMINO_ACIDS)))
        for i, a in enumerate(AMINO_ACIDS):
            for j, b in enumerate(AMINO_ACIDS):
                first, second = sorted((residue_class[a], residue_class[b]), key=names.index)
                table[i, j] = (
                    getattr(self, f'r_{first}_{second}')
                    + self.o_CYS_CYS * (a == b == 'C')
                    + self.o_GLY_GLY * (a == b == 'G')
                    + self.o_ARO_ARO * (a in AROMATIC and b in AROMATIC)
                )
        return table

    def __post_init__(self):
        if isinstance(self.N, bool) or not isinstance(self.N, Integral) or self.N < 1:
            raise ValueError('N must be an integer >= 1')
        if self.pair_mutation_approach not in APPROACHES:
            raise ValueError(f'Unknown pair_mutation_approach: {self.pair_mutation_approach}')
        if self.substitution_background not in BACKGROUNDS:
            raise ValueError(f'Unknown substitution_background: {self.substitution_background}')
        for field in fields(self):
            value = getattr(self, field.name)
            if is_table_parameter(field.name):
                if (isinstance(value, bool) or not isinstance(value, Real)
                        or not math.isfinite(value) or abs(value) > TABLE_BOUND):
                    raise ValueError(f'{field.name} must be a finite number between '
                                     f'{-TABLE_BOUND} and {TABLE_BOUND}')
            elif field.name not in ('N', 'pair_mutation_approach', 'substitution_background'):
                probability(value, field.name)


def search_space():
    """HPO search entries for the new protein parameters (NePS YAML format).

    JSON is valid YAML, so the printed entries paste under
    generators.protein.search next to the existing parameters.
    """
    space = {
        'pair_mutation_approach': {'type': 'categorical', 'choices': list(APPROACHES)},
        'substitution_background': {'type': 'categorical', 'choices': list(BACKGROUNDS)},
    }
    for field in fields(MutationParameters):
        if is_table_parameter(field.name):
            space[field.name] = {'type': 'float', 'lower': -TABLE_BOUND, 'upper': TABLE_BOUND}
    return space


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


if __name__ == '__main__':
    print(json.dumps(search_space(), indent=2))
