"""Protein generator: residue-class pair table, substitution background, regression.

The protein and RNA generators use the same module names, so each test imports
the protein modules freshly and restores whatever was imported before.
"""
from dataclasses import fields
import json
import math
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
PROTEIN_DIR = ROOT / 'SHS-Generator-Protein'
MODULES = ('generator_config', 'pair_map', 'shs_generator', 'json_generator')
CASES = json.loads((Path(__file__).with_name('fixtures') / 'protein_generator_regression.json').read_text())
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
CLASSES = {'POS': 'KRH', 'NEG': 'DE', 'POL': 'STNQC', 'SMALL': 'GAVLIMP', 'BIG': 'FWY'}
R_NAMES = [f'r_{a}_{b}' for i, a in enumerate(CLASSES) for b in list(CLASSES)[i:]]
TABLE_NAMES = R_NAMES + ['o_CYS_CYS', 'o_ARO_ARO', 'o_GLY_GLY']
EMPTY = (np.array([], str), np.array([], str), np.array([]), np.array([]))


@pytest.fixture
def protein(monkeypatch):
    saved = {name: sys.modules.pop(name) for name in MODULES if name in sys.modules}
    monkeypatch.syspath_prepend(str(PROTEIN_DIR))
    try:
        import generator_config
        import shs_generator
        yield SimpleNamespace(config=generator_config, generator=shs_generator,
                              Parameters=generator_config.MutationParameters)
    finally:
        for name in MODULES:
            sys.modules.pop(name, None)
        sys.modules.update(saved)


def index(residue):
    return AMINO_ACIDS.index(residue)


def block(value, first, second):
    expected = np.zeros((20, 20))
    for a in first:
        for b in second:
            expected[index(a), index(b)] = expected[index(b), index(a)] = value
    return expected


def ready_generator(protein, query='KKDE', **parameters):
    """A generator after generate() set its background; the rng only records draws."""
    generator = protein.generator.MsaGenerator(protein.Parameters(N=1, **parameters), seed=1)
    generator.generate(query, protein.config.build_pair_map(query, []))
    generator.rng = Mock()
    generator.rng.random.return_value = 0.0
    generator.rng.choices.return_value = ['A']
    return generator


def drawn_weights(generator):
    (population,), kwargs = generator.rng.choices.call_args
    return dict(zip(population, kwargs['weights']))


def table_draw(generator, nt, partners, strengths):
    generator.mutate_residue_class(nt, np.array(partners, str), np.array(partners, str),
                                   np.array(strengths, float), np.full(len(partners), .5), .5, False)
    return drawn_weights(generator)


# Table structure

def test_classes_partition_the_amino_acids_in_the_fixed_order(protein):
    assert protein.config.RESIDUE_CLASSES == CLASSES
    assert list(protein.config.RESIDUE_CLASSES) == list(CLASSES)
    assert sorted(''.join(CLASSES.values())) == sorted(AMINO_ACIDS) == sorted(protein.config.AMINO_ACIDS)


def test_the_eighteen_table_parameters_default_to_zero(protein):
    names = [f.name for f in fields(protein.Parameters) if protein.config.is_table_parameter(f.name)]
    assert names == TABLE_NAMES
    assert not protein.Parameters().residue_table().any()


@pytest.mark.parametrize('name', R_NAMES)
def test_each_class_entry_sets_exactly_its_block(protein, name):
    _, first, second = name.split('_')
    table = protein.Parameters(**{name: 1.5}).residue_table()
    assert np.array_equal(table, block(1.5, CLASSES[first], CLASSES[second]))


@pytest.mark.parametrize('name, residues', [('o_CYS_CYS', 'C'), ('o_GLY_GLY', 'G'), ('o_ARO_ARO', 'HFWY')])
def test_each_override_sets_exactly_its_cells(protein, name, residues):
    table = protein.Parameters(**{name: -2.0}).residue_table()
    assert np.array_equal(table, block(-2.0, residues, residues))


def test_overrides_add_to_the_class_value(protein):
    table = protein.Parameters(r_BIG_BIG=1.0, r_POS_BIG=0.25, o_ARO_ARO=-0.5,
                               r_POL_POL=0.5, o_CYS_CYS=-3.0).residue_table()
    assert table[index('F'), index('W')] == 0.5      # BIG-BIG + ARO
    assert table[index('H'), index('Y')] == -0.25    # POS-BIG + ARO
    assert table[index('H'), index('H')] == -0.5     # POS-POS (0) + ARO
    assert table[index('K'), index('F')] == 0.25     # POS-BIG only
    assert table[index('C'), index('C')] == -2.5     # POL-POL + CYS
    assert table[index('C'), index('S')] == 0.5      # POL-POL only


def test_table_is_symmetric(protein):
    values = {name: round(-3 + 6 * i / len(TABLE_NAMES), 3) for i, name in enumerate(TABLE_NAMES)}
    table = protein.Parameters(**values).residue_table()
    assert np.array_equal(table, table.T)


# Validation

@pytest.mark.parametrize('value', [-3, 3, 0, 2.5, -0.1])
def test_table_values_within_bounds_are_accepted(protein, value):
    protein.Parameters(r_POS_NEG=value, o_GLY_GLY=value)


@pytest.mark.parametrize('value', [3.01, -3.5, math.nan, math.inf, True, '1'])
def test_invalid_table_values_are_rejected(protein, value):
    with pytest.raises(ValueError, match='o_CYS_CYS'):
        protein.Parameters(o_CYS_CYS=value)


def test_other_settings_are_still_validated(protein):
    with pytest.raises(ValueError, match='substitution_background'):
        protein.Parameters(substitution_background='blosum')
    with pytest.raises(ValueError, match='loop_single_insertion_prob'):
        protein.Parameters(loop_single_insertion_prob=1.5)
    with pytest.raises(ValueError, match='r_FOO'):
        protein.config.split_parameters({'r_FOO': 0.0})


# Replacement weights

def test_uniform_background_without_partners_is_flat_and_excludes_the_current_residue(protein):
    weights = table_draw(ready_generator(protein), 'K', [], [])
    assert weights.pop('K') == 0
    assert weights == pytest.approx({a: 1 / 20 for a in AMINO_ACIDS if a != 'K'})


def test_query_background_is_the_smoothed_query_composition(protein):
    # Query KKDE: f_a = (count + 1) / (4 + 20).
    expected = {a: 1 / 24 for a in AMINO_ACIDS} | {'K': 0, 'D': 2 / 24, 'E': 2 / 24}
    generator = ready_generator(protein, substitution_background='query')
    assert table_draw(generator, 'K', [], []) == pytest.approx(expected)
    # Non-table substitutions (unpaired, none, covariance) use the same background.
    generator.mutate_random('K', 1.0)
    assert drawn_weights(generator) == pytest.approx(expected)


def test_favourable_entry_raises_weights_by_the_strength_weighted_sum(protein):
    generator = ready_generator(protein, r_POS_NEG=-3.0)
    weights = table_draw(generator, 'S', ['K'], [0.5])
    assert weights['D'] / weights['A'] == pytest.approx(math.exp(1.5))
    assert weights['E'] / weights['A'] == pytest.approx(math.exp(1.5))
    assert weights['R'] == pytest.approx(weights['A'])
    # Contributions of several partners are summed: K (1.0) + R (0.5).
    weights = table_draw(generator, 'S', ['K', 'R'], [1.0, 0.5])
    assert weights['D'] / weights['A'] == pytest.approx(math.exp(4.5))


@pytest.mark.parametrize('partners, strengths', [(['-'], [1.0]), (['K'], [0.0]), (['-', 'K'], [1.0, 0.0])])
def test_deleted_or_zero_strength_partners_change_nothing(protein, partners, strengths):
    generator = ready_generator(protein, r_POS_NEG=-3.0, r_NEG_SMALL=2.0)
    assert table_draw(generator, 'S', partners, strengths) == table_draw(generator, 'S', [], [])


def test_decision_comes_from_the_rate_or_covariance_not_the_table(protein):
    generator = ready_generator(protein, o_CYS_CYS=-3.0)
    generator.rng.random.return_value = 0.9
    assert generator.mutate_residue_class('C', *EMPTY, 0.5, False) == 'C'
    generator.rng.choices.assert_not_called()
    generator.mutate_cov = Mock(return_value='C')
    assert generator.mutate_residue_class('C', *EMPTY, 0.5, True) == 'C'
    generator.mutate_cov.return_value = 'A'
    assert generator.mutate_residue_class('C', *EMPTY, 0.5, True) == 'A'
    assert drawn_weights(generator)['C'] == 0


# Whole sequences

SEQUENCE = 'MKCAYICKQRQGSFVKSHFDREWEERLGLIEYQ'
PAIRS = [[1, 19], [2, 7], [4, 30], [12, 27], [6, 25], [9, 21]]
NO_INDELS = {name: 0.0 for name in (
    'stem_single_insertion_prob', 'stem_long_insertion_prob', 'stem_single_deletion_prob',
    'stem_pair_deletion_prob', 'loop_single_insertion_prob', 'loop_single_deletion_prob',
    'loop_long_insertion_prob', 'loop_long_deletion_prob')}
ALL_COMBINATIONS = [(a, b) for a in ('none', 'covariance', 'residue_class_pref', 'covariance_residue_class_pref')
                    for b in ('uniform', 'query')]


def generate(protein, approach, background, rate, **table):
    parameters = protein.Parameters(N=40, pair_mutation_approach=approach, substitution_background=background,
                                    **NO_INDELS, **table)
    pair_map = protein.config.build_pair_map(SEQUENCE, PAIRS, rate, rate)
    return protein.generator.MsaGenerator(parameters, seed=3).generate(SEQUENCE, pair_map)


@pytest.mark.parametrize('approach, background', ALL_COMBINATIONS)
def test_zero_mutation_rate_changes_nothing(protein, approach, background):
    assert set(generate(protein, approach, background, 0.0, o_CYS_CYS=-3.0, r_POS_NEG=-3.0)) == {SEQUENCE}


@pytest.mark.parametrize('background', ['uniform', 'query'])
def test_full_rate_changes_every_position_even_against_a_strong_table(protein, background):
    msa = generate(protein, 'residue_class_pref', background, 1.0, o_CYS_CYS=-3.0, r_POL_POL=-3.0)
    assert all(a != b for row in msa[1:] for a, b in zip(row, SEQUENCE))


@pytest.mark.parametrize('case', CASES, ids=[case['id'] for case in CASES])
def test_seeded_output_matches_regression_fixture(protein, case):
    pair_map = protein.config.build_pair_map(case['sequence'], case['pairs'],
                                             case['mutation_rate_paired'], case['mutation_rate_unpaired'])
    parameters = protein.Parameters(**case['parameters'])
    assert protein.generator.MsaGenerator(parameters, seed=case['seed']).generate(case['sequence'], pair_map) == case['msa']


# Request and HPO interface

def test_request_records_every_new_parameter(protein, tmp_path):
    table = {name: 0.5 for name in TABLE_NAMES}
    request = tmp_path / 'request.json'
    request.write_text(json.dumps(dict(
        task_name='p', sequence=SEQUENCE, pairs=PAIRS, shs_seed=1, af3_seed=1,
        parameters=dict(N=5, pair_mutation_approach='covariance_residue_class_pref',
                        substitution_background='query', **table))))
    output = tmp_path / 'out.json'
    assert protein.generator.main(['--request-json', str(request), '--output-json', str(output)]) == 0
    effective = json.loads((tmp_path / 'shs_effective_parameters.json').read_text())['parameters']
    assert effective['pair_mutation_approach'] == 'covariance_residue_class_pref'
    assert effective['substitution_background'] == 'query'
    assert {name: effective[name] for name in TABLE_NAMES} == table
    assert json.loads(output.read_text())['sequences'][0]['protein']['sequence'] == SEQUENCE


def test_spec_dump_lists_every_new_search_entry():
    result = subprocess.run([sys.executable, str(PROTEIN_DIR / 'generator_config.py')],
                            capture_output=True, text=True, check=True)
    space = json.loads(result.stdout)
    assert list(space) == ['pair_mutation_approach', 'substitution_background'] + TABLE_NAMES
    assert space['pair_mutation_approach'] == {
        'type': 'categorical',
        'choices': ['covariance', 'none', 'residue_class_pref', 'covariance_residue_class_pref']}
    assert space['substitution_background'] == {'type': 'categorical', 'choices': ['uniform', 'query']}
    assert all(space[name] == {'type': 'float', 'lower': -3.0, 'upper': 3.0} for name in TABLE_NAMES)
