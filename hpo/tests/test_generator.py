"""SHS generator internals, imported in-process: requests, pair map, weights, sampling.

The generator's command-line contract is in ``test_generator_cli.py``.
"""
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'SHS-Generator-RNA'))
from generator_config import MutationParameters, split_parameters  # noqa: E402
from pair_map import build_pair_map  # noqa: E402
import shs_generator  # noqa: E402
from shs_generator import MsaGenerator, load_request, parse_args  # noqa: E402

ROW = dict(id='example', sequence='ACGU', pairs=[[0, 3]])
# Deliberately distinct rows catch transposition, row reuse and global normalization.
RAW_WEIGHTS = {'A': [1, 0, .5, .5], 'C': [.25, .5, .75, 1], 'G': [0, 1, 0, 0], 'U': [.6, .2, 0, .2]}
WEIGHTS = {f'pair_weight_{a}_{b}': value for a, row in RAW_WEIGHTS.items() for b, value in zip('ACGU', row)}


def request(tmp_path, row, parameters=None):
    """Write a generator request as the trial does and load it."""
    path = tmp_path / 'request.json'
    path.write_text(json.dumps(dict(**row, task_name='test', parameters=parameters or {},
                                    shs_seed=42, af3_seed=42)))
    return load_request(path)


def mock_rng(generator, random=0, choice='A', choices=('C',)):
    """Fix only the random draws; the real mutation decisions and weights run."""
    generator.rng = Mock()
    generator.rng.random.return_value = random
    generator.rng.choice.return_value = choice
    generator.rng.choices.return_value = list(choices)
    return generator.rng


# ---- Requests and pair map --------------------------------------------------

def test_pairs_only_request_uses_full_interactions_and_default_rates(tmp_path):
    _, _, pairs, parameters = request(tmp_path, ROW)
    assert pairs.interaction(0, 3) == 1
    assert pairs.mutation_rate(0) == .2
    matrix = parameters.pair_mutation_probabilities()
    assert matrix['A'] == dict(A=0, C=0, G=0, U=1)
    assert matrix['G'] == dict(A=0, C=.75, G=0, U=.25)
    assert matrix['C'] == dict(A=0, C=0, G=1, U=0)
    assert matrix['U'] == dict(A=.75, C=0, G=.25, U=0)


def test_weighted_interactions_and_absolute_rates_override_parameters(tmp_path):
    row = dict(ROW, pairs=[[0, 3]], interactions=[[0, 3, .2], [1, 2, 0]], mutation_rates=[0, .3, .5, 1])
    _, _, pairs, _ = request(tmp_path, row, dict(mutation_rate_paired=.8, mutation_rate_unpaired=.9))
    assert (pairs.interaction(0, 3), pairs.interaction(1, 2)) == (.2, 0)
    assert [pairs.mutation_rate(i) for i in range(4)] == [0, .3, .5, 1]


def test_empty_interactions_mean_no_mutation_pairs(tmp_path):
    _, _, pairs, _ = request(tmp_path, dict(ROW, interactions=[]))
    assert pairs.unique_pairs == []


@pytest.mark.parametrize('overrides', [
    dict(interactions=[[0, 4, .2]]), dict(interactions=[[0, 3, -1]]),
    dict(mutation_rates=[]), dict(mutation_rates=[.1, .2, .3, 2]),
    dict(pairs=[[0, 5]], interactions=[]),
])
def test_request_is_validated_independently_of_the_csv_reader(tmp_path, overrides):
    with pytest.raises(ValueError):
        request(tmp_path, dict(ROW, **overrides))


def test_request_execution_records_effective_parameters(tmp_path, monkeypatch):
    request(tmp_path, dict(ROW, interactions=[[0, 3, .2]], mutation_rates=[0, 1, 0, 1]),
            dict(N=3, pair_weight_A_C=1))
    monkeypatch.setitem(sys.modules, 'json_generator', SimpleNamespace(
        build_input_json=lambda sequence, msa, **kwargs: dict(sequence=sequence, msa=msa, **kwargs)))
    data = shs_generator.prepare_af3_json(parse_args(['--request-json', str(tmp_path / 'request.json')]))
    assert len(data['msa']) == 3
    effective = json.loads((tmp_path / 'shs_effective_parameters.json').read_text())
    assert effective['mutation_rates'] == [0, 1, 0, 1]
    assert effective['pair_mutation_probabilities']['A']['C'] == .5


# ---- Pair weights -----------------------------------------------------------

def test_weight_matrix_normalizes_each_partner_row_independently():
    expected = {
        'A': dict(A=.5, C=0, G=.25, U=.25),
        'C': dict(A=.1, C=.2, G=.3, U=.4),
        'G': dict(A=0, C=1, G=0, U=0),
        'U': dict(A=.6, C=.2, G=0, U=.2),
    }
    parameters, _, _ = split_parameters(WEIGHTS)
    matrix = parameters.pair_mutation_probabilities()
    assert set(matrix) == set('ACGU')
    for base in 'ACGU':
        assert matrix[base] == pytest.approx(expected[base])


@pytest.mark.parametrize('base', 'ACGU')
def test_all_zero_weight_row_is_rejected(base):
    with pytest.raises(ValueError, match='positive sum'):
        MutationParameters(**{f'pair_weight_{base}_{other}': 0 for other in 'ACGU'})


@pytest.mark.parametrize('value', [-.1, 1.1, True, float('nan'), float('inf'), '0.5'])
def test_invalid_weight_is_rejected(value):
    with pytest.raises(ValueError):
        MutationParameters(pair_weight_U_A=value)


def test_cli_accepts_every_weight():
    for base in 'ACGU':
        for candidate in 'ACGU':
            args = parse_args(['--rna-seq', 'ACGU', '--structure', '....',
                               f'--pair-weight-{base}-{candidate}', '.7'])
            assert getattr(args, f'pair_weight_{base}_{candidate}') == .7


# ---- Sampling ---------------------------------------------------------------

@pytest.mark.parametrize('partner,expected', [
    ('A', dict(A=.50001, C=.00001, U=.25001)),
    ('C', dict(A=.10001, C=.20001, U=.40001)),
    ('G', dict(A=.00001, C=1.00001, U=.00001)),
    ('U', dict(A=.60001, C=.20001, U=.20001)),
])
def test_sampling_uses_the_partner_row_and_excludes_the_current_base(partner, expected):
    generator = MsaGenerator(MutationParameters(**WEIGHTS))
    rng = mock_rng(generator)
    result = generator.mutate_wc('G', np.array(['U']), np.array([partner]),
                                 np.array([1.]), np.array([1.]), 1., False)
    candidates, sampling_weights = rng.choices.call_args.args
    assert dict(zip(candidates, sampling_weights)) == pytest.approx(expected)
    assert 'G' not in candidates
    assert result == 'C'


def test_sampling_combines_partner_strengths():
    generator = MsaGenerator(MutationParameters())
    rng = mock_rng(generator, choices=('U',))
    generator.mutate_wc('C', np.array(['G', 'U']), np.array(['A', 'G']),
                        np.array([.2, .8]), np.array([1., 1.]), 1., False)
    candidates, weights = rng.choices.call_args.args
    # A contributes .2 to U; G contributes .8*.25 to U. C is excluded.
    assert dict(zip(candidates, weights)) == pytest.approx(dict(A=.00001, G=.00001, U=.40001))


def test_zero_mutation_rate_samples_no_replacement():
    generator = MsaGenerator(MutationParameters())
    rng = mock_rng(generator, random=.5)
    assert generator.mutate_wc('G', np.array(['A']), np.array(['A']),
                               np.array([1.]), np.array([1.]), 0., False) == 'G'
    rng.choices.assert_not_called()


@pytest.mark.parametrize('approach', ['watson_crick', 'watson_crick_cov'])
def test_generate_uses_the_configured_pair_weights(approach):
    no_indels = {f'{region}_{kind}_prob': 0 for region in ('stem', 'loop')
                 for kind in ('single_insertion', 'long_insertion', 'single_deletion')}
    parameters = MutationParameters(
        N=2, pair_mutation_approach=approach, pair_weight_A_U=0, pair_weight_A_C=1,
        stem_pair_deletion_prob=0, loop_long_deletion_prob=0, **no_indels,
    )
    generator = MsaGenerator(parameters)
    rng = mock_rng(generator, random=.25)
    rng.choice.side_effect = lambda candidates: candidates[0]
    rng.choices.side_effect = lambda candidates, weights: [candidates[int(np.argmax(weights))]]
    pairs = build_pair_map('AG', [[0, 1]], mutation_rates=[0, 1])
    assert generator.generate('AG', pairs) == ['AG', 'AC']
    rng.choices.assert_called_once()
