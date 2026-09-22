"""Run from repository root: python -m pytest hpo/tests/test_extended_inputs.py."""
import csv
import importlib.util
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'SHS-Generator'))
sys.path.insert(0, str(ROOT))
from generator_config import MutationParameters, split_parameters
from pair_map import build_pair_map
from shs_generator import MsaGenerator, load_request, parse_args
from hpo.pipeline.run import create_dataset, fingerprint
from hpo.pipeline.scoring import score_row


def dataset(tmp_path, **values):
    path = tmp_path / 'input.csv'
    row = dict(id='example', sequence='ACGU', **values)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader()
        writer.writerow({k: json.dumps(v) if isinstance(v, list) else v for k, v in row.items()})
    return create_dataset(path)


def request(tmp_path, row, parameters=None):
    path = tmp_path / 'request.json'
    path.write_text(json.dumps(dict(**row, task_name='test', parameters=parameters or {}, shs_seed=42, af3_seed=42)))
    return load_request(path)


def test_old_dataset_and_default_matrix(tmp_path):
    row = dataset(tmp_path, pairs=[[3, 0]])[0]
    assert row == dict(id='example', sequence='ACGU', pairs=[(0, 3)])
    _, _, pm, parameters = request(tmp_path, row)
    assert pm.interaction(0, 3) == 1
    assert pm.mutation_rate(0) == .2
    matrix = parameters.pair_mutation_probabilities()
    assert matrix['A'] == dict(A=0, C=0, G=0, U=1)
    assert matrix['G'] == dict(A=0, C=.75, G=0, U=.25)
    assert matrix['C'] == dict(A=0, C=0, G=1, U=0)
    assert matrix['U'] == dict(A=.75, C=0, G=.25, U=0)


def test_weighted_input_and_absolute_rates(tmp_path):
    row = dataset(tmp_path, interactions=[[3, 0, .2], [1, 2, 0]], mutation_rates=[0, .3, .5, 1])[0]
    assert row['pairs'] == [(0, 3)]
    assert row['interactions'] == [[0, 3, .2], [1, 2, 0]]
    _, _, pm, _ = request(tmp_path, row, dict(mutation_rate_paired=.8, mutation_rate_unpaired=.9))
    assert pm.interaction(0, 3) == .2
    assert pm.interaction(1, 2) == 0
    assert [pm.mutation_rate(i) for i in range(4)] == [0, .3, .5, 1]
    assert score_row(row, {(0, 3)})['loss'] == 0


def test_explicit_pairs_and_empty_interactions(tmp_path):
    row = dataset(tmp_path, pairs=[[0, 3]], interactions=[])[0]
    _, _, pm, _ = request(tmp_path, row)
    assert pm.unique_pairs == []
    assert score_row(row, set())['loss'] == 1


def test_blank_fields_fall_back(tmp_path):
    row = dataset(tmp_path, pairs=[[0, 3]], interactions='', mutation_rates='')[0]
    assert set(row) == {'id', 'sequence', 'pairs'}


@pytest.mark.parametrize('values', [
    dict(interactions=[[0, 4, 1]]), dict(interactions=[[0, 0, 1]]),
    dict(interactions=[[0, 3, -1]]), dict(interactions=[[0, 3, float('nan')]]),
    dict(interactions=[[0, 3, True]]), dict(interactions=[[False, 3, 1]]),
    dict(interactions=[[0, 3, 1.1]]), dict(interactions=[[0, 3, float('inf')]]), dict(interactions=[[0, 3, .1], [3, 0, .2]]),
    dict(pairs=[], mutation_rates=[]), dict(pairs=[], mutation_rates=[0, 0, 0, 2]),
    dict(pairs=[], mutation_rates=[False, 0, 0, 0]), dict(pairs='null'),
])
def test_invalid_dataset(tmp_path, values):
    with pytest.raises(ValueError):
        dataset(tmp_path, **values)


def test_matrix_normalizes_each_partner_row_independently():
    # Deliberately distinct rows catch transposition, row reuse and global normalization.
    raw = {
        'A': [1, 0, .5, .5],
        'C': [.25, .5, .75, 1],
        'G': [0, 1, 0, 0],
        'U': [.6, .2, 0, .2],
    }
    expected = {
        'A': dict(A=.5, C=0, G=.25, U=.25),
        'C': dict(A=.1, C=.2, G=.3, U=.4),
        'G': dict(A=0, C=1, G=0, U=0),
        'U': dict(A=.6, C=.2, G=0, U=.2),
    }
    weights = {f'pair_weight_{a}_{b}': value
               for a, row in raw.items() for b, value in zip('ACGU', row)}
    parameters, _, _ = split_parameters(weights)
    matrix = parameters.pair_mutation_probabilities()
    assert set(matrix) == set('ACGU')
    for base in 'ACGU':
        assert matrix[base] == pytest.approx(expected[base])


@pytest.mark.parametrize('base', 'ACGU')
def test_zero_weight_rows_are_rejected(base):
    with pytest.raises(ValueError, match='positive sum'):
        MutationParameters(**{f'pair_weight_{base}_{other}': 0 for other in 'ACGU'})


@pytest.mark.parametrize('value', [-.1, 1.1, True, float('nan'), float('inf'), '0.5'])
def test_invalid_weights_are_rejected(value):
    with pytest.raises(ValueError):
        MutationParameters(pair_weight_U_A=value)


@pytest.mark.parametrize('base', 'ACGU')
@pytest.mark.parametrize('candidate', 'ACGU')
def test_cli_accepts_each_weight(base, candidate):
    args = parse_args(['--rna-seq', 'ACGU', '--structure', '....',
                       f'--pair-weight-{base}-{candidate}', '.7'])
    assert getattr(args, f'pair_weight_{base}_{candidate}') == .7


@pytest.mark.parametrize('partner,expected', [
    ('A', dict(A=.50001, C=.00001, U=.25001)),
    ('C', dict(A=.10001, C=.20001, U=.40001)),
    ('G', dict(A=.00001, C=1.00001, U=.00001)),
    ('U', dict(A=.60001, C=.20001, U=.20001)),
])
def test_sampling_uses_partner_row_and_excludes_current_base(partner, expected):
    import numpy as np
    from unittest.mock import Mock
    raw = {'A': [1, 0, .5, .5], 'C': [.25, .5, .75, 1],
           'G': [0, 1, 0, 0], 'U': [.6, .2, 0, .2]}
    weights = {f'pair_weight_{a}_{b}': v for a, row in raw.items()
               for b, v in zip('ACGU', row)}
    generator = MsaGenerator(MutationParameters(**weights))
    # Fix only the random draws; the real mutation decision and weight calculation run.
    generator.rng = Mock()
    generator.rng.random.return_value = 0
    generator.rng.choice.return_value = 'A'
    generator.rng.choices.return_value = ['C']
    result = generator.mutate_wc('G', np.array(['U']), np.array([partner]),
                                 np.array([1.]), np.array([1.]), 1., False)
    generator.rng.choices.assert_called_once()
    candidates, sampling_weights = generator.rng.choices.call_args.args
    assert dict(zip(candidates, sampling_weights)) == pytest.approx(expected)
    assert 'G' not in candidates
    assert result == 'C'


def test_sampling_combines_partner_strengths():
    import numpy as np
    from unittest.mock import Mock
    generator = MsaGenerator(MutationParameters())
    generator.rng = Mock()
    generator.rng.random.return_value = 0
    generator.rng.choice.return_value = 'A'
    generator.rng.choices.return_value = ['U']
    generator.mutate_wc('C', np.array(['G', 'U']), np.array(['A', 'G']),
                        np.array([.2, .8]), np.array([1., 1.]), 1., False)
    candidates, weights = generator.rng.choices.call_args.args
    # A contributes .2 to U; G contributes .8*.25 to U. C is excluded.
    assert dict(zip(candidates, weights)) == pytest.approx(dict(A=.00001, G=.00001, U=.40001))


def test_zero_mutation_rate_does_not_sample_replacement():
    import numpy as np
    from unittest.mock import Mock
    generator = MsaGenerator(MutationParameters())
    generator.rng = Mock()
    generator.rng.random.return_value = .5
    assert generator.mutate_wc('G', np.array(['A']), np.array(['A']),
                               np.array([1.]), np.array([1.]), 0., False) == 'G'
    generator.rng.choices.assert_not_called()


def test_fingerprint_includes_new_inputs(tmp_path):
    first = dataset(tmp_path, interactions=[[0, 3, .2]], mutation_rates=[.1]*4)
    second = dataset(tmp_path, interactions=[[0, 3, .3]], mutation_rates=[.1]*4)
    third = dataset(tmp_path, interactions=[[0, 3, .2]], mutation_rates=[.2]*4)
    assert len({fingerprint({}, rows) for rows in (first, second, third)}) == 3


def test_request_execution_records_effective_parameters(tmp_path, monkeypatch):
    from types import SimpleNamespace
    import shs_generator
    row = dataset(tmp_path, interactions=[[0, 3, .2]], mutation_rates=[0, 1, 0, 1])[0]
    request(tmp_path, row, dict(N=3, pair_weight_A_C=1))
    monkeypatch.setitem(sys.modules, 'json_generator', SimpleNamespace(
        build_input_json=lambda sequence, msa, **kwargs: dict(sequence=sequence, msa=msa, **kwargs)))
    data = shs_generator.prepare_af3_json(parse_args(['--request-json', str(tmp_path/'request.json')]))
    assert len(data['msa']) == 3
    effective = json.loads((tmp_path/'shs_effective_parameters.json').read_text())
    assert effective['mutation_rates'] == [0, 1, 0, 1]
    assert effective['pair_mutation_probabilities']['A']['C'] == .5


def test_trial_forwards_and_records_complete_row(tmp_path, monkeypatch):
    from types import SimpleNamespace
    # External DSSR and AF3 execution are mocked; real request parsing and scoring run.
    import hpo.pipeline
    dssr_stub = SimpleNamespace(run_dssr=lambda *args: [[0, 3]])
    monkeypatch.setitem(sys.modules, 'hpo.pipeline.dssr', dssr_stub)
    # Also isolate the package attribute if another test already imported DSSR.
    monkeypatch.setattr(hpo.pipeline, 'dssr', dssr_stub, raising=False)
    spec = importlib.util.spec_from_file_location('extended_trial', ROOT/'hpo/pipeline/trial.py')
    trial = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(trial)
    monkeypatch.setattr(trial, 'run_af3_tasks', lambda *args: None)
    monkeypatch.setattr(trial, 'completed_af3_model', lambda *args: tmp_path/'model.cif')
    seen = []
    def generate(args, **kwargs):
        req_path = Path(args[args.index('--request-json')+1])
        req, seq, pm, params = load_request(req_path)
        seen.append((req, params.pair_mutation_probabilities()))
        assert pm.interaction(0, 3) == .3
        assert [pm.mutation_rate(i) for i in range(4)] == [.2]*4
        MsaGenerator(params, seed=req['shs_seed']).generate(seq, pm)
        Path(args[args.index('--output-json')+1]).write_text('{}')
    monkeypatch.setattr(trial.subprocess, 'run', generate)
    rows = dataset(tmp_path, interactions=[[0, 3, .3]], mutation_rates=[.2]*4)
    config = dict(fixed={'pair_weight_A_C': .5}, shs_seed=42, af3_seed=42,
                  shs_python='python', generator_timeout_seconds=10,
                  model_dir=str(tmp_path), dssr='dssr', dssr_timeout_seconds=10)
    loss = trial.run_pipeline(config, rows, tmp_path/'trial', 'identity', {'pair_weight_A_U': .5})
    assert loss == 0
    assert len(seen) == 1
    assert seen[0][0]['mutation_rates'] == [.2]*4
    assert seen[0][0]['interactions'] == [[0, 3, .3]]
    assert seen[0][1]['A']['C'] == .5
    evaluation = json.loads((tmp_path/'trial/rna_00000/base_pair_evaluation.json').read_text())
    assert evaluation['interactions'] == rows[0]['interactions']
    assert evaluation['mutation_rates'] == rows[0]['mutation_rates']


def test_scoring_uses_explicit_pairs_not_generator_interactions(tmp_path):
    row = dataset(tmp_path, pairs=[[0, 3], [1, 2]], interactions=[[0, 1, .2]])[0]
    _, _, pm, _ = request(tmp_path, row)
    assert pm.interaction(0, 1) == .2
    assert pm.interaction(0, 3) == 0
    result = score_row(row, {(0, 3), (0, 1)})
    assert result == dict(tp=1, fp=1, fn=1, f1=.5, loss=.5)


@pytest.mark.parametrize('overrides', [
    dict(interactions=[[0, 4, .2]]), dict(interactions=[[0, 3, -1]]),
    dict(mutation_rates=[]), dict(mutation_rates=[.1, .2, .3, 2]),
    dict(pairs=[[0, 5]], interactions=[]),
])
def test_request_validation_cannot_be_bypassed_by_skipping_csv(tmp_path, overrides):
    row = dict(id='example', sequence='ACGU', pairs=[[0, 3]])
    row.update(overrides)
    with pytest.raises(ValueError):
        request(tmp_path, row)


@pytest.mark.parametrize('approach', ['watson_crick', 'watson_crick_cov'])
def test_generate_dispatches_to_configured_preferences(approach):
    from unittest.mock import Mock
    parameters = MutationParameters(
        N=2, pair_mutation_approach=approach,
        pair_weight_A_U=0, pair_weight_A_C=1,
        stem_single_insertion_prob=0, stem_long_insertion_prob=0,
        stem_single_deletion_prob=0, stem_pair_deletion_prob=0,
        loop_single_insertion_prob=0, loop_long_insertion_prob=0,
        loop_single_deletion_prob=0, loop_long_deletion_prob=0,
    )
    generator = MsaGenerator(parameters)
    generator.rng = Mock()
    generator.rng.random.return_value = .25
    generator.rng.choice.side_effect = lambda candidates: candidates[0]
    generator.rng.choices.side_effect = lambda candidates, weights: [candidates[max(range(len(weights)), key=weights.__getitem__)]]
    pm = build_pair_map('AG', [[0, 1]], mutation_rates=[0, 1])
    assert generator.generate('AG', pm) == ['AG', 'AC']
    generator.rng.choices.assert_called_once()