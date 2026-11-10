"""The trial comparison reads trial parameters saved per chain type."""
from analysis.modules import trial_comparison


def result(parameters):
    return dict(parameters=parameters, scores=[dict(id='x', predicted_pairs=[[0, 3]])])


def test_per_type_parameters_are_compared_as_type_prefixed_names():
    first, _ = trial_comparison._parse_trial(result({'rna': {'N': 8}, 'protein': {'N': 5}}))
    assert first == {'rna.N': 8, 'protein.N': 5}


def test_flat_parameters_are_left_unchanged():
    flat, _ = trial_comparison._parse_trial(result({'N': 8}))
    assert flat == {'N': 8}


def test_differing_generator_parameters_are_reported_per_type():
    trials = [trial_comparison._parse_trial(result({'rna': {'N': n}, 'protein': {'N': 5}})) for n in (8, 9)]
    report, identical = trial_comparison.compare(['config_1', 'config_2'], trials)
    text = '\n'.join(report)
    assert identical
    assert 'rna.N:' in text and 'protein.N = 5' in text
