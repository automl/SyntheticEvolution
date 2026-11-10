"""analyze_checks on a tiny check export with hand-computed values."""
import csv
import json
import math

import pytest

from analysis.analyze_checks import analyze

TRAIN = dict(id='a', sequence='ACGU', pairs=[[0, 3]])
HELDOUT = dict(id='b', sequence='GGCC', pairs=[[0, 3], [1, 2]])


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.fixture
def exports(tmp_path):
    hpo = tmp_path / 'results'
    write(hpo / 'run_metadata.json', dict(dataset=[TRAIN]))
    (hpo / 'summary').mkdir()
    (hpo / 'summary' / 'full.csv').write_text(
        'id,objective_to_minimize,state\n7,0.2,success\n8,1.0,crashed\n')

    checks = tmp_path / 'checks' / 'results'
    write(checks / 'run_metadata.json', dict(dataset=[TRAIN, HELDOUT]))
    trials = [
        dict(name='c7_hpo', check='reproduction', base=7),
        dict(name='c7_s1', check='seed', base=7),
        dict(name='c7_s2', check='seed', base=7),
        dict(name='nomsa_s1', check='baseline'),
        dict(name='nomsa_s2', check='baseline'),  # failed: no result
    ]
    write(checks / 'manifest.json', dict(trials=trials))
    # (loss of a, prediction of a, loss of b, prediction of b) per trial
    values = {
        'c7_hpo': (0.25, [[0, 3]], 0.0, [[0, 3], [1, 2]]),
        'c7_s1': (0.0, [[0, 3]], 0.5, [[0, 3]]),
        'c7_s2': (0.5, [], 0.25, [[1, 2]]),
        'nomsa_s1': (1.0, [], 1.0, []),
    }
    for name, (loss_a, pred_a, loss_b, pred_b) in values.items():
        write(checks / 'trials' / name / 'trial_result.json', dict(
            loss=(loss_a + loss_b) / 2, shs_seed=1, af3_seed=1,
            scores=[dict(id='a', loss=loss_a, predicted_pairs=pred_a),
                    dict(id='b', loss=loss_b, predicted_pairs=pred_b)]))
    return checks, hpo


def read(path):
    with path.open(newline='') as stream:
        return {row.get('group', '') + row.get('rna_id', ''): row for row in csv.DictReader(stream)}


def test_summary_splits_train_and_heldout_and_compares_with_hpo(exports, tmp_path, caplog):
    checks, hpo = exports
    output = tmp_path / 'checks' / 'analysis'
    output.mkdir()
    summary = {row['group']: row for row in analyze(checks, hpo, output)}
    c7 = summary['c7']
    assert (c7['hpo_loss'], c7['reproduction_train'], c7['seeds']) == (0.2, 0.25, 2)
    assert c7['reproduction_diff'] == pytest.approx(0.05)
    assert (c7['train_mean'], c7['heldout_mean']) == (0.25, 0.375)
    assert c7['train_sd'] == pytest.approx(math.sqrt(0.125))
    baseline = summary['no MSA']
    assert math.isnan(baseline['hpo_loss']) and math.isnan(baseline['train_sd'])
    assert (baseline['seeds'], baseline['train_mean']) == (1, 1.0)
    assert 'nomsa_s2 has no result' in caplog.text

    stability = read(output / 'stability.csv')
    # Seeds only: a: {(0,3)} vs {} -> 0; b: {(0,3)} vs {(1,2)} -> 0.
    assert (stability['c7a']['split'], float(stability['c7a']['pairwise_f1'])) == ('train', 0.0)
    assert stability['c7b']['split'] == 'heldout'
    assert stability['no MSAa']['pairwise_f1'] == 'nan'  # one seed: nothing to compare
    assert sorted(path.name for path in output.iterdir()) == [
        'seed_spread.png', 'stability.csv', 'stability.png', 'summary.csv', 'trials.csv']


def test_an_rna_id_with_another_sequence_than_in_the_hpo_dataset_is_an_error(exports, tmp_path):
    checks, hpo = exports
    write(hpo / 'run_metadata.json', dict(dataset=[dict(TRAIN, sequence='AAAA')]))
    with pytest.raises(ValueError, match='another sequence'):
        analyze(checks, hpo, tmp_path)
