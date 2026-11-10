"""make_checks: a check run config and dataset from a small fake export."""
import csv
import json
from pathlib import Path

import pytest
import yaml

import hpo.pipeline.run as run
from hpo.make_checks import make_checks
from hpo.pipeline.run import create_dataset, load_config

FIXTURE = Path(__file__).parent / 'fixtures' / 'valid_standalone_config.yaml'
# As create_dataset returns them: one RNA and one protein chain.
ROWS = [
    {'id': 'a', 'chains': [{'id': 'A', 'type': 'rna', 'sequence': 'ACGU'}], 'type': 'rna',
     'sequence': 'ACGU', 'pairs': [[0, 3]], 'interactions': [[0, 3, 0.5]]},
    {'id': 'b', 'chains': [{'id': 'A', 'type': 'protein', 'sequence': 'MKTA'}], 'type': 'protein',
     'sequence': 'MKTA', 'pairs': []},
]


@pytest.fixture
def export(tmp_path):
    search = {'mutation_rate_paired': {'type': 'float', 'lower': 0, 'upper': 1},
              'N_extra': {'type': 'integer', 'lower': 1, 'upper': 9},
              'pair_mutation_approach': {'type': 'categorical', 'choices': ['covariance', 'none']}}
    hpo = yaml.safe_load(FIXTURE.read_text())
    hpo['generators']['rna'].update(fixed={'N': 64}, search=search)
    hpo['generators']['protein'] = {'python': '/example/protein/bin/python', 'search': {'N': search['N_extra']}}
    hpo.update(
        input='/helix/exp/v1/dataset.csv', shs_seed=2, af3_seed=2,
        evaluations=10, failures={'max_rate': 0.2},
        controller=dict(count=2, partition='gpu-single', cpus=1, memory='20G', gres='gpu:A100:1', time='04:00:00'),
    )
    results = tmp_path / 'v1' / 'results'
    (results / 'summary').mkdir(parents=True)
    (results / 'run_metadata.json').write_text(json.dumps(dict(config=hpo, dataset=ROWS, fingerprint='x')))
    with (results / 'summary' / 'full.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['id', 'config.mutation_rate_paired', 'config.N_extra', 'config.pair_mutation_approach',
                         'objective_to_minimize', 'state', 'evaluation_duration'])
        writer.writerows([
            [1, 0.1, 3.0, 'none', 0.5, 'success', 10],
            [2, 0.2, 4.0, 'covariance', 0.3, 'success', 10],
            [3, 0.3, 5.0, 'none', 0.3, 'success', 10],   # tie with 2: lower id first
            [4, 0.4, 6.0, 'none', 0.0, 'crashed', 10],    # never a top config
        ])
    return results


def test_trials_reproduce_reseed_and_compare_with_the_baseline(export):
    config = make_checks(export, export.parent / 'checks', top_n=2, seeds=2)
    trials = {trial['name']: trial for trial in config['trials']}
    # Seeds 1 and 3: the HPO seed 2 is kept for the reproduction trials.
    assert list(trials) == ['c2_hpo', 'c2_s1', 'c2_s3', 'c3_hpo', 'c3_s1', 'c3_s3', 'nomsa_s1', 'nomsa_s3']
    assert trials['c2_hpo'] == dict(name='c2_hpo', check='reproduction', base=2, shs_seed=2, af3_seed=2,
                                    parameters={'mutation_rate_paired': 0.2, 'N_extra': 4,
                                                'pair_mutation_approach': 'covariance'})
    assert (trials['c3_s3']['check'], trials['c3_s3']['shs_seed'], trials['c3_s3']['af3_seed']) == ('seed', 3, 3)
    assert trials['nomsa_s1'] == dict(name='nomsa_s1', check='baseline', shs_seed=1, af3_seed=1,
                                      parameters={'N': 1})
    assert config['input'] == '/helix/exp/v1/checks/dataset'
    assert config['generators']['rna']['fixed'] == {'N': 64}
    assert config['optimize'] == 'rna'
    # Search spaces of every generator are NePS-only.
    assert not any('search' in generator for generator in config['generators'].values())
    assert not {'search', 'evaluations', 'failures'} & config.keys()
    assert 'count' not in config['controller']


def test_written_config_is_a_valid_standalone_run(export, monkeypatch):
    monkeypatch.setattr(run, 'repo_path', lambda value: Path(value))
    make_checks(export, export.parent / 'checks', top_n=1, seeds=1, input_path='/helix/checks/dataset')
    loaded = load_config(export.parent / 'checks' / 'config.yaml', mode='standalone', submitting_controller=True)
    assert loaded['run_name'] == 'test-run-checks'
    assert loaded['input'] == '/helix/checks/dataset'
    assert [trial['name'] for trial in loaded['trials']] == ['c2_hpo', 'c2_s1', 'nomsa_s1']


def test_dataset_holds_the_training_rows_and_the_heldout_rnas(export, tmp_path):
    heldout = tmp_path / 'heldout.csv'
    heldout.write_text('id,sequence,pairs\nc,GGGAAACCC,"[[0,8]]"\n')
    make_checks(export, export.parent / 'checks', heldout=heldout)
    dataset = export.parent / 'checks' / 'dataset'
    assert sorted(path.name for path in dataset.iterdir()) == ['1_train.csv', '2_heldout.csv']
    # The chain types survive the CSV round trip (a protein must not become an RNA).
    assert 'protein' in (dataset / '1_train.csv').read_text()
    rows = json.loads(json.dumps(create_dataset(dataset)))
    assert rows[:2] == ROWS
    assert rows[2]['id'] == 'c'


def test_heldout_ids_must_differ_from_the_training_ids(export, tmp_path):
    heldout = tmp_path / 'heldout.csv'
    heldout.write_text('id,sequence,pairs\na,GGGAAACCC,[]\n')
    with pytest.raises(ValueError, match='(?i)duplicate'):
        make_checks(export, export.parent / 'checks', heldout=heldout)


def test_an_existing_check_config_is_not_overwritten(export):
    make_checks(export, export.parent / 'checks')
    with pytest.raises(SystemExit, match='exists'):
        make_checks(export, export.parent / 'checks')
