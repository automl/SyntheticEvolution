"""Export of a finished run, including error files; no cluster or NePS needed."""
import json

import pytest

from hpo import export_results as export


@pytest.fixture
def run(tmp_path):
    """A minimal finished run with one trajectory config and one successful trial."""
    workspace = tmp_path / 'workspace'
    root = workspace / 'shs_hpo' / 'test-run'
    neps = root / 'neps'
    (neps / 'summary').mkdir(parents=True)
    for name in ('best_config.txt', 'short.csv'):
        (neps / 'summary' / name).write_text('summary\n')
    (neps / 'summary' / 'best_config_trajectory.txt').write_text('Config ID: 1\n')
    (root / 'run_metadata.json').write_text('{}')
    artifacts = neps / 'configs' / 'config_1' / 'artifacts'
    (artifacts / 'rna_00000').mkdir(parents=True)
    (artifacts / 'trial_result.json').write_text(json.dumps(dict(loss=0.5)))
    for name in ('af3_input.json', 'shs_generation_request.json'):
        (artifacts / 'rna_00000' / name).write_text('{}')
    (artifacts.parent / 'report.yaml').write_text('reported_as: success\n')
    config = tmp_path / 'config' / 'run.yaml'
    config.parent.mkdir()
    config.write_text(f'run_name: test-run\nworkspace_name: ws\nworkspace_path: {workspace}\n')
    return config, neps


def manifest(destination):
    return json.loads((destination / 'manifest.json').read_text())


def test_failures_block_is_written_without_failures(run, capsys):
    config, _ = run
    destination = export.export_results(config)
    assert manifest(destination)['failures'] == dict(
        status='none', count=0, finished=1, rate=0.0, groups=[], files=[],
    )
    assert not (destination / 'errors').exists()
    assert 'WARNING' not in capsys.readouterr().err


def test_run_without_successful_trials_exports_errors(run, capsys):
    config, neps = run
    (neps / 'summary' / 'best_config_trajectory.txt').write_text('')
    failed = neps / 'configs' / 'config_1'
    (failed / 'artifacts' / 'failure_report.json').write_text(json.dumps(dict(
        error_type='CalledProcessError', stage='generation', rna_id='rna_a')))
    destination = export.export_results(config)
    data = manifest(destination)
    assert data['config_ids'] == []
    assert data['failures']['count'] == 1
    assert (destination / 'errors' / 'reports' / 'config_1.json').is_file()
    assert not (destination / 'configs').exists()
    assert 'no trial succeeded' in capsys.readouterr().err


def test_failures_are_copied_and_reported(run, capsys):
    config, neps = run
    failed = neps / 'configs' / 'config_2'
    (failed / 'artifacts').mkdir(parents=True)
    (failed / 'report.yaml').write_text('reported_as: crashed\n')
    (failed / 'artifacts' / 'failure_report.json').write_text(json.dumps(dict(
        error_type='CalledProcessError', signal='SIGSEGV', stage='dssr', rna_id='rna_b',
    )))
    # Trial 3 has no failure report (older run); trial 2 is also in NePS's log.
    # The pasted shell prompt mirrors the saved file from the measured run.
    (neps / 'shared_errors.jsonl').write_text('\n'.join([
        '(neps) user@helix:~$ cat shared_errors.jsonl',
        json.dumps(dict(trial_id='2', err_type='CalledProcessError', err='x')),
        json.dumps(dict(trial_id='3', err_type='CalledProcessError',
                        err="Command '['dssr']' died with <Signals.SIGSEGV: 11>.")),
    ]) + '\n')
    (neps / 'configs' / 'config_3').mkdir()
    (neps / 'configs' / 'config_3' / 'report.yaml').write_text('reported_as: crashed\n')

    destination = export.export_results(config)

    block = manifest(destination)['failures']
    assert (block['status'], block['count'], block['finished'], block['rate']) == ('found', 2, 3, 0.6667)
    assert [(g['count'], g['description']) for g in block['groups']] == [
        (1, 'SIGSEGV in dssr on rna_b'), (1, 'SIGSEGV in unknown stage'),
    ]
    assert block['files'] == ['errors/reports/config_2.json', 'errors/shared_errors.jsonl']
    for name in block['files']:
        assert (destination / name).is_file()
    err = capsys.readouterr().err
    assert 'Skipped 1 non-JSON line(s)' in err
    assert 'WARNING: 2 of 3 trials failed (66.7%)' in err
    assert 'See results/errors/' in err
