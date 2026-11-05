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
    for name in ('best_config.txt', 'short.csv', 'full.csv'):
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


# ---- Failures ---------------------------------------------------------------

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


# ---- SHS files from the per-RNA archives -------------------------------------

@pytest.fixture
def archived_run(run):
    """The same run in the archive layout: no rna_<index>/ folders."""
    from hpo.pipeline import archive
    config, neps = run
    root = neps.parent
    rows = [dict(id='ex_1', sequence='ACGU', pairs=[[0, 3]]),
            dict(id='ex_2', sequence='AAAA', pairs=[])]
    (root / 'run_metadata.json').write_text(json.dumps(dict(
        fingerprint='f', config=dict(shs_seed=11, af3_seed=2), dataset=rows)))
    artifacts = neps / 'configs' / 'config_1' / 'artifacts'
    import shutil
    shutil.rmtree(artifacts / 'rna_00000')
    (artifacts / 'trial_result.json').write_text(json.dumps(dict(loss=0.5, parameters={'N': 20})))
    for index in range(2):
        content = root.parent / f'content_{index}'
        content.mkdir()
        (content / f'rna_{index:05d}_model.cif').write_text('model')
        (content / f'rna_{index:05d}_data.json').write_text(json.dumps({'msa': index}))
        archive.publish(content, artifacts / f'rna_{index:05d}.tar.gz')
    return config, neps, rows


def test_export_derives_shs_files_from_archives(archived_run):
    config, _, rows = archived_run
    destination = export.export_results(config)
    shs = destination / 'configs' / 'config_1' / 'shs'
    for index in range(2):
        task = shs / f'rna_{index:05d}'
        assert json.loads((task / 'af3_input.json').read_text()) == {'msa': index}
        # Rebuilt exactly as the trial built it for the generator.
        assert json.loads((task / 'shs_generation_request.json').read_text()) == dict(
            rows[index], task_name=f'rna_{index:05d}', parameters={'N': 20}, shs_seed=11, af3_seed=2)
    assert not (destination / 'configs' / 'config_1' / 'archives').exists()
    assert 'configs/config_1/shs/rna_00001/af3_input.json' in manifest(destination)['files']


def test_export_can_include_archives(archived_run):
    config, neps, _ = archived_run
    destination = export.export_results(config, include_archives=True)
    copied = destination / 'configs' / 'config_1' / 'archives'
    assert sorted(path.name for path in copied.iterdir()) == ['rna_00000.tar.gz', 'rna_00001.tar.gz']


def test_export_rejects_missing_archives(archived_run):
    config, neps, _ = archived_run
    (neps / 'configs' / 'config_1' / 'artifacts' / 'rna_00001.tar.gz').unlink()
    with pytest.raises(FileNotFoundError, match='Expected 2 RNA archives'):
        export.export_results(config)


# ---- Pruning after the run --------------------------------------------------

@pytest.fixture
def finished_run(tmp_path):
    """Five configs in the archive layout; 4 crashed, the others succeeded."""
    import csv
    workspace = tmp_path / 'workspace'
    neps = workspace / 'shs_hpo' / 'test-run' / 'neps'
    (neps / 'summary').mkdir(parents=True)
    losses = {'1': .5, '2': .3, '3': .6, '4': None, '5': .4}
    with (neps / 'summary' / 'full.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['id', 'state', 'time_sampled', 'objective_to_minimize'])
        for trial_id, loss in losses.items():
            writer.writerow([trial_id, 'crashed' if loss is None else 'success', trial_id, '' if loss is None else loss])
    (neps / 'summary' / 'best_config_trajectory.txt').write_text('Config ID: 1\nConfig ID: 2\n')
    for trial_id, loss in losses.items():
        artifacts = neps / 'configs' / f'config_{trial_id}' / 'artifacts'
        artifacts.mkdir(parents=True)
        (artifacts / 'rna_00000.tar.gz').write_bytes(b'x' * 1000)
        if loss is not None:
            (artifacts / 'trial_result.json').write_text(json.dumps(dict(loss=loss)))
    config = tmp_path / 'run.yaml'
    config.write_text(f'run_name: test-run\nworkspace_name: ws\nworkspace_path: {workspace}\n'
                      'storage: {keep_top_k: 1}\n')
    return config, neps


def configs_with_archives(neps):
    return sorted(path.parent.parent.name for path in neps.glob('configs/*/artifacts/*.tar.gz'))


def test_prune_keeps_trajectory_top_k_and_failed_trials(finished_run, capsys):
    config, neps = finished_run
    totals = export.prune_run(config)
    # Trajectory 1, 2 and top 1 (config 2) stay; config 4 crashed and stays.
    assert configs_with_archives(neps) == ['config_1', 'config_2', 'config_4']
    assert totals == dict(configs=2, files=2, bytes=2000)
    result = json.loads((neps / 'configs/config_3/artifacts/trial_result.json').read_text())
    assert result['retention'] == dict(kept=False, reasons=[], pruned_after='export_results --prune')
    assert 'Removed 2 archives' in capsys.readouterr().out


def test_prune_dry_run_changes_nothing(finished_run, capsys):
    config, neps = finished_run
    assert export.prune_run(config, dry_run=True)['files'] == 2
    assert len(configs_with_archives(neps)) == 5
    assert 'Would remove 2 archives' in capsys.readouterr().out
