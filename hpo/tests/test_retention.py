"""Retention rules: which configs keep their archives in ``trajectory`` mode."""
import json
import random

import pytest

from hpo.pipeline import retention


def trials(*losses):
    """Trial IDs '1', '2', ... sampled in this order with these losses."""
    return {str(index): (float(index), loss) for index, loss in enumerate(losses, start=1)}


# ---- Which trials are kept --------------------------------------------------

def test_trajectory_follows_sampling_order_and_keeps_ties():
    kept = retention.kept_trials(trials(0.5, 0.6, 0.5, 0.4, 0.45), top_k=0)
    assert kept == {'1': ['trajectory'], '3': ['trajectory'], '4': ['trajectory']}


def test_late_finishing_earlier_sampled_trial_pushes_later_ones_off():
    # Trial 2 (sampled second) finishes last, with the best loss.
    early = {'1': (1.0, 0.5), '3': (3.0, 0.4)}
    assert set(retention.kept_trials(early, top_k=0)) == {'1', '3'}
    assert set(retention.kept_trials(dict(early, **{'2': (2.0, 0.3)}), top_k=0)) == {'1', '2'}


def test_top_k_keeps_the_lowest_losses_in_addition():
    kept = retention.kept_trials(trials(0.5, 0.6, 0.7, 0.55), top_k=2)
    assert kept == {'1': ['trajectory', 'top_k'], '4': ['top_k']}


def test_kept_trials_only_ever_drop_out():
    """Pruning a config as soon as it drops out is safe: it never comes back."""
    rng = random.Random(3)
    seen, kept_before = {}, None
    for trial_id in rng.sample(range(1, 60), 59):  # finishing order differs from sampling order
        seen[str(trial_id)] = (float(trial_id), rng.random())
        kept_now = set(retention.kept_trials(seen, top_k=3))
        if kept_before is not None:
            assert kept_now & set(seen) - {str(trial_id)} <= kept_before
        kept_before = kept_now


# ---- Pruning ----------------------------------------------------------------

def trial_directory(tmp_path, result=True):
    directory = tmp_path / 'config_1' / 'artifacts'
    directory.mkdir(parents=True)
    for index in range(2):
        (directory / f'rna_{index:05d}.tar.gz').write_bytes(b'x' * 10)
    if result:
        (directory / 'trial_result.json').write_text(json.dumps(dict(loss=0.5, storage=dict(files=2, bytes=20))))
    return directory


def test_prune_removes_only_the_archives_and_records_it(tmp_path):
    directory = trial_directory(tmp_path)
    assert retention.prune(directory) == dict(files=2, bytes=20)
    retention.record(directory, None, pruned_by='7')
    assert [path.name for path in directory.iterdir()] == ['trial_result.json']
    result = json.loads((directory / 'trial_result.json').read_text())
    assert result['retention'] == dict(kept=False, reasons=[], pruned_after='7')
    assert result['storage'] == dict(files=0, bytes=0)


def test_dry_run_and_failed_trials_keep_their_archives(tmp_path):
    directory = trial_directory(tmp_path)
    assert retention.prune(directory, dry_run=True) == dict(files=2, bytes=20)
    assert len(list(directory.glob('*.tar.gz'))) == 2
    failed = trial_directory(tmp_path / 'failed', result=False)
    assert retention.prune(failed) == dict(files=0, bytes=0)
    assert len(list(failed.glob('*.tar.gz'))) == 2


def test_kept_trial_records_its_reasons(tmp_path):
    directory = trial_directory(tmp_path)
    retention.record(directory, ['trajectory', 'top_k'])
    result = json.loads((directory / 'trial_result.json').read_text())
    assert result['retention'] == dict(kept=True, reasons=['trajectory', 'top_k'])
    assert result['storage'] == dict(files=2, bytes=20)
