"""shs_best_config: which standalone trials are plotted."""
import json
from types import SimpleNamespace

from analysis.modules.shs_best_config import _selected_config_ids


def trial(tmp_path, name, rows):
    """A trial whose one chain has an MSA of ``rows`` sequences (the query first)."""
    msa = ''.join(f'>row_{i}\nMKTA\n' for i in range(rows))
    path = tmp_path / name / 'shs' / 'rna_00000' / 'af3_input.json'
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(dict(sequences=[dict(protein=dict(id='A', unpairedMsa=msa))])))
    return path.parents[2]


def test_trials_without_msa_samples_are_skipped(tmp_path):
    # nomsa (N = 1) has the lowest loss but only the query; the best trial with samples is plotted.
    losses = dict(nomsa=0.3, default=0.5, deep=0.7)
    rows = dict(nomsa=1, default=20, deep=256)
    context = SimpleNamespace(
        analysis_config={}, standalone=True, config_ids=list(losses),
        trial_dirs={name: trial(tmp_path, name, rows[name]) for name in losses},
        trial_results={name: dict(loss=loss) for name, loss in losses.items()},
    )
    assert _selected_config_ids(context) == ['default']
