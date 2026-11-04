"""Retention: which NePS configs keep their per-RNA archives.

With ``storage.keep_artifacts: trajectory``, a successful config keeps its
archives while it is on NePS's best-config trajectory or among the
``keep_top_k`` lowest losses. Other successful configs keep only their
``trial_result.json`` (and NePS's own files). Failed trials are never pruned.

Kept configs can only drop out, never come back: more finished trials only
make the trajectory and the top k stricter. Pruning a config as soon as it
drops out is therefore final and safe.
"""

import json
import math
from pathlib import Path
from typing import Optional

from hpo.pipeline import archive
from hpo.pipeline.run import write_json

TRIAL_RESULT = 'trial_result.json'


def kept_trials(trials: dict[str, tuple[float, float]], top_k: int) -> dict[str, list[str]]:
    """Reasons why each successful trial keeps its archives.

    ``trials`` maps trial IDs to ``(time_sampled, loss)``. Trials not in the
    result are pruned. The trajectory follows NePS: in order of sampling, a
    trial is on it when no earlier-sampled trial has a lower loss (ties stay).
    """
    reasons: dict[str, list[str]] = {}
    best = math.inf
    for trial_id, (_, loss) in sorted(trials.items(), key=lambda item: item[1][0]):
        if loss <= best:
            reasons.setdefault(trial_id, []).append('trajectory')
            best = loss
    for trial_id in sorted(trials, key=lambda trial_id: (trials[trial_id][1], trials[trial_id][0]))[:top_k]:
        reasons.setdefault(trial_id, []).append('top_k')
    return reasons


def record(trial_directory: Path, reasons: Optional[list[str]], pruned_by: Optional[str] = None) -> None:
    """Note the retention decision in the trial's result, for transparency."""
    path = trial_directory / TRIAL_RESULT
    result = json.loads(path.read_text())
    result['retention'] = dict(kept=bool(reasons), reasons=reasons or [])
    if pruned_by is not None:
        result['retention']['pruned_after'] = pruned_by
    if not reasons:
        result['storage'] = dict(files=0, bytes=0)
    write_json(path, result)


def prune(trial_directory: Path, dry_run: bool = False) -> dict[str, int]:
    """Delete a successful trial's archives; return the files and bytes removed.

    A trial without ``trial_result.json`` (failed or unfinished) is never pruned.
    """
    if not (trial_directory / TRIAL_RESULT).is_file():
        return dict(files=0, bytes=0)
    archives = sorted(trial_directory.glob(f'rna_*{archive.SUFFIX}'))
    removed = dict(files=len(archives), bytes=sum(path.stat().st_size for path in archives))
    if not dry_run:
        for path in archives:
            path.unlink()
    return removed
