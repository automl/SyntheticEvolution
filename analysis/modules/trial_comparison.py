"""Compare parameters and predicted pairs across trajectory configurations.

Parameters use exact equality, without rounding or floating-point tolerance.
Predictions are compared per RNA ID as sets of undirected pairs: ordering,
reversed endpoints, and duplicate pairs do not affect equality. Indices retain
their original numbering. Missing parameters/IDs are reported as differences.
"""

import json
import logging
from pathlib import Path
from typing import Any

from analysis.results_loader import ResultsContext

LOGGER = logging.getLogger("analysis.trial_comparison")
MISSING = object()


def _flatten_parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    """Per-chain-type parameters ({'rna': {'N': 8}}) become {'rna.N': 8}.

    Generator parameters are scalars, so a dict value marks the per-type format.
    """
    flat = {}
    for key, value in parameters.items():
        if isinstance(value, dict):
            flat.update({f"{key}.{name}": item for name, item in value.items()})
        else:
            flat[key] = value
    return flat


def _parse_trial(result: dict[str, Any]) -> tuple[dict[str, Any], dict[str, frozenset]]:
    if not isinstance(result.get("parameters"), dict):
        raise ValueError("Expected a 'parameters' object")
    if not isinstance(result.get("scores"), list) or not result["scores"]:
        raise ValueError("Expected a nonempty 'scores' list")

    predictions = {}
    for score in result["scores"]:
        if not isinstance(score, dict):
            raise ValueError("Each score must be an object")
        rna_id = score.get("id")
        if not isinstance(rna_id, str) or not rna_id:
            raise ValueError("Each score needs a nonempty string 'id'")
        if rna_id in predictions:
            raise ValueError(f"Duplicate RNA ID: {rna_id}")
        raw_pairs = score.get("predicted_pairs")
        if not isinstance(raw_pairs, list):
            raise ValueError(f"{rna_id}: missing/invalid 'predicted_pairs'")
        pairs = set()
        for pair in raw_pairs:
            if (not isinstance(pair, list) or len(pair) != 2
                    or any(type(i) is not int or i < 0 for i in pair)
                    or pair[0] == pair[1]):
                raise ValueError(f"{rna_id}: invalid pair {pair!r}")
            pairs.add(tuple(sorted(pair)))
        predictions[rna_id] = frozenset(pairs)
    return _flatten_parameters(result["parameters"]), predictions


def _display(value: Any) -> str:
    return "<MISSING>" if value is MISSING else json.dumps(value, sort_keys=True)


def _format_pairs(pairs: set) -> str:
    return json.dumps(sorted(pairs))


def compare(names: list[str], trials: list[tuple[dict, dict]]) -> tuple[list[str], bool]:
    """Return the report lines and whether all predictions are identical."""
    lines = [f"Comparing {len(names)} configurations"]
    lines += [f"  {name}" for name in names]

    parameters = [trial[0] for trial in trials]
    keys = sorted(set().union(*(set(p) for p in parameters)))
    same, different = [], []
    for key in keys:
        values = [p.get(key, MISSING) for p in parameters]
        if all(v is not MISSING and v == values[0] for v in values):
            same.append((key, values[0]))
        else:
            different.append((key, values))

    lines.append(
        f"\nPARAMETERS IDENTICAL IN ALL {len(names)} CONFIGURATIONS "
        f"({len(same)}/{len(keys)})"
    )
    lines += [f"  {key} = {_display(value)}" for key, value in same]
    if not same:
        lines.append("  None")
    lines.append(f"\nDIFFERING OR MISSING PARAMETERS ({len(different)})")
    for key, values in different:
        lines.append(f"  {key}:")
        lines += [f"    {name}: {_display(value)}" for name, value in zip(names, values)]
    if not different:
        lines.append("  None")

    predictions = [trial[1] for trial in trials]
    ids = sorted(set().union(*(set(p) for p in predictions)))
    all_identical = True
    lines.append("\nPREDICTED PAIRS BY RNA ID (set comparison)")
    for rna_id in ids:
        values = [p.get(rna_id, MISSING) for p in predictions]
        identical = all(v is not MISSING and v == values[0] for v in values)
        all_identical &= identical
        if identical:
            lines.append(
                f"  {rna_id}: IDENTICAL in all configurations ({len(values[0])} pairs)"
            )
            continue
        lines.append(f"  {rna_id}: DIFFERENT or MISSING")
        present = [v for v in values if v is not MISSING]
        common = set.intersection(*(set(v) for v in present))
        union = set.union(*(set(v) for v in present))
        scope = (
            "all configurations" if len(present) == len(values)
            else "configurations containing this RNA ID"
        )
        lines.append(f"    Common to {scope}: {_format_pairs(common)}")
        for name, pairs in zip(names, values):
            if pairs is MISSING:
                lines.append(f"    {name}: RNA ID MISSING")
            else:
                lines.append(f"    {name}: {len(pairs)} pairs")
                lines.append(f"      Beyond common: {_format_pairs(pairs - common)}")
                lines.append(f"      Absent from union: {_format_pairs(union - pairs)}")

    lines.append(
        "\nAre predicted pairs identical across all configurations for every RNA ID? "
        + ("YES" if all_identical else "NO")
    )
    return lines, all_identical


def run(context: ResultsContext, output_dir: Path) -> None:
    """Write a parameter and prediction comparison of all trajectory configurations."""
    if len(context.config_ids) < 2:
        LOGGER.warning(
            "Skipping trial comparison: trajectory has only %d configuration",
            len(context.config_ids),
        )
        return
    names = [f"config_{config_id}" for config_id in context.config_ids]
    trials = []
    for name, config_id in zip(names, context.config_ids):
        try:
            trials.append(_parse_trial(context.trial_results[config_id]))
        except ValueError as exc:
            raise ValueError(f"{name}/trial_result.json: {exc}") from exc

    lines, all_identical = compare(names, trials)
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "trial_comparison.txt"
    report_path.write_text("\n".join(lines) + "\n")
    LOGGER.info(
        "Wrote trial comparison to %s (predictions identical: %s)",
        report_path,
        "yes" if all_identical else "no",
    )
