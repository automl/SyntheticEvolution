"""Write a check run (standalone config with `trials`) for the best configs of an HPO export.

Runs locally on an export (`results/`). Writes `<out>/config.yaml` and
`<out>/dataset/` (the HPO training rows and an optional held-out CSV). Trials:

- `c<id>_hpo`: each of the top n configs with the HPO seeds (reproduction),
- `c<id>_s<seed>`: each of the top n configs with k new seed pairs,
- `nomsa_s<seed>`: baseline `N: 1` (only the query: AF3 without MSA), same seeds.
"""
import argparse
import csv
import json
import shutil
import statistics
from itertools import count, islice
from pathlib import Path, PurePosixPath

import yaml

from hpo.pipeline.run import create_dataset

# HPO settings a standalone run does not use.
NEPS_ONLY = ('search', 'evaluations', 'optimizer', 'ignore_errors', 'failures', 'neps_python')


def top_configs(full_csv: Path, search: dict, n: int) -> list[tuple[int, float, dict]]:
    """The n successful configs with the lowest loss (ties: lower id first).

    Values are converted to the search space types, so they equal what NePS
    passed to the pipeline.
    """
    with full_csv.open(newline='') as stream:
        rows = [row for row in csv.DictReader(stream) if row['state'] == 'success']
    rows.sort(key=lambda row: (float(row['objective_to_minimize']), int(row['id'])))
    configs = []
    for row in rows[:n]:
        parameters = {}
        for name, spec in search.items():
            value = row[f'config.{name}']
            if spec['type'] == 'integer':
                parameters[name] = int(float(value))
            elif spec['type'] == 'float':
                parameters[name] = float(value)
            else:
                parameters[name] = next(choice for choice in spec['choices'] if str(choice) == value)
        configs.append((int(row['id']), float(row['objective_to_minimize']), parameters))
    return configs


def write_rows(rows: list[dict], path: Path) -> None:
    """Write validated dataset rows back as an input CSV."""
    columns = ['id', 'sequence', 'pairs'] + [
        key for key in ('interactions', 'mutation_rates') if any(key in row for row in rows)
    ]
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: row[key] if key in ('id', 'sequence') else json.dumps(row[key])
                for key in columns if key in row
            })


def make_checks(results: Path, out: Path, heldout: Path = None, top_n: int = 3,
                seeds: int = 5, input_path: str = None) -> dict:
    metadata = json.loads((results / 'run_metadata.json').read_text())
    hpo = metadata['config']
    if (out / 'config.yaml').exists():
        raise SystemExit(f'{out / "config.yaml"} exists; delete it to write a new check run.')

    configs = top_configs(results / 'summary' / 'full.csv', hpo['search'], top_n)
    hpo_seeds = (hpo['shs_seed'], hpo['af3_seed'])
    new_seeds = list(islice((seed for seed in count(1) if seed not in hpo_seeds), seeds))
    trials = []
    for config_id, _, parameters in configs:
        trials.append(dict(name=f'c{config_id}_hpo', check='reproduction', base=config_id,
                           shs_seed=hpo_seeds[0], af3_seed=hpo_seeds[1], parameters=parameters))
        trials += [dict(name=f'c{config_id}_s{seed}', check='seed', base=config_id,
                        shs_seed=seed, af3_seed=seed, parameters=dict(parameters)) for seed in new_seeds]
    trials += [dict(name=f'nomsa_s{seed}', check='baseline', shs_seed=seed, af3_seed=seed,
                    parameters={'N': 1}) for seed in new_seeds]

    # Training rows exactly as the HPO run used them; held-out RNAs next to them.
    dataset = out / 'dataset'
    dataset.mkdir(parents=True, exist_ok=True)
    write_rows(metadata['dataset'], dataset / '1_train.csv')
    if heldout is not None:
        shutil.copyfile(heldout, dataset / '2_heldout.csv')
    rows = create_dataset(dataset)  # fails on invalid rows or ids used twice
    if json.loads(json.dumps(rows[:len(metadata['dataset'])])) != metadata['dataset']:
        raise ValueError('The written training rows differ from the HPO dataset')

    config = dict(
        run_name=f"{hpo['run_name']}-checks",
        # The config runs on Helix: by default the dataset lies in checks/ next to the HPO dataset.
        input=input_path or str(PurePosixPath(hpo['input']).parent / 'checks' / 'dataset'),
        controller={k: v for k, v in hpo['controller'].items() if k != 'count'},
        **{k: v for k, v in hpo.items() if k not in (*NEPS_ONLY, 'run_name', 'input', 'controller')},
        trials=trials,
    )
    (out / 'config.yaml').write_text(
        f'# Check run for {hpo["run_name"]}, written by hpo/make_checks.py.\n'
        + yaml.safe_dump(config, sort_keys=False, default_flow_style=None, width=120)
    )

    # Rough run time: HPO trial time scaled to the larger dataset.
    with (results / 'summary' / 'full.csv').open(newline='') as stream:
        durations = [float(row['evaluation_duration']) for row in csv.DictReader(stream)
                     if row['state'] == 'success']
    estimate = statistics.median(durations) * len(trials) * len(rows) / len(metadata['dataset'])
    print(f"Wrote {out / 'config.yaml'}: {len(trials)} trials on {len(rows)} RNAs "
          f"(~{estimate / 3600:.1f} h from the HPO trial times; controller.time is "
          f"{config['controller']['time']}).")
    print(f"Top configs (HPO loss): {', '.join(f'{i} ({loss:.4f})' for i, loss, _ in configs)}")
    print(f"input on Helix: {config['input']}")
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('results', type=Path, help='Exported HPO results directory')
    parser.add_argument('--heldout', type=Path, help='CSV of held-out RNAs (dataset format)')
    parser.add_argument('--top-n', type=int, default=3)
    parser.add_argument('--seeds', type=int, default=5, help='New seed pairs per config')
    parser.add_argument('--out', type=Path, help='Check directory; default: checks/ next to results/')
    parser.add_argument('--input', help="Dataset path written into the config (Helix); "
                                         "default: checks/dataset next to the HPO input")
    args = parser.parse_args()
    make_checks(args.results.resolve(), (args.out or args.results.resolve().parent / 'checks'),
                args.heldout, args.top_n, args.seeds, args.input)


if __name__ == '__main__':
    main()
