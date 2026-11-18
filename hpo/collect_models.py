"""Copy the top-ranked AF3 model of every RNA out of its archive, for evaluation/.

Usage: python -m hpo.collect_models <archive dir> <output dir> [--seed S]

<archive dir> is either
- a trial directory of a run (``<run>/trials/<name>/`` or
  ``<run>/standalone_trial/artifacts/``) with ``rna_<index>.tar.gz``; the IDs come
  from ``<run>/run_metadata.json``, or
- a directory of ``<ID>.tar.gz``, each a packed AF3 output folder (e.g. written by
  a plain ``run_alphafold.py`` run).

Writes ``<output dir>/<ID>.cif``: the sample with the highest ranking score in
``ranking_scores.csv``, among the samples of seed S if given. For a trial (one
seed) this is AF3's top-ranked model, the one DSSR scored.
Runs with Python 3.9 and the standard library only.
"""
import argparse
import csv
import io
import json
import re
import tarfile
from pathlib import Path

SUFFIX = '.tar.gz'


def member_name(member):
    return member.name[2:] if member.name.startswith('./') else member.name


def best_model(archive, seed):
    with tarfile.open(archive, 'r:gz') as tar:
        members = {member_name(m): m for m in tar.getmembers() if m.isfile()}
        with tar.extractfile(members['ranking_scores.csv']) as stream:
            rows = list(csv.DictReader(io.TextIOWrapper(stream)))
        if seed is not None:
            rows = [row for row in rows if int(row['seed']) == seed]
        if not rows:
            raise SystemExit(f'{archive}: no sample of seed {seed}')
        best = max(rows, key=lambda row: float(row['ranking_score']))
        with tar.extractfile(members[f"seed-{best['seed']}_sample-{best['sample']}/model.cif"]) as stream:
            return stream.read()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('archives', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--seed', type=int)
    args = parser.parse_args()

    archives = sorted(args.archives.glob(f'*{SUFFIX}'))
    if not archives:
        raise SystemExit(f'No archives in {args.archives}')
    names = [path.name[:-len(SUFFIX)] for path in archives]
    if all(re.fullmatch(r'rna_\d+', name) for name in names):
        run = args.archives.parent.parent
        metadata = run / 'run_metadata.json'
        if not metadata.is_file():
            metadata = run.parent / 'run_metadata.json'  # standalone_trial/artifacts/
        ids = [row['id'] for row in json.loads(metadata.read_text())['dataset']]
        if len(archives) != len(ids):
            raise SystemExit(f'{len(archives)} archives for {len(ids)} RNAs in {args.archives}')
        names = [ids[int(name[len('rna_'):])] for name in names]

    args.output.mkdir(parents=True, exist_ok=True)
    for archive, name in zip(archives, names):
        (args.output / f'{name}.cif').write_bytes(best_model(archive, args.seed))
    print(f'{len(archives)} models -> {args.output}')


if __name__ == '__main__':
    main()
