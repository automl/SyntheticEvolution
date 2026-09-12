"""One shared configuration across every dataset row; synchronous Slurm controller.

Only one controller may use a run directory. Job IDs and submission intent are
persisted before waiting. On ambiguous submission, stop rather than duplicate work.
"""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import time
from pathlib import Path
import yaml
from dssr import run_dssr
from scoring import normalize_pairs, score_pairs, aggregate

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True))
    temporary.replace(path)


def repo_path(value):
    path = Path(value).expanduser()
    return (ROOT / path).resolve() if not path.is_absolute() else path.resolve()


def command(args):
    return subprocess.run(args, check=True, capture_output=True, text=True, timeout=60).stdout.strip()


def load_config(path):
    config = yaml.safe_load(repo_path(path).read_text())
    if not re.fullmatch(r'[A-Za-z0-9_-]+', config['run_name']):
        raise ValueError('run_name must contain only letters, digits, underscore or hyphen')
    return config


def dataset(path):
    path = repo_path(path)
    files = sorted(path.glob('*.csv')) if path.is_dir() else [path]
    rows, seen = [], set()
    for file in files:
        with file.open(newline='') as stream:
            for row in csv.DictReader(stream):
                name, seq = row['id'].strip(), row['sequence'].strip().upper()
                if not name or name in seen or not seq or set(seq) - set('ACGU'):
                    raise ValueError(f'Invalid/duplicate ID or RNA sequence in {file}: {name}')
                seen.add(name)
                pairs = sorted(normalize_pairs(json.loads(row['pairs']), len(seq)))
                rows.append(dict(id=name, sequence=seq, pairs=pairs))
    if not rows:
        raise ValueError('Dataset is empty')
    return rows


def workspace(config):
    raw = config.get('workspace_path') or command(['ws_find', config['workspace_name']])
    path = Path(raw).expanduser()
    if not path.is_absolute() or not path.is_dir():
        raise ValueError(f'Workspace must exist and have an absolute path: {raw}')
    return path.resolve() / 'shs_hpo' / config['run_name']


def fingerprint(config, rows):
    # Guard against reusing outcomes after input, generator or pipeline code changes.
    code = {}
    for file in sorted(HERE.glob('*.py')) + sorted((HERE / 'slurm').glob('*.slurm')):
        code[str(file.relative_to(ROOT))] = hashlib.sha256(file.read_bytes()).hexdigest()
    for file in sorted(repo_path(config['generator_dir']).glob('*.py')):
        code[str(file)] = hashlib.sha256(file.read_bytes()).hexdigest()
    scientific = {k:v for k,v in config.items() if k not in ('evaluations', 'controller', 'poll_seconds', 'wait_timeout_seconds')}
    return hashlib.sha256(json.dumps([scientific, rows, code], sort_keys=True).encode()).hexdigest()


def prepare(config):
    rows = dataset(config['input'])
    for name in ('shs_generator.py', 'pair_map.py', 'json_generator.py'):
        if not (repo_path(config['generator_dir']) / name).is_file():
            raise FileNotFoundError(repo_path(config['generator_dir']) / name)
    root = workspace(config)
    root.mkdir(parents=True, exist_ok=True)
    identity = fingerprint(config, rows)
    guard = root / 'run.json'
    if guard.exists() and json.loads(guard.read_text())['fingerprint'] != identity:
        raise ValueError('Run inputs/code/config changed. Choose a new run_name.')
    write_json(guard, dict(fingerprint=identity, config=config, dataset=rows))
    return root, rows, identity


def wait_job(job, count, config):
    started = time.monotonic()
    while True:
        # sacct must list every array element with successful state AND exit code.
        accounting = command(['sacct', '-j', job, '--noheader', '--parsable2', '--format=JobID%64,State%40,ExitCode'])
        records = {}
        for line in accounting.splitlines():
            fields = line.split('|')
            if len(fields) >= 3:
                records[fields[0].strip()] = (fields[1].split()[0].rstrip('+'), fields[2].strip())
        expected = [records.get(f'{job}_{i}') for i in range(count)]
        failed = [x for x in expected if x and x[0] in ('FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE')]
        if any(x and x[0] == 'COMPLETED' and x[1] != '0:0' for x in expected):
            raise RuntimeError(f'Array {job} has a nonzero completed-task exit code')
        if failed:
            raise RuntimeError(f'Array {job} failed: {failed}; inspect logs. No partial score returned.')
        if all(x == ('COMPLETED','0:0') for x in expected):
            return
        if time.monotonic() - started > config['wait_timeout_seconds']:
            raise TimeoutError(f'Waiting for array {job} timed out; job may still run. Saved ID permits reattachment.')
        print(f'Waiting for AF3 array {job}: {sum(x == ("COMPLETED","0:0") for x in expected)}/{count}', flush=True)
        time.sleep(config['poll_seconds'])


def run_pipeline(parameters, config, rows, directory, identity):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    merged = {**config.get('fixed', {}), **parameters}
    signature = dict(fingerprint=identity, parameters=merged)
    state_path = directory / 'trial.json'
    if state_path.exists() and json.loads(state_path.read_text()) != signature:
        raise ValueError('Trial parameters changed; choose a fresh trial directory')
    write_json(state_path, signature)
    result = directory / 'result.json'
    if result.exists():
        return json.loads(result.read_text())['loss']
    marker = directory / 'job.json'
    if not marker.exists():
        tasks = []
        for index, row in enumerate(rows):
            task_dir = directory / f'rna_{index:05d}'
            task_dir.mkdir(exist_ok=True)
            request = dict(**row, name=f'rna_{index:05d}', parameters=merged,
                           shs_seed=config['shs_seed'], af3_seed=config['af3_seed'])
            request_file, input_file = task_dir / 'request.json', task_dir / 'input.json'
            write_json(request_file, request)
            with (task_dir / 'generator.log').open('w') as log:
                subprocess.run([config['shs_python'], str(HERE / 'generate_input.py'),
                                str(repo_path(config['generator_dir'])), str(request_file), str(input_file)],
                               stdout=log, stderr=subprocess.STDOUT, check=True,
                               timeout=config['generator_timeout_seconds'])
            tasks.append(dict(input=str(input_file), output=str(task_dir / 'af3'),
                              model_dir=str(Path(config['model_dir']).expanduser())))
        manifest = directory / 'manifest.json'
        write_json(manifest, tasks)
        gpu = config['gpu']
        logs = directory / 'logs'
        logs.mkdir(exist_ok=True)
        args = ['sbatch', '--parsable', '--job-name=shs-' + directory.name,
                '--array=0-' + str(len(rows)-1) + '%' + str(gpu['concurrency']),
                '--partition=' + gpu['partition'], '--cpus-per-task=' + str(gpu['cpus']),
                '--mem=' + gpu['memory'], '--gres=' + gpu['gres'], '--time=' + gpu['time'],
                '--output=' + str(logs / '%A_%a.out'), '--error=' + str(logs / '%A_%a.err'),
                str(HERE / 'slurm/af3-inference-array.slurm'), str(HERE / 'af3_task.py'),
                str(manifest), config['module']]
        intent = directory / 'submission-intent.json'
        if intent.exists():
            raise RuntimeError(f'Ambiguous earlier submission. Inspect {intent} and Slurm before recovery; do not resubmit blindly.')
        write_json(intent, dict(command=args))
        job = command(args).split(';')[0]
        if not job.isdigit():
            raise RuntimeError('Unrecognized sbatch job ID; inspect submission intent')
        write_json(marker, dict(job_id=job))
    job = json.loads(marker.read_text())['job_id']
    wait_job(job, len(rows), config)
    scores = []
    for index, row in enumerate(rows):
        task_dir = directory / f'rna_{index:05d}'
        model = json.loads((task_dir / 'af3/success.json').read_text())['model']
        if not Path(model).is_file():
            raise FileNotFoundError(model)
        predicted = run_dssr(model, task_dir / 'dssr', row['sequence'], config['dssr'], config['dssr_timeout_seconds'])
        metrics = score_pairs(map(tuple,row['pairs']), predicted)
        item = dict(id=row['id'], **metrics, target_pairs=row['pairs'], predicted_pairs=sorted(predicted), model=model)
        write_json(task_dir / 'score.json', item)
        scores.append(item)
    loss = aggregate(scores)
    write_json(result, dict(loss=loss, scores=scores, parameters=merged, job_id=job))
    return loss


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='hpo/configs/helix.example.yaml')
    parser.add_argument('--resume-trial', help='Absolute artifacts/ or standalone/ directory to reattach without resubmission')
    parser.add_argument('--validate-only', action='store_true')
    parser.add_argument('--parameters', help='Repo-relative JSON parameter file for one standalone trial')
    args = parser.parse_args()
    config = load_config(args.config)
    if args.validate_only:
        print(json.dumps(dataset(config['input']), indent=2))
        return
    root, rows, identity = prepare(config)
    if args.resume_trial:
        directory = Path(args.resume_trial).resolve()
        if not directory.is_relative_to(root):
            raise ValueError('Recovery directory must be within this run')
        saved = json.loads((directory / 'trial.json').read_text())
        print(run_pipeline(saved['parameters'], config, rows, directory, identity))
        return
    parameters = json.loads(repo_path(args.parameters).read_text()) if args.parameters else dict(N=20, mutation_rate_paired=.2, mutation_rate_unpaired=.2, pair_mutation_approach='watson_crick_cov')
    print(run_pipeline(parameters, config, rows, root / 'standalone', identity))


if __name__ == '__main__':
    main()
