"""Submit a standalone trial with YAML resources and workspace logs; run from any cwd."""
import argparse
from hpo.pipeline.run import (
    ROOT,
    repo_path,
    load_config,
    get_workspace_path,
    command,
    create_dataset,
    fingerprint,
    check_existing_run,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='hpo/examples/helix.example.yaml')
    args = parser.parse_args()
    config = load_config(
        args.config,
        mode='standalone',
        submitting_controller=True,
    )
    ws_root = get_workspace_path(config)
    rows = create_dataset(config['input'])
    identity = fingerprint(config, rows)

    try:
        existing_run = check_existing_run(ws_root, identity)
    except ValueError as exc:
        raise SystemExit(f'Submission cancelled: {exc}') from exc

    if existing_run:
        answer = input(
            "A run with matching inputs/code/config already exists. "
            "Resume it? [y/N] "
        )
        if answer.strip().lower() not in {'y', 'yes'}:
            raise SystemExit('Submission cancelled.')

    logs = ws_root / 'controller_logs'
    logs.mkdir(parents=True, exist_ok=True)
    resources = config['controller']

    print(command([
        'sbatch',
        '--parsable',
        '--job-name=shs-standalone',
        '--chdir=' + str(ROOT),
        '--nodes=1',
        '--ntasks=1',
        '--partition=' + resources['partition'],
        '--time=' + resources['time'],
        '--cpus-per-task=' + str(resources['cpus']),
        '--mem=' + resources['memory'],
        '--gres=' + resources['gres'],
        '--output=' + str(logs / '%j.out'),
        '--error=' + str(logs / '%j.err'),
        str(ROOT / 'hpo/slurm/controller.slurm'),
        config['module'],
        config['controller_python'],
        'hpo.run_standalone_trial',
        str(repo_path(args.config)),
    ]))


if __name__ == '__main__':
    main()