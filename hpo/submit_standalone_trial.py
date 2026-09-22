"""Submit a standalone trial with YAML resources and workspace logs; run from any cwd."""
import argparse
from hpo.pipeline.run import ROOT, repo_path, load_config, get_workspace_path, command


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='hpo/examples/helix.example.yaml')
    args = parser.parse_args()
    config = load_config(
        args.config,
        mode='standalone',
        submitting_controller=True,
    )
    logs = get_workspace_path(config) / 'controller_logs'
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