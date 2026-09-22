"""Submit controller with YAML resources and workspace logs; run from any cwd."""
import argparse
from hpo.pipeline.run import ROOT, repo_path, load_config, get_workspace_path, command


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='hpo/examples/helix.example.yaml')
    args = parser.parse_args()
    config = load_config(
        args.config,
        mode='neps',
        submitting_controller=True,
    )
    logs = get_workspace_path(config) / 'controller_logs'
    logs.mkdir(parents=True, exist_ok=True)
    resources = config['controller']
    count = resources.get('count', 1)
    log_name = '%A_%a' if count > 1 else '%j'

    submission = [
        'sbatch',
        '--parsable',
        '--job-name=shs-hpo',
        '--chdir=' + str(ROOT),
        '--nodes=1',
        '--ntasks=1',
        '--partition=' + resources['partition'],
        '--time=' + resources['time'],
        '--cpus-per-task=' + str(resources['cpus']),
        '--mem=' + resources['memory'],
        '--gres=' + resources['gres'],
        '--output=' + str(logs / (log_name + '.out')),
        '--error=' + str(logs / (log_name + '.err')),
    ]
    if count > 1:
        submission.append(f'--array=0-{count - 1}')

    submission.extend([
        str(ROOT / 'hpo/slurm/controller.slurm'),
        config['module'],
        config['neps_python'],
        'hpo.run_neps_hpo',
        str(repo_path(args.config)),
    ])
    print(command(submission))


if __name__ == '__main__':
    main()