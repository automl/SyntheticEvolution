"""Submit controller with YAML resources and workspace logs; run from any cwd."""
import argparse
import shlex
from hpo.pipeline.trial import ROOT, repo_path, load_config, get_workspace_path, command


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

    # Command executed on the allocated CPU node.
    # shlex.join quotes each argument safely for the shell used by --wrap.
    controller_command = shlex.join([
        config['neps_python'],          # Python to use
        '-u',                           # Write log output without buffering
        '-m',
        'hpo.run_neps_hpo',             # NePS controller module
        '--config',
        str(repo_path(args.config)),    # Absolute path to config file
    ])

    print(command([
        'sbatch',
        '--parsable',
        '--job-name=shs-hpo',
        '--chdir=' + str(ROOT),
        '--nodes=1',                    # Previously set in controller.slurm
        '--ntasks=1',                   # Previously set in controller.slurm
        '--partition=' + resources['partition'],
        '--time=' + resources['time'],
        '--cpus-per-task=' + str(resources['cpus']),
        '--mem=' + resources['memory'],
        '--output=' + str(logs / '%j.out'),
        '--error=' + str(logs / '%j.err'),
        # Slurm creates the shell wrapper; exec replaces that shell with Python.
        '--wrap=exec ' + controller_command,
    ]))


if __name__ == '__main__':
    main()