"""Submit controller with YAML resources and workspace logs; run from any cwd."""
import argparse
from run_shs_dssr_pipeline import HERE, repo_path, load_config, workspace, command


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='hpo/configs/helix.example.yaml')
    args = parser.parse_args()
    config = load_config(args.config)
    logs = workspace(config) / 'controller_logs'
    logs.mkdir(parents=True, exist_ok=True)
    resources = config['controller']
    print(command(['sbatch', '--parsable', '--job-name=shs-hpo',
        '--partition=' + resources['partition'], '--time=' + resources['time'],
        '--cpus-per-task=' + str(resources['cpus']), '--mem=' + resources['memory'],
        '--output=' + str(logs / '%j.out'), '--error=' + str(logs / '%j.err'),
        str(HERE / 'slurm/controller.slurm'), config['neps_python'],
        str(HERE / 'run_neps_hpo.py'), str(repo_path(args.config))]))


if __name__ == '__main__':
    main()
