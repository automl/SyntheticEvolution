"""Run the NePS hyperparameter search for one shared RNA dataset."""
import argparse
import fcntl
from importlib.metadata import version

from hpo.pipeline.trial import load_config, prepare_run, run_pipeline


NEPS_PACKAGE = 'neural-pipeline-search'
NEPS_VERSION = '0.16.0'


def main():
    """Load the search configuration and run the NePS controller."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--config',
        default='hpo/examples/helix.example.yaml',
        help='YAML configuration file for the dataset and NePS search',
    )
    args = parser.parse_args()

    # The search-space and neps.run API are tied to this tested package version.
    installed_neps_version = version(NEPS_PACKAGE)
    if installed_neps_version != NEPS_VERSION:
        raise RuntimeError(
            f'This implementation requires {NEPS_PACKAGE} {NEPS_VERSION}; '
            f'found {installed_neps_version}'
        )
    import neps

    config = load_config(args.config, mode='neps')

    run_root, dataset, fingerprint = prepare_run(config)

    # Only one controller may operate on a run directory at a time.
    with (run_root / 'controller.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

        # Convert YAML search definitions into NePS parameter objects.
        space = {}
        for name, raw in config['search'].items():
            spec = dict(raw)
            kind = spec.pop('type')
            parameter_class = {
                'integer': neps.Integer,
                'float': neps.Float,
                'categorical': neps.Categorical,
            }[kind]
            space[name] = parameter_class(**spec)

        def evaluate_pipeline(pipeline_directory, **parameters):
            return run_pipeline(
                config,
                dataset,
                pipeline_directory / 'artifacts',
                fingerprint,
                parameters,
            )

        neps_root = run_root / "neps"
        if (neps_root / "pipeline_space.pkl").is_file():
            space_kwargs = {}
        else:
            space_kwargs = {"pipeline_space": space}

        optimizer_kwargs = {}
        if config.get('optimizer') is not None:
            optimizer_kwargs['optimizer'] = config['optimizer']

        neps.run(
            evaluate_pipeline= evaluate_pipeline,
            root_directory= run_root / 'neps',
            evaluations_to_spend= config['evaluations'],
            ignore_errors= config['ignore_errors'],
            **space_kwargs,
            **optimizer_kwargs,
        )


if __name__ == '__main__':
    main()
