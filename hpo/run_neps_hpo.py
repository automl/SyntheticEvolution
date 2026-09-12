"""NePS 0.16.0 entry point. One controller, shared configuration across all RNAs."""
import argparse
import fcntl
from importlib.metadata import version
from run_shs_dssr_pipeline import load_config, prepare, run_pipeline


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='hpo/configs/helix.example.yaml')
    args = parser.parse_args()
    if version('neural-pipeline-search') != '0.16.0':
        raise RuntimeError('This implementation targets neural-pipeline-search 0.16.0')
    import neps
    config = load_config(args.config)
    root, rows, identity = prepare(config)
    with (root / 'controller.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        space = {}
        for name, raw in config['search'].items():
            spec = dict(raw)
            kind = spec.pop('type')
            cls = {'integer': neps.Integer, 'float': neps.Float, 'categorical': neps.Categorical}[kind]
            space[name] = cls(**spec)
        def evaluate_pipeline(pipeline_directory, **parameters):
            return run_pipeline(parameters, config, rows, pipeline_directory / 'artifacts', identity)
        neps.run(evaluate_pipeline=evaluate_pipeline, pipeline_space=space,
                 root_directory=root / 'neps', evaluations_to_spend=config['evaluations'],
                 optimizer=config['optimizer'], ignore_errors=False)


if __name__ == '__main__':
    main()
