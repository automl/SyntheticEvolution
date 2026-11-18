import argparse
import json
import logging
logger = logging.getLogger(__name__)

from hpo.pipeline import signals
from hpo.pipeline.run import load_config, create_dataset, prepare_run
from hpo.pipeline.trial import run_pipeline

def main():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='hpo/examples/helix.example.yaml')
    parser.add_argument('--validate-only', action='store_true',
                        help='Validate and print the input dataset, then exit without running the pipeline')
    args = parser.parse_args()

    config = load_config(args.config, mode='standalone')

    if args.validate_only:
        print(json.dumps(create_dataset(config['input']), indent=2))
        return

    ws_root, dataset, fingerprint = prepare_run(config)

    # Without `trials`, one trial of the `fixed` parameters that keeps its logs.
    # Each listed trial runs in trials/<name>/ with its own parameters and
    # seeds; its logs are not kept (many files), a failure report holds them.
    if 'trials' in config:
        trials = [
            (
                trial['name'], ws_root / 'trials' / trial['name'],
                {**config, **{key: trial[key] for key in ('shs_seed', 'af3_seed') if key in trial}},
                trial.get('parameters'), False,
            )
            for trial in config['trials']
        ]
    else:
        trials = [('standalone_trial', ws_root / 'standalone_trial', config, None, True)]

    # All trials share this controller's AF3 worker. A failed trial does not
    # stop the others; submitting again reruns only unfinished trials.
    signals.install()
    failed = []
    for name, directory, trial_config, parameters, keep_logs in trials:
        if signals.requested():
            logger.error('Stopped before trial %s: %s received', name, signals.requested())
            raise SystemExit(signals.EXIT_CODE)
        try:
            with signals.trial_running():
                loss = run_pipeline(
                    trial_config, dataset, directory, fingerprint, parameters,
                    keep_logs_in_trial=keep_logs,
                )
        except signals.ControllerStop as stop:
            # run_pipeline has written the failure report; AF3 outputs that are
            # complete are reused when the run is submitted again.
            logger.error('Standalone trial %s stopped: %s', name, stop)
            raise SystemExit(signals.EXIT_CODE) from stop
        except Exception:
            logger.exception('Standalone trial %s failed', name)
            failed.append(name)
            continue
        print(name, loss, flush=True)
    if failed:
        raise SystemExit(f"Failed trials: {', '.join(failed)}. Submit again to rerun them.")


if __name__ == '__main__':
    main()
