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
    if 'search' in config:
        logger.warning(
            "You are running the pipeline in Standalone mode." 
            "Using 'fixed' parameters from settings and ignoring 'search' parameters."
            "For missing parameters the generators defaults are used."
        )

    signals.install()
    try:
        with signals.trial_running():
            loss = run_pipeline(
                config, dataset, ws_root / 'standalone_trial', fingerprint,
                keep_logs_in_trial=True,
            )
    except signals.ControllerStop as stop:
        # run_pipeline has written the failure report; AF3 outputs that are
        # complete are reused when the trial is submitted again.
        logger.error('Standalone trial stopped: %s', stop)
        raise SystemExit(signals.EXIT_CODE) from stop
    print(loss)


if __name__ == '__main__':
    main()