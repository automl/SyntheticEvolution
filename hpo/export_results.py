import argparse
import csv
import json
from hpo.pipeline.trial import ROOT, load_config, get_workspace_path, write_json


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            'Export completed NePS trial summaries for a configured run.\n\n'
            'The following files are written to hpo/reports/<run_name>/:\n'
            '  scores.csv        Per-RNA scores from every completed trial\n'
            '  best.json         Parameters and scores from the best trial\n'
            '  run_metadata.json Fingerprint, dataset, and run configuration'
        ),
    )
    parser.add_argument(
        '--config',
        required=True,
        help='YAML configuration file identifying the run workspace to export',
    )
    args = parser.parse_args()

    ### Make folder 
    config = load_config(args.config)
    root = get_workspace_path(config)
    destination = ROOT / 'hpo/reports' / config['run_name']
    destination.mkdir(parents=True, exist_ok=True)

    ### read results
    results = [
        (path, json.loads(path.read_text()))
        for path in (root / 'neps').rglob('artifacts/trial_result.json')
    ]
    if not results:
        raise ValueError('No completed NePS trials to export')

    ############################## scores.csv ##############################
    scores_path = destination / 'scores.csv'
    score_fields = ['trial', 'id', 'loss', 'f1', 'tp', 'fp', 'fn']
    with scores_path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=score_fields)
        writer.writeheader()
        for path, result in results:
            for score in result['scores']:
                writer.writerow({
                    'trial': str(path.relative_to(root)),
                    **{key: score[key] for key in score_fields[1:]},
                })

    ############################## best.json ##############################
    best_path, best_result = min(results, key=lambda item: item[1]['loss'])
    write_json(destination / 'best.json', {
        'source': str(best_path),
        **best_result,
    })

    ########################## run_metadata.json ##########################
    run_metadata = json.loads((root / 'run_metadata.json').read_text())
    write_json(destination / 'run_metadata.json', run_metadata)
    print(destination)


if __name__ == '__main__':
    main()
