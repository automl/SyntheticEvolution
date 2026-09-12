"""Export compact completed NePS trial summaries; no structures or large arrays copied."""
import argparse
import csv
import json
from run_shs_dssr_pipeline import ROOT, load_config, workspace, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='hpo/configs/helix.example.yaml')
    args = parser.parse_args()
    config = load_config(args.config)
    root = workspace(config)
    results = [(path,json.loads(path.read_text())) for path in (root/'neps').rglob('artifacts/result.json')]
    if not results:
        raise ValueError('No completed NePS trials to export')
    destination = ROOT / 'hpo/reports' / config['run_name']
    destination.mkdir(parents=True, exist_ok=True)
    with (destination/'scores.csv').open('w',newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['trial','id','loss','f1','tp','fp','fn'])
        writer.writeheader()
        for path,result in results:
            for score in result['scores']:
                writer.writerow(dict(trial=str(path.relative_to(root)), **{k:score[k] for k in ['id','loss','f1','tp','fp','fn']}))
    path,best = min(results,key=lambda item:item[1]['loss'])
    write_json(destination/'best.json',dict(source=str(path),**best))
    write_json(destination/'run.json',json.loads((root/'run.json').read_text()))
    print(destination)


if __name__ == '__main__':
    main()
