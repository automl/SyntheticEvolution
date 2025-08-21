from pathlib import Path

import pandas as pd
import json
from tqdm import tqdm

def read_json(path):
    with open(path, 'r') as f:
        return json.load(f)

def get_msa_info(json_path):
    data = read_json(json_path)
    name = data['name']
    samples = []
    for s in data['sequences']:
        for k, v in s.items():
            seq = s[k]['sequence']
            sample = {
                    'sequence': seq,
                    'has_msa': not s[k]['unpairedMsa'] == f'>query\n{seq}\n',
                    'chain_id': s[k]['id'],
                    'mol_type': k,
                    'pdb_id': name[:4],
                    'json_id': name,
            }
            samples.append(sample)
    return pd.DataFrame(samples)

def main(json_dir):
    json_dir = Path(json_dir)
    # json_dir.mkdir(exist_ok=True, parents=True)

    data = []

    print(f'Start collecting json data from {json_dir.resolve()}')

    jsons = list(Path(json_dir).rglob('*_data.json'))

    print('done.')

    print('Found', len(jsons), 'jsons.')
    print('Start Processing')

    pbar = tqdm(total=len(jsons))

    for p in jsons:
        # print('Extract MSA info from:', p.stem)
        try:
            df = get_msa_info(p)
            data.append(df)
            pbar.update(1)
        except Exception as e:
            print(e)
            pbar.update(1)
            continue
    results = pd.concat(data)
    return results


if __name__ == '__main__':
    # json_path = Path("/work/dlclarge1/runget-af3_evaluation/data/nonx_data_af3_preds_for_iris/7ojn_s1/7ojn_s1_data.json")
    # results = get_msa_info(json_path)
    # for i, row in results.iterrows():
    #     print(row)
    json_dir = "/work/dlclarge1/runget-af3_evaluation/data/nonx_data_af3_preds_for_iris"
    results = main(json_dir)

    # print(results)
    print(results.groupby(['has_msa', 'mol_type']).count())

    results.to_csv("msa_info_nonx_data_af3_preds_for_iris.csv", index=False)  #  sep="\t", index=False)

