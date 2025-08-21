import subprocess
import json

import pandas as pd

from typing import Optional, Dict, Tuple
from pathlib import Path


def get_json(cif_file: str, dssr_dir: str, out_dir: Optional[str]=None):
    """
    Run DSSR (x3dna) on a provided cif file.
    This will create a json file of the output of DSSR.

    params:
    cif_file <string>: path to the cif file.
    dssr_dir <string>: directory of the DSSR binaries.
    out_dir (optional) <string>: output directory for the json file. If not provided, output goes to dssr_dir.
    """
    cif_file = Path(cif_file)

    print(f"Processing {cif_file.resolve()}")

    if out_dir is None:
        out_dir = ''
    else:
        Path(out_dir).mkdir(exist_ok=True, parents=True)

    # json_path = Path(out_dir, f"{cif_file.stem}.json")
    json_path = Path('/'.join(str(cif_file).split('/')[:-1]), cif_file.name.replace('.cif', '_dssr.json'))

    if json_path.is_file():
        print(f"File {json_path.resolve()} already exists. Skipping...")
        return

    args = [
                "./x3dna-dssr",
                f"-i={cif_file.resolve()}",
                f"-o={json_path.resolve()}",
                "--json",
                # "--pair-only",
            ]
    subprocess.call(args, cwd=dssr_dir)

def process_files_to_json(cif_dir: str, dssr_dir: str, out_dir: Optional[str]=None):
    """
    Parse all cif files in a directory and run DSSR on them to obtain json outputs.

    params:
    cif_dir <string>: directory containing cif files.
    dssr_dir <string>: directory of the DSSR binaries.
    out_dir (optional) <string>: output directory for the json files. If not provided, output goes to dssr_dir.
    """
    file_list = Path(cif_dir).rglob("*.cif")
    
    for f in file_list:
        get_json(f, dssr_dir=dssr_dir, out_dir=out_dir)


def parse_single_dssr_json(json_file: str) -> Tuple[str, str, Dict[str, Dict]]:
    """
    Get sequence and secondary structure information from a single DSSR json file.

    params:
    json_file <string>: path to the json file.
    """
    with open(json_file, 'r') as file:
        data = json.load(file)

    results = {}
    
    if 'chains' in data:
        for chain_key, chain_data in data['chains'].items():
            sequence = chain_data['bseq'].replace('&', '')  # do we want to replace t with u?
            dot_bracket = chain_data['sstr'].replace('&', '')
            # TODO: parse more of the json file to get more information
            
            results[chain_key] = {
                'sequence': sequence,
                'dot_bracket': dot_bracket
            }
            
    else:
        raise UserWarning("No 'chains' key found in json file")

    return results


def parse_dssr_json(json_dir: Optional[str] = None) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Get sequence and secondary structure information from DSSR json files.

    params:
    json_dir <string>: directory containing json files.
    """
    result = {}


    if json_dir is None:
        json_dir = '.'
    

    jsons = list(Path(json_dir).rglob("*_dssr.json"))
    
    for file_path in jsons:
        print('Extracting data from', file_path.resolve())
        try:
            with open(file_path, 'r') as file:
                data = json.load(file)
    
            file_stem = file_path.stem.upper()
            
            result[file_stem] = {}
    
            if 'chains' in data:
                for chain_key, chain_data in data['chains'].items():
                    sequence = chain_data['bseq'].replace('&', '')
                    dot_bracket = chain_data['sstr'].replace('&', '')
    
                    #TODO: parse base pairs to get exact information about pairs
                    # id_map, pair_df = get_pairs(data)
    
                    # pos1 = [(id_map[p][0].upper(), int(id_map[p][2])) 
                    #         if not '^' in p.split('.')[1][1:] 
                    #         else (id_map[p][0].upper(), int(id_map[p][2])) 
                    #         for p in pair_df['nt1']]
                    # 
                    # pos2 = [(id_map[p][0].upper(), int(id_map[p][2])) 
                    #         if not '^' in p.split('.')[1][1:] 
                    #         else (id_map[p][0].upper(), int(id_map[p][2])) 
                    #         for p in pair_df['nt2']]
                    # from collections import Counter
                    # c1 = Counter([p[1] for p in pos1])
                    # c2 = Counter([p[1] for p in pos2])
#     
                    # multiplets = sorted([p for p in c1.keys() if c1[p] > 1] + [p for p in c2.keys() if c2[p] > 1])
                    # 
                    # pairs = [[p1[1]-1, p2[1]-1] for p1, p2 in zip(pos1, pos2)]
                   
                    assert len(sequence) == len(dot_bracket), f"Sequence and dot-bracket length mismatch for {file_stem} chain {chain_key}"
                   
                    result[file_stem][chain_key] = {
                        'sequence': sequence,
                        'dot_bracket': dot_bracket,
                        # 'base_pairs': pairs,
                        # 'multiplet_positions': multiplets,
                        # 'pos1': pos1,
                        # 'pos2': pos2,
                    }
            else:
                raise UserWarning(f"No 'chains' key found in {file_stem}")   
        except Exception as e:
            RED = "\033[31m"
            RESET = "\033[0m"
            print(f"{RED}Error processing {file_path}: {e}{RESET}")
            continue
    return result


def get_pairs(json_data):
    if 'nts' in json_data.keys():
        id_2_index_map = {nts['nt_id']: (nts['nt_name'], nts['nt_code'], nts['index_chain']) for nts in json_data['nts']}
    else:
        id_2_index_map = None
    if 'pairs' in json_data.keys():
        df = pd.DataFrame(json_data['pairs'])
    else:
        df = None
    return id_2_index_map, df




if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Evaluate secondary structure prediction')

    parser.add_argument('--cif_dir',
                        '-i',
                        type=str,
                        # nargs="+",
                        required=True,
                        help='Path to directory containing cif files',
                        )
    parser.add_argument('--dssr_dir',
                        '-d',
                        type=str,
                        required=True,
                        help='Path to DSSR directory',
                        )
    parser.add_argument('--out_dir',
                        '-o',
                        type=str,
                        required=False,
                        help='Output directory for json files',
                        )

    args = parser.parse_args()

    process_files_to_json(args.cif_dir, dssr_dir=args.dssr_dir, out_dir=args.out_dir)
    result = parse_dssr_json(args.out_dir)
    print(result)

    print()

    json_path = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/done/a2021_s1/dssr/8ssw_af_dssr.json'

    result = parse_single_dssr_json(json_path)
    print(result)

    
