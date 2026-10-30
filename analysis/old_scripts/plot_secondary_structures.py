import pandas as pd
from pathlib import Path
import json
from typing import Dict, Optional
import collections

from matplotlib import pyplot as plt
import numpy as np

import ast

import pickle
from Bio.PDB import MMCIFParser

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


def parse_dssr_json(json_path: str) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Get sequence and secondary structure information from DSSR json files.

    params:
    json_dir <string>: directory containing json files.
    """
    result = {}

    json_path = Path(json_path)
    
    print('Extracting data from', json_path.resolve())
    try:
        with open(json_path, 'r') as file:
            data = json.load(file)

        file_stem = json_path.stem.upper()
        
        result[file_stem] = {}

        if 'chains' in data:
            for chain_key, chain_data in data['chains'].items():
                sequence = chain_data['bseq'].replace('&', '')
                dot_bracket = chain_data['sstr'].replace('&', '')

                #TODO: parse base pairs to get exact information about pairs
                id_map, pair_df = get_pairs(data)
                pos1 = [(id_map[p][0].upper(), int(id_map[p][2])) 
                        if not '^' in p.split('.')[1][1:] 
                        else (id_map[p][0].upper(), int(id_map[p][2])) 
                        for p in pair_df['nt1']]
                
                pos2 = [(id_map[p][0].upper(), int(id_map[p][2])) 
                        if not '^' in p.split('.')[1][1:] 
                        else (id_map[p][0].upper(), int(id_map[p][2])) 
                        for p in pair_df['nt2']]
                from collections import Counter
                c1 = Counter([p[1] for p in pos1])
                c2 = Counter([p[1] for p in pos2])
#  
                multiplets = sorted([p for p in c1.keys() if c1[p] > 1] + [p for p in c2.keys() if c2[p] > 1])
                
                pairs = [[p1[1]-1, p2[1]-1] for p1, p2 in zip(pos1, pos2)]
               
                assert len(sequence) == len(dot_bracket), f"Sequence and dot-bracket length mismatch for {file_stem} chain {chain_key}"
               
                result[file_stem][chain_key] = {
                    'sequence': sequence,
                    'dot_bracket': dot_bracket,
                    'base_pairs': pairs,
                    # 'multiplet_positions': multiplets,
                    # 'pos1': pos1,
                    # 'pos2': pos2,
                }
        else:
            raise UserWarning(f"No 'chains' key found in {file_stem}")   
    except Exception as e:
        RED = "\033[31m"
        RESET = "\033[0m"
        print(f"{RED}Error processing {json_path}: {e}{RESET}")
    return result

def db2pairs(structure, start_index=0):
    """
    Converts dot-bracket string into a list of pairs.

    Input:
      structure <string>: A sequence in dot-bracket format.
      start_index <int>: Starting index of first nucleotide (default zero-indexing).

    Returns:
      pairs <list>: A list of tuples of (index1, index2, pk_level).

    """
    level_stacks = collections.defaultdict(list)
    closing_partners = {')': '(', ']': '[', '}': '{', '>': '<'}
    levels = {')': 0, ']': 1, '}': 2, '>': 3}

    pairs = []

    for i, sym in enumerate(structure, start_index):
        if sym == '.':
            continue
        # high order pks are alphabetical characters
        if sym.isalpha():
            if sym.isupper():
                level_stacks[sym].append(i)
            else:
                try:  # in case we have invalid preditions, we continue with next bracket
                    op = level_stacks[sym.upper()].pop()
                    pairs.append((op, i,
                                  ord(sym.upper()) - 61))  # use asci code if letter is used to asign PKs, start with level 4 (A has asci code 65)
                except:
                    continue
        else:
            if sym in closing_partners.values():
                level_stacks[sym].append(i)
            else:
                try:  # in case we have invalid preditions, we continue with next bracket
                    op = level_stacks[closing_partners[sym]].pop()
                    pairs.append([op, i, levels[sym]])
                except:
                    continue
    return sorted(pairs, key=lambda x: x[0])

def pairs2mat(pairs, length, symmetric=True, no_pk=False):
    """
    Convert list of pairs to matrix representation of structure.
    """
    # print(pairs)
    mat = np.zeros((length, length))
    if no_pk:
        for p1, p2 in pairs:
            mat[p1, p2] = 1
            if symmetric:
                mat[p2, p1] = 1
        return mat
    for p1, p2, _ in pairs:
        mat[p1, p2] = 1
        if symmetric:
            mat[p2, p1] = 1
    return mat

def get_sequence_from_json(data_dir: Path, pdb_id: str, chain_id: str) -> Optional[str]:
    """
    Get sequence from json file in data_dir for given pdb_id.
    """
    json_path = data_dir / f"{pdb_id.lower()}_s1_data.json"
    if not json_path.is_file():
        json_path = data_dir / f"{pdb_id.upper()}_data.json"
    if not json_path.is_file():
        print(f"JSON file {json_path} does not exist.")
        return None

    try:
        with open(json_path, 'r') as file:
            data = json.load(file)
        for key, value in data.items():
            # print(key)
            if key == "sequences":
                print('Found sequences')
                for chain in value:
                    # print(chain, 'rna' in chain)
                    if "rna" in chain:
                        print(f"Found RNA chain in {json_path}")
                        if chain["rna"]["id"] == chain_id:
                            print(f"Chain {chain_id} found in {json_path}")
                            return chain["rna"]["sequence"]
                        else:
                            print(f"Chain {chain_id} not found in {json_path}")
                            return None
                    else:
                        print(f"No RNA sequence found for chain {chain_id} in {json_path}")
                        return None
            
        print(f"No sequences found in {json_path}")
        return None
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return None


def adjust_base_pairs(bp, sequence):
    """
    Adjust base pairs to match the sequence length.
    """
    max_index = max([p[1] for p in bp])
    if len(sequence) <= max_index:
        print(f"Adjusting base pairs for sequence length {len(sequence)} and max index {max_index}.")
        return [[p[0], p[1]] for p in bp if p[1] < len(sequence)]
    else:
        print(f"No adjustment needed for sequence length {len(sequence)} and max index {max_index}.")
        return bp


def choose_sequence(bp, sequences):
    """
    Choose the best sequence based on base pairs and available sequences.
    """
    print('Selecting best sequence...')
    max_index = max([p[1] for p in bp])
    print('Max base pair index:', max_index)
    # print('sequence lengths', len(sequence), len(cif_sequence), len(json_sequence))

    longest = max([s for s in sequences if s is not None], key=len)
    if all([s is None for s in sequences]):
        return 'A' * (max_index + 1)
    
    return longest
    

if __name__ == "__main__":
    from collections import defaultdict
    write_data = True
    multichain = True
    pdb_ids = [str(p.stem)[:4].upper() for p in Path('evaluation/dssr/gt').glob('*.json')]

    base_pairs = collections.defaultdict(list)

    rnaformer_preds = pd.read_pickle('data/rnaformer_predictions.pkl')
    print(rnaformer_preds.columns)

    data_dir = Path('data', 'datafiles', 'alphafold3')

    plotting_dir = Path('plots', 'secondary_structure')
    Path(plotting_dir).mkdir(parents=True, exist_ok=True)

    pickle_dir = Path('results')
    pickle_dir.mkdir(parents=True, exist_ok=True)
    pickle_path = pickle_dir / f'secondary_structure_information_multichain_{multichain}.pkl'
    
    for pdb_id in pdb_ids:
        if not pdb_id.upper() == '1SJ4' and not pdb_id.upper() == '4P8Z':
            continue
        jsons = [
            ('RNAformer_dssr', f'evaluation/dssr/rnaformer/{pdb_id.upper()}_rnaformer_dssr.json'),
            ('AlphaFold_dssr', f'evaluation/dssr/alphafold3/{pdb_id.upper()}_alphafold_dssr.json'),
            ('GroundTruth', f'evaluation/dssr/gt/{pdb_id.upper()}_gt_dssr.json'),
            ('SpotRNA_dssr', f'evaluation/dssr/spotrna/{pdb_id.upper()}_spotrna_dssr.json'),
            ('RNAfold_dssr', f'evaluation/dssr/rnafold/{pdb_id.upper()}_rnafold_dssr.json'),
            ('DSSR-3D', f'evaluation/dssr/dssr/{pdb_id.upper()}_dssr_dssr.json'),
        ]
        sequences = defaultdict(list)
        for model, json_path in jsons:
            try:
                res = parse_dssr_json(json_path)
                if not res:
                    print(f"No data found in {json_path}, skipping...")
                    continue
                for key, value in res.items():
                    for chain in value:
                        seq = value[chain]['sequence']
                        seq = ''.join([s for s in seq if not s.islower()])
                        sequence = get_sequence_from_json(data_dir, pdb_id, str(chain).split('_')[-1])
                        print(f"Sequence in json not found for {pdb_id} {chain}.")
                        print('Try reading from cif file.')
                        
                        parser = MMCIFParser()
                        name = model.split('_')[0].lower() if not model == 'GroundTruth' else 'gt'
                        name2 = name + '3' if model.split('_')[0].lower() == 'alphafold' else name
                        print(name)
                        print(name2)
                        cif_path = Path(json_path).parent.parent.parent / f"predictions/{name2}/{pdb_id}_{name}.cif"
                        print(cif_path)
                        if not cif_path.is_file():
                            print(f"CIF file {cif_path} does not exist. Using sequence from DSSR")
                            cif_sequence = None
                        else:
                            structure = parser.get_structure(pdb_id, cif_path)
                            for m in structure:
                                for chain_obj in m:
                                    print(chain_obj.id, str(chain).split('_')[-1])
                                    
                                    if chain_obj.id == str(chain).split('_')[-1]:
                                        # for res in chain_obj:
                                        #     print(res.get_resname())
                                        cif_sequence = ''.join([res.get_resname() for res in chain_obj if res.get_resname() in ['A', 'C', 'G', 'U']])
                                        
                                        print('Sequence found in cif file.')
                                        break
                            print(f"Sequence from cif: {cif_sequence}")
                        sequences[chain+'_'+model].append(seq)
                        sequences[chain+'_'+model].append(cif_sequence)
                        sequences[chain+'_'+model].append(sequence)
            except Exception as e:
                print(f"Error processing {json_path} during sequence selection: {e}")
                continue
        
        # First, keep longest per model+chain
        intermediate_seqs = {}
        for chain, seqs in sequences.items():
            if seqs:
                longest = max((s for s in seqs if s), key=len)
                intermediate_seqs[chain] = longest
        
        # Then, keep longest per chain ID across all models
        final_seqs = {}
        for chain_with_model, seq in intermediate_seqs.items():
            chain_id = '_'.join(chain_with_model.split('_')[:3])  # pure chain name
            if chain_id not in final_seqs or len(seq) > len(final_seqs[chain_id]):
                final_seqs[chain_id] = seq
        print(final_seqs)
        ids_with_sequence_mismatch = []
        for model, json_path in jsons:
            
            try:
                res = parse_dssr_json(json_path)
                if not res:
                    print(f"No data found in {json_path}, skipping...")
                    continue
                bp = []
                for key, value in res.items():
                    for chain in value:
                        if multichain:
                            print(f"Chain: {chain}")
                            # print(f"Sequence: {value[chain]['sequence']}")
                            # print(f"Dot-bracket: {value[chain]['dot_bracket']}")
                            print(f"Base pairs: {value[chain]['base_pairs']}")
                            bp = value[chain]['base_pairs']
                            # print(type(bp))
                            # print(bp)
                        else:
                            bp += value[chain]['base_pairs']
                        # seq = value[chain]['sequence']
                        # seq = ''.join([s for s in seq if not s.islower()])
                        # sequence = get_sequence_from_json(data_dir, pdb_id, str(chain).split('_')[-1])
                        # if sequence is None:
                        #     print(f"Sequence in json not found for {pdb_id} {chain}.")
                        #     print('Try reading from cif file.')
                        #     
                        #     parser = MMCIFParser()
                        #     name = model.split('_')[0].lower() if not model == 'GroundTruth' else 'gt'
                        #     name2 = name + '3' if model.split('_')[0].lower() == 'alphafold' else name
                        #     print(name)
                        #     print(name2)
                        #     cif_path = Path(json_path).parent.parent.parent / f"predictions/{name2}/{pdb_id}_{name}.cif"
                        #     print(cif_path)
                        #     if not cif_path.is_file():
                        #         print(f"CIF file {cif_path} does not exist. Using sequence from DSSR")
                        #         sequence = seq
                        #         cif_sequence = None
                        #     else:
                        #         structure = parser.get_structure(pdb_id, cif_path)
                        #         for m in structure:
                        #             for chain_obj in m:
                        #                 print(chain_obj.id, str(chain).split('_')[-1])
                        #                 
                        #                 if chain_obj.id == str(chain).split('_')[-1]:
                        #                     # for res in chain_obj:
                        #                     #     print(res.get_resname())
                        #                     cif_sequence = ''.join([res.get_resname() for res in chain_obj if res.get_resname() in ['A', 'C', 'G', 'U']])
                        #                     
                        #                     print('Sequence found in cif file.')
                        #                     break
                        #         print(f"Sequence from cif: {cif_sequence}")
                        #         sequence = cif_sequence
                        # if seq == sequence:
                        #     sequence = seq
                        # else:
                        #     if len(seq) >= len(sequence):
                        #         sequence = seq
                        #     else:
                        #         print(f"Sequence mismatch for {pdb_id} {chain}: {seq} vs {sequence}. Using sequence from json.")
                        #         print(f"Max pair index is: {max([p[1] for p in bp])} vs CIF: {len(cif_sequence)} vs DSSR: {len(seq)} vs JSON: {len(sequence)}")
                        sequence = final_seqs.get(chain, None)
                        if sequence is None:
                            print(f"No sequence found for {pdb_id} {chain}. Using A only sequence.")
                            sequence = 'A' * (max([p[1] for p in bp]) + 1)
                        if len(sequence) < max([p[1] for p in bp]) + 1:
                            print(f"Longest sequence found for {pdb_id} {chain} is too short. Start removing base pairs.")
                            bp_original = bp
                            bp = adjust_base_pairs(bp, sequence)
                            ids_with_sequence_mismatch.append((model, pdb_id, chain, sequence, bp_original, bp))
                            
                        
                        print(f"Using sequence: {sequence}")
                        print('Selected sequence:', sequence)
                        print('Sequence length:', len(sequence))
                        print('Max base pair index:', max([p[1] for p in bp]))
                        if multichain:
                            base_pairs[pdb_id.upper()].append((model, chain, sequence, bp))
                if not multichain:
                    base_pairs[pdb_id.upper()].append((model, chain, sequence, bp))

            except Exception as e:
                print(f"Error processing {json_path}: {e}")
                continue
    
    data = []

    # print(base_pairs)

    for pdb_id, plotting_data in base_pairs.items():
        for model, chain, sequence, pairs in plotting_data:
            # print(sequence)
            print(f"Processing {pdb_id} {chain} {model}...")
            if model == 'RNAformer_dssr':
                try:
                    print('Try getting RNAformer prediction.')
                    if multichain:
                        c = chain.split('_')[-1]
                        mask = (
                            (rnaformer_preds['pdb_id'].str.upper() == pdb_id.upper())
                            & (rnaformer_preds['chain_id'].str.upper() == c.upper())
                        )                        
                        pred_pairs = rnaformer_preds.loc[mask, 'pairs'].iloc[0]
                        seq = rnaformer_preds.loc[mask, 'sequence'].iloc[0]
                        
                        print(pdb_id, c)
                        print(pred_pairs)
                        print(seq)
                        print(len(seq), len(sequence))

                        rnaformer_pred_mat = pairs2mat(pred_pairs, len(seq), no_pk=True)
                        rnaformer_sample = {
                            'pdb_id': pdb_id,
                            'chain': chain,
                            'model': 'RNAformer_pred',
                            'sequence': seq,
                            'pairs': pred_pairs,
                        }

                    else:
                        pred_pairs = rnaformer_preds[rnaformer_preds['pdb_id'].str.upper() == pdb_id.upper()]['pairs'].values[0]
    
                        rnaformer_pred_mat = pairs2mat(pred_pairs, len(sequence), no_pk=True)
                        rnaformer_sample = {
                            'pdb_id': pdb_id,
                            'chain': chain,
                            'model': 'RNAformer_pred',
                            'sequence': sequence,
                            'pairs': pred_pairs,
                        }
                    data.append(rnaformer_sample)
                    if not Path(f"{plotting_dir}", f"{pdb_id}_{chain}_RNAformer_pred.png").is_file():
                        plt.imshow(rnaformer_pred_mat, cmap='gray', interpolation='nearest')
                        plt.tight_layout()
                        plt.savefig(f"{plotting_dir}/{pdb_id}_{chain}_RNAformer_pred.png", dpi=300)
                        plt.close()
                except Exception as e:
                    print(f"Error with RNAformer pred: {e}")
                
        # print(pairs)
            print(pairs)
            print(len(sequence))
            print(sequence)

            try:
                mat = pairs2mat(pairs, len(sequence), no_pk=True)
                sample = {
                    'pdb_id': pdb_id,
                    'chain': chain,
                    'model': model,
                    'sequence': sequence,
                    'pairs': pairs,
                }
                data.append(sample)
                
                if Path(f"{plotting_dir}", f"{pdb_id}_{chain}_{model}.png").is_file():
                    print(f"{plotting_dir}/{pdb_id}_{chain}_{model}.png exists, skipping...")
                    continue
                else:
                    print(f"Plotting {pdb_id} {chain} {model}...")
                    plt.imshow(mat, cmap='gray', interpolation='nearest')
                    # plt.title(f"{pdb_id} - {model}")
                    # plt.colorbar()
                    # plt.xticks(range(len(sequence)), sequence, rotation=90)
                    # plt.yticks(range(len(sequence)), sequence)
                    plt.tight_layout()
                    plt.savefig(f"{plotting_dir}/{pdb_id}_{chain}_{model}.png", dpi=300)
                    plt.close()
                    # plt.show()
            except Exception as e:
                print(f"Error plotting {pdb_id} {chain} {model}: {e}")
                continue
    if write_data:
        df = pd.DataFrame(data)
        with open(pickle_path, 'wb') as f:
            pickle.dump(df, f)