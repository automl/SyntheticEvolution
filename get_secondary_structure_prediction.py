import argparse
import random
import json
import logging
from pathlib import Path
from typing import Dict, List, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import subprocess

import RnaBench

from RnaBench.lib.visualization import RnaVisualizer
from RnaBench.lib.rna_folding_algorithms.rnafold import RNAFold
# from RnaBench.lib.rna_folding_algorithms.contrafold import ContraFold
# from RnaBench.lib.rna_folding_algorithms.ipknot import IpKnot
# from RnaBench.lib.rna_folding_algorithms.pkiss import PKiss
# from RnaBench.lib.rna_folding_algorithms.linearfold import LinearFoldC, LinearFoldV
# from RnaBench.lib.rna_folding_algorithms.rnastructure import Fold
from RnaBench.lib.rna_folding_algorithms.DL.spotrna import SpotRna
# from RnaBench.lib.rna_folding_algorithms.DL.mxfold2 import MxFold2
# from RnaBench.lib.rna_folding_algorithms.DL.RNAformer.rnaformer import RNAformer
# from RnaBench.lib.rna_folding_algorithms.DL.ProbTransformer.probtransformer import ProbabilisticTransformer
# from RnaBench.lib.rna_folding_algorithms.DL.ufold.ufold import UFold
from RnaBench.lib.utils import pairs2db, pairs2mat

def convert_pred_pairs(pred_pairs: Any) -> Dict[int, int]:
    """
    Convert predicted base pairs into a dictionary.
    If `pred_pairs` is already a dict, it is returned unchanged.
    If it's a list of lists/tuples (each with at least two elements),
    then the first two elements are interpreted as the paired positions.
    """
    if isinstance(pred_pairs, dict):
        return pred_pairs
    elif isinstance(pred_pairs, list):
        new_dict = {}
        for item in pred_pairs:
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                i, j = int(item[0]), int(item[1])
                new_dict[i] = j
                new_dict[j] = i
        return new_dict
    else:
        raise ValueError("Unknown format for predicted pairs: {}".format(type(pred_pairs)))


def get_seq_from_json(pdb_id: str) -> str:
    pdb_id = pdb_id.lower()[:4]
    try:
        with open(f'/home/fred/current_projects/github/RnaBench/single_chain_data_files_for_custom_msa_runs/{pdb_id}_s1_data.json', 'r') as f:
            data = json.load(f)
    except Exception as e:
        try:
            with open(f'/home/fred/current_projects/github/RnaBench/a2021_rnaonly_data_files/{pdb_id}_s1_data.json', 'rb') as f:
                data = json.load(f)
        except Exception as e:
            try:
                with open(f'/home/fred/current_projects/github/RnaBench/a2021_rnaonly_data_files/{pdb_id}_s1_data.json', 'rb') as f:
                    data = json.load(f)
            except Exception as e:
                raise FileNotFoundError(f"No JSON file found with pdb_id {pdb_id}.")
    
    for key, value in data.items():
        if key == "sequences":
            for chain in value:
                if "rna" in chain:
                    return chain["rna"]["sequence"]
    raise ValueError(f"Sequence for PDB ID {pdb_id} not found in JSON data.")


def get_structure(predictor: str, pdb_id: str = None, sequence: str = None) -> Dict[int, int]:
    # Use provided structure predictor if available.
    pdb_id = pdb_id.lower()[:4]
    seq = sequence.upper() if sequence else get_seq_from_json(pdb_id).upper()
    if predictor == "rnafold":
        # seq = get_seq_from_json(pdb_id).upper()
        model = RNAFold()
        pred_pairs, _ = model(seq)
    # elif predictor == "contrafold":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = ContraFold()
    #     pred_pairs = model(seq)
    # elif predictor == "ipknot":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = IpKnot()
    #     pred_pairs = model(seq)
    # elif predictor == "pkiss":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = PKiss()
    #     pred_pairs = model(seq)
    # elif predictor == "linearfoldc":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = LinearFoldC()
    #     pred_pairs = model(seq)
    # elif predictor == "linearfoldv":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = LinearFoldV()
    #     pred_pairs = model(seq)
    elif predictor == "spotrna":
        # seq = get_seq_from_json(pdb_id).upper()
        model = SpotRna()
        pred_pairs = model(seq)
    elif predictor == "rnaformer":
        # seq = get_seq_from_json(pdb_id).upper()
        print(seq)
        known = pd.read_pickle('synthetic_msa_selection_max_rna_len_200_min_rna_len_15_max_protein_len_200_rnaformer2025_pdb_id2pairs_mapping.pkl')

        if seq in known['sequence'].unique():
            # print('Known sequence, using cached RNAformer prediction.')
            # print(known[known['sequence'] == seq]['pairs'].values)
            pred_pairs = known[known['sequence'] == seq]['pairs'].values[0]
            # print(pred_pairs)
            pred_pairs = sorted(pred_pairs, key=lambda x: x[0])
        else:
            print('Unknown sequence. Start RNAformer prediction.')
            sample = {'Id': pdb_id.upper(), 'sequence': seq}
            df = pd.DataFrame([sample])
            
            import pickle
            
            with open(Path(f'RNAformer/datasets/{pdb_id.upper()}.pkl'), 'wb') as f:
                pickle.dump(df, f)
            
            subprocess.run(['./RNAformer/run_cpu.sh', 'infer_RNAformer.py', '-c', '1', '-m', 'models/af3_like_finetune', '-p', f'datasets/{pdb_id.upper()}.pkl'])
            
            pred = pd.read_pickle(f'RNAformer/datasets/{pdb_id.upper()}_processed.pkl')
    
            pred_pairs = sorted(pred['pairs'].values[0], key=lambda x: x[0])
    # elif predictor == "mxfold2":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = MxFold2()
    #     pred_pairs = model(seq)
    # elif predictor == "rnaformer2025rnaprotein":
    #     df = pd.read_pickle('synthetic_msa_nature_methods/from_helix/pickles/synthetic_msa_rnaformer2025_pdb_id2pairs_mapping.pkl')
    #     pred_pairs = df[df['pdb_id'].str.lower() == pdb_id]['pairs'].values[0]
    #     pred_pairs = [[p1, p2, 0] for p1, p2 in pred_pairs]
    # elif predictor == "rnaformer2025rnaonlyA2021":
    #     df = pd.read_pickle('synthetic_msa_nature_methods/from_helix/pickles/a2021_rnaformer2025_pdb_id2pairs_mapping.pkl')
    #     pred_pairs = df[df['pdb_id'].str.lower() == pdb_id]['pairs'].values[0]
    #     pred_pairs = [[p1, p2, 0] for p1, p2 in pred_pairs]
    # elif predictor == "rnaformer2025rnaonlyB2021":
    #     df = pd.read_pickle('synthetic_msa_nature_methods/from_helix/pickles/b2021_rnaformer2025_pdb_id2pairs_mapping.pkl')
    #     pred_pairs = df[df['pdb_id'].str.lower() == pdb_id]['pairs'].values[0]
    #     pred_pairs = [[p1, p2, 0] for p1, p2 in pred_pairs]
    # elif predictor == 'dssrRnaprotein':
    #     pred_pairs = None
    #     with open('synthetic_msa_nature_methods/from_helix/jsons/synthetic_msa_sequences.json', 'r') as f:
    #         data = json.load(f)
    #     # print(data)
    #     # Iterate through the sequences list and return the matching record.
    #     for record in data.get("sequences", []):
    #         if record.get("pdb_id", "").lower() == pdb_id:
    #             pred_pairs = record.get("base-pairs", [])
    #     if pred_pairs is None:
    #         raise ValueError(f"Base pairs not found for PDB ID: {pdb_id}")
    #     pred_pairs = [[p1, p2, 0] for p1, p2 in pred_pairs] 
    # elif predictor == 'dssrRnaonlyA2021':
    #     pred_pairs = None
    #     with open('synthetic_msa_nature_methods/from_helix/jsons/a2021_secondary_structure_data_dssr.json', 'r') as f:
    #         data = json.load(f)
    #     # print(data)
    #     # Iterate through the sequences list and return the matching record.
    #     for record in data.get("sequences", []):
    #         if record.get("pdb_id", "").lower() == pdb_id:
    #             pred_pairs = record.get("base-pairs", [])
    #     if pred_pairs is None:
    #         raise ValueError(f"Base pairs not found for PDB ID: {pdb_id}")
    #     pred_pairs = [[p1, p2, 0] for p1, p2 in pred_pairs] 
    # elif predictor == 'dssrRnaonlyB2021':
    #     pred_pairs = None
    #     with open('synthetic_msa_nature_methods/from_helix/jsons/b2021_secondary_structure_data_dssr.json', 'r') as f:
    #         data = json.load(f)
    #     # print(data)
    #     # Iterate through the sequences list and return the matching record.
    #     for record in data.get("sequences", []):
    #         if record.get("pdb_id", "").lower() == pdb_id:
    #             pred_pairs = record.get("base-pairs", [])
    #     if pred_pairs is None:
    #         raise ValueError(f"Base pairs not found for PDB ID: {pdb_id}")
    #     pred_pairs = [[p1, p2, 0] for p1, p2 in pred_pairs] 
    # elif predictor == "rnaformer64":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = RNAformer(dim=64, cycling=False)
    #     pred_pairs = model(seq)
    # elif predictor == "rnaformer128":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = RNAformer(dim=128, cycling=False)
    #     pred_pairs = model(seq)
    # elif predictor == "rnaformer256":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = RNAformer(dim=256, cycling=False)
    #     pred_pairs = model(seq)
    # elif predictor == "rnaformer256cycling":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = RNAformer(dim=256, cycling=True)
    #     pred_pairs = model(seq)
    # elif predictor == "probtransformer":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = ProbabilisticTransformer()
    #     pred_pairs = model(seq)
    # elif predictor == "ufold_train":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = UFold(model='ufold_train.pt', nc=False)
    #     pred_pairs = model(seq)
    # elif predictor == "ufold_train_nc":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = UFold(model='ufold_train.pt', nc=True)
    #     pred_pairs = model(seq)
    # elif predictor == "ufold_train_all":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = UFold(model='ufold_train_alldata.pt', nc=False)
    #     pred_pairs = model(seq)
    # elif predictor == "ufold_train_all_nc":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = UFold(model='ufold_train_alldata.pt', nc=True)
    #     pred_pairs = model(seq)
    # elif predictor == "ufold_pdb":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = UFold(model='ufold_train_pdbfinetune.pt', nc=False)
    #     pred_pairs = model(seq)
    # elif predictor == "ufold_pdb_nc":
    #     seq = get_seq_from_json(pdb_id).upper()
    #     model = UFold(model='ufold_train_pdbfinetune.pt', nc=True)
    #     pred_pairs = model(seq)
    else:
        raise ValueError(f"Unknown structure predictor: {predictor}")
    # Convert the predicted pairs to a dictionary, if necessary.
    pred_pairs = convert_pred_pairs(pred_pairs)
    logging.info("%s predicted %d unique base pairs.",
                 predictor,
                 len({(min(i, j), max(i, j)) for i, j in pred_pairs.items() if i < j}))
    return pred_pairs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyse MSA data.")
    parser.add_argument('--pdb_id', type=str, required=False, help='PDB ID of the RNA structure.')
    parser.add_argument('--sequence', type=str, required=False, default=None, help='Input sequence. If provided, it will be used instead of fetching from PDB ID.')
    parser.add_argument('--predictor', type=str, default='rnafold',
                        choices=['rnafold', 'spotrna', 'rnaformer'],
                        help='Structure prediction algorithm to use.')
    args = parser.parse_args()

    if args.sequence is None and args.pdb_id is None:
        raise ValueError("Either --pdb_id or --sequence must be provided.")
    elif args.sequence is not None and args.pdb_id is not None:
        logging.warning("Both --pdb_id and --sequence provided. Using --sequence.")
        pdb_id = args.pdb_id
    elif args.sequence is not None and args.pdb_id is None:
        pdb_id = 'custom_prediction'
        logging.info("Using custom sequence prediction with ID: %s", pdb_id)
    else:
        pdb_id = args.pdb_id.upper()[:4]
        logging.info("Using PDB ID: %s", pdb_id)
    # Get the predicted base pairs for the given PDB ID.
    pred_pairs = get_structure(args.predictor, pdb_id, args.sequence)
    
    # Print the predicted base pairs.
    # print(f"Predicted base pairs for {args.pdb_id} using {args.predictor}:")
    # for i, j in sorted(pred_pairs.items()):
    #     print(f"{i} - {j}")
    
    seq = get_seq_from_json(pdb_id) if args.sequence is None else args.sequence
    logging.info(f"Sequence for {pdb_id}: {seq}")
    mat = pairs2mat(sorted(pred_pairs.items()), length=len(seq), no_pk=True, symmetric=True)

    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    
    plotting_dir = Path('secondary_structure_predictions')
    plotting_dir.mkdir(exist_ok=True, parents=True)

    # cmap = ListedColormap(["white", "#9bd0ff"])  # light blue as other plots
    # cmap = ListedColormap(["white", "#D81B60"])  # white + nice magenta
    cmap = ListedColormap(["white", "#b20000"])  # white + red from dot plot
    #cmap='gray'


    plt.imshow(mat, cmap=cmap, interpolation='nearest')
    # plt.colorbar()
    # plt.title(f"Base Pair Matrix for {pdb_id} using {args.predictor}")
    # plt.xlabel("Position")
    # plt.ylabel("Position")
    # plt.show()
    plt.tight_layout()
    plt.savefig(f"{plotting_dir}/{pdb_id}_{args.predictor}_white_red.png", dpi=300)
    plt.close()

    vis = RnaVisualizer()
    pred_pairs = sorted(set((min(p1, p2), max(p1, p2)) for p1, p2 in pred_pairs.items()), key=lambda x: x[0])
    vis.visualize_rna(pred_pairs, seq, f'{pdb_id}_{args.predictor}', algorithm=args.predictor, plot_dir=plotting_dir, plot_type='radiate', resolution='8.0')
    vis.visualize_rna(pred_pairs, seq, f'{pdb_id}_{args.predictor}', algorithm=args.predictor, plot_dir=plotting_dir, plot_type='line', resolution='8.0')

