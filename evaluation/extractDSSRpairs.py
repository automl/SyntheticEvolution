import torch
import numpy as np
import pandas as pd
import json
import sys
import os

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from matplotlib import pyplot as plt
from scipy import signal
from grakel import Graph
from grakel.kernels import WeisfeilerLehman, VertexHistogram, WeisfeilerLehmanOptimalAssignment, ShortestPath
# from pathlib import Path
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig

def parse_json(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    try:
        chain_id = list(data['chains'].keys())[0]
        sequence = data['chains'][chain_id]['bseq']
    except:
        sequence = None
    #
    if 'multiplets' in data.keys():
        has_multiplet = True
    else:
        has_multiplet = False
    if 'nts' in data.keys():
        id_2_index_map = {nts['nt_id']: (nts['nt_name'], nts['nt_code'], nts['index_chain']) for nts in data['nts']}
    else:
        id_2_index_map = None
    if 'pairs' in data.keys():
        df = pd.DataFrame(data['pairs'])
    else:
        df = None

    return sequence, id_2_index_map, df, has_multiplet  # return data

def get_pdb_residue_mapping(pdb_id, aligned_folder):
    # Create mapping between PDB residue numbers and sequential numbering (1-based)
    pdb_file = os.path.join(aligned_folder, f"{pdb_id}_pdb_aligned_rna.pdb")
    if not os.path.exists(pdb_file):
        print(f"Warning: No aligned PDB file found for {pdb_id}")
        return None
        
    residue_map = {}
    current_seq_num = 1  # Start from 1
    prev_res_num = None
    first_res = None
    last_res = None
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    # Handle potential formatting issues by looking at fixed columns
                    res_num_str = line[22:26].strip()
                    # Handle negative numbers
                    if res_num_str.startswith('-'):
                        res_num = -int(res_num_str[1:])  # Convert to negative number
                    else:
                        # Skip if residue number is malformed (not a number or negative sign)
                        if not (res_num_str.isdigit() or (res_num_str.startswith('-') and res_num_str[1:].isdigit())):
                            continue
                        res_num = int(res_num_str)
                    
                    # Track first and last residue numbers
                    if first_res is None:
                        first_res = res_num
                    last_res = res_num
                    
                    # Create sequential mapping if this is a new residue
                    if res_num != prev_res_num:
                        residue_map[res_num] = current_seq_num
                        current_seq_num += 1
                        prev_res_num = res_num
                        
                except ValueError as e:
                    print(f"Warning: Could not parse residue number in line: {line.strip()}")
                    continue
    
    if residue_map:
        print(f"Created mapping for {len(residue_map)} residues")
        print(f"Residue range in PDB: {first_res} -> {last_res}")
        print(f"Sequential range: 1 -> {len(residue_map)}")
            
    return residue_map

def extract_base_pairs(dssr_json_file, is_af=False, residue_map=None, target_length=None):
    with open(dssr_json_file, 'r') as f:
        data = json.load(f)
    true_pairs = []

    if "pairs" in data:
        # Get sequence length
        if "chains" in data:
            first_chain = next(iter(data["chains"].values()))
            if "bseq" in first_chain:
                sequence_length = len(first_chain["bseq"])
                if not is_af and residue_map:  # For PDB files
                    sequence_length = len(residue_map)  # Use mapping length instead
                print(f"{'AF' if is_af else 'PDB'} sequence length: {sequence_length}")

        # Extract residue numbers from nt identifiers
        pairs_info = []
        
        for pair in data["pairs"]:
            try:
                nt1_id = pair["nt1"]
                nt2_id = pair["nt2"]

                # Handle different formats
                def extract_residue_number(nt_id):
                    if '.' in nt_id:  # Format: "A.U1"
                        _, res = nt_id.split('.')
                        return int(''.join(filter(str.isdigit, res)))
                    elif '/' in nt_id:  # Format: "///D/2"
                        parts = nt_id.split('/')
                        # Get last non-empty part
                        res = next(p for p in reversed(parts) if p)
                        return int(res)
                    else:
                        return int(''.join(filter(str.isdigit, nt_id)))

                nt1 = extract_residue_number(nt1_id)
                nt2 = extract_residue_number(nt2_id)

                # Map PDB residue numbers to sequential if needed
                if not is_af and residue_map:
                    nt1 = residue_map.get(nt1, nt1)
                    nt2 = residue_map.get(nt2, nt2)

                # Only add pairs that are within sequence bounds
                if 0 < nt1 <= sequence_length and 0 < nt2 <= sequence_length:
                    true_pairs.append([nt1, nt2])
                else:
                    print(f"!Warning: Residue numbers ({nt1}, {nt2}) out of bounds for sequence length {sequence_length}")
                    
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse pair: {pair}, Error: {e}")
                continue

        print(f"Total valid pairs found: {len(true_pairs)}")

    return true_pairs

# def get_json_data(json_path):
#     NUCS = {
#         'T': 'U',
#         'P': 'U',
#         'R': 'A',  # or 'G'
#         'Y': 'C',  # or 'T'
#         'M': 'C',  # or 'A'
#         'K': 'U',  # or 'G'
#         'S': 'C',  # or 'G'
#         'W': 'U',  # or 'A'
#         'H': 'C',  # or 'A' or 'U'
#         'B': 'U',  # or 'G' or 'C'
#         'V': 'C',  # or 'G' or 'A'
#         'D': 'A',  # or 'G' or 'U'
#         'N': 'N',
#         'A': 'A',
#         'U': 'U',
#         'C': 'C',
#         'G': 'G',
#     }
#
#     json_path = Path(json_path)
#
#     sequence_from_json, id_map, pair_df, has_multiplet = parse_json(json_path)
#
#     if pair_df is None:
#         return []
#
#     sequence = sequence_from_json.replace('T', 'U')
#     sequence = sequence.replace('&', '')
#     json_stem = json_path.stem
#
#     pos1 = [(NUCS[id_map[p][0]], int(id_map[p][2])) if not '^' in p.split('.')[1][1:] else (
#     NUCS[id_map[p][0]], int(id_map[p][2])) for p in pair_df['nt1']]
#
#     pos2 = [(NUCS[id_map[p][0]], int(id_map[p][2])) if not '^' in p.split('.')[1][1:] else (
#     NUCS[id_map[p][0]], int(id_map[p][2])) for p in pair_df['nt2']]
#
#     pairs = []
#
#     sequence = sequence.upper()
#     try:
#         for i in range(len(pos1)):
#             assert str(sequence[pos1[i][1] - 1]) == str(pos1[i][0]) or str(sequence[pos1[i][1] - 1]) == 'N' or str(
#                 pos1[i][0]) == 'N', f"{i}: {pos1[i]};\n{sequence[pos1[i][1] - 1]} != {pos1[i][0]}"
#         for i in range(len(pos2)):
#             assert str(sequence[pos2[i][1] - 1]) == str(pos2[i][0]) or str(sequence[pos2[i][1] - 1]) == 'N' or str(
#                 pos2[i][0]) == 'N', f"{i}: {pos2[i]};\n{sequence[pos2[i][1] - 1]} != {pos2[i][0]}"
#         # assert all(str(sequence[pos1[i][1]-1]) == str(pos1[i][0]) for i in range(len(pos1))), f"{sequence[pos1[0][1]-1]} != {pos1[0][0]}"
#         # assert all(str(sequence[pos2[i][1]-1]) == str(pos2[i][0]) for i in range(len(pos2))), f"{sequence[pos2[0][1]-1]} != {pos2[0][0]}"
#         pairs = [[p1[1] - 1, p2[1] - 1] for p1, p2 in zip(pos1, pos2)]
#     except:
#         try:
#             for i in range(len(pos1)):
#                 assert str(sequence[pos1[i][1]]) == str(pos1[i][0]) or str(sequence[pos1[i][1]]) == 'N' or str(
#                     pos1[i][0]) == 'N', f"{i}: {pos1[i]};\n{sequence[pos1[i][1]]} != {pos1[i][0]}"
#             for i in range(len(pos2)):
#                 assert str(sequence[pos2[i][1]]) == str(pos2[i][0]) or str(sequence[pos2[i][1]]) == 'N' or str(
#                     pos2[i][0]) == 'N', f"{i}: {pos2[i]};\n{sequence[pos2[i][1]]} != {pos2[i][0]}"
#             # assert all(str(sequence[pos1[i][1]]) == str(pos1[i][0]) for i in range(len(pos1))), f"{sequence[pos1[0][1]]} != {pos1[0][0]}"
#             # assert all(str(sequence[pos2[i][1]]) == str(pos2[i][0]) for i in range(len(pos2))), f"{sequence[pos2[0][1]]} != {pos2[0][0]}"
#             pairs = [[p1[1], p2[1]] for p1, p2 in zip(pos1, pos2)]
#         except:
#             try:
#                 for i in range(len(pos1)):
#                     assert str(sequence[pos1[i][1] + 1]) == str(pos1[i][0]) or str(
#                         sequence[pos1[i][1] + 1]) == 'N' or str(
#                         pos1[i][0]) == 'N', f"{i}: {pos1[i]};\n{sequence[pos1[i][1] + 1]} != {pos1[i][0]}"
#                 for i in range(len(pos2)):
#                     assert str(sequence[pos2[i][1] + 1]) == str(pos2[i][0]) or str(
#                         sequence[pos2[i][1] + 1]) == 'N' or str(
#                         pos2[i][0]) == 'N', f"{i}: {pos2[i]};\n{sequence[pos2[i][1] + 1]} != {pos2[i][0]}"
#                 # assert all(str(sequence[pos1[i][1]+1]) == str(pos1[i][0]) for i in range(len(pos1))), f"{sequence[pos1[0][1]+1]} != {pos1[0][0]}"
#                 # assert all(str(sequence[pos2[i][1]+1]) == str(pos2[i][0]) for i in range(len(pos2))), f"{sequence[pos2[0][1]+1]} != {pos2[0][0]}"
#                 pairs = [[p1[1] + 1, p2[1] + 1] for p1, p2 in zip(pos1, pos2)]
#             except:
#                 pass
#
#     if not pairs:
#         return []
#
#     return sorted(pairs, key=lambda x: x[0]), sequence

def pairs2mat(pairs, length, symmetric=True):
    mat = np.zeros((length, length))
    # print(pairs)

    for p1, p2 in pairs:
        # Convert 1-based to 0-based indices!!!
        p1 -= 1
        p2 -= 1
        # Ensure indices are within bounds
        if p1 < 0 or p2 < 0:
            print(f"Warning: Invalid residue numbers ({p1}, {p2}), must be positive")
            continue

        if 0 <= p1 < length and 0 <= p2 < length:
            mat[p1, p2] = 1
            if symmetric:
                mat[p2, p1] = 1
        else:
            print(f"Warning: Residue numbers ({p1}, {p2}) out of bounds for sequence length {length}")

    return torch.from_numpy(mat)

# def pairs2mat(pairs, length, symmetric=True):
#     """
#     Convert list of pairs to matrix representation of structure.
#     """
#     # print(pairs)
#     mat = np.zeros((length, length))
#
#     for p1, p2 in pairs:
#         mat[p1, p2] = 1
#         if symmetric:
#             mat[p2, p1] = 1
#     return torch.from_numpy(mat)

def get_graph_kernel(kernel, n_iter=5, normalize=True):
    if kernel == 'WeisfeilerLehman':
        return WeisfeilerLehman(n_iter=n_iter,
                                normalize=normalize,
                                base_graph_kernel=VertexHistogram)
    elif kernel == 'WeisfeilerLehmanOptimalAssignment':
        return WeisfeilerLehmanOptimalAssignment(n_iter=n_iter,
                                                 normalize=normalize)
    elif kernel == 'ShortestPath':
        return ShortestPath(normalize=normalize)


def mat2graph(matrix, node_labels=None):
    if node_labels is not None:
        graph = Graph(initialization_object=matrix.astype(int),
                      node_labels=node_labels)  # TODO: Think about if we need to label the nodes differenty
    else:
        graph = Graph(initialization_object=matrix.astype(int),
                      node_labels={s: str(s) for s in
                                   range(
                                       matrix.shape[0])})  # TODO: Think about if we need to label the nodes differenty

    return graph


def evaluate_shifted(pred_a, true_a):
    kernel = np.array([[0.0, 1.0, 0.0],
                       [1.0, 1.0, 1.0],
                       [0.0, 1.0, 0.0]])
    pred_a_filtered = signal.convolve2d(pred_a, kernel, 'same')
    fn = len(torch.where((true_a - torch.Tensor(pred_a_filtered)) == 1)[0])
    pred_p = torch.sign(torch.Tensor(pred_a)).sum()
    true_p = true_a.sum()
    tp = true_p - fn
    fp = pred_p - tp
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)
    f1_score = 2 * tp / (2 * tp + fp + fn)
    return precision.item(), recall.item(), f1_score.item()


def graph_distance_score_from_matrices(pred, true, kernel, true_labels=None, pred_labels=None):
    pred_graph = mat2graph(pred, node_labels=pred_labels)
    true_graph = mat2graph(true, node_labels=true_labels)
    kernel = get_graph_kernel(kernel=kernel)
    kernel.fit_transform([true_graph])
    distance_score = kernel.transform([pred_graph])  # TODO: Check output, might be list or list of lists

    return distance_score[0][0]


def eval_pairs(sequence, pred_pairs, true_pairs, plot_matrix=False):
    sequence = sequence.upper()
    length = len(sequence)

    pred_mat = pairs2mat(pred_pairs, length)
    true_mat = pairs2mat(true_pairs, length)

    wc_pairs = ['AU', 'UA', 'CG', 'GC']
    wobble_pairs = ['GU', 'UG']

    # Initialise masks
    wc = torch.zeros_like(true_mat)
    wobble = torch.zeros_like(true_mat)
    nc = torch.zeros_like(true_mat)
    canonicals = torch.zeros_like(true_mat)
    pk = torch.zeros_like(true_mat)
    multi = torch.zeros_like(true_mat)

    # Get positions of true base pairs
    triu_true_mat = torch.triu(torch.ones_like(true_mat), diagonal=1).bool()
    filtered_true_mat = true_mat * triu_true_mat
    true_pair_positions = torch.nonzero(filtered_true_mat, as_tuple=True)

    # print(true_pair_positions)
    multi_pred = []
    for i in range(len(sequence)):
        multi_pred.append(pred_mat[i, :].sum() > 1)
    # get specific masks for each base pair type
    for i, (p1, p2) in enumerate(zip(true_pair_positions[0], true_pair_positions[1])):

        if ''.join(sequence[p1] + sequence[p2]) in wc_pairs:
            wc[p1, p2] = 1
            wc[p2, p1] = 1
            canonicals[p1, p2] = 1
            canonicals[p2, p1] = 1
        elif ''.join(sequence[p1] + sequence[p2]) in wobble_pairs:
            wobble[p1, p2] = 1
            wobble[p2, p1] = 1
            canonicals[p1, p2] = 1
            canonicals[p2, p1] = 1
        else:
            nc[p1, p2] = 1
            nc[p2, p1] = 1
            # print('Non-canonical pair:', sequence[p1], sequence[p2])

        if true_mat[p1, :].sum() > 1:
            multi[p1, p2] = 1
            multi[p2, p1] = 1
            # print('Multi pair:', true_mat[p1, :])

    if plot_matrix:
        fig, axs = plt.subplots(1, 7, figsize=(12, 4))
        axs[0].imshow(true_mat.cpu().numpy())
        axs[0].set_title("GT")
        axs[1].imshow(wc.cpu().numpy())
        axs[1].set_title("WC")
        axs[2].imshow(wobble.cpu().numpy())
        axs[2].set_title("Wobble")
        axs[3].imshow(nc.cpu().numpy())
        axs[3].set_title("NC")
        axs[4].imshow(canonicals.cpu().numpy())
        axs[4].set_title("Canonicals")
        axs[5].imshow(pk.cpu().numpy())
        axs[5].set_title("PK")
        axs[6].imshow(multi.cpu().numpy())
        axs[6].set_title("Multi")
        plt.show()

    assert (wc + wobble + nc).sum() == true_mat.sum()
    assert (canonicals + nc).sum() == true_mat.sum()

    # compute metrics
    solved = torch.equal(true_mat, pred_mat).__int__()

    num_pred_pairs = torch.sum(torch.triu(pred_mat, diagonal=1))
    num_gt_pairs = torch.sum(torch.triu(true_mat, diagonal=1))

    metrics = {'solved': int(solved),
               'num_pred_pairs': int(num_pred_pairs.cpu().item()),
               'num_gt_pairs': int(num_gt_pairs.cpu().item()),
               'length': true_mat.size(0),
               'has_pk': pk.sum() > 0,
               'has_multi': multi.sum() > 0,
               'multi_pred': any(multi_pred),
               }

    tp = torch.logical_and(pred_mat, true_mat).sum()
    tn = torch.logical_and(torch.logical_not(pred_mat), torch.logical_not(true_mat)).sum()
    fp = torch.logical_and(pred_mat, torch.logical_not(true_mat)).sum()
    fn = torch.logical_and(torch.logical_not(pred_mat), true_mat).sum()

    assert pred_mat.size().numel() == tp + tn + fp + fn

    metrics['tp'] = tp.cpu().item()
    metrics['tn'] = tn.cpu().item()
    metrics['fp'] = fp.cpu().item()
    metrics['fn'] = fn.cpu().item()

    accuracy = tp / pred_mat.size().numel()
    precision = tp / (1e-4 + tp + fp)
    recall = tp / (1e-4 + tp + fn)
    f1_score = 2 * tp / (1e-4 + (2 * tp + fp + fn))
    mcc = (tp * tn - fp * fn) / (1e-4 + torch.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))

    # print(sequence)
    true_labels = pred_labels = {i: sequence[i] for i in range(len(sequence))}

    wl = graph_distance_score_from_matrices(pred_mat.cpu().numpy(), true_mat.cpu().numpy(), 'WeisfeilerLehman')
    wl_with_seq = graph_distance_score_from_matrices(pred_mat.cpu().numpy(), true_mat.cpu().numpy(), 'WeisfeilerLehman',
                                                     true_labels=true_labels, pred_labels=pred_labels)

    prec_shift, rec_shift, f1_shift = evaluate_shifted(pred_mat.cpu(), true_mat.cpu())

    metrics['accuracy'] = accuracy.cpu().item()
    metrics['precision'] = precision.cpu().item()
    metrics['recall'] = recall.cpu().item()
    metrics['f1_score'] = f1_score.cpu().item()
    metrics['mcc'] = mcc.cpu().item()
    metrics['wl'] = wl
    metrics['wl_with_seq'] = wl_with_seq
    metrics['prec_shift'] = prec_shift
    metrics['rec_shift'] = rec_shift
    metrics['f1_shift'] = f1_shift

    for (name, mask) in [('wc', wc), ('wobble', wobble), ('nc', nc), ('canonicals', canonicals), ('pk', pk),
                         ('multi', multi)]:
        solved = torch.equal(mask, pred_mat).__int__()
        metrics[f'{name}_solved'] = int(solved)

        tp = torch.logical_and(pred_mat, mask).sum()
        tn = torch.logical_and(torch.logical_not(pred_mat), torch.logical_not(mask)).sum()
        fp = torch.logical_and(pred_mat, torch.logical_not(mask)).sum()
        fn = torch.logical_and(torch.logical_not(pred_mat), mask).sum()
        assert pred_mat.size().numel() == tp + tn + fp + fn

        accuracy = tp / pred_mat.size().numel()
        precision = tp / (1e-4 + tp + fp)
        recall = tp / (1e-4 + tp + fn)
        f1_score = 2 * tp / (1e-4 + (2 * tp + fp + fn))
        mcc = (tp * tn - fp * fn) / (1e-4 + torch.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))

        metrics[f'{name}_accuracy'] = accuracy.cpu().item()
        metrics[f'{name}_precision'] = precision.cpu().item()
        metrics[f'{name}_recall'] = recall.cpu().item()
        metrics[f'{name}_f1_score'] = f1_score.cpu().item()
        metrics[f'{name}_mcc'] = mcc.cpu().item()

        if name == 'pk':
            metrics['pk_hit'] = tp > 0
            metrics['num_pk_hits'] = tp.cpu().item()
            metrics['num_pk_gt'] = mask.sum().cpu().item()
        if name == 'multi':
            metrics['multi_hit'] = tp > 0
            metrics['num_multi_hits'] = tp.cpu().item()
            metrics['num_multi_gt'] = mask.sum().cpu().item()

    pred_pair_indices = torch.where(pred_mat == 1)
    list_of_indices = list(zip(pred_pair_indices[0].numpy(), pred_pair_indices[1].numpy()))

    pair_types = [sequence[p1] + sequence[p2] for p1, p2 in list_of_indices]

    metrics['pred_has_wc'] = any([pair in wc_pairs for pair in pair_types])
    metrics['pred_has_wobble'] = any([pair in wobble_pairs for pair in pair_types])
    metrics['pred_has_nc'] = any([pair not in wc_pairs + wobble_pairs for pair in pair_types])
    metrics['pred_has_canonicals'] = any([pair in wc_pairs + wobble_pairs for pair in pair_types])

    true_pair_indices = torch.where(true_mat == 1)
    list_of_indices = list(zip(true_pair_indices[0].numpy(), true_pair_indices[1].numpy()))

    pair_types = [sequence[p1] + sequence[p2] for p1, p2 in list_of_indices]

    metrics['gt_has_wc'] = any([pair in wc_pairs for pair in pair_types])
    metrics['gt_has_wobble'] = any([pair in wobble_pairs for pair in pair_types])
    metrics['gt_has_nc'] = any([pair not in wc_pairs + wobble_pairs for pair in pair_types])
    metrics['gt_has_canoncical'] = any([pair in wc_pairs + wobble_pairs for pair in pair_types])

    return metrics


# if __name__ == '__main__':
#     json_path_pred = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/dssr/7d8o_af_dssr.json'
#     json_path_exp = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/dssr/7d8o_pdb_dssr.json'
#     plot_matrix = False
#     all_metrics = True
#
#     sequence, _, _, _ = parse_json(json_path_pred)
#     print("pred sequence", sequence)
#     pred_pairs = extract_base_pairs(json_path_pred)
#     print("pred_pairs", pred_pairs)
#     true_pairs = extract_base_pairs(json_path_exp)
#     print("true_pairs", true_pairs)
#     sequence, _, _, _ = parse_json(json_path_exp)
#     print("exp sequence", sequence)
#
#     metrics = eval_pairs(sequence, pred_pairs, true_pairs, plot_matrix=plot_matrix)
#     if all_metrics:
#         print('Results:')
#         for k, v in metrics.items():
#             if isinstance(v, torch.Tensor):
#                 print(f'{k}: {v.cpu().item()}')
#             else:
#                 print(f'{k}: {v}')
#     else:
#         print('F1-Score:', metrics['f1_score'])
#         print('Precision:', metrics['precision'])
#         print('Recall:', metrics['recall'])
#         print('MCC:', metrics['mcc'])

def get_sequence_length_from_pdb(pdb_file):
    # Get sequence length from PDB file by counting unique residue numbers
    if not os.path.exists(pdb_file):
        print(f"Warning: File not found: {pdb_file}")
        return None
        
    residue_numbers = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    res_num_str = line[22:26].strip()
                    if res_num_str.startswith('-'):
                        res_num = -int(res_num_str[1:])
                    else:
                        if not (res_num_str.isdigit() or (res_num_str.startswith('-') and res_num_str[1:].isdigit())):
                            continue
                        res_num = int(res_num_str)
                    residue_numbers.add(res_num)
                except ValueError:
                    continue
    
    return len(residue_numbers)

def process_dssr_folder(dssr_folder):
    log_file = os.path.join(os.path.dirname(dssr_folder), 'dssr_processing.log')
    original_stdout = sys.stdout
    with open(log_file, 'w') as f:
        sys.stdout = f  # Redirect stdout to file
        
        try:
            database = DatabaseMethods()
            config = StartConfig()

            pdb_groups = {}
            for filename in os.listdir(dssr_folder):
                if not filename.endswith('.json'):
                    continue

                if '_af_' in filename:
                    pdb_id = filename.split('_')[0].lower()
                    if pdb_id not in pdb_groups:
                        pdb_groups[pdb_id] = {'af': None, 'pdb': None}
                    pdb_groups[pdb_id]['af'] = filename
                    print(f"Found AF file: {filename} for {pdb_id}")
                elif '_pdb_' in filename:
                    pdb_id = filename.split('_')[0].lower()
                    if pdb_id not in pdb_groups:
                        pdb_groups[pdb_id] = {'af': None, 'pdb': None}
                    pdb_groups[pdb_id]['pdb'] = filename

            updates = []
            for pdb_id, files in pdb_groups.items():
                if not (files['af'] and files['pdb']):
                    continue
                    
                try:
                    print(f"\nProcessing {pdb_id}:")

                    af_path = os.path.join(dssr_folder, files['af'])
                    pdb_path = os.path.join(dssr_folder, files['pdb'])
                    parent_folder = os.path.dirname(dssr_folder)
                    aligned_folder = os.path.join(parent_folder, 'aligned')
                    
                    # Get AF length from aligned PDB file
                    af_aligned_file = os.path.join(aligned_folder, f"{pdb_id}_af_aligned_rna.pdb")
                    af_length_from_pdb = get_sequence_length_from_pdb(af_aligned_file)
                    
                    # Get sequence lengths from DSSR
                    af_seq, af_id_map, _, _ = parse_json(af_path)
                    pdb_seq, pdb_id_map, _, _ = parse_json(pdb_path)
                    af_length_from_dssr = len(af_seq)
                    
                    print(f"AF sequence length from DSSR: {af_length_from_dssr}")
                    print(f"AF sequence length from PDB: {af_length_from_pdb}")
                    
                    # Use PDB length if available, otherwise fall back to DSSR
                    if af_length_from_pdb is not None:
                        af_length = af_length_from_pdb
                    else:
                        af_length = af_length_from_dssr
                        
                    print(f"Using AF length: {af_length}")

                    residue_map = get_pdb_residue_mapping(pdb_id, aligned_folder)

                    # Validate mapping against AF length
                    if residue_map:
                        mapped_length = max(residue_map.values())
                        print(f"PDB mapped length: {mapped_length}")
                        if mapped_length != af_length:
                            print(f"Warning: Length mismatch - AF: {af_length}, PDB mapped: {mapped_length}")
                            print(f"Skipping {pdb_id}")
                            continue

                    pred_pairs = extract_base_pairs(af_path, is_af=True, residue_map=residue_map)
                    true_pairs = extract_base_pairs(pdb_path, is_af=False, residue_map=residue_map, target_length=af_length)

                    print("\nAF pairs:")
                    print(pred_pairs)
                    print("\nPDB pairs:")
                    print(true_pairs)

                    metrics = eval_pairs(af_seq, pred_pairs, true_pairs, plot_matrix=False)
                    values = tuple(round(float(metrics[m]), 2) for m in 
                                 ['f1_score', 'precision', 'recall', 'mcc', 'wl'])
                    
                    updates.append((pdb_id.upper(), values))

                    print(f"Results for {pdb_id}:")
                    print(f"F1-Score: {metrics['f1_score']:.2f}")
                    print(f"Precision: {metrics['precision']:.2f}")
                    print(f"Recall: {metrics['recall']:.2f}")
                    print(f"MCC: {metrics['mcc']:.2f}")
                    print(f"WL: {metrics['wl']:.2f}")
                    print(f"WL with seq: {metrics['wl_with_seq']:.2f}")
                    print("=" * 40)

                except Exception as e:
                    print(f"Error processing {pdb_id}: {e}")
                    continue

            if updates:
                try:
                    update_query = f"""
                        UPDATE {config.pred_table}
                        SET RNA_f1_score = ?,
                            RNA_precision = ?,
                            RNA_recall = ?,
                            RNA_mcc = ?,
                            RNA_wl = ?
                        WHERE exp_db_id = ?
                    """

                    with database.connection:
                        for pdb_id, values in updates:
                            database.connection.execute(update_query, values + (pdb_id,))
                            
                    print(f"\nSuccessfully updated {len(updates)} records in database")
                    
                except Exception as e:
                    print(f"Error updating database: {e}")
                    
            database.close_connection()
            
        finally:
            sys.stdout = original_stdout  # Restore stdout
            print(f"Processing complete. Log saved to: {log_file}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python extractDSSRpairs.py <dssr_folder>")
        sys.exit(1)

    dssr_folder = sys.argv[1]
    process_dssr_folder(dssr_folder)