import torch
import numpy as np
import pandas as pd
import json

from matplotlib import pyplot as plt
from scipy import signal
from grakel import Graph
from grakel.kernels import WeisfeilerLehman, VertexHistogram, WeisfeilerLehmanOptimalAssignment, ShortestPath
from pathlib import Path


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

    return sequence, id_2_index_map, df, has_multiplet    # return data


def get_json_data(json_path):
    NUCS = {
        'T': 'U',
        'P': 'U',
        'R': 'A',  # or 'G'
        'Y': 'C',  # or 'T'
        'M': 'C',  # or 'A'
        'K': 'U',  # or 'G'
        'S': 'C',  # or 'G'
        'W': 'U',  # or 'A'
        'H': 'C',  # or 'A' or 'U'
        'B': 'U',  # or 'G' or 'C'
        'V': 'C',  # or 'G' or 'A'
        'D': 'A',  # or 'G' or 'U'
        'N': 'N',
        'A': 'A',
        'U': 'U',
        'C': 'C',
        'G': 'G',
    }

    json_path = Path(json_path)

    sequence_from_json, id_map, pair_df, has_multiplet = parse_json(json_path)
        
    if pair_df is None:
        return []
    
    sequence = sequence_from_json.replace('T', 'U')
    sequence = sequence.replace('&', '')
    json_stem = json_path.stem

    pos1 = [(NUCS[id_map[p][0]], int(id_map[p][2])) if not '^' in p.split('.')[1][1:] else (NUCS[id_map[p][0]], int(id_map[p][2])) for p in pair_df['nt1']]

    pos2 = [(NUCS[id_map[p][0]], int(id_map[p][2])) if not '^' in p.split('.')[1][1:] else (NUCS[id_map[p][0]], int(id_map[p][2])) for p in pair_df['nt2']]

    pairs = []
    
    sequence = sequence.upper()
    try:
        for i in range(len(pos1)):
            assert str(sequence[pos1[i][1]-1]) == str(pos1[i][0]) or str(sequence[pos1[i][1]-1]) == 'N' or str(pos1[i][0]) == 'N', f"{i}: {pos1[i]};\n{sequence[pos1[i][1]-1]} != {pos1[i][0]}"
        for i in range(len(pos2)):
            assert str(sequence[pos2[i][1]-1]) == str(pos2[i][0]) or str(sequence[pos2[i][1]-1]) == 'N' or str(pos2[i][0]) == 'N', f"{i}: {pos2[i]};\n{sequence[pos2[i][1]-1]} != {pos2[i][0]}"
        # assert all(str(sequence[pos1[i][1]-1]) == str(pos1[i][0]) for i in range(len(pos1))), f"{sequence[pos1[0][1]-1]} != {pos1[0][0]}"
        # assert all(str(sequence[pos2[i][1]-1]) == str(pos2[i][0]) for i in range(len(pos2))), f"{sequence[pos2[0][1]-1]} != {pos2[0][0]}"
        pairs = [[p1[1]-1, p2[1]-1] for p1, p2 in zip(pos1, pos2)]
    except:
        try:
            for i in range(len(pos1)):
                assert str(sequence[pos1[i][1]]) == str(pos1[i][0]) or str(sequence[pos1[i][1]]) == 'N' or str(pos1[i][0]) == 'N', f"{i}: {pos1[i]};\n{sequence[pos1[i][1]]} != {pos1[i][0]}"
            for i in range(len(pos2)):
                assert str(sequence[pos2[i][1]]) == str(pos2[i][0]) or str(sequence[pos2[i][1]]) == 'N' or str(pos2[i][0]) == 'N', f"{i}: {pos2[i]};\n{sequence[pos2[i][1]]} != {pos2[i][0]}"
            # assert all(str(sequence[pos1[i][1]]) == str(pos1[i][0]) for i in range(len(pos1))), f"{sequence[pos1[0][1]]} != {pos1[0][0]}"
            # assert all(str(sequence[pos2[i][1]]) == str(pos2[i][0]) for i in range(len(pos2))), f"{sequence[pos2[0][1]]} != {pos2[0][0]}"
            pairs = [[p1[1], p2[1]] for p1, p2 in zip(pos1, pos2)]
        except:
            try:
                for i in range(len(pos1)):
                    assert str(sequence[pos1[i][1]+1]) == str(pos1[i][0]) or str(sequence[pos1[i][1]+1]) == 'N' or str(pos1[i][0]) == 'N', f"{i}: {pos1[i]};\n{sequence[pos1[i][1]+1]} != {pos1[i][0]}"
                for i in range(len(pos2)):
                    assert str(sequence[pos2[i][1]+1]) == str(pos2[i][0]) or str(sequence[pos2[i][1]+1]) == 'N' or str(pos2[i][0]) == 'N', f"{i}: {pos2[i]};\n{sequence[pos2[i][1]+1]} != {pos2[i][0]}"
                # assert all(str(sequence[pos1[i][1]+1]) == str(pos1[i][0]) for i in range(len(pos1))), f"{sequence[pos1[0][1]+1]} != {pos1[0][0]}"
                # assert all(str(sequence[pos2[i][1]+1]) == str(pos2[i][0]) for i in range(len(pos2))), f"{sequence[pos2[0][1]+1]} != {pos2[0][0]}"
                pairs = [[p1[1]+1, p2[1]+1] for p1, p2 in zip(pos1, pos2)]
            except:
                pass
    
    if not pairs:
        return []    
    
    return sorted(pairs, key=lambda x: x[0]), sequence

def pairs2mat(pairs, length, symmetric=True):
    """
    Convert list of pairs to matrix representation of structure.
    """
    # print(pairs)
    mat = np.zeros((length, length))

    for p1, p2 in pairs:
        mat[p1, p2] = 1
        if symmetric:
            mat[p2, p1] = 1
    return torch.from_numpy(mat)


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
                                   range(matrix.shape[0])})  # TODO: Think about if we need to label the nodes differenty

    return graph


def evaluate_shifted(pred_a, true_a):
    kernel = np.array([[0.0,1.0,0.0],
                        [1.0,1.0,1.0],
                        [0.0,1.0,0.0]])
    pred_a_filtered = signal.convolve2d(pred_a, kernel, 'same')
    fn = len(torch.where((true_a - torch.Tensor(pred_a_filtered))==1)[0])
    pred_p = torch.sign(torch.Tensor(pred_a)).sum()
    true_p = true_a.sum()
    tp = true_p - fn
    fp = pred_p - tp
    recall = tp/(tp+fn)
    precision = tp/(tp+fp)
    f1_score = 2*tp/(2*tp + fp + fn)
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

    # Watson-Crick pairs: AU, UA, CG, GC
    wc_pairs = ['AU', 'UA', 'CG', 'GC']
    # Wobble pairs: GU, UG
    wobble_pairs = ['GU', 'UG']
    
    # Initialize masks
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
    for i, (p1,p2) in enumerate(zip(true_pair_positions[0], true_pair_positions[1])):

        if ''.join(sequence[p1]+sequence[p2]) in wc_pairs:
            wc[p1, p2] = 1
            wc[p2, p1] = 1
            canonicals[p1, p2] = 1
            canonicals[p2, p1] = 1
        elif ''.join(sequence[p1]+sequence[p2]) in wobble_pairs:
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

    assert (wc+wobble+nc).sum() == true_mat.sum()
    assert (canonicals+nc).sum() == true_mat.sum()

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
    wl_with_seq = graph_distance_score_from_matrices(pred_mat.cpu().numpy(), true_mat.cpu().numpy(), 'WeisfeilerLehman', true_labels=true_labels, pred_labels=pred_labels)

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
    
    for (name, mask) in [('wc', wc), ('wobble', wobble), ('nc', nc), ('canonicals', canonicals), ('pk', pk), ('multi', multi)]:
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
            metrics['num_pk_hits']  = tp.cpu().item()
            metrics['num_pk_gt'] = mask.sum().cpu().item()
        if name == 'multi':
            metrics['multi_hit'] = tp > 0
            metrics['num_multi_hits']  = tp.cpu().item()
            metrics['num_multi_gt'] = mask.sum().cpu().item()
    

    pred_pair_indices = torch.where(pred_mat == 1)
    list_of_indices = list(zip(pred_pair_indices[0].numpy(), pred_pair_indices[1].numpy()))

    pair_types = [sequence[p1]+sequence[p2] for p1, p2 in list_of_indices]

    metrics['pred_has_wc'] = any([pair in wc_pairs for pair in pair_types])
    metrics['pred_has_wobble'] = any([pair in wobble_pairs for pair in pair_types])
    metrics['pred_has_nc'] = any([pair not in wc_pairs+wobble_pairs for pair in pair_types])
    metrics['pred_has_canonicals'] = any([pair in wc_pairs+wobble_pairs for pair in pair_types])

    true_pair_indices = torch.where(true_mat == 1)
    list_of_indices = list(zip(true_pair_indices[0].numpy(), true_pair_indices[1].numpy()))

    pair_types = [sequence[p1]+sequence[p2] for p1, p2 in list_of_indices]

    metrics['gt_has_wc'] = any([pair in wc_pairs for pair in pair_types])
    metrics['gt_has_wobble'] = any([pair in wobble_pairs for pair in pair_types])
    metrics['gt_has_nc'] = any([pair not in wc_pairs+wobble_pairs for pair in pair_types])
    metrics['gt_has_canoncical'] = any([pair in wc_pairs+wobble_pairs for pair in pair_types])


    return metrics


if __name__ == '__main__':
    # json_path = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/json'
    # plot_matrix = True
    # all_metrics = False
    #
    # pred_pairs, sequence = get_json_data(json_path)
    #
    # # true pairs match prediction exactly. remove / add pairs to test different metrics
    # true_pairs = [[0, 26], [1, 25], [2, 24], [3, 23], [4, 22], [5, 9], [8, 21], [9, 20], [10, 19], [11, 18], [12, 16]]

    json_path_pred = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/dssr/6fql_af_dssr.json'
    json_path_exp = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/dssr/6fql_pdb_dssr.json'
    plot_matrix = True
    all_metrics = True

    # Get predicted pairs and sequence
    pred_pairs, sequence = get_json_data(json_path_pred)
    print("pred_pairs", pred_pairs)

    # Get experimental pairs and use them as true pairs
    true_pairs, _ = get_json_data(json_path_exp)
    print("true_pairs", true_pairs)

    metrics = eval_pairs(sequence, pred_pairs, true_pairs, plot_matrix=plot_matrix)
    if all_metrics:
        print('Results:')
        for k, v in metrics.items():
            if isinstance(v, torch.Tensor):
                print(f'{k}: {v.cpu().item()}')
            else:
                print(f'{k}: {v}')
    else:
        print('F1-Score:', metrics['f1_score'])
        print('Precision:', metrics['precision'])
        print('Recall:', metrics['recall'])
        print('MCC:', metrics['mcc'])
        

