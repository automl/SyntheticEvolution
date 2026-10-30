import pandas as pd

from pathlib import Path
import numpy as np

from grakel import Graph
from grakel.kernels import WeisfeilerLehman, VertexHistogram, WeisfeilerLehmanOptimalAssignment, ShortestPath

import torch
import pickle

def pairs2mat(pairs, length, symmetric=True, no_pk=False):
    """
    Convert list of pairs to matrix representation of structure.
    """
    try:
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
    except Exception as e:
        print(f"Error converting pairs to matrix: {e}")
        return np.nan

#####METRICS#####
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


def graph_distance_score_from_matrices(pred, true, kernel, node_labels=None):
    pred_graph = mat2graph(pred, node_labels=node_labels)
    true_graph = mat2graph(true, node_labels=node_labels)
    kernel = get_graph_kernel(kernel=kernel)
    kernel.fit_transform([true_graph])
    distance_score = kernel.transform([pred_graph])  # TODO: Check output, might be list or list of lists

    return distance_score[0][0]

# is tested
def f1(tp, fp, tn, fn):
    f1_score = 2 * tp / (2 * tp + fp + fn + 1e-8)
    return f1_score

# is tested
def recall(tp, fp, tn, fn):
    recall = tp / (tp + fn + 1e-8)
    return recall

# is tested
def specificity(tp, fp, tn, fn):
    specificity = tn / (tn + fp + 1e-8)
    return specificity

# is tested
def precision(tp, fp, tn, fn):
    precision = tp / (tp + fp + 1e-8)
    return precision

# is tested
def mcc(tp, fp, tn, fn):
    mcc = (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) + 1e-8)
    return mcc

# is tested
def tp_from_matrices(pred, true):
    tp = np.logical_and(pred, true).sum()
    return tp

# is tested
def tn_from_matrices(pred, true):
    tn = np.logical_and(np.logical_not(pred), np.logical_not(true)).sum()
    return tn

# is tested
def get_fp(pred, tp):
    fp = pred.sum() - tp
    return fp

# is tested
def get_fn(true, tp):
    fn = true.sum() - tp
    return fn


# mat_data = pd.read_pickle('secondary_structure_data.pkl')
# mat_data = pd.read_pickle('a2021_secondary_structure_data.pkl')
# mat_data = pd.read_pickle('b2021_secondary_structure_data.pkl')
# mat_data = pd.read_pickle('secondary_structure_data_synthetic_msa_selection_max_rna_len_200_min_rna_len_15_max_protein_len_200.pkl')
# mat_data = pd.read_pickle('secondary_structure_data_synthetic_msa_selection_max_rna_len_200_min_rna_len_15_max_protein_len_200_2.pkl')
# mat_data = pd.read_pickle('rnaprotein_rnaformerafc23_msa_size100_secondary_structure_data.pkl')
# mat_data = pd.read_pickle('rnaprotein_rnaformerafc23_msa_size500_secondary_structure_data.pkl')
# mat_data = pd.read_pickle('rnaprotein_rnaformerafc23_af3_dssr_gt_secondary_structure_data.pkl')
# mat_data = pd.read_pickle('a2021_rnaonly_rnaformerafc23_af3_dssr_gt_secondary_structure_data.pkl')
mat_data = pd.read_pickle('results/secondary_structure_information_multichain_True.pkl')

# outpath = 'complex_secondary_structure_evaluation_results2.pkl'
# outpath = 'a2021_secondary_structure_evaluation_results_no_empty.pkl'
# outpath = 'b2021_secondary_structure_evaluation_results_no_empty.pkl'
# outpath = 'rnaprotein_secondary_structure_evaluation_results_no_empty.pkl'
# outpath = 'a2021_rnaonly_rnaformerafc23_af3_dssr_gt_secondary_structure_evaluation_results.pkl'
outpath = 'results/secondary_structure_eval_results_multichain_True.pkl'

ignore_empty_mat = False

mat_data = mat_data.drop_duplicates(subset=['pdb_id', 'chain', 'model'], keep='first')
mat_data.loc[:, 'length'] = mat_data['sequence'].apply(lambda x: len(x))

mat_data.loc[:, 'pair_mat'] = mat_data.apply(lambda x: pairs2mat(x.pairs, len(x.sequence), symmetric=True, no_pk=True), axis=1)
if ignore_empty_mat:
    print('Remove empty gt mats...')
    print(len(mat_data))
    mat_data.loc[:, 'valid_gt_mat'] = True
    mat_data.loc[:, 'valid_gt_mat'] = mat_data.apply(lambda x: False if x.model == 'GroundTruth' and not any(x.pair_mat[i, j] for i in range(x.pair_mat.shape[0]) for j in range(x.pair_mat.shape[1])) else True, axis=1)
    mat_data = mat_data[mat_data['valid_gt_mat']]
    print(len(mat_data))

# print(mat_data)

gt_data = mat_data[mat_data['model'] == 'GroundTruth']
mat_data = mat_data[mat_data['model'] != 'GroundTruth']
# gt_data = gt_data.drop_duplicates(subset=['pdb_id'], keep='first')

# print(gt_data)
print(mat_data)
print(gt_data)

results = []

for pdb_id, group in mat_data.groupby('pdb_id'):
    # print(pdb_id)
    # print(group)
    try:
        gt_mat = gt_data[gt_data['pdb_id'] == pdb_id]['pair_mat'].values[0]
        # print(gt_mat)
        gt_sequence = gt_data[gt_data['pdb_id'] == pdb_id]['sequence'].values[0]
    
        for i, row in group.iterrows():
            pred_sequence = row['sequence']
            pred_length = len(pred_sequence)
            pred_mat = row['pair_mat']
            model = row['model']
            if not pred_mat.shape[0] == gt_mat.shape[0]:
                print(f"Length mismatch for {pdb_id} model {model}: {pred_length} vs {gt_mat.shape[0]}")
                continue
    
            sample = {
                'pdb_id': pdb_id,
                'model': model,
                'length': pred_length,
            }
    
            sample['tp'] = tp_from_matrices(pred_mat, gt_mat)
            sample['tn'] = tn_from_matrices(pred_mat, gt_mat)
            sample['fp'] = get_fp(pred_mat, sample['tp'])
            sample['fn'] = get_fn(gt_mat, sample['tp'])
            sample['precision'] = precision(sample['tp'], sample['fp'], sample['tn'], sample['fn'])
            sample['recall'] = recall(sample['tp'], sample['fp'], sample['tn'], sample['fn'])
            sample['specificity'] = specificity(sample['tp'], sample['fp'], sample['tn'], sample['fn'])
            sample['f1'] = f1(sample['tp'], sample['fp'], sample['tn'], sample['fn'])
            sample['mcc'] = mcc(sample['tp'], sample['fp'], sample['tn'], sample['fn'])
            sample['wl'] = graph_distance_score_from_matrices(pred_mat, gt_mat, kernel='WeisfeilerLehman')
    
            results.append(sample)
    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")
        continue

results_df = pd.DataFrame(results)

with open(outpath, 'wb') as f:
    pickle.dump(results_df, f)




        # print(model)
        # print(row)
        # if row['pair_mat'].sum() == 0:
        #     continue
        # print(row['pair_mat'])
        # print(row['sequence'])
        # print(row['dot_bracket'])
        # print(row['pairs'])
        # print(row['multiplet_positions'])
        # print(row['pos1'])
        # print(row['pos2'])



