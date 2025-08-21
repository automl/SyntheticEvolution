import pandas as pd
import numpy as np

from pathlib import Path

import pickle

def mat2pairs(matrix, symmetric=True):  # TODO: Verify that this one is working correctly.
    """
    Convert matrix representation of structure to list of pairs.
    """
    if symmetric:
        return list(set(tuple(sorted(pair)) for pair in np.argwhere(matrix == 1)))
    else:
        return list(tuple(pair) for pair in np.argwhere(matrix == 1))



data_dir = Path('datasets')
df_path = Path(data_dir, 'b2021_rnaonly_l500_predictions.pkl')
# data_path = Path(data_dir, 'a2021_rnaonly.pkl')

df = pd.read_pickle(df_path)
# data_df = pd.read_pickle(data_path)

seq_vocab = ['A', 'C', 'G', 'U', 'N']
seq_stoi = dict(zip(seq_vocab, range(len(seq_vocab))))
seq_itos = {v: k for k, v in seq_stoi.items()}

df.loc[:, 'Id'] = df['Id'].apply(lambda x: ''.join([chr(i.item()) for i in x]))

df.loc[:, 'sequence'] = df['sequence'].apply(lambda x: ''.join([seq_itos[s] for s in x[0]]))

df.loc[:, 'pairs'] = df['pred_mat'].apply(lambda x: mat2pairs(x, symmetric=True))

# for i, row in df.iterrows():
#     print('Sample', row['Id'])
#     for i in range(row['pred_mat'].shape[0]):
#         for j in range(row['pred_mat'].shape[0]):
#             if row['pred_mat'][i,j]:
#                 assert (i, j) in row['pairs'] or (j, i) in row['pairs']
#                 print('found pair', (i, j))

# print(df)
# print(data_df)

outpath = Path(data_dir, f"{df_path.stem}_processed.pkl")

with open(outpath, 'wb') as f:
    pickle.dump(df, f)

print('processed Dataframe written to', outpath.resolve())
df = pd.read_pickle(outpath)
print(df)

