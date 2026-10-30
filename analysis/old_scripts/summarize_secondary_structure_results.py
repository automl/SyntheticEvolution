import pandas as pd
import numpy as np

res_df = pd.read_pickle('results/secondary_structure_eval_results_multichain_True.pkl')


data = []
algos = ['AlphaFold_dssr', 'RNAformer_dssr', 'RNAformer_pred', 'SpotRNA_dssr', 'RNAfold_dssr']
for pdb_id, group in res_df.groupby('pdb_id'):
    if len(group) == len(algos):
        data.append(group)

eval_df = pd.concat(data)

short = [] # <20
medium = [] # 20-100
medium_2 = [] # 100-200
long = [] # 200-500
very_long = [] # >500

for model, group in eval_df.groupby('model'):
    # print(group['length'].max())
    # elif group['length'].max() >= 500:
    #     very_long.append(group)
    print(model, len(group), np.round(np.mean(group['f1']), 3), np.round(np.mean(group['mcc']), 3), np.round(np.mean(group['wl']), 3), np.round(np.mean(group['precision']), 3), np.round(np.mean(group['recall']), 3), np.round(np.mean(group['specificity']), 3))

for pdb_id, group in eval_df.groupby('pdb_id'):
    # print(len(group), pdb_id, group['model'].unique())
    if len(group) == len(algos):
        max_length = group['length'].max()
        if max_length < 20:
            short.append(group)
        elif 20 <= max_length < 100:
            medium.append(group)
        elif 100 <= max_length < 200:
            medium_2.append(group)
        elif 200 <= max_length < 500:
            long.append(group)
        
        elif max_length >= 500:
            very_long.append(group)

print()
if short:
    short_df = pd.concat(short)
    print('Short sequences (<20):')
    for model, group in short_df.groupby('model'):
        print(model, len(group), np.round(np.mean(group['f1']), 3), np.round(np.mean(group['mcc']), 3), np.round(np.mean(group['wl']), 3), np.round(np.mean(group['precision']), 3), np.round(np.mean(group['recall']), 3), np.round(np.mean(group['specificity']), 3))
    print()
    # print("Short sequences (<20):", len(short_df), np.round(np.mean(short_df['f1']), 3), np.round(np.mean(short_df['mcc']), 3), np.round(np.mean(short_df['wl']), 3), np.round(np.mean(short_df['precision']), 3), np.round(np.mean(short_df['recall']), 3), np.round(np.mean(short_df['specificity']), 3))

if medium:
    medium_df = pd.concat(medium)
    print('Medium sequences (20-100):')
    for model, group in medium_df.groupby('model'):
        print(model, len(group), np.round(np.mean(group['f1']), 3), np.round(np.mean(group['mcc']), 3), np.round(np.mean(group['wl']), 3), np.round(np.mean(group['precision']), 3), np.round(np.mean(group['recall']), 3), np.round(np.mean(group['specificity']), 3))
    print()
    # print("Medium sequences (20-100):", len(medium_df), np.round(np.mean(medium_df['f1']), 3), np.round(np.mean(medium_df['mcc']), 3), np.round(np.mean(medium_df['wl']), 3), np.round(np.mean(medium_df['precision']), 3), np.round(np.mean(medium_df['recall']), 3), np.round(np.mean(medium_df['specificity']), 3))
if medium_2:
    print('Medium sequences (100-200):')
    medium_2_df = pd.concat(medium_2)
    for model, group in medium_2_df.groupby('model'):
        print(model, len(group), np.round(np.mean(group['f1']), 3), np.round(np.mean(group['mcc']), 3), np.round(np.mean(group['wl']), 3), np.round(np.mean(group['precision']), 3), np.round(np.mean(group['recall']), 3), np.round(np.mean(group['specificity']), 3))
    print()
    # print("Medium sequences (100-200):", len(medium_2_df), np.round(np.mean(medium_2_df['f1']), 3), np.round(np.mean(medium_2_df['mcc']), 3), np.round(np.mean(medium_2_df['wl']), 3), np.round(np.mean(medium_2_df['precision']), 3), np.round(np.mean(medium_2_df['recall']), 3), np.round(np.mean(medium_2_df['specificity']), 3))
if long:
    long_df = pd.concat(long)
    print('Long sequences (200-500):')
    for model, group in long_df.groupby('model'):
        print(model, len(group), np.round(np.mean(group['f1']), 3), np.round(np.mean(group['mcc']), 3), np.round(np.mean(group['wl']), 3), np.round(np.mean(group['precision']), 3), np.round(np.mean(group['recall']), 3), np.round(np.mean(group['specificity']), 3))
    print()
    # print("Long sequences (200-500):", len(long_df), np.round(np.mean(long_df['f1']), 3), np.round(np.mean(long_df['mcc']), 3), np.round(np.mean(long_df['wl']), 3), np.round(np.mean(long_df['precision']), 3), np.round(np.mean(long_df['recall']), 3), np.round(np.mean(long_df['specificity']), 3))

if very_long:
    print('Very Long sequences (>500):')
    very_long_df = pd.concat(very_long)
    for model, group in very_long_df.groupby('model'):
        print(model, len(group), np.round(np.mean(group['f1']), 3), np.round(np.mean(group['mcc']), 3), np.round(np.mean(group['wl']), 3), np.round(np.mean(group['precision']), 3), np.round(np.mean(group['recall']), 3), np.round(np.mean(group['specificity']), 3))
    print()

per_id_diffs = []

for pdb_id, group in eval_df.groupby('pdb_id'):
    if len(group) > len(algos):
        continue
    diff = group[group['model'] == 'AlphaFold_dssr']['mcc'].values[0] - group[group['model'] == 'RNAformer_dssr']['mcc'].values[0]
    per_id_diffs.append((pdb_id, diff))

per_id_diffs = sorted(per_id_diffs, key=lambda x: x[1])
for i, mcc in per_id_diffs:
    print(i, mcc)
