import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import pandas.api.types as ptypes

from scipy.stats import wilcoxon, ttest_rel

def plot_algorithm_distributions(
    df: pd.DataFrame,
    algorithm_col: str,
    value_col: str,
    kind: str = 'violin',
    palette: list = None,
    fig_size: tuple = (10, 6),
    title: str = None,
    title_size: int = 14,
    xlabel: str = None,
    xlabel_size: int = 12,
    ylabel: str = None,
    ylabel_size: int = 12,
    xtick_label_size: int = 10,
    ytick_label_size: int = 10,
    xtick_label_rotation: float = 45,
    remove_spines: bool = True,
    annotate: bool = True,
    ref_algorithm: str = None,
    test_method: str = 'wilcoxon',  # or 'ttest'
    sig_levels: list = [(0.001, '***'), (0.01, '**'), (0.05, '*')],
    num_algos_for_statistics: int = 2,
    box_width: float = 0.6,
    spacing: float = 1.0,
    annot_offset: float = 0.05,
    annot_step: float = 0.05,
    step_factor: float = 1.0,
    outpath: Path = Path('algorithm_comparison.svg')
):
    """
    Plot violin/box for each algorithm, with full control over box_width, spacing,
    annotation offsets, and detailed Wilcoxon/t-test output.
    """
    def wilcoxon_stats(x, y):
        diffs = x - y
        nz = diffs != 0
        d = diffs[nz]
        n = len(d)
        if n == 0:
            return None
        W, p = wilcoxon(x, y)
        n_pos = np.sum(diffs > 0)
        n_neg = np.sum(diffs < 0)
        r_rb = (n_pos - n_neg) / n
        mean_W = n*(n+1)/4
        sd_W   = np.sqrt(n*(n+1)*(2*n+1)/24)
        Z      = (W - mean_W) / sd_W
        r      = Z / np.sqrt(n)
        return W, p, r_rb, Z, r

    def ttest_stats(x, y):
        t, p = ttest_rel(x, y)
        d = np.mean(x - y) / np.std(x - y, ddof=1) if len(x)>1 else np.nan
        return t, p, d

    data = df.copy()
    # pick algorithms
    if pd.api.types.is_categorical_dtype(data[algorithm_col]):
        algos = list(data[algorithm_col].cat.categories)
    else:
        algos = sorted(data[algorithm_col].unique())
    n = len(algos)

    # map & print means
    data_map = {}
    for alg in algos:
        vals = data.loc[data[algorithm_col]==alg, value_col].dropna().values
        data_map[alg] = vals
        print(f"{alg}: mean={np.mean(vals):.3f}, n={len(vals)}")

    data_to_plot = [data_map[alg] for alg in algos]
    positions    = np.arange(n) * spacing

    # pick colors
    if palette is None:
        # cmap = plt.get_cmap('Blues')
        # palette = [cmap(0.4 + 0.6 * i/(n-1)) for i in range(n)]
        palette = ['skyblue' for _ in range(n)]

    fig, ax = plt.subplots(figsize=fig_size)

    # draw
    if kind=='violin':
        parts = ax.violinplot(data_to_plot, positions=positions, widths=box_width,
                              showmedians=True)
        for body,color in zip(parts['bodies'], palette):
            body.set_facecolor(color)
            body.set_edgecolor('black')
            body.set_alpha(0.8)
        if 'cmedians' in parts:
            parts['cmedians'].set(color='firebrick', linewidth=2)
    else:
        parts = ax.boxplot(data_to_plot, positions=positions, widths=box_width,
                           patch_artist=True, showfliers=False)
        for patch,color in zip(parts['boxes'], palette):
            patch.set_facecolor(color)
            patch.set_edgecolor('black')
            patch.set_linewidth(1)
        for med in parts['medians']:
            med.set(color='firebrick', linewidth=2)

    if remove_spines:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # clip at 1.0
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, 1.0)
    y_min, y_max = ymin, 1.0
    y_span = y_max - y_min

    # labels
    ax.set_xticks(positions)
    ax.set_xticklabels(algos, rotation=xtick_label_rotation,
                       ha='right', fontsize=xtick_label_size)
    if xlabel: ax.set_xlabel(xlabel, fontsize=xlabel_size)
    if ylabel: ax.set_ylabel(ylabel, fontsize=ylabel_size)
    ax.tick_params(axis='y', labelsize=ytick_label_size)
    if title: ax.set_title(title, loc='left', fontsize=title_size, pad=10)

    # annotate
    if annotate and n>1:
        ref = ref_algorithm if ref_algorithm in algos else algos[0]
        ref_idx = algos.index(ref)
        ref_vals = data_map[ref]

        base_offset = annot_offset * y_span
        step = annot_step * y_span * step_factor

        for idx,alg in enumerate(algos[:num_algos_for_statistics]):
            if alg==ref: continue
            vals = data_map[alg]
            if test_method=='ttest':
                t,p,d = ttest_stats(ref_vals, vals)
                print(f"T-test {ref} vs {alg}: t={t:.3f}, p={p:.3g}, d={d:.3f}")
            else:
                res = wilcoxon_stats(ref_vals, vals)
                if res:
                    W,p,r_rb,Z,r = res
                    print(f"Wilcoxon {ref} vs {alg}: W={W:.1f}, p={p:.3g}, r_rb={r_rb:.3f}, Z={Z:.3f}, r={r:.3f}")
            # significance?
            stars = next((s for thr,s in sig_levels if p<thr), '')
            if not stars:
                continue

            # level above max box value
            ymax_box = max(np.max(ref_vals), np.max(vals))
            y = ymax_box + base_offset

            # if you want to keep inside plot, optionally uncomment:
            # y = min(y, y_max - 0.01*y_span)

            x1 = positions[ref_idx]
            x2 = positions[idx]
            ax.plot([x1,x1,x2,x2], [y, y+step/2, y+step/2, y],
                    color='black', linewidth=1, clip_on=False)
            ax.text((x1+x2)/2, y+step/2, stars,
                    ha='center', va='bottom',
                    fontsize=ytick_label_size, clip_on=False)

            base_offset += step

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.show()
    plt.close(fig)


#############################################################################################################
# Load and process secondary structure results for RNAformer and other algorithms

dataset = 'pdb_ts123'
datadir = Path('evaluation/predictions/secondary_structure_rnaformer_eval/final_large_table')
Path('plots', 'performance').mkdir(exist_ok=True, parents=True)
outpath = Path('plots', 'performance', 'secondary_structure_ts_pdb_comparison_box.svg')
score_of_interest = 'mcc'  # 'f1_score'

fl = list(datadir.glob(f"*{dataset}*.pkl"))

data = []
rnaformer_data = []

for f in fl:
    algo = str(f.stem).split('_')[0]
    if algo == 'finetune':
        df = pd.read_pickle(f)
        # df['algorithm'] = 'RNAformerFT'
        # # print(df.columns)
        if '_s1_' in str(f.stem):
            df['algorithm'] = 'RNAformerFT_s1'
        elif '_s2_' in str(f.stem):
            df['algorithm'] = 'RNAformerFT_s2'
        elif '_s3_' in str(f.stem):
            df['algorithm'] = 'RNAformerFT_s3'
        rnaformer_data.append(df)
    elif algo == 'pretrain':
        df = pd.read_pickle(f)
        df['algorithm'] = 'RNAformerPT'
        data.append(df)
        # rnaformer_data.append(df)
    elif 'SPOT-RNA2' in algo:
        continue
    elif 'spot2' in algo:
        continue
    # elif 'RNAStructure' in algo:
    #     df['algorithm'] = 'RNAStructure'
    #     data.append(df)
    else:
        df = pd.read_pickle(f)
        df['algorithm'] = algo
        data.append(df)
data = pd.concat(data, ignore_index=True)



rnaformer_mean = []

cols = ['f1_score', 'mcc', 'wl', 'precision', 'recall']
data = data[cols + ['algorithm']]

for i, row in rnaformer_data[0].iterrows():
    r2 = rnaformer_data[1].loc[i]
    r3 = rnaformer_data[2].loc[i]
    print(row['f1_score'], r2['f1_score'], r3['f1_score'])
    sample = {}
    for c in cols:
        mean_value = np.mean([row[c], r2[c], r3[c]])
        sample[c] = mean_value
    rnaformer_mean.append(sample)
rnaformer_mean_df = pd.DataFrame(rnaformer_mean)
rnaformer_mean_df['algorithm'] = 'RNAformer'
rnaformer_mean_df = rnaformer_mean_df[cols + ['algorithm']]

data = pd.concat([rnaformer_mean_df, data], ignore_index=True)

order = (
    data
    .groupby('algorithm')[score_of_interest]
    .mean()
    .sort_values(ascending=False)
    .index
    .tolist()
)

# 2. Make `algorithm` a categorical with that specific order
data['algorithm'] = pd.Categorical(data['algorithm'], categories=order, ordered=True)

# for algo, group in data.groupby('algorithm'):
#     print(algo, len(group), np.round(np.mean(group['f1_score']), 3), np.round(np.mean(group['mcc']), 3), np.round(np.mean(group['wl']), 3), np.round(np.mean(group['precision']), 3), np.round(np.mean(group['recall']), 3))
# 
# print(data)
plot_algorithm_distributions(
    df=data,
    algorithm_col='algorithm',
    value_col=score_of_interest,
    kind='box',  # 'violin', 'box',
    # palette=[],  # ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'],
    fig_size=(12, 6),
    # title='Secondary Structure F1 Score Comparison',
    # xlabel='Algorithm',
    ylabel='Matthews Correlation Coefficient',  # 'F1 Score',
    xtick_label_rotation=45,
    remove_spines=True,
    xlabel_size=18,
    ylabel_size=18,
    xtick_label_size=18,
    ytick_label_size=16,
    ref_algorithm='RNAformer',
    num_algos_for_statistics=2,
    test_method='wilcoxon',  # 'wilcoxon', 'ttest'
    annotate=True,
    # annot_offset=1.0,
    annot_step=0.07,
    step_factor=1.0,
    box_width=0.5,
    outpath=Path(outpath)
)


############################################################################################################
# RNAformer vs all others on ts-hard
dataset = 'pdb_ts_hard'
datadir = Path('evaluation/predictions/secondary_structure_rnaformer_eval/final_large_table')
outpath = Path('plots', 'performance', 'secondary_structure_ts_hard_comparison_box.svg')
score_of_interest = 'mcc'  # 'f1_score'

fl = list(datadir.glob(f"*{dataset}*.pkl"))

data = []
rnaformer_data = []

for f in fl:
    algo = str(f.stem).split('_')[0]
    if algo == 'finetune':
        df = pd.read_pickle(f)
        # df['algorithm'] = 'RNAformerFT'
        # # print(df.columns)
        if '_s1_' in str(f.stem):
            df['algorithm'] = 'RNAformerFT_s1'
        elif '_s2_' in str(f.stem):
            df['algorithm'] = 'RNAformerFT_s2'
        elif '_s3_' in str(f.stem):
            df['algorithm'] = 'RNAformerFT_s3'
        rnaformer_data.append(df)
    elif algo == 'pretrain':
        df = pd.read_pickle(f)
        df['algorithm'] = 'RNAformerPT'
        data.append(df)
        # rnaformer_data.append(df)
    elif 'SPOT-RNA2' in algo:
        continue
    elif 'spot2' in algo:
        continue
    # elif 'RNAStructure' in algo:
    #     df['algorithm'] = 'RNAStructure'
    #     data.append(df)
    else:
        df = pd.read_pickle(f)
        df['algorithm'] = algo
        data.append(df)
data = pd.concat(data, ignore_index=True)


rnaformer_mean = []

cols = ['f1_score', 'mcc', 'wl', 'precision', 'recall']
data = data[cols + ['algorithm']]

for i, row in rnaformer_data[0].iterrows():
    r2 = rnaformer_data[1].loc[i]
    r3 = rnaformer_data[2].loc[i]
    print(row['f1_score'], r2['f1_score'], r3['f1_score'])
    sample = {}
    for c in cols:
        mean_value = np.mean([row[c], r2[c], r3[c]])
        sample[c] = mean_value
    rnaformer_mean.append(sample)
rnaformer_mean_df = pd.DataFrame(rnaformer_mean)
rnaformer_mean_df['algorithm'] = 'RNAformer'
rnaformer_mean_df = rnaformer_mean_df[cols + ['algorithm']]

data = pd.concat([rnaformer_mean_df, data], ignore_index=True)

order = (
    data
    .groupby('algorithm')[score_of_interest]
    .mean()
    .sort_values(ascending=False)
    .index
    .tolist()
)

# 2. Make `algorithm` a categorical with that specific order
data['algorithm'] = pd.Categorical(data['algorithm'], categories=order, ordered=True)

# for algo, group in data.groupby('algorithm'):
#     print(algo, len(group), np.round(np.mean(group['f1_score']), 3), np.round(np.mean(group['mcc']), 3), np.round(np.mean(group['wl']), 3), np.round(np.mean(group['precision']), 3), np.round(np.mean(group['recall']), 3))
# 
# print(data)
plot_algorithm_distributions(
    df=data,
    algorithm_col='algorithm',
    value_col=score_of_interest,
    kind='box',  # 'violin', 'box',
    # palette=[],  # ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'],
    fig_size=(12, 6),
    # title='Secondary Structure F1 Score Comparison',
    # xlabel='Algorithm',
    ylabel='Matthews Correlation Coefficient',  # 'F1 Score',
    xtick_label_rotation=45,
    remove_spines=True,
    xlabel_size=18,
    ylabel_size=18,
    xtick_label_size=18,
    ytick_label_size=16,
    ref_algorithm='RNAformer',
    num_algos_for_statistics=2,
    test_method='wilcoxon',  # 'wilcoxon', 'ttest'
    annotate=True,
    # annot_offset=1.0,
    annot_step=0.07,
    step_factor=1.0,
    box_width=0.5,
    outpath=Path(outpath)
)

############################################################################################################
# AF3 vs RNAformer comparison for TS-Hard dataset
rf_preds_path = Path('evaluation/predictions/secondary_structure_rnaformer_eval', 'af3_eval', 'results_TS-Hard_all', 'finetune_model_long_test_samples_pdb_ts_hard_results.pkl')
af3_preds_path = Path('evaluation/predictions/secondary_structure_rnaformer_eval', 'af3_eval', 'results_TS-Hard_all', 'results_pdb_ts_hard.npy')
outpath = Path('plots', 'performance', 'secondary_structure_ts_hard_comparison_af3_rnaformer_box.svg')

af3_preds = np.load(af3_preds_path, allow_pickle=True)
af3_data = []
for pred in af3_preds:
    sample = {
        # 'pdb_id': pred['pdb_id'],
        'algorithm': 'AlphaFold 3',
        'f1_score': pred['f1_score'],
        'mcc': pred['mcc'],
        # 'wl': pred['wl'],
        'precision': pred['precision'],
        'recall': pred['recall'],
    }
    af3_data.append(sample)
af3_df = pd.DataFrame(af3_data)
af3_df.loc[:, 'algorithm'] = 'AlphaFold 3'

rf_df = pd.read_pickle(rf_preds_path)
rf_df.loc[:, 'algorithm'] = 'RNAformer'

# print(f"AF3 preds: {len(af3_preds)}, RF preds: {len(rf_preds)}")
print(rf_df)
print(af3_df)

data = pd.concat([rf_df, af3_df], ignore_index=True)

order = (
    data
    .groupby('algorithm')[score_of_interest]
    .mean()
    .sort_values(ascending=False)
    .index
    .tolist()
)

# 2. Make `algorithm` a categorical with that specific order
data['algorithm'] = pd.Categorical(data['algorithm'], categories=order, ordered=True)


# print(len(af3_preds), len(rf_preds))

plot_algorithm_distributions(
    df=data,
    algorithm_col='algorithm',
    value_col=score_of_interest,
    kind='box',  # 'violin', 'box',
    # palette=[],  # ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'],
    fig_size=(4, 6),
    # title='Secondary Structure F1 Score Comparison',
    # xlabel='Algorithm',
    ylabel='Matthews Correlation Coefficient',  # 'F1 Score',
    xtick_label_rotation=45,
    remove_spines=True,
    xlabel_size=24,
    ylabel_size=18,
    xtick_label_size=18,
    ytick_label_size=16,
    ref_algorithm='RNAformer',
    test_method='wilcoxon',  # 'wilcoxon', 'ttest'
    annotate=True,
    box_width=0.3,
    outpath=Path(outpath)
)


############################################################################################################
# AF3 vs RNAformer comparison for TS-PDB dataset

rf_preds_path = Path('evaluation/predictions/secondary_structure_rnaformer_eval', 'af3_eval', 'results_TS-PDB_all', 'finetune_model_long_test_samples_pdb_ts123_results.pkl')
af3_preds_path = Path('evaluation/predictions/secondary_structure_rnaformer_eval', 'af3_eval', 'results_TS-PDB_all', 'results_pdb_ts_pdb.npy')
outpath = Path('plots', 'secondary_structure_ts_pdb_comparison_af3_rnaformer_box.svg')

af3_preds = np.load(af3_preds_path, allow_pickle=True)
af3_data = []
for pred in af3_preds:
    sample = {
        # 'pdb_id': pred['pdb_id'],
        'algorithm': 'AlphaFold 3',
        'f1_score': pred['f1_score'],
        'mcc': pred['mcc'],
        # 'wl': pred['wl'],
        'precision': pred['precision'],
        'recall': pred['recall'],
    }
    af3_data.append(sample)
af3_df = pd.DataFrame(af3_data)
af3_df.loc[:, 'algorithm'] = 'AlphaFold 3'

rf_df = pd.read_pickle(rf_preds_path)
rf_df.loc[:, 'algorithm'] = 'RNAformer'

# print(f"AF3 preds: {len(af3_preds)}, RF preds: {len(rf_preds)}")
print(rf_df)
print(af3_df)

data = pd.concat([rf_df, af3_df], ignore_index=True)

order = (
    data
    .groupby('algorithm')[score_of_interest]
    .mean()
    .sort_values(ascending=False)
    .index
    .tolist()
)

# 2. Make `algorithm` a categorical with that specific order
data['algorithm'] = pd.Categorical(data['algorithm'], categories=order, ordered=True)


# print(len(af3_preds), len(rf_preds))

plot_algorithm_distributions(
    df=data,
    algorithm_col='algorithm',
    value_col=score_of_interest,
    kind='box',  # 'violin', 'box',
    # palette=[],  # ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'],
    fig_size=(4, 6),
    # title='Secondary Structure F1 Score Comparison',
    # xlabel='Algorithm',
    ylabel='Matthews Correlation Coefficient',  # 'F1 Score',
    xtick_label_rotation=45,
    remove_spines=True,
    xlabel_size=24,
    ylabel_size=18,
    xtick_label_size=18,
    ytick_label_size=16,
    ref_algorithm='RNAformer',
    test_method='wilcoxon',  # 'wilcoxon', 'ttest'
    annotate=True,
    box_width=0.3,
    outpath=Path(outpath)
)
