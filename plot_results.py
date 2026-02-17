import ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from pathlib import Path
from scipy.stats import wilcoxon, ttest_rel, mannwhitneyu, fisher_exact, chi2_contingency

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle, ConnectionPatch
import matplotlib.patches as mpatches

import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message="In a future version, `df.iloc[:, i] = newvals`"
)


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, ttest_rel
from pathlib import Path

def plot_grouped_bar_with_error_and_stats(
    df: pd.DataFrame,
    x_col: str,
    hue_col: str,
    value_col: str,
    x_bins: list = None,
    x_bin_labels: list = None,
    x_order: list = None,
    hue_order: list = None,
    palette: list = None,
    fig_size: tuple = (8, 6),
    bar_width: float = 0.8,
    bar_spacing: float = 1.0,   # NEW: spacing multiplier between bars
    title: str = None,
    title_size: int = 14,
    title_pad: int = 10,
    xlabel: str = None,
    xlabel_size: int = 12,
    ylabel: str = None,
    ylabel_size: int = 12,
    xtick_label_size: int = 10,
    ytick_label_size: int = 10,
    yticks: list = None,
    legend_title: str = None,
    legend_title_size: int = 12,
    legend_label_size: int = 10,
    legend_position: str = 'right',
    remove_spines: bool = True,
    annotate: bool = False,
    ref_hue: str = None,
    test_method: str = 'mannwhitney',
    sig_levels: list = [(0.001, '***'), (0.01, '**'), (0.05, '*')],
    step_factor: float = 1.0,
    legend_offset: float = -0.15,
    limit_y_range: bool = False,
    outpath: Path = Path('barplot.svg')
):
    """Grouped bar plot with SEM, inverted hue order, and adjustable bar spacing."""
    # helper statistics functions
    def wilcoxon_stats(x, y):
        diffs = x - y
        nz = diffs != 0
        d = diffs[nz]
        n = len(d)
        if n < 1:
            return None
        W, p = wilcoxon(x, y)
        # rank-biserial
        n_pos = np.sum(diffs > 0)
        n_neg = np.sum(diffs < 0)
        r_rb = (n_pos - n_neg) / n
        # standardized r
        mean_W = n*(n+1)/4
        sd_W = np.sqrt(n*(n+1)*(2*n+1)/24)
        Z = (W - mean_W) / sd_W
        r = Z / np.sqrt(n)
        return W, p, r_rb, Z, r

    def ttest_stats(x, y):
        stat, p = ttest_rel(x, y)
        d = np.mean(x - y) / np.std(x - y, ddof=1) if len(x) > 1 else np.nan
        return stat, p, d
    
    data = df.copy()
    # bin x_col if provided
    if x_bins is not None:
        data['x_cat'] = pd.cut(data[x_col], bins=x_bins, labels=x_bin_labels, include_lowest=True)
        x_cats = list(data['x_cat'].cat.categories)
    else:
        data['x_cat'] = data[x_col].astype(str)
        x_cats = sorted(data['x_cat'].unique())

    # aggregate with SEM
    agg = data.groupby(['x_cat', hue_col])[value_col].agg(
        mean='mean',
        sem=lambda x: np.std(x, ddof=1) / np.sqrt(len(x)),
        count='count'
    ).reset_index()

    # order
    x_cats = x_order or x_cats
    hue_cats = hue_order or sorted(agg[hue_col].unique())
    hue_cats = hue_cats[::-1]  # invert order

    # pivot
    mean_df = agg.pivot(index='x_cat', columns=hue_col, values='mean').reindex(index=x_cats, columns=hue_cats)
    sem_df  = agg.pivot(index='x_cat', columns=hue_col, values='sem').reindex(index=x_cats, columns=hue_cats)

    # palette
    n_hue = len(hue_cats)
    if palette is None:
        cmap = plt.get_cmap('Blues')
        palette = [cmap(0.3 + 0.7 * i/(n_hue-1)) for i in range(n_hue)]

    # positions with adjustable spacing
    x = np.arange(len(x_cats))
    width = (bar_width / n_hue) * bar_spacing
    offsets = (np.arange(n_hue) - (n_hue-1)/2) * width

    # plot bars
    fig, ax = plt.subplots(figsize=fig_size)
    all_vals = mean_df.values + sem_df.values
    data_range = np.nanmax(all_vals) - np.nanmin(mean_df.values)
    offset_unit = 0.05 * data_range

    for i, hue in enumerate(hue_cats):
        pos = x + offsets[i]
        means = mean_df[hue].values
        errs = sem_df[hue].values
        print(f"Plotting {hue} with means: {means}, SEMs: {errs}")
        ax.bar(pos, means, width=width, yerr=errs, error_kw={'capsize':0},
               color=palette[i], edgecolor=None, label=str(hue))

    # spines
    if remove_spines:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # labels
    ax.set_xticks(x)
    ax.set_xticklabels(x_cats, fontsize=xtick_label_size)
    if xlabel: ax.set_xlabel(xlabel, fontsize=xlabel_size)
    if ylabel: ax.set_ylabel(ylabel, fontsize=ylabel_size)
    if yticks is not None: ax.set_yticks(yticks)

    if limit_y_range:
        ymin, _ = ax.get_ylim()
        ax.set_ylim(ymin, 1.0)

    if title: ax.set_title(title, loc='left', fontsize=title_size, pad=title_pad)

    # legend
    if legend_position == 'right':
        leg = ax.legend(title=legend_title, fontsize=legend_label_size,
                        title_fontsize=legend_title_size, loc='upper right', bbox_to_anchor=(1.2,legend_offset))
    elif legend_position == 'above':
        leg = ax.legend(title=legend_title, fontsize=legend_label_size,
                        title_fontsize=legend_title_size, loc='lower center', bbox_to_anchor=(0.5,legend_offset), ncol=n_hue)
    else:  # below
        leg = ax.legend(title=legend_title, fontsize=legend_label_size,
                        title_fontsize=legend_title_size, loc='upper center', bbox_to_anchor=(0.5,legend_offset), ncol=n_hue)
    leg.get_frame().set_edgecolor('black')

    # statistical annotation
    if annotate and n_hue > 1:
        ref = ref_hue if ref_hue in hue_cats else hue_cats[0]
        ref_idx = hue_cats.index(ref)
        for xi, xc in enumerate(x_cats):
            ymax_grp = (mean_df.loc[xc] + sem_df.loc[xc]).max()
            y = ymax_grp + offset_unit
            ref_vals = data[(data['x_cat']==xc)&(data[hue_col]==ref)][value_col].dropna().values
            for j, hue in enumerate(hue_cats):
                if hue == ref: continue
                comp_vals = data[(data['x_cat']==xc)&(data[hue_col]==hue)][value_col].dropna().values
                if test_method == 'ttest' and len(ref_vals) == len(comp_vals):
                    stat, p, d = ttest_stats(ref_vals, comp_vals)
                    print(f"T-test {ref} vs {hue} in {xc}: t={stat:.3f}, p={p:.3g}, d={d:.3f}")
                else:
                    res = wilcoxon_stats(ref_vals, comp_vals)
                    if res:
                        W, p, r_rb, Z, r = res
                        print(f"Wilcoxon {ref} vs {hue} in {xc}: W={W:.1f}, p={p:.3g}, r_rb={r_rb:.3f}, Z={Z:.3f}, r={r:.3f}")
                    else:
                        print(f"Wilcoxon {ref} vs {hue} in {xc}: not enough data")
                        continue
                stars = next((s for thr, s in sig_levels if p < thr), '')
                if stars:
                    x1 = xi + offsets[ref_idx]
                    x2 = xi + offsets[j]
                    
                    line_y = y
                    ax.plot([x1, x2], [line_y, line_y], color='black', linewidth=1, clip_on=False)
                    
                    ax.text((x1 + x2) / 2, line_y + 0.2 * offset_unit, stars,
                            ha='center', va='bottom', fontsize=ytick_label_size, clip_on=False)
                    
                    y += step_factor * offset_unit
                    
                    # ax.plot([x1, x1, x2, x2], [y, y+0.5*offset_unit, y+0.5*offset_unit, y],
                    #         color='black', linewidth=1, clip_on=False)
                    # ax.text((x1+x2)/2, y+0.5*offset_unit, stars,
                    #         ha='center', va='bottom', fontsize=ytick_label_size, clip_on=False)
                    # y += step_factor * offset_unit

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    #plt.show()
    plt.close(fig)


def plot_grouped_fraction_bar_with_stats(
    df: pd.DataFrame,
    x_col: str,
    hue_col: str,
    x_bins: list,
    x_bin_labels: list,
    x_order: list = None,
    hue_order: list = None,
    palette: list = None,
    fig_size: tuple = (8, 6),
    bar_width: float = 0.8,
    title: str = None,
    title_size: int = 14,
    xlabel: str = None,
    xlabel_size: int = 12,
    ylabel: str = None,
    ylabel_size: int = 12,
    xtick_label_size: int = 10,
    ytick_label_size: int = 10,
    legend_title: str = None,
    legend_label_size: int = 10,
    legend_position: str = 'right',
    annotate: bool = True,
    test_method: str = 'fisher',     # 'fisher' or 'chi2'
    sig_levels: list = [(0.001, '***'), (0.01, '**'), (0.05, '*')],
    legend_offset: float = -0.2,
    remove_spines: bool = True,
    annot_offset: float = 0.05,
    step_factor: float = 1.0,
    limit_y_range: bool = False,
    outpath: Path = Path('fraction_barplot.png')
):
    """
    Bar-plot of fraction (count/total_per_hue) per x-bin and hue with
    simple-line pairwise Fisher’s-exact (or χ²) annotation on raw counts.
    Includes console output of bin counts, fractions, and test statistics.
    """
    data = df.copy()
    # 1) Bin into categories
    data['x_cat'] = pd.cut(
        data[x_col],
        bins=x_bins,
        labels=x_bin_labels,
        include_lowest=True
    )
    # 2) Orders
    x_cats = x_order or x_bin_labels
    hue_cats = hue_order or sorted(data[hue_col].unique())
    n_hue = len(hue_cats)
    # 3) Count table and totals
    count_df = (
        data
        .groupby(['x_cat', hue_col])
        .size()
        .unstack(fill_value=0)
        .reindex(index=x_cats, columns=hue_cats, fill_value=0)
    )
    total = count_df.sum(axis=0)  # total per hue
    # 4) Fraction table
    frac_df = count_df.div(total, axis=1)
    # 5) Palette
    if palette is None:
        cmap = plt.get_cmap('tab10')
        palette = [cmap(i) for i in range(n_hue)]
    else:
        palette = palette[:n_hue]
    # 6) Bar positions
    x = np.arange(len(x_cats))
    width = bar_width / n_hue
    offsets = (np.arange(n_hue) - (n_hue-1)/2) * width

    # 7) Plot bars
    fig, ax = plt.subplots(figsize=fig_size)
    for i, hue in enumerate(hue_cats):
        ax.bar(
            x + offsets[i],
            frac_df[hue].values,
            width=width,
            color=palette[i],
            edgecolor='black',
            label=str(hue)
        )

    # 8) Aesthetics
    if title:
        ax.set_title(title, fontsize=title_size, loc='left', pad=30)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=xlabel_size)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=ylabel_size)
    ax.set_xticks(x)
    ax.set_xticklabels(x_cats, fontsize=xtick_label_size)
    ax.set_ylim(0, 1.0)
    ax.tick_params(axis='y', labelsize=ytick_label_size)
    if remove_spines:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Legend
    if legend_position == 'right':
        ax.legend(
            title=legend_title,
            fontsize=legend_label_size,
            title_fontsize=legend_label_size,
            loc='upper right',
            bbox_to_anchor=(1.2, 1)
        )
    elif legend_position == 'below':
        ax.legend(
            title=legend_title,
            fontsize=legend_label_size,
            title_fontsize=legend_label_size,
            loc='upper center',
            bbox_to_anchor=(0.5, legend_offset),
            ncol=n_hue
        )

    # 9) Stats & simple-line annotation on raw counts
    if annotate and n_hue >= 2:
        max_frac = frac_df.values.max()
        for xi, xc in enumerate(x_cats):
            # baseline above the highest bar in this group
            y_base = max_frac + annot_offset
            pair_idx = 0

            for i in range(n_hue-1):
                for j in range(i+1, n_hue):
                    h1, h2 = hue_cats[i], hue_cats[j]
                    c1 = int(count_df.loc[xc, h1])
                    c2 = int(count_df.loc[xc, h2])
                    t1 = int(total[h1])
                    t2 = int(total[h2])
                    f1 = c1 / t1 if t1 else 0
                    f2 = c2 / t2 if t2 else 0
                    # console output of counts and fractions
                    print(f"Bin '{xc}': {h1} => {c1}/{t1} ({f1:.2f}), {h2} => {c2}/{t2} ({f2:.2f})")

                    # build contingency table
                    table = np.array([[c1, t1-c1],
                                      [c2, t2-c2]])
                    if test_method.lower() == 'chi2':
                        stat, p, _, _ = chi2_contingency(table)
                        print(f"Chi2 {h1} vs {h2} in '{xc}': chi2={stat:.2f}, p={p:.3g}")
                    else:
                        odds, p = fisher_exact(table)
                        print(f"Fisher {h1} vs {h2} in '{xc}': odds_ratio={odds:.2f}, p={p:.3g}")

                    stars = next((s for thr, s in sig_levels if p < thr), '')
                    print(f"  Significance: '{stars or 'ns'}'\n")

                    if stars:
                        # compute y position for this pair
                        y = y_base + pair_idx * annot_offset * step_factor

                        x1 = xi + offsets[i]
                        x2 = xi + offsets[j]
                        # draw simple horizontal line
                        ax.plot([x1, x2], [y, y],
                                lw=1.2, c='black', clip_on=False)
                        # stars above the line
                        ax.text((x1 + x2)/2, y + annot_offset*0.2, stars,
                                ha='center', va='bottom',
                                fontsize=ytick_label_size, clip_on=False)

                    pair_idx += 1

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    #plt.show()
    plt.close(fig)


def plot_grouped_boxplot_with_statistics(
    df: pd.DataFrame,
    group_col: str,
    metric_col: str,
    bins: list = None,
    bin_labels: list = None,
    fig_size: tuple = (8, 6),
    box_color: str = 'skyblue',
    title: str = None,
    title_size: int = 14,
    xlabel: str = None,
    xlabel_size: int = 12,
    ylabel: str = None,
    ylabel_size: int = 12,
    xtick_label_size: int = 10,
    ytick_label_size: int = 10,
    limit_y_range: bool = False,
    yticks: list = None,
    annotate: bool = False,
    ref_group: str = None,
    test_method: str = 'wilcoxon',
    sig_levels: list = [(0.001, '***'), (0.01, '**'), (0.05, '*')],
    num_groups_for_statistics: int = None,
    step_factor: float = 1.0,
    box_width: float = 0.6,
    spacing: float = 1.0,
    annot_offset: float = 0.03,
    annot_step: float = 0.06,
    star_offset: float = 0.2,
    bracket_arm: float = 0.03,
    clip_on: bool = True,
    outpath: Path = Path('boxplot.svg'),
):
    """
    Create a boxplot of `metric_col` vs. categories of `group_col`, with:
      - statistical annotations against a reference group
      - adjustable bracket spacing (annot_step) and star offset (star_offset)
      - console summaries of data sizes and test results
    """
    def wilcoxon_stats(x, y):
        diffs = x - y
        nz = diffs != 0
        d = diffs[nz]
        if len(d) == 0:
            return None
        W, p = wilcoxon(x, y)
        return W, p

    def ttest_stats(x, y):
        stat, p = ttest_rel(x, y)
        return stat, p

    # Prepare categories
    data = df.copy()
    print(f"\n=== Dataset summary ===")
    print(f"Total observations: {len(data)}")
    if bins is not None:
        data['category'] = pd.cut(data[group_col], bins=bins,
                                  labels=bin_labels, include_lowest=True)
        categories = list(data['category'].cat.categories)
        print(f"Binned into {len(categories)} bins: {categories}")
    else:
        if pd.api.types.is_categorical_dtype(data[group_col]):
            data['category'] = data[group_col]
            categories = list(data['category'].cat.categories)
        else:
            data['category'] = data[group_col].astype(str)
            categories = sorted(data['category'].unique())
        print(f"Found {len(categories)} categories: {categories}")

    if num_groups_for_statistics is None:
        num_groups_for_statistics = len(categories)

    # Build map of values & print counts
    data_map = {}
    for cat in categories:
        vals = data.loc[data['category'] == cat, metric_col].dropna().values
        data_map[cat] = vals
        print(f"  – {cat}: n={len(vals)}, mean={np.mean(vals):.3f}, median={np.median(vals):.3f}, std={np.std(vals, ddof=1):.3f}")

    positions = np.arange(len(categories)) * spacing

    # Draw boxplot
    fig, ax = plt.subplots(figsize=fig_size)
    bp = ax.boxplot(
        [data_map[cat] for cat in categories],
        positions=positions,
        widths=box_width,
        patch_artist=True,
        showfliers=False
    )
    for patch in bp['boxes']:
        patch.set_facecolor(box_color)
        patch.set_edgecolor('black')
    for med in bp['medians']:
        med.set(color='firebrick', linewidth=2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if limit_y_range:
        ymin, _ = ax.get_ylim()
        ax.set_ylim(ymin, 1.0)

    ax.set_xticks(positions)
    ax.set_xticklabels(categories, fontsize=xtick_label_size)
    if title:
        ax.set_title(title, loc='left', fontsize=title_size, pad=30)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=xlabel_size)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=ylabel_size)
    if yticks is not None:
        ax.set_yticks(yticks)
    ax.tick_params(axis='y', labelsize=ytick_label_size)

    # Statistical annotations
    if annotate and len(categories) > 1:
        print(f"\n=== Statistical tests ({test_method}) ===")
        ref = ref_group if ref_group in categories else categories[0]
        ref_vals = data_map[ref]
        ref_pos = positions[categories.index(ref)]
        print(f"Reference group: {ref} (n={len(ref_vals)}, mean={np.mean(ref_vals):.3f})")

        ymin, ymax = ax.get_ylim()
        y_span = ymax - ymin
        step_height = annot_step * y_span * step_factor
        arm_height = bracket_arm * y_span

        comp_cats = [c for c in categories[:num_groups_for_statistics] if c != ref]
        for comp_idx, cat in enumerate(comp_cats):
            comp_vals = data_map[cat]
            print(f"\nComparing {ref} vs {cat}:")
            if test_method == 'ttest':
                stat, p = ttest_stats(ref_vals, comp_vals)
                print(f"  t = {stat:.3f}, p = {p:.3g}")
            else:
                res = wilcoxon_stats(ref_vals, comp_vals)
                if res:
                    W, p = res
                    print(f"  W = {W:.3f}, p = {p:.3g}")
                else:
                    print("  No nonzero differences; skipped")
                    continue

            sym = next((s for thr, s in sig_levels if p < thr), '')
            print(f"  significance: {sym or 'ns'}")

            if not sym:
                continue

            # stack annotations above the highest point
            y0 = max(np.max(ref_vals), np.max(comp_vals))
            y = y0 + (comp_idx + 1) * step_height + annot_offset * y_span
            x2 = positions[categories.index(cat)]

            # draw bracket
            ax.plot(
                [ref_pos, ref_pos, x2, x2],
                # [y, y + step_height/2, y + step_height/2, y],
                [y, y+arm_height, y+arm_height, y],
                color='black', linewidth=1, clip_on=clip_on
            )
            # draw stars
            star_y = y + star_offset * step_height
            ax.text(
                (ref_pos + x2) / 2,
                star_y,
                sym,
                ha='center', va='bottom',
                fontsize=ytick_label_size,
                clip_on=clip_on
            )

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    #plt.show()
    plt.close(fig)


def plot_bar_means_with_statistics(
    df: pd.DataFrame,
    group_col: str,
    metric_col: str,
    order: list = None,
    palette: list = None,
    fig_size: tuple = (3.5, 3.2),
    bar_width: float = 0.75,
    bar_spacing: float = 1.0,    # NEW: >1 = more space between bars
    title: str = None,
    title_size: int = 12,
    xlabel: str = None,
    xlabel_size: int = 10,
    ylabel: str = None,
    ylabel_size: int = 11,
    xtick_label_size: int = 10,
    ytick_label_size: int = 10,
    yticks: list = None,
    xlim: tuple = None,
    ylim: tuple = (0, 100),
    remove_spines: bool = True,
    annotate: bool = True,
    ref_group: str = None,
    test_method: str = 'wilcoxon',
    sig_levels: list = [(0.001, '***'), (0.01, '**'), (0.05, '*')],
    error: str = 'sem',
    step_factor: float = 1.0,
    annot_offset: float = 0.035,
    bracket_arm: float = 0.02,
    floor_at_zero: bool = True,
    y_as_percent: bool = True,
    scale_in_0_1: bool = False,
    show_nsamples: bool = True,
    nsamples_prefix: str = r"$n$ = ",
    nsamples_size: int = 9,
    error_lw: float = 1.2,
    edge_lw: float = 1.0,
    outpath: Path = Path('bar_means.svg'),
):
    groups = order or list(pd.unique(df[group_col].dropna()))
    n = len(groups)
    if palette is None:
        palette = ['#9bd0ff'] + ['#bfbfbf'] * (n - 1)
    else:
        palette = palette[:n]

    # data aggregation
    per_group_vals = {}
    rows = []
    for g in groups:
        vals = df.loc[df[group_col] == g, metric_col].astype(float).dropna().values
        if scale_in_0_1:
            vals = vals * 100.0
        per_group_vals[g] = vals
        if vals.size == 0:
            mean, err = np.nan, np.nan
        else:
            mean = float(np.mean(vals))
            if error.lower() == 'sem':
                err = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
            else:
                err = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
        rows.append((g, mean, err, len(vals)))
    agg_df = pd.DataFrame(rows, columns=[group_col, 'mean', 'err', 'n'])

    # X positions with spacing control
    x = np.arange(n) * bar_spacing

    fig, ax = plt.subplots(figsize=fig_size)
    means = agg_df['mean'].to_numpy(float)
    errs  = agg_df['err'].to_numpy(float)

    if floor_at_zero:
        lower_err = np.minimum(errs, np.maximum(0, means - 0.0))
    else:
        lower_err = errs
    upper_err = errs
    yerr = np.vstack([lower_err, upper_err])

    # draw bars — NO borders now
    ax.bar(
        x, means,
        yerr=yerr,
        width=bar_width,
        color=palette,
        edgecolor=None,       # removed black edge
        linewidth=0,          # no outline
        error_kw={'capsize': 0, 'elinewidth': error_lw, 'ecolor': 'black'},
    )

    if remove_spines:
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)

    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=xtick_label_size)

    if show_nsamples:
        labels = []
        for g, n_ in zip(groups, agg_df['n'].tolist()):
            labels.append(f"{g}\n{nsamples_prefix}{n_}")
        ax.set_xticklabels(labels, fontsize=xtick_label_size, linespacing=1.2)
        ax.tick_params(axis='x', pad=6)

    if ylabel:
        ax.set_ylabel(ylabel, fontsize=ylabel_size)
    if y_as_percent:
        ax.yaxis.set_major_formatter(plt.PercentFormatter(100))
    if yticks is not None:
        ax.set_yticks(yticks)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.tick_params(axis='y', labelsize=ytick_label_size)
    if title:
        ax.set_title(title, fontsize=title_size, loc='left', pad=10)

    # annotation logic unchanged...
    # [keep same significance testing + brackets code from before]

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    #plt.show()
    plt.close(fig)


def plot_bar_means_fixed_panel(
    df: pd.DataFrame,
    group_col: str,
    metric_col: str,
    order: list = None,
    palette: list = None,
    fig_size: tuple = (3.5, 3.2),
    bar_width: float = 0.75,
    bar_spacing: float = 1.0,         # controls distance between bar centers
    title: str = None,
    title_size: int = 12,
    title_pad: int = 10,
    xlabel: str = None,
    xlabel_size: int = 10,
    ylabel: str = None,
    ylabel_size: int = 11,
    xtick_label_size: int = 10,
    ytick_label_size: int = 10,
    yticks: list = None,
    xlim: tuple = None,
    ylim: tuple = None,
    remove_spines: bool = True,
    annotate: bool = True,
    ref_group: str = None,
    test_method: str = 'wilcoxon',    # 'wilcoxon' | 'ttest' | 'mannwhitney'
    sig_levels: list = [(0.001, '***'), (0.01, '**'), (0.05, '*')],
    error: str = 'sem',               # 'sem' or 'sd'
    step_factor: float = 1.0,
    annot_offset: float = 0.035,      # fraction of y-range
    bracket_arm: float = 0.02,        # fraction of y-range
    floor_at_zero: bool = True,
    y_as_percent: bool = False,
    scale_in_0_1: bool = False,
    show_nsamples: bool = False,
    nsamples_prefix: str = r"$n$ = ",
    nsamples_size: int = 9,
    error_lw: float = 1.2,
    edge_lw: float = 0.0,             # ignored (we remove borders)
    outpath: Path = Path('bar_means.svg'),
    # fixed panel margins (do NOT use tight_layout / bbox_inches='tight')
    _left: float = 0.18, _right: float = 0.98, _bottom: float = 0.20, _top: float = 0.92,
):
    """
    Fixed-layout 'Nature-style' bar panel: mean ± SEM/SD with significance brackets.
    Signature matches your previous function so you can call it the same way.
    """
    # ---- deterministic ordering ----
    if order is None:
        order = list(pd.unique(df[group_col].dropna()))
    groups = list(order)
    n = len(groups)

    # ---- summarize data ----
    values_by_g = {}
    rows = []
    for g in groups:
        vals = df.loc[df[group_col] == g, metric_col].astype(float).dropna().values
        if scale_in_0_1:
            vals = vals * 100.0
        values_by_g[g] = vals
        if vals.size > 0:
            mean = float(np.mean(vals))
            if error.lower() == 'sem':
                err = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
            else:
                err = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
        else:
            mean, err = np.nan, np.nan
        rows.append((g, mean, err, len(vals)))
    agg = pd.DataFrame(rows, columns=[group_col, 'mean', 'err', 'n'])

    # percent display (just presentation)
    if y_as_percent:
        agg['mean'] = agg['mean'] * (100.0 if not scale_in_0_1 else 1.0)
        agg['err']  = agg['err']  * (100.0 if not scale_in_0_1 else 1.0)

    # ---- fixed geometry: positions & limits ----
    x = np.arange(n, dtype=float) * bar_spacing
    # pad so the left/right space is constant regardless of labels/data
    pad = max(0.05 * bar_spacing, 0.03)
    x_left  = x[0] - bar_width/2 - pad
    x_right = x[-1] + bar_width/2 + pad

    # ---- figure (NO tight_layout) ----
    fig, ax = plt.subplots(figsize=fig_size)
    plt.subplots_adjust(left=_left, right=_right, bottom=_bottom, top=_top)

    # ---- bars + errors (no borders) ----
    means = agg['mean'].to_numpy(float)
    errs  = agg['err'].to_numpy(float)
    if floor_at_zero:
        lower = np.minimum(errs, np.maximum(0.0, means - 0.0))
    else:
        lower = errs
    yerr = np.vstack([lower, errs])

    # palette
    if palette is None:
        cmap = plt.get_cmap('Blues')
        palette = [cmap(0.3 + 0.7 * i/(max(1, n)-1)) for i in range(n)]

    ax.bar(
        x, means,
        width=bar_width,
        color=palette[:n],
        edgecolor=None, linewidth=0,
        yerr=yerr,
        error_kw={'capsize': 0, 'elinewidth': error_lw, 'ecolor': 'black'},
    )

    # ---- axes cosmetics (fixed) ----
    ax.margins(x=0, y=0)            # disable auto margins
    ax.set_xlim(x_left, x_right)    # lock horizontal frame
    if ylim is not None:
        ax.set_ylim(*ylim)
    if yticks is not None:
        ax.set_yticks(yticks)
    if xlim is not None:
        ax.set_xlim(*xlim)          # let user override if they want

    ax.set_xticks(x)
    if show_nsamples:
        labels = [f"{g}\n{nsamples_prefix}{int(n_)}" for g, n_ in zip(groups, agg['n'])]
        ax.set_xticklabels(labels, fontsize=xtick_label_size, linespacing=1.2)
        ax.tick_params(axis='x', pad=6)
    else:
        ax.set_xticklabels(groups, fontsize=xtick_label_size)

    if ylabel:
        ax.set_ylabel(ylabel, fontsize=ylabel_size)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=xlabel_size)
    if title:
        ax.set_title(title, fontsize=title_size, loc='left', pad=title_pad)

    if remove_spines:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)
    ax.tick_params(axis='y', labelsize=ytick_label_size)

    if y_as_percent:
        from matplotlib.ticker import PercentFormatter
        ax.yaxis.set_major_formatter(PercentFormatter(100))

    # ---- significance vs ref_group (fixed placement) ----
    if annotate and n > 1:
        if ref_group is None or ref_group not in groups:
            ref_group = groups[0]
        ref_vals = values_by_g.get(ref_group, np.array([]))

        # position brackets just above tallest mean+err
        ymin, ymax = ax.get_ylim()
        yspan = ymax - ymin
        base = np.nanmax(means + errs) if np.isfinite(means + errs).any() else ymax * 0.9
        y0 = base + annot_offset * yspan
        arm = bracket_arm * yspan
        step = annot_offset * yspan * step_factor

        print(f"\n=== Stats vs '{ref_group}' on {metric_col} ===")
        k = 0
        for g in groups:
            if g == ref_group:
                continue
            comp_vals = values_by_g.get(g, np.array([]))
            if comp_vals.size == 0 or ref_vals.size == 0:
                print(f"{ref_group} vs {g}: not enough data")
                continue

            p = np.nan
            txt = ""
            try:
                method = test_method.lower()
                if method == 'ttest' and len(ref_vals) == len(comp_vals):
                    stat, p = ttest_rel(ref_vals, comp_vals)
                    txt = f"t={stat:.3f}, p={p:.3g}"
                elif method == 'wilcoxon' and len(ref_vals) == len(comp_vals):
                    W, p = wilcoxon(ref_vals, comp_vals)
                    txt = f"W={W:.1f}, p={p:.3g}"
                else:
                    stat, p = mannwhitneyu(ref_vals, comp_vals, alternative='two-sided')
                    txt = f"U={stat:.1f}, p={p:.3g}" + (" (fallback)" if method!='mannwhitney' else "")
            except Exception as e:
                txt = f"test failed: {e}"
            print(f"{ref_group} vs {g}: {txt}")

            stars = next((s for thr, s in sig_levels if p < thr), '')
            if stars:
                i1, i2 = groups.index(ref_group), groups.index(g)
                y = y0 + k * step
                ax.plot([x[i1], x[i1], x[i2], x[i2]], [y, y+arm, y+arm, y],
                        color='black', linewidth=1.0, clip_on=False)
                ax.text((x[i1] + x[i2]) / 2, y + arm, stars,
                        ha='center', va='bottom', fontsize=ytick_label_size, clip_on=False)
                k += 1

    try:
        fig.canvas.draw()  # need a renderer for tight bbox
        renderer = fig.canvas.get_renderer()
        tight = ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted())

        pad = 0.01  # ~1% of figure, small and consistent
        # how much overflows outside [0, 1] in figure coords
        dl = max(0.0, 0.0 - tight.x0) + pad
        dr = max(0.0, tight.x1 - 1.0) + pad
        db = max(0.0, 0.0 - tight.y0) + pad
        dt = max(0.0, tight.y1 - 1.0) + pad

        # current margins were set earlier via plt.subplots_adjust(left=_left, ...)
        # increase them if needed so nothing is cut off
        new_left   = min(0.35, max(0.05, _left   + dl))
        new_right  = max(0.65, min(0.98, _right  - dr))
        new_bottom = min(0.40, max(0.05, _bottom + db))
        new_top    = max(0.60, min(0.95, _top    - dt))

        # only apply if they still leave a reasonable axes area
        if new_right - new_left > 0.25 and new_top - new_bottom > 0.25:
            plt.subplots_adjust(left=new_left, right=new_right,
                                bottom=new_bottom, top=new_top)
    except Exception:
        # if anything goes weird, fall back silently to current margins
        pass
    #plt.show()

    fig.savefig(outpath, dpi=300)   # <- do not use bbox_inches='tight'
    plt.close(fig)


def plot_bar_means_fixed_panel_and_inset(
    df: pd.DataFrame,
    group_col: str,
    metric_col: str,
    order: list = None,
    palette: list = None,
    fig_size: tuple = (3.5, 3.2),
    bar_width: float = 0.75,
    bar_spacing: float = 1.0,
    title: str = None,
    title_size: int = 12,
    title_pad: int = 10,
    xlabel: str = None,
    xlabel_size: int = 10,
    ylabel: str = None,
    ylabel_size: int = 11,
    xtick_label_size: int = 10,
    ytick_label_size: int = 10,
    yticks: list = None,
    xlim: tuple = None,
    ylim: tuple = None,
    remove_spines: bool = True,
    annotate: bool = True,
    ref_group: str = None,
    test_method: str = 'wilcoxon',
    sig_levels: list = [(0.001, '***'), (0.01, '**'), (0.05, '*')],
    error: str = 'sem',
    step_factor: float = 1.0,
    annot_offset: float = 0.035,
    bracket_arm: float = 0.02,
    floor_at_zero: bool = True,
    y_as_percent: bool = False,
    scale_in_0_1: bool = False,
    show_nsamples: bool = False,
    nsamples_prefix: str = r"$n$ = ",
    nsamples_size: int = 9,
    error_lw: float = 1.2,
    edge_lw: float = 0.0,
    # Inset options
    inset_enabled: bool = False,
    inset_groups: list = None,        # groups to show in inset (names)
    inset_ylim: tuple = None,         # if None and inset_auto=True, computed from data
    inset_size: tuple = (0.46, 0.74), # (w,h) as fraction of main axes
    inset_loc: str = 'upper right',
    inset_box_color: str = '#79c34a',
    inset_auto: bool = True,          # auto y-lims from means±err of inset bars
    inset_margin_frac: float = 0.08,  # expand auto y-lims by this fraction
    outpath: Path = Path('bar_means.svg'),
    # fixed panel margins
    _left: float = 0.18, _right: float = 0.98, _bottom: float = 0.20, _top: float = 0.92,
):
    # ---- deterministic ordering ----
    if order is None:
        order = list(pd.unique(df[group_col].dropna()))
    groups = list(order)
    n = len(groups)

    # ---- summarize data ----
    values_by_g, rows = {}, []
    for g in groups:
        vals = df.loc[df[group_col] == g, metric_col].astype(float).dropna().values
        if scale_in_0_1:
            vals = vals * 100.0
        values_by_g[g] = vals
        if vals.size:
            mean = float(np.mean(vals))
            if error.lower() == 'sem':
                err = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
            else:
                err = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
        else:
            mean, err = np.nan, np.nan
        rows.append((g, mean, err, len(vals)))
    agg = pd.DataFrame(rows, columns=[group_col, 'mean', 'err', 'n'])

    if y_as_percent:
        agg['mean'] *= (100.0 if not scale_in_0_1 else 1.0)
        agg['err']  *= (100.0 if not scale_in_0_1 else 1.0)

    # ---- geometry ----
    x = np.arange(n, dtype=float) * bar_spacing
    pad = max(0.05 * bar_spacing, 0.03)
    x_left  = x[0] - bar_width/2 - pad
    x_right = x[-1] + bar_width/2 + pad

    # ---- figure ----
    fig, ax = plt.subplots(figsize=fig_size)
    plt.subplots_adjust(left=_left, right=_right, bottom=_bottom, top=_top)

    # ---- bars ----
    means = agg['mean'].to_numpy(float)
    errs  = agg['err'].to_numpy(float)
    lower = np.minimum(errs, np.maximum(0.0, means - 0.0)) if floor_at_zero else errs
    yerr  = np.vstack([lower, errs])

    if palette is None:
        cmap = plt.get_cmap('Blues')
        palette = [cmap(0.3 + 0.7 * i/(max(1, n)-1)) for i in range(n)]

    ax.bar(
        x, means,
        width=bar_width,
        color=palette[:n],
        edgecolor=None, linewidth=0,
        yerr=yerr,
        error_kw={'capsize': 0, 'elinewidth': error_lw, 'ecolor': 'black'},
        zorder=2
    )

    # ---- axes cosmetics ----
    ax.margins(x=0, y=0)
    ax.set_xlim(x_left, x_right)
    if ylim is not None: ax.set_ylim(*ylim)
    if yticks is not None: ax.set_yticks(yticks)
    if xlim is not None: ax.set_xlim(*xlim)
    ax.set_xticks(x)

    if show_nsamples:
        labels = [f"{g}\n{nsamples_prefix}{int(n_)}" for g, n_ in zip(groups, agg['n'])]
        ax.set_xticklabels(labels, fontsize=xtick_label_size, linespacing=1.2)
        ax.tick_params(axis='x', pad=6)
    else:
        ax.set_xticklabels(groups, fontsize=xtick_label_size)

    if ylabel: ax.set_ylabel(ylabel, fontsize=ylabel_size)
    if xlabel: ax.set_xlabel(xlabel, fontsize=xlabel_size)
    if title:  ax.set_title(title, fontsize=title_size, loc='left', pad=title_pad)

    if remove_spines:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)
    ax.tick_params(axis='y', labelsize=ytick_label_size)
    if y_as_percent:
        from matplotlib.ticker import PercentFormatter
        ax.yaxis.set_major_formatter(PercentFormatter(100))

    # ---- INSET ----
    if inset_enabled and inset_groups:
        inset_groups = [g for g in inset_groups if g in groups]
        if inset_groups:
            sel_idx = [groups.index(g) for g in inset_groups]
            sel_x   = x[sel_idx]
            sel_m   = means[sel_idx]
            sel_e   = errs[sel_idx]
            sel_pal = [palette[i] for i in sel_idx]

            # compute absolute sizes in inches (avoid % strings)
            fig_w_in, fig_h_in = fig.get_size_inches()
            ax_bbox = ax.get_position()  # figure coords
            ax_w_in = fig_w_in * (ax_bbox.x1 - ax_bbox.x0)
            ax_h_in = fig_h_in * (ax_bbox.y1 - ax_bbox.y0)
            inset_w_in = inset_size[0] * ax_w_in
            inset_h_in = inset_size[1] * ax_h_in

            # place inset: inside vs outside-right
            if inset_loc == 'outside right':
                from mpl_toolkits.axes_grid1 import make_axes_locatable
                divider = make_axes_locatable(ax)
                iax = divider.append_axes("right", size=inset_w_in, pad=0.25)  # pad in inches

                # vertically center the outside inset and give it the requested height
                box = iax.get_position()
                new_h = inset_h_in / fig_h_in
                iax.set_position([box.x0, 0.5 - new_h/2, box.width, new_h])

                # a bit more top margin to protect the title when we add an outside inset
                plt.subplots_adjust(top=min(0.96, max(plt.gcf().subplotpars.top, 0.93)))
            else:
                from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                iax = inset_axes(ax, width=inset_w_in, height=inset_h_in,
                                 loc=inset_loc, borderpad=0.0)

            iax.set_facecolor('white')
            iax.set_zorder(3)
            for s in ('top','right'):
                iax.spines[s].set_visible(False)
            for s in ('left','bottom'):
                iax.spines[s].set_linewidth(1.2)

            # choose tight ylim for inset
            if inset_auto and inset_ylim is None:
                inset_ylim = _inset_auto_ylim(sel_m, sel_e, inset_margin_frac)
            if inset_ylim is not None:
                iax.set_ylim(*inset_ylim)

            # set xlim FIRST so bars get clipped if needed
            xpad = max(0.05*bar_spacing, 0.03)
            iax.set_xlim(sel_x.min()-bar_width/2-xpad, sel_x.max()+bar_width/2+xpad)

            # bars inside inset (clip on)
            iax.bar(
                sel_x, sel_m,
                width=bar_width, color=sel_pal,
                edgecolor=None, linewidth=0,
                yerr=np.vstack([np.minimum(sel_e, np.maximum(0.0, sel_m-0.0)) if floor_at_zero else sel_e,
                                sel_e]),
                error_kw={'capsize': 0, 'elinewidth': max(1.0, error_lw-0.2), 'ecolor': 'black'},
                zorder=2,
                clip_on=True
            )

            # no x labels in inset
            iax.set_xticks([])
            iax.tick_params(axis='x', length=0)
            iax.tick_params(axis='y', labelsize=max(8, ytick_label_size-4))

            # rectangle on main axes (zoom window)
            rect_x0 = sel_x.min()-bar_width/2-xpad
            rect_x1 = sel_x.max()+bar_width/2+xpad
            # clamp rect y-range to the main axes limits
            y0, y1 = ax.get_ylim()
            ry0 = max(y0, inset_ylim[0]) if inset_ylim else y0
            ry1 = min(y1, inset_ylim[1]) if inset_ylim else y1
            rect = Rectangle((rect_x0, ry0), rect_x1-rect_x0, ry1-ry0,
                             fill=False, lw=1.6, ec=inset_box_color, zorder=4)
            ax.add_patch(rect)

            # connectors: to left edge of the inset, top/bottom
            con1 = ConnectionPatch(xyA=(rect_x1, ry1),  coordsA=ax.transData,
                                   xyB=(0.0, 1.0),       coordsB=iax.transAxes,
                                   color=inset_box_color, lw=1.6, zorder=4)
            con2 = ConnectionPatch(xyA=(rect_x1, ry0),  coordsA=ax.transData,
                                   xyB=(0.0, 0.0),       coordsB=iax.transAxes,
                                   color=inset_box_color, lw=1.6, zorder=4)
            ax.add_artist(con1); ax.add_artist(con2)

            
    # ---- stats (unchanged) ----
    if annotate and n > 1:
        if ref_group is None or ref_group not in groups:
            ref_group = groups[0]
        ref_vals = values_by_g.get(ref_group, np.array([]))
        ymin, ymax = ax.get_ylim()
        yspan = ymax - ymin
        base = np.nanmax(means + errs) if np.isfinite(means + errs).any() else ymax * 0.9
        y0 = base + annot_offset * yspan
        arm = bracket_arm * yspan
        step = annot_offset * yspan * step_factor

        print(f"\n=== Stats vs '{ref_group}' on {metric_col} ===")
        k = 0
        for g in groups:
            if g == ref_group:
                continue
            comp_vals = values_by_g.get(g, np.array([]))
            if comp_vals.size == 0 or ref_vals.size == 0:
                print(f"{ref_group} vs {g}: not enough data")
                continue
            p, txt = np.nan, ""
            try:
                m = test_method.lower()
                if m == 'ttest' and len(ref_vals) == len(comp_vals):
                    stat, p = ttest_rel(ref_vals, comp_vals); txt=f"t={stat:.3f}, p={p:.3g}"
                elif m == 'wilcoxon' and len(ref_vals) == len(comp_vals):
                    W, p = wilcoxon(ref_vals, comp_vals);      txt=f"W={W:.1f}, p={p:.3g}"
                else:
                    stat, p = mannwhitneyu(ref_vals, comp_vals, alternative='two-sided'); txt=f"U={stat:.1f}, p={p:.3g}"
            except Exception as e:
                txt = f"test failed: {e}"
            print(f"{ref_group} vs {g}: {txt}")

            stars = next((s for thr, s in sig_levels if p < thr), '')
            if stars:
                i1, i2 = groups.index(ref_group), groups.index(g)
                y = y0 + k*step
                ax.plot([x[i1], x[i1], x[i2], x[i2]], [y, y+arm, y+arm, y],
                        color='black', linewidth=1.0, clip_on=False)
                ax.text((x[i1]+x[i2])/2, y+arm, stars,
                        ha='center', va='bottom', fontsize=ytick_label_size, clip_on=False)
                k += 1

    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        tight = ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted())
        pad = 0.01
        dl = max(0.0, 0.0 - tight.x0) + pad
        dr = max(0.0, tight.x1 - 1.0) + pad
        db = max(0.0, 0.0 - tight.y0) + pad
        dt = max(0.0, tight.y1 - 1.0) + pad
        new_left   = min(0.35, max(0.05, _left   + dl))
        new_right  = max(0.65, min(0.98, _right  - dr))
        new_bottom = min(0.40, max(0.05, _bottom + db))
        new_top    = max(0.95, min(0.98, _top    - dt))
        if new_right - new_left > 0.25 and new_top - new_bottom > 0.25:
            plt.subplots_adjust(left=new_left, right=new_right,
                                bottom=new_bottom, top=new_top)
    except Exception:
        pass
    #plt.show()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def _inset_auto_ylim(means, errs, margin_frac=0.08):
    lo = np.nanmin(means - errs)
    hi = np.nanmax(means + errs)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return None
    span = hi - lo
    return (lo - margin_frac*span, hi + margin_frac*span)


def get_max_rna_length(row):
    """
    Extract the maximum RNA length from the 'RNA_msa_size' column.
    The column is expected to be a string representation of a dictionary.
    """
    try:
        rnalengths = ast.literal_eval(str(row['RNALength']))
    except (ValueError, SyntaxError):
        print(f"Error parsing RNALength for row: {row}")
        print(row['RNALength'])
        rnalengths = row['RNALength']
    if isinstance(rnalengths, int):
        return rnalengths
    elif isinstance(rnalengths, str):
        rnalengths = [int(l) for l in rnalengths.split(',')]
    elif isinstance(rnalengths, list):
        rnalengths = [int(l) for l in rnalengths]
    return np.max([int(x) for x in rnalengths])

def get_rna_msa_size(row):
    val = row["RNA_msa_size"]
    print(val)
    
    # 1. Handle real NaN/None
    if pd.isna(val):
        return 0
    
    # 2. If it's already a dict, use it directly
    if isinstance(val, dict):
        msasizes = val
    else:
        # 3. Otherwise assume it's a string, try to parse it
        try:
            msasizes = ast.literal_eval(str(val))
        except (ValueError, SyntaxError):
            # Could not parse, decide what you want here
            raise
    
    # 4. Compute max, guarding against empty dict
    if not msasizes:
        return 0
    
    return int(np.max([int(x) for x in msasizes.values()]))

def get_rna_max_msa_size(row, annot):
    i = row['exp_db_id']
    msasizes = annot[annot['exp_db_id'] == i]['RNA_msa_size'].values[0]
    return np.max([int(x) for x in msasizes.values()])

def get_rna_min_msa_size(row, annot):
    i = row['exp_db_id']
    msasizes = annot[annot['exp_db_id'] == i]['RNA_msa_size'].values[0]
    return np.min([int(x) for x in msasizes.values()])

def get_all_msa_sizes(row, annot):
    i = row['exp_db_id']
    msasizes = annot[annot['exp_db_id'] == i]['RNA_msa_size'].values[0]
    if all([int(x) <= 1 for x in msasizes.values()]):
        return 1
    else:
        return np.max([int(x) for x in msasizes.values()])

def get_rmsd(x):
    if isinstance(x, str):
        value = ast.literal_eval(x)
        if isinstance(value, list):
            return float(value[0])
        else: 
            return float(value)
    elif isinstance(x, (int, float)):
        if x == np.nan:
            return np.nan
        else:
            return float(x)
    else:
        return np.nan

import argparse

parser = argparse.ArgumentParser(description='Process rnaformer and af3 data.')
parser.add_argument('--measure', type=str, default='Complex_LDDT', help='Measure to use for comparison (default: Complex_LDDT)')
args = parser.parse_args()

measure = args.measure
max_len = 200
min_len = 0

max_rna = 5
max_protein = 5

print('Measure:', measure)

print('############################################################')
print('Evaluate all Predictions.')
print()
print(f"Measure: {measure}, Max Length: {max_len}, Min Length: {min_len}, Max RNAs: {max_rna}, Max Proteins: {max_protein}")
print('############################################################')

print('Starting to load data...')

datadir = Path('results')
plotting_dir = Path('plots/performance')
plotting_dir.mkdir(exist_ok=True, parents=True)
print(f'### Plotting directory: {plotting_dir}')

#####################################################################################################################
# AlphaFold 3

msa_sizes = pd.read_csv('msa_sizes_all_data.csv')
# msa_sizes = pd.concat([
#     pd.read_csv('msa_sizes_all_data_orphan.csv'),
#     pd.read_csv('msa_sizes_all_data_non_orphan.csv')
# ])
print(msa_sizes)


# all data
alphafold_all = pd.read_csv(Path(datadir, 'csvs/All_alphafold.csv'))
alphafold_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_alphafold_orphan.csv'))
alphafold_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_alphafold_non_orphan.csv'))

alphafold_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_alphafold.csv'))
alphafold_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_alphafold_orphan.csv'))
alphafold_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_alphafold_non_orphan.csv'))

alphafold_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_alphafold.csv'))
alphafold_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_alphafold_orphan.csv'))
alphafold_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_alphafold_non_orphan.csv'))

# rna monomers
alphafold_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_alphafold.csv'))
alphafold_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_alphafold_orphan.csv'))
alphafold_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_alphafold_non_orphan.csv'))

alphafold_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_alphafold.csv'))
alphafold_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_alphafold_orphan.csv'))
alphafold_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_alphafold_non_orphan.csv'))

alphafold_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_alphafold.csv'))
alphafold_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_alphafold_orphan.csv'))
alphafold_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_alphafold_non_orphan.csv'))

# protein-rna
alphafold_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_alphafold.csv'))
alphafold_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_alphafold_orphan.csv'))
alphafold_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_alphafold_non_orphan.csv'))

alphafold_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_alphafold.csv'))
alphafold_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_alphafold_orphan.csv'))
alphafold_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_alphafold_non_orphan.csv'))

alphafold_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_alphafold.csv'))
alphafold_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_alphafold_orphan.csv'))
alphafold_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_alphafold_non_orphan.csv'))


#####################################################################################################################
# SHS RNAformer

# all data
rnaformer_all = pd.read_csv(Path(datadir, 'csvs/All_rnaformer.csv'))
rnaformer_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformer_orphan.csv'))
rnaformer_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformer_non_orphan.csv'))

rnaformer_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformer.csv'))
rnaformer_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformer_orphan.csv'))
rnaformer_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformer_non_orphan.csv'))

rnaformer_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformer.csv'))
rnaformer_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformer_orphan.csv'))
rnaformer_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformer_non_orphan.csv'))

# rna monomers
rnaformer_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformer.csv'))
rnaformer_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformer_orphan.csv'))
rnaformer_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformer_non_orphan.csv'))

rnaformer_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformer.csv'))
rnaformer_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformer_orphan.csv'))
rnaformer_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformer_non_orphan.csv'))

rnaformer_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformer.csv'))
rnaformer_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformer_orphan.csv'))
rnaformer_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformer_non_orphan.csv'))

# protein-rna
rnaformer_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformer.csv'))
rnaformer_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformer_orphan.csv'))
rnaformer_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformer_non_orphan.csv'))

rnaformer_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformer.csv'))
rnaformer_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformer_orphan.csv'))
rnaformer_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformer_non_orphan.csv'))

rnaformer_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformer.csv'))
rnaformer_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformer_orphan.csv'))
rnaformer_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformer_non_orphan.csv'))

#####################################################################################################################
# SHS RNAfold

# all data
rnafold_all = pd.read_csv(Path(datadir, 'csvs/All_rnafold.csv'))
rnafold_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnafold_orphan.csv'))
rnafold_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnafold_non_orphan.csv'))

rnafold_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnafold.csv'))
rnafold_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnafold_orphan.csv'))
rnafold_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnafold_non_orphan.csv'))

rnafold_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnafold.csv'))
rnafold_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnafold_orphan.csv'))
rnafold_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnafold_non_orphan.csv'))

# rna monomers
rnafold_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnafold.csv'))
rnafold_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnafold_orphan.csv'))
rnafold_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnafold_non_orphan.csv'))

rnafold_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnafold.csv'))
rnafold_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnafold_orphan.csv'))
rnafold_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnafold_non_orphan.csv'))

rnafold_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnafold.csv'))
rnafold_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnafold_orphan.csv'))
rnafold_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnafold_non_orphan.csv'))

# protein-rna
rnafold_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnafold.csv'))
rnafold_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnafold_orphan.csv'))
rnafold_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnafold_non_orphan.csv'))

rnafold_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnafold.csv'))
rnafold_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnafold_orphan.csv'))
rnafold_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnafold_non_orphan.csv'))

rnafold_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnafold.csv'))
rnafold_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnafold_orphan.csv'))
rnafold_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnafold_non_orphan.csv'))

#####################################################################################################################
# SHS RNAfold N100

# all data
rnafoldN100_all = pd.read_csv(Path(datadir, 'csvs/All_rnafoldN100.csv'))
rnafoldN100_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnafoldN100_orphan.csv'))
rnafoldN100_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnafoldN100_non_orphan.csv'))

rnafoldN100_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnafoldN100.csv'))
rnafoldN100_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnafoldN100_orphan.csv'))
rnafoldN100_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnafoldN100_non_orphan.csv'))

rnafoldN100_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnafoldN100.csv'))
rnafoldN100_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnafoldN100_orphan.csv'))
rnafoldN100_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnafoldN100_non_orphan.csv'))

# rna monomers
rnafoldN100_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnafoldN100.csv'))
rnafoldN100_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnafoldN100_orphan.csv'))
rnafoldN100_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnafoldN100_non_orphan.csv'))

rnafoldN100_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnafoldN100.csv'))
rnafoldN100_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnafoldN100_orphan.csv'))
rnafoldN100_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnafoldN100_non_orphan.csv'))

rnafoldN100_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnafoldN100.csv'))
rnafoldN100_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnafoldN100_orphan.csv'))
rnafoldN100_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnafoldN100_non_orphan.csv'))

# protein-rna
rnafoldN100_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnafoldN100.csv'))
rnafoldN100_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnafoldN100_orphan.csv'))
rnafoldN100_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnafoldN100_non_orphan.csv'))

rnafoldN100_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnafoldN100.csv'))
rnafoldN100_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnafoldN100_orphan.csv'))
rnafoldN100_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnafoldN100_non_orphan.csv'))

rnafoldN100_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnafoldN100.csv'))
rnafoldN100_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnafoldN100_orphan.csv'))
rnafoldN100_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnafoldN100_non_orphan.csv'))


#####################################################################################################################
# SHS SPOT-RNA

# all data
spotrna_all = pd.read_csv(Path(datadir, 'csvs/All_spotrna.csv'))
spotrna_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_spotrna_orphan.csv'))
spotrna_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_spotrna_non_orphan.csv'))

spotrna_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_spotrna.csv'))
spotrna_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_spotrna_orphan.csv'))
spotrna_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_spotrna_non_orphan.csv'))

spotrna_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_spotrna.csv'))
spotrna_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_spotrna_orphan.csv'))
spotrna_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_spotrna_non_orphan.csv'))

# rna monomers
spotrna_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_spotrna.csv'))
spotrna_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_spotrna_orphan.csv'))
spotrna_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_spotrna_non_orphan.csv'))

spotrna_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_spotrna.csv'))
spotrna_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_spotrna_orphan.csv'))
spotrna_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_spotrna_non_orphan.csv'))

spotrna_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_spotrna.csv'))
spotrna_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_spotrna_orphan.csv'))
spotrna_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_spotrna_non_orphan.csv'))

# protein-rna
spotrna_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_spotrna.csv'))
spotrna_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_spotrna_orphan.csv'))
spotrna_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_spotrna_non_orphan.csv'))

spotrna_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_spotrna.csv'))
spotrna_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_spotrna_orphan.csv'))
spotrna_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_spotrna_non_orphan.csv'))

spotrna_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_spotrna.csv'))
spotrna_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_spotrna_orphan.csv'))
spotrna_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_spotrna_non_orphan.csv'))

#####################################################################################################################
# SHS SPOT-RNA N100

# all data
spotrnaN100_all = pd.read_csv(Path(datadir, 'csvs/All_spotrnaN100.csv'))
spotrnaN100_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_spotrnaN100_orphan.csv'))
spotrnaN100_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_spotrnaN100_non_orphan.csv'))

spotrnaN100_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_spotrnaN100.csv'))
spotrnaN100_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_spotrnaN100_orphan.csv'))
spotrnaN100_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_spotrnaN100_non_orphan.csv'))

spotrnaN100_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_spotrnaN100.csv'))
spotrnaN100_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_spotrnaN100_orphan.csv'))
spotrnaN100_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_spotrnaN100_non_orphan.csv'))

# rna monomers
spotrnaN100_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_spotrnaN100.csv'))
spotrnaN100_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_spotrnaN100_orphan.csv'))
spotrnaN100_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_spotrnaN100_non_orphan.csv'))

spotrnaN100_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_spotrnaN100.csv'))
spotrnaN100_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_spotrnaN100_orphan.csv'))
spotrnaN100_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_spotrnaN100_non_orphan.csv'))

spotrnaN100_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_spotrnaN100.csv'))
spotrnaN100_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_spotrnaN100_orphan.csv'))
spotrnaN100_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_spotrnaN100_non_orphan.csv'))

# protein-rna
spotrnaN100_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_spotrnaN100.csv'))
spotrnaN100_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_spotrnaN100_orphan.csv'))
spotrnaN100_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_spotrnaN100_non_orphan.csv'))

spotrnaN100_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_spotrnaN100.csv'))
spotrnaN100_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_spotrnaN100_orphan.csv'))
spotrnaN100_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_spotrnaN100_non_orphan.csv'))

spotrnaN100_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_spotrnaN100.csv'))
spotrnaN100_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_spotrnaN100_orphan.csv'))
spotrnaN100_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_spotrnaN100_non_orphan.csv'))



#####################################################################################################################
# SHS RNAformerN100

# all data
rnaformerN100_all = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN100.csv'))
rnaformerN100_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN100_orphan.csv'))
rnaformerN100_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN100_non_orphan.csv'))

rnaformerN100_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN100.csv'))
rnaformerN100_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN100_orphan.csv'))
rnaformerN100_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN100_non_orphan.csv'))

rnaformerN100_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN100.csv'))
rnaformerN100_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN100_orphan.csv'))
rnaformerN100_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN100_non_orphan.csv'))

# rna monomers
rnaformerN100_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN100.csv'))
rnaformerN100_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN100_orphan.csv'))
rnaformerN100_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN100_non_orphan.csv'))

rnaformerN100_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN100.csv'))
rnaformerN100_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN100_orphan.csv'))
rnaformerN100_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN100_non_orphan.csv'))

rnaformerN100_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN100.csv'))
rnaformerN100_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN100_orphan.csv'))
rnaformerN100_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN100_non_orphan.csv'))

# protein-rna
rnaformerN100_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN100.csv'))
rnaformerN100_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN100_orphan.csv'))
rnaformerN100_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN100_non_orphan.csv'))

rnaformerN100_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN100.csv'))
rnaformerN100_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN100_orphan.csv'))
rnaformerN100_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN100_non_orphan.csv'))

rnaformerN100_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN100.csv'))
rnaformerN100_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN100_orphan.csv'))
rnaformerN100_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN100_non_orphan.csv'))


#####################################################################################################################
# SHS RNAformerN5000

# all data
rnaformerN5000_all = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN5000.csv'))
rnaformerN5000_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN5000_orphan.csv'))
rnaformerN5000_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN5000_non_orphan.csv'))

rnaformerN5000_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN5000.csv'))
rnaformerN5000_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN5000_orphan.csv'))
rnaformerN5000_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN5000_non_orphan.csv'))

rnaformerN5000_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN5000.csv'))
rnaformerN5000_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN5000_orphan.csv'))
rnaformerN5000_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN5000_non_orphan.csv'))

# rna monomers
rnaformerN5000_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN5000.csv'))
rnaformerN5000_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN5000_orphan.csv'))
rnaformerN5000_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN5000_non_orphan.csv'))

rnaformerN5000_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN5000.csv'))
rnaformerN5000_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN5000_orphan.csv'))
rnaformerN5000_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN5000_non_orphan.csv'))

rnaformerN5000_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN5000.csv'))
rnaformerN5000_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN5000_orphan.csv'))
rnaformerN5000_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN5000_non_orphan.csv'))

# protein-rna
rnaformerN5000_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN5000.csv'))
rnaformerN5000_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN5000_orphan.csv'))
rnaformerN5000_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN5000_non_orphan.csv'))

rnaformerN5000_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN5000.csv'))
rnaformerN5000_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN5000_orphan.csv'))
rnaformerN5000_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN5000_non_orphan.csv'))

rnaformerN5000_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN5000.csv'))
rnaformerN5000_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN5000_orphan.csv'))
rnaformerN5000_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN5000_non_orphan.csv'))

#####################################################################################################################
# SHS RNAformerN10000

# all data
rnaformerN10000_all = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN10000.csv'))
rnaformerN10000_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN10000_orphan.csv'))
rnaformerN10000_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN10000_non_orphan.csv'))

rnaformerN10000_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN10000.csv'))
rnaformerN10000_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN10000_orphan.csv'))
rnaformerN10000_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN10000_non_orphan.csv'))

rnaformerN10000_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN10000.csv'))
rnaformerN10000_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN10000_orphan.csv'))
rnaformerN10000_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN10000_non_orphan.csv'))

# rna monomers
rnaformerN10000_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN10000.csv'))
rnaformerN10000_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN10000_orphan.csv'))
rnaformerN10000_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN10000_non_orphan.csv'))

rnaformerN10000_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN10000.csv'))
rnaformerN10000_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN10000_orphan.csv'))
rnaformerN10000_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN10000_non_orphan.csv'))

rnaformerN10000_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN10000.csv'))
rnaformerN10000_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN10000_orphan.csv'))
rnaformerN10000_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN10000_non_orphan.csv'))

# protein-rna
rnaformerN10000_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN10000.csv'))
rnaformerN10000_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN10000_orphan.csv'))
rnaformerN10000_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN10000_non_orphan.csv'))

rnaformerN10000_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN10000.csv'))
rnaformerN10000_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN10000_orphan.csv'))
rnaformerN10000_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN10000_non_orphan.csv'))

rnaformerN10000_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN10000.csv'))
rnaformerN10000_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN10000_orphan.csv'))
rnaformerN10000_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN10000_non_orphan.csv'))


#####################################################################################################################
# SHS RNAformerN20000

# all data
rnaformerN20000_all = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN20000.csv'))
rnaformerN20000_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN20000_orphan.csv'))
rnaformerN20000_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_rnaformerN20000_non_orphan.csv'))

rnaformerN20000_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN20000.csv'))
rnaformerN20000_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN20000_orphan.csv'))
rnaformerN20000_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_rnaformerN20000_non_orphan.csv'))

rnaformerN20000_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN20000.csv'))
rnaformerN20000_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN20000_orphan.csv'))
rnaformerN20000_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_rnaformerN20000_non_orphan.csv'))

# rna monomers
rnaformerN20000_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN20000.csv'))
rnaformerN20000_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN20000_orphan.csv'))
rnaformerN20000_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_rnaformerN20000_non_orphan.csv'))

rnaformerN20000_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN20000.csv'))
rnaformerN20000_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN20000_orphan.csv'))
rnaformerN20000_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_rnaformerN20000_non_orphan.csv'))

rnaformerN20000_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN20000.csv'))
rnaformerN20000_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN20000_orphan.csv'))
rnaformerN20000_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_rnaformerN20000_non_orphan.csv'))

# protein-rna
rnaformerN20000_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN20000.csv'))
rnaformerN20000_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN20000_orphan.csv'))
rnaformerN20000_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_rnaformerN20000_non_orphan.csv'))

rnaformerN20000_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN20000.csv'))
rnaformerN20000_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN20000_orphan.csv'))
rnaformerN20000_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_rnaformerN20000_non_orphan.csv'))

rnaformerN20000_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN20000.csv'))
rnaformerN20000_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN20000_orphan.csv'))
rnaformerN20000_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_rnaformerN20000_non_orphan.csv'))

#####################################################################################################################
# SHS RNAformerN20000

# all data
noMSA_all = pd.read_csv(Path(datadir, 'csvs/All_noMSA.csv'))
noMSA_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_noMSA_orphan.csv'))
noMSA_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_noMSA_non_orphan.csv'))

noMSA_a2021_all = pd.read_csv(Path(datadir, 'csvs/All_a2021_noMSA.csv'))
noMSA_a2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_noMSA_orphan.csv'))
noMSA_a2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_a2021_noMSA_non_orphan.csv'))

noMSA_b2021_all = pd.read_csv(Path(datadir, 'csvs/All_b2021_noMSA.csv'))
noMSA_b2021_all_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_noMSA_orphan.csv'))
noMSA_b2021_all_non_orphan = pd.read_csv(Path(datadir, 'csvs/All_b2021_noMSA_non_orphan.csv'))

# rna monomers
noMSA_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_noMSA.csv'))
noMSA_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_noMSA_orphan.csv'))
noMSA_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_noMSA_non_orphan.csv'))

noMSA_a2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_noMSA.csv'))
noMSA_a2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_noMSA_orphan.csv'))
noMSA_a2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_a2021_noMSA_non_orphan.csv'))

noMSA_b2021_rnamonomers = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_noMSA.csv'))
noMSA_b2021_rnamonomers_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_noMSA_orphan.csv'))
noMSA_b2021_rnamonomers_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA_Monomers_b2021_noMSA_non_orphan.csv'))

# protein-rna
noMSA_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_noMSA.csv'))
noMSA_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_noMSA_orphan.csv'))
noMSA_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_noMSA_non_orphan.csv'))

noMSA_a2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_noMSA.csv'))
noMSA_a2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_noMSA_orphan.csv'))
noMSA_a2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_a2021_noMSA_non_orphan.csv'))

noMSA_b2021_protein_rna = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_noMSA.csv'))
noMSA_b2021_protein_rna_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_noMSA_orphan.csv'))
noMSA_b2021_protein_rna_non_orphan = pd.read_csv(Path(datadir, 'csvs/RNA-Protein_b2021_noMSA_non_orphan.csv'))


def suggest_figsize(
        n_bars: int,
        height: float = 3.0,      # keep panels short for grids
        bar_width: float = 0.7,
        bar_spacing: float = 1.0,  # 1.0 = default center-to-center
        has_title: bool = True,
        has_nsamples: bool = True
    ):
    """
    Returns (width, height) in inches for compact 'Nature-style' bar panels.
    """
    # base margins: y-axis labels + a bit of breathing room
    base = 1.6
    if has_title:
        base += 0.2
    if has_nsamples:
        # second line under ticks needs a bit more width
        base += 0.1

    # scale with bars and gaps
    per_bar  = 0.50 * (bar_width / 0.7)        # each bar’s visual footprint
    per_gap  = 0.28 * bar_spacing              # space between bar centers

    # total width
    width = base + n_bars * per_bar + max(0, n_bars - 1) * per_gap
    # clamp to sane bounds
    width = max(2.6, min(width, 6.5))
    return (round(width, 2), height)


def suggest_figsize_with_inset(
    n_bars: int,
    height: float = 5.0,
    bar_width: float = 0.3,
    bar_spacing: float = 0.4,
    has_title: bool = True,
    has_nsamples: bool = False,
    inset: bool = False,
    inset_outside: bool = False,
    inset_rel_width: float = 0.38,   # fraction of main axes width to reserve
    inset_gap: float = 0.08          # fraction gap between main axes and inset
):
    base = 1.6 + (0.2 if has_title else 0.0) + (0.1 if has_nsamples else 0.0)
    per_bar = 0.50 * (bar_width / 0.7)
    per_gap = 0.28 * bar_spacing

    # width of the *main* axes
    main_w = base + n_bars * per_bar + max(0, n_bars - 1) * per_gap
    main_w = max(2.8, main_w)   # small floor

    # if we park an inset outside, add room (no artificial small cap)
    if inset and inset_outside:
        extra = main_w * (inset_rel_width + inset_gap)
        width = main_w + extra
        # give a generous ceiling: panels with outside insets are wide by design
        width = min(width, 12.0)
    else:
        width = min(main_w, 7.5)

    return (round(width, 2), height)

#######################################################################################################################################



dfs = [
    ['AD', (alphafold_all, rnaformer_all, rnafold_all, spotrna_all, rnaformerN100_all, rnaformerN5000_all, rnaformerN10000_all, rnaformerN20000_all, rnafoldN100_all, spotrnaN100_all, noMSA_all)], 
    ['ADO', (alphafold_all_orphan, rnaformer_all_orphan, rnafold_all_orphan, spotrna_all_orphan, rnaformerN100_all_orphan, rnaformerN5000_all_orphan, rnaformerN10000_all_orphan, rnaformerN20000_all_orphan, rnafoldN100_all_orphan, spotrnaN100_all_orphan, noMSA_all_orphan)],
    ['ADNO', (alphafold_all_non_orphan, rnaformer_all_non_orphan, rnafold_all_non_orphan, spotrna_all_non_orphan, rnaformerN100_all_non_orphan, rnaformerN5000_all_non_orphan, rnaformerN10000_all_non_orphan, rnaformerN20000_all_non_orphan, rnafoldN100_all_non_orphan, spotrnaN100_all_non_orphan, noMSA_all_non_orphan)],
    ['ADb2021', (alphafold_b2021_all, rnaformer_b2021_all, rnafold_b2021_all, spotrna_b2021_all, rnaformerN100_b2021_all, rnaformerN5000_b2021_all, rnaformerN10000_b2021_all, rnaformerN20000_b2021_all, rnafoldN100_b2021_all, spotrnaN100_b2021_all, noMSA_b2021_all)],
    ['ADOb2021', (alphafold_b2021_all_orphan, rnaformer_b2021_all_orphan, rnafold_b2021_all_orphan, spotrna_b2021_all_orphan, rnaformerN100_b2021_all_orphan, rnaformerN5000_b2021_all_orphan, rnaformerN10000_b2021_all_orphan, rnaformerN20000_b2021_all_orphan, rnafoldN100_b2021_all_orphan, spotrnaN100_b2021_all_orphan, noMSA_b2021_all_orphan)],
    ['ADNOb2021', (alphafold_b2021_all_non_orphan, rnaformer_b2021_all_non_orphan, rnafold_b2021_all_non_orphan, spotrna_b2021_all_non_orphan, rnaformerN100_b2021_all_non_orphan, rnaformerN5000_b2021_all_non_orphan, rnaformerN10000_b2021_all_non_orphan, rnaformerN20000_b2021_all_non_orphan, rnafoldN100_b2021_all_non_orphan, spotrnaN100_b2021_all_non_orphan, noMSA_b2021_all_non_orphan)],
    ['ADa2021', (alphafold_a2021_all, rnaformer_a2021_all, rnafold_a2021_all, spotrna_a2021_all, rnaformerN100_a2021_all, rnaformerN5000_a2021_all, rnaformerN10000_a2021_all, rnaformerN20000_a2021_all, rnafoldN100_a2021_all, spotrnaN100_a2021_all, noMSA_a2021_all)],
    ['ADOa2021', (alphafold_a2021_all_orphan, rnaformer_a2021_all_orphan, rnafold_a2021_all_orphan, spotrna_a2021_all_orphan, rnaformerN100_a2021_all_orphan, rnaformerN5000_a2021_all_orphan, rnaformerN10000_a2021_all_orphan, rnaformerN20000_a2021_all_orphan, rnafoldN100_a2021_all_orphan, spotrnaN100_a2021_all_orphan, noMSA_a2021_all_orphan)],
    ['ADNOa2021', (alphafold_a2021_all_non_orphan, rnaformer_a2021_all_non_orphan, rnafold_a2021_all_non_orphan, spotrna_a2021_all_non_orphan, rnaformerN100_a2021_all_non_orphan, rnaformerN5000_a2021_all_non_orphan, rnaformerN10000_a2021_all_non_orphan, rnaformerN20000_a2021_all_non_orphan, rnafoldN100_a2021_all_non_orphan, spotrnaN100_a2021_all_non_orphan, noMSA_a2021_all_non_orphan)],
    ['RM', (alphafold_rnamonomers, rnaformer_rnamonomers, rnafold_rnamonomers, spotrna_rnamonomers, rnaformerN100_rnamonomers, rnaformerN5000_rnamonomers, rnaformerN10000_rnamonomers, rnaformerN20000_rnamonomers, rnafoldN100_rnamonomers, spotrnaN100_rnamonomers, noMSA_rnamonomers)], 
    ['RMO', (alphafold_rnamonomers_orphan, rnaformer_rnamonomers_orphan, rnafold_rnamonomers_orphan, spotrna_rnamonomers_orphan, rnaformerN100_rnamonomers_orphan, rnaformerN5000_rnamonomers_orphan, rnaformerN10000_rnamonomers_orphan, rnaformerN20000_rnamonomers_orphan, rnafoldN100_rnamonomers_orphan, spotrnaN100_rnamonomers_orphan, noMSA_rnamonomers_orphan)],
    ['RMNO', (alphafold_rnamonomers_non_orphan, rnaformer_rnamonomers_non_orphan, rnafold_rnamonomers_non_orphan, spotrna_rnamonomers_non_orphan, rnaformerN100_rnamonomers_non_orphan, rnaformerN5000_rnamonomers_non_orphan, rnaformerN10000_rnamonomers_non_orphan, rnaformerN20000_rnamonomers_non_orphan, rnafoldN100_rnamonomers_non_orphan, spotrnaN100_rnamonomers_non_orphan, noMSA_rnamonomers_non_orphan)],
    ['RMb2021', (alphafold_b2021_rnamonomers, rnaformer_b2021_rnamonomers, rnafold_b2021_rnamonomers, spotrna_b2021_rnamonomers, rnaformerN100_b2021_rnamonomers, rnaformerN5000_b2021_rnamonomers, rnaformerN10000_b2021_rnamonomers, rnaformerN20000_b2021_rnamonomers, rnafoldN100_b2021_rnamonomers, spotrnaN100_b2021_rnamonomers, noMSA_b2021_rnamonomers)],
    ['RMOb2021', (alphafold_b2021_rnamonomers_orphan, rnaformer_b2021_rnamonomers_orphan, rnafold_b2021_rnamonomers_orphan, spotrna_b2021_rnamonomers_orphan, rnaformerN100_b2021_rnamonomers_orphan, rnaformerN5000_b2021_rnamonomers_orphan, rnaformerN10000_b2021_rnamonomers_orphan, rnaformerN20000_b2021_rnamonomers_orphan, rnafoldN100_b2021_rnamonomers_orphan, spotrnaN100_b2021_rnamonomers_orphan, noMSA_b2021_rnamonomers_orphan)],
    ['RMNOb2021', (alphafold_b2021_rnamonomers_non_orphan, rnaformer_b2021_rnamonomers_non_orphan, rnafold_b2021_rnamonomers_non_orphan, spotrna_b2021_rnamonomers_non_orphan, rnaformerN100_b2021_rnamonomers_non_orphan, rnaformerN5000_b2021_rnamonomers_non_orphan, rnaformerN10000_b2021_rnamonomers_non_orphan, rnaformerN20000_b2021_rnamonomers_non_orphan, rnafoldN100_b2021_rnamonomers_non_orphan, spotrnaN100_b2021_rnamonomers_non_orphan, noMSA_b2021_rnamonomers_non_orphan)],
    ['RMa2021', (alphafold_a2021_rnamonomers, rnaformer_a2021_rnamonomers, rnafold_a2021_rnamonomers, spotrna_a2021_rnamonomers, rnaformerN100_a2021_rnamonomers, rnaformerN5000_a2021_rnamonomers, rnaformerN10000_a2021_rnamonomers, rnaformerN20000_a2021_rnamonomers, rnafoldN100_a2021_rnamonomers, spotrnaN100_a2021_rnamonomers, noMSA_a2021_rnamonomers)],
    ['RMOa2021', (alphafold_a2021_rnamonomers_orphan, rnaformer_a2021_rnamonomers_orphan, rnafold_a2021_rnamonomers_orphan, spotrna_a2021_rnamonomers_orphan, rnaformerN100_a2021_rnamonomers_orphan, rnaformerN5000_a2021_rnamonomers_orphan, rnaformerN10000_a2021_rnamonomers_orphan, rnaformerN20000_a2021_rnamonomers_orphan, rnafoldN100_a2021_rnamonomers_orphan, spotrnaN100_a2021_rnamonomers_orphan, noMSA_a2021_rnamonomers_orphan)],
    ['RMNOa2021', (alphafold_a2021_rnamonomers_non_orphan, rnaformer_a2021_rnamonomers_non_orphan, rnafold_a2021_rnamonomers_non_orphan, spotrna_a2021_rnamonomers_non_orphan, rnaformerN100_a2021_rnamonomers_non_orphan, rnaformerN5000_a2021_rnamonomers_non_orphan, rnaformerN10000_a2021_rnamonomers_non_orphan, rnaformerN20000_a2021_rnamonomers_non_orphan, rnafoldN100_a2021_rnamonomers_non_orphan, spotrnaN100_a2021_rnamonomers_non_orphan, noMSA_a2021_rnamonomers_non_orphan)],
    ['RP', (alphafold_protein_rna, rnaformer_protein_rna, rnafold_protein_rna, spotrna_protein_rna, rnaformerN100_protein_rna, rnaformerN5000_protein_rna, rnaformerN10000_protein_rna, rnaformerN20000_protein_rna, rnafoldN100_protein_rna, spotrnaN100_protein_rna, noMSA_protein_rna)], 
    ['RPO', (alphafold_protein_rna_orphan, rnaformer_protein_rna_orphan, rnafold_protein_rna_orphan, spotrna_protein_rna_orphan, rnaformerN100_protein_rna_orphan, rnaformerN5000_protein_rna_orphan, rnaformerN10000_protein_rna_orphan, rnaformerN20000_protein_rna_orphan, rnafoldN100_protein_rna_orphan, spotrnaN100_protein_rna_orphan, noMSA_protein_rna_orphan)],
    ['RPNO', (alphafold_protein_rna_non_orphan, rnaformer_protein_rna_non_orphan, rnafold_protein_rna_non_orphan, spotrna_protein_rna_non_orphan, rnaformerN100_protein_rna_non_orphan, rnaformerN5000_protein_rna_non_orphan, rnaformerN10000_protein_rna_non_orphan, rnaformerN20000_protein_rna_non_orphan, rnafoldN100_protein_rna_non_orphan, spotrnaN100_protein_rna_non_orphan, noMSA_protein_rna_non_orphan)],
    ['RPb2021', (alphafold_b2021_protein_rna, rnaformer_b2021_protein_rna, rnafold_b2021_protein_rna, spotrna_b2021_protein_rna, rnaformerN100_b2021_protein_rna, rnaformerN5000_b2021_protein_rna, rnaformerN10000_b2021_protein_rna, rnaformerN20000_b2021_protein_rna, rnafoldN100_b2021_protein_rna, spotrnaN100_b2021_protein_rna, noMSA_b2021_protein_rna)],
    ['RPOb2021', (alphafold_b2021_protein_rna_orphan, rnaformer_b2021_protein_rna_orphan, rnafold_b2021_protein_rna_orphan, spotrna_b2021_protein_rna_orphan, rnaformerN100_b2021_protein_rna_orphan, rnaformerN5000_b2021_protein_rna_orphan, rnaformerN10000_b2021_protein_rna_orphan, rnaformerN20000_b2021_protein_rna_orphan, rnafoldN100_b2021_protein_rna_orphan, spotrnaN100_b2021_protein_rna_orphan, noMSA_b2021_protein_rna_orphan)],
    ['RPNOb2021', (alphafold_b2021_protein_rna_non_orphan, rnaformer_b2021_protein_rna_non_orphan, rnafold_b2021_protein_rna_non_orphan, spotrna_b2021_protein_rna_non_orphan, rnaformerN100_b2021_protein_rna_non_orphan, rnaformerN5000_b2021_protein_rna_non_orphan, rnaformerN10000_b2021_protein_rna_non_orphan, rnaformerN20000_b2021_protein_rna_non_orphan, rnafoldN100_b2021_protein_rna_non_orphan, spotrnaN100_b2021_protein_rna_non_orphan, noMSA_b2021_protein_rna_non_orphan)],
    ['RPa2021', (alphafold_a2021_protein_rna, rnaformer_a2021_protein_rna, rnafold_a2021_protein_rna, spotrna_a2021_protein_rna, rnaformerN100_a2021_protein_rna, rnaformerN5000_a2021_protein_rna, rnaformerN10000_a2021_protein_rna, rnaformerN20000_a2021_protein_rna, rnafoldN100_a2021_protein_rna, spotrnaN100_a2021_protein_rna, noMSA_a2021_protein_rna)],
    ['RPOa2021', (alphafold_a2021_protein_rna_orphan, rnaformer_a2021_protein_rna_orphan, rnafold_a2021_protein_rna_orphan, spotrna_a2021_protein_rna_orphan, rnaformerN100_a2021_protein_rna_orphan, rnaformerN5000_a2021_protein_rna_orphan, rnaformerN10000_a2021_protein_rna_orphan, rnaformerN20000_a2021_protein_rna_orphan, rnafoldN100_a2021_protein_rna_orphan, spotrnaN100_a2021_protein_rna_orphan, noMSA_a2021_protein_rna_orphan)],
    ['RPNOa2021', (alphafold_a2021_protein_rna_non_orphan, rnaformer_a2021_protein_rna_non_orphan, rnafold_a2021_protein_rna_non_orphan, spotrna_a2021_protein_rna_non_orphan, rnaformerN100_a2021_protein_rna_non_orphan, rnaformerN5000_a2021_protein_rna_non_orphan, rnaformerN10000_a2021_protein_rna_non_orphan, rnaformerN20000_a2021_protein_rna_non_orphan, rnafoldN100_a2021_protein_rna_non_orphan, spotrnaN100_a2021_protein_rna_non_orphan, noMSA_a2021_protein_rna_non_orphan)],
    ]

print('Evaluate', measure)

val_ids = ['3BWP', '255D', '6E80', '4FAQ', '3Q50', '2ZY6', '4E8V', '7KD1', '3AM1', '2Q1R', '3GCA', '3ND3', '3DHS', '3NPN', '7M5O', '4E8P', '4E8Q', '4RBQ', '1U9S', '4WJ4', '6UES', '4EN5', '4E8M', '4C40', '6TF3', '5C5W', '4CS1', '4E8N', '5DA6', '6TB7', '4P8Z', '2A2E', '6IV9', '2A64', '5HSW', '413D', '3R4F', '2DVI', '4GMA', '6TFE', '3D0M', '4DS6', '387D', '7D7W', '6TF1', '6UET', '6T3S', '6DTD', '6PQ7', '4AOB']

msa_sizes = msa_sizes[~msa_sizes['exp_db_id'].isin(val_ids)]

print(msa_sizes)

for dset, (af, rf, rnafold, spot, n100, n5000, n10000, n20000, rnafoldN100, spotN100, noMSA) in dfs:
    af = af[~af['exp_db_id'].isin(val_ids)]
    rf = rf[~rf['exp_db_id'].isin(val_ids)]
    rnafold = rnafold[~rnafold['exp_db_id'].isin(val_ids)]
    spot = spot[~spot['exp_db_id'].isin(val_ids)]
    n100 = n100[~n100['exp_db_id'].isin(val_ids)]
    n5000 = n5000[~n5000['exp_db_id'].isin(val_ids)]
    n10000 = n10000[~n10000['exp_db_id'].isin(val_ids)]
    n20000 = n20000[~n20000['exp_db_id'].isin(val_ids)]
    rnafoldN100 = rnafoldN100[~rnafoldN100['exp_db_id'].isin(val_ids)]
    spotN100 = spotN100[~spotN100['exp_db_id'].isin(val_ids)]
    noMSA = noMSA[~noMSA['exp_db_id'].isin(val_ids)]

    af = af[af['NumberRNAs'] <= max_rna]
    af = af[af['NumberProteins'] <= max_protein]
    rf = rf[rf['exp_db_id'].isin(af['exp_db_id'].unique())]
    rnafold = rnafold[rnafold['exp_db_id'].isin(af['exp_db_id'].unique())]
    spot = spot[spot['exp_db_id'].isin(af['exp_db_id'].unique())]
    n100 = n100[n100['exp_db_id'].isin(af['exp_db_id'].unique())]
    n5000 = n5000[n5000['exp_db_id'].isin(af['exp_db_id'].unique())]
    n20000 = n20000[n20000['exp_db_id'].isin(af['exp_db_id'].unique())]
    rnafoldN100 = rnafoldN100[rnafoldN100['exp_db_id'].isin(af['exp_db_id'].unique())]
    spotN100 = spotN100[spotN100['exp_db_id'].isin(af['exp_db_id'].unique())]
    noMSA = noMSA[noMSA['exp_db_id'].isin(af['exp_db_id'].unique())]

    # af['msa_size'] = af.apply(get_rna_msa_size, axis=1)
    # rf['msa_size'] = rf.apply(get_rna_msa_size, axis=1)
    # rnafold['msa_size'] = rnafold.apply(get_rna_msa_size, axis=1)
    # spot['msa_size'] = spot.apply(get_rna_msa_size, axis=1)
    # n100['msa_size'] = n100.apply(get_rna_msa_size, axis=1)
    # n5000['msa_size'] = n5000.apply(get_rna_msa_size, axis=1)
    # n10000['msa_size'] = n10000.apply(get_rna_msa_size, axis=1)
    # n20000['msa_size'] = n20000.apply(get_rna_msa_size, axis=1)
    # rnafoldN100['msa_size'] = rnafoldN100.apply(get_rna_msa_size, axis=1)
    # spotN100['msa_size'] = spotN100.apply(get_rna_msa_size, axis=1)
    # noMSA['msa_size'] = noMSA.apply(get_rna_msa_size, axis=1)
    # # rf = rf[rf['NumberRNAs'] <= max_rna]
    # # rf = rf[rf['NumberProteins'] <= max_protein]
    # rnafold = rnafold[rnafold['NumberRNAs'] <= max_rna]
    # rnafold = rnafold[rnafold['NumberProteins'] <= max_protein]
    # spot = spot[spot['NumberRNAs'] <= max_rna]
    # spot = spot[spot['NumberProteins'] <= max_protein]
    # n100 = n100[n100['NumberRNAs'] <= max_rna]
    # n100 = n100[n100['NumberProteins'] <= max_protein]
    # n5000 = n5000[n5000['NumberRNAs'] <= max_rna]
    # n5000 = n5000[n5000['NumberProteins'] <= max_protein]
    # n20000 = n20000[n20000['NumberRNAs'] <= max_rna]
    # n20000 = n20000[n20000['NumberProteins'] <= max_protein]
    # rnafoldN100 = rnafoldN100[rnafoldN100['NumberRNAs'] <= max_rna]
    # rnafoldN100 = rnafoldN100[rnafoldN100['NumberProteins'] <= max_protein]
    # spotN100 = spotN100[spotN100['NumberRNAs'] <= max_rna]
    # spotN100 = spotN100[spotN100['NumberProteins'] <= max_protein]


    print()
    print(dset, f'(n={len(af)})')
    print(f'Alphafold 3 (n={len(af)}):', np.round(af['Complex_RMSD'].mean(), 3), np.round(af['Complex_RMSD'].std(), 3), np.round(af['Complex_LDDT'].mean(), 3), np.round(af['Complex_LDDT'].std(), 3), np.round(af['Complex_TM'].mean(), 3), np.round(af['Complex_TM'].std(), 3))
    
    print(f'SHS RNAformer (n={len(rf)}):', np.round(rf['Complex_RMSD'].mean(), 3), np.round(rf['Complex_RMSD'].std(), 3), np.round(rf['Complex_LDDT'].mean(), 3), np.round(rf['Complex_LDDT'].std(), 3), np.round(rf['Complex_TM'].mean(), 3), np.round(rf['Complex_TM'].std(), 3))
    
    print(f'SHS RNAfold (n={len(rnafold)}):', np.round(rnafold['Complex_RMSD'].mean(), 3), np.round(rnafold['Complex_RMSD'].std(), 3), np.round(rnafold['Complex_LDDT'].mean(), 3), np.round(rnafold['Complex_LDDT'].std(), 3), np.round(rnafold['Complex_TM'].mean(), 3), np.round(rnafold['Complex_TM'].std(), 3))
    
    print(f'SHS SPOT-RNA (n={len(spot)}):', np.round(spot['Complex_RMSD'].mean(), 3), np.round(spot['Complex_RMSD'].std(), 3), np.round(spot['Complex_LDDT'].mean(), 3), np.round(spot['Complex_LDDT'].std(), 3), np.round(spot['Complex_TM'].mean(), 3), np.round(spot['Complex_TM'].std(), 3))
    
    print(f"SHS {n100['algorithm'].unique()} (n={len(n100)}):", np.round(n100['Complex_RMSD'].mean(), 3), np.round(n100['Complex_RMSD'].std(), 3), np.round(n100['Complex_LDDT'].mean(), 3), np.round(n100['Complex_LDDT'].std(), 3), np.round(n100['Complex_TM'].mean(), 3), np.round(n100['Complex_TM'].std(), 3))
    
    print(f'SHS RNAformerN5000 (n={len(n5000)}):', np.round(n5000['Complex_RMSD'].mean(), 3), np.round(n5000['Complex_RMSD'].std(), 3), np.round(n5000['Complex_LDDT'].mean(), 3), np.round(n5000['Complex_LDDT'].std(), 3), np.round(n5000['Complex_TM'].mean(), 3), np.round(n5000['Complex_TM'].std(), 3))
    
    print(f'SHS RNAformerN10000 (n={len(n10000)}):', np.round(n10000['Complex_RMSD'].mean(), 3), np.round(n10000['Complex_RMSD'].std(), 3), np.round(n10000['Complex_LDDT'].mean(), 3), np.round(n10000['Complex_LDDT'].std(), 3), np.round(n10000['Complex_TM'].mean(), 3), np.round(n10000['Complex_TM'].std(), 3))
    
    print(f'SHS RNAformerN20000 (n={len(n20000)}):', np.round(n20000['Complex_RMSD'].mean(), 3), np.round(n20000['Complex_RMSD'].std(), 3), np.round(n20000['Complex_LDDT'].mean(), 3), np.round(n20000['Complex_LDDT'].std(), 3), np.round(n20000['Complex_TM'].mean(), 3), np.round(n20000['Complex_TM'].std(), 3))
    
    print(f"SHS RNAfoldN100 (n={len(rnafoldN100)}):", np.round(rnafoldN100['Complex_RMSD'].mean(), 3), np.round(rnafoldN100['Complex_RMSD'].std(), 3), np.round(rnafoldN100['Complex_LDDT'].mean(), 3), np.round(rnafoldN100['Complex_LDDT'].std(), 3), np.round(rnafoldN100['Complex_TM'].mean(), 3), np.round(rnafoldN100['Complex_TM'].std(), 3))
    
    print(f"SHS SPOT-RNA N100 (n={len(spotN100)}):", np.round(spotN100['Complex_RMSD'].mean(), 3), np.round(spotN100['Complex_RMSD'].std(), 3), np.round(spotN100['Complex_LDDT'].mean(), 3), np.round(spotN100['Complex_LDDT'].std(), 3), np.round(spotN100['Complex_TM'].mean(), 3), np.round(spotN100['Complex_TM'].std(), 3))
    
    print(f"Alphafold noMSA (n={len(noMSA)}):", np.round(noMSA['Complex_RMSD'].mean(), 3), np.round(noMSA['Complex_RMSD'].std(), 3), np.round(noMSA['Complex_LDDT'].mean(), 3), np.round(noMSA['Complex_LDDT'].std(), 3), np.round(noMSA['Complex_TM'].mean(), 3), np.round(noMSA['Complex_TM'].std(), 3))


    if dset == 'AD':
        dataset = f'All Data (n={len(af)})'
    elif dset == 'ADO':
        dataset = f'All Data w/o Homologs (n={len(af)})'
    elif dset == 'ADNO':
        dataset = f'All Data w/ Homologs (n={len(af)})'
    elif dset == 'ADb2021':
        dataset = f'All Data (n={len(af)}; b2021)'
    elif dset == 'ADOb2021':
        dataset = f'All Data w/o Homologs (n={len(af)}; b2021)'
    elif dset == 'ADNOb2021':
        dataset = f'All Data w/ Homologs (n={len(af)}; b2021)'
    elif dset == 'ADa2021':
        dataset = f'All Data (n={len(af)}; a2021)'
    elif dset == 'ADOa2021':
        dataset = f'All Data w/o Homologs (n={len(af)}; a2021)'
    elif dset == 'ADNOa2021':
        dataset = f'All Data w/ Homologs (n={len(af)}; a2021)'
    elif dset == 'RM':
        dataset = f'RNA (n={len(af)})'
    elif dset == 'RMO':
        dataset = f'RNA w/o Homologs (n={len(af)})'
    elif dset == 'RMNO':
        dataset = f'RNA w/ Homologs (n={len(af)})'
    elif dset == 'RMb2021':
        dataset = f'RNA (n={len(af)}; b2021)'
    elif dset == 'RMOb2021':
        dataset = f'RNA w/o Homologs (n={len(af)}; b2021)'
    elif dset == 'RMNOb2021':
        dataset = f'RNA w/ Homologs (n={len(af)}; b2021)'
    elif dset == 'RMa2021':
        dataset = f'RNA (n={len(af)}; a2021)'        
    elif dset == 'RMOa2021':
        dataset = f'RNA w/o Homologs (n={len(af)}; a2021)'
    elif dset == 'RMNOa2021':
        dataset = f'RNA w/ Homologs (n={len(af)}; a2021)'
    elif dset == 'RP':
        dataset = f'RNA-Protein (n={len(af)})'
    elif dset == 'RPO':
        dataset = f'RNA-Protein w/o Homologs (n={len(af)})'
    elif dset == 'RPNO':
        dataset = f'RNA-Protein w/ Homologs (n={len(af)})'
    elif dset == 'RPb2021':
        dataset = f'RNA-Protein (n={len(af)}; b2021)'
    elif dset == 'RPOb2021':
        dataset = f'RNA-Protein w/o Homologs (n={len(af)}; b2021)'
    elif dset == 'RPNOb2021':
        dataset = f'RNA-Protein w/ Homologs (n={len(af)}; b2021)'
    elif dset == 'RPa2021':
        dataset = f'RNA-Protein (n={len(af)}; a2021)'
    elif dset == 'RPOa2021':
        dataset = f'RNA-Protein w/o Homologs (n={len(af)}; a2021)'
    elif dset == 'RPNOa2021':
        dataset = f'RNA-Protein w/ Homologs (n={len(af)}; a2021)'
    else:
        raise UserWarning(f'Unknown dataset {dset} during plotting of results.')
    

###########################################################################################################################################
# Plotting


###########################################################################################################################################
# overview

    print('Plotting', dataset)

    # combined_df = pd.concat([af, rf, spot, rnafold])
    combined_df = pd.concat([af, n100, spotN100, rnafoldN100])
    print(len(combined_df), combined_df['algorithm'].unique())
    
    plotting_data = []
    
    for model, group in combined_df.groupby('algorithm'):
        g = group.copy()
        print(model)
        print(model, f"n={len(g)}", np.round(np.mean(g['Complex_RMSD']), 3), np.round(np.std(g['Complex_RMSD']), 3), np.round(np.mean(g['Complex_LDDT']), 3), np.round(np.std(g['Complex_LDDT']), 3), np.round(np.median(g['Complex_LDDT']), 3), np.round(np.mean(g['Complex_TM']), 3), np.round(np.std(g['Complex_TM']), 3), np.round(np.median(g['Complex_TM']), 3), g[g['Complex_TM'] >= 0.6].shape[0], np.round(np.mean(g['iLDDT']), 3), np.round(np.std(g['iLDDT']), 3), np.round(np.median(g['iLDDT']), 3))
        if model == 'alphafold':
            g.loc[:, 'algorithm'] = 'MSA'
        elif model == 'spotrnaN100':
            g.loc[:, 'algorithm'] = 'SHS\n(SPOT-RNA)'
            # continue
        elif model == 'rnaformerN100':
            g.loc[:, 'algorithm'] = 'SHS\n(RNAformer)'  # 'Synthetic MSA\n(RNAformer)'
            # continue
        elif model == 'rnafoldN100':
            g.loc[:, 'algorithm'] = 'SHS\n(RNAfold)'
            # continue
        else:
            raise UserWarning(f'Unknown algorithm {model} in overview')
        plotting_data.append(g)
    
    plot_df = pd.concat(plotting_data)

    print(len(plot_df['algorithm'].unique()))

    n_bars = len(plot_df['algorithm'].unique())
    height=5.0 
    bar_width=2.5 
    bar_spacing=3.5
    width_add = 1.5
    title_pad = 10  # 55
    annot_offset = 0.04
    bracket_arm = 0.0
    step_factor = 2.8

    fig_size = suggest_figsize(n_bars, height=height, bar_width=bar_width, bar_spacing=bar_spacing)

    fig_size = (fig_size[0]+width_add, fig_size[1])

    print(fig_size)

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_RMSD',
        order=['SHS\n(RNAformer)', 'MSA', 'SHS\n(SPOT-RNA)', 'SHS\n(RNAfold)'],
        palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex RMSD',
        title=dataset.replace('; a2021', '').replace('; b2021', ''),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 15),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(RNAformer)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_bar_means_Complex_RMSD.svg"),
    )

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_LDDT',
        order=['SHS\n(RNAformer)', 'MSA', 'SHS\n(SPOT-RNA)', 'SHS\n(RNAfold)'],
        palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex LDDT',
        title=dataset.replace('; a2021', '').replace('; b2021', ''),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        # yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 1.0),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(RNAformer)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        # scale_in_0_1=True,
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_bar_means_Complex_LDDT.svg"),
    )

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_TM',
        order=['SHS\n(RNAformer)', 'MSA', 'SHS\n(SPOT-RNA)', 'SHS\n(RNAfold)'],
        palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex TM',
        title=dataset.replace('; a2021', '').replace('; b2021', ''),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        # yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 1.0),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(RNAformer)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        # scale_in_0_1=True,
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_bar_means_Complex_TM.svg"),
    )


###########################################################################################################################################
# RNAformer vs AF3

    print('Plotting', dataset)

    # combined_df = pd.concat([af, rf, spot, rnafold])
    combined_df = pd.concat([af, n100, spotN100, rnafoldN100])
    print(len(combined_df), combined_df['algorithm'].unique())
    
    plotting_data = []
    
    for model, group in combined_df.groupby('algorithm'):
        g = group.copy()
        print(model)
        print(model, f"n={len(g)}", np.round(np.mean(g['Complex_RMSD']), 3), np.round(np.mean(g['Complex_LDDT']), 3), np.round(np.std(g['Complex_LDDT']), 3), np.round(np.median(g['Complex_LDDT']), 3), np.round(np.mean(g['Complex_TM']), 3), np.round(np.std(g['Complex_TM']), 3), np.round(np.median(g['Complex_TM']), 3), g[g['Complex_TM'] >= 0.6].shape[0], np.round(np.mean(g['iLDDT']), 3), np.round(np.std(g['iLDDT']), 3), np.round(np.median(g['iLDDT']), 3))
        if model == 'alphafold':
            g.loc[:, 'algorithm'] = 'MSA'
        elif model == 'spotrnaN100':
            # g.loc[:, 'algorithm'] = 'AF3shs\n(SPOT-RNA)'
            continue
        elif model == 'rnaformerN100':
            g.loc[:, 'algorithm'] = 'SHS\n(RNAformer)'  # 'Synthetic MSA\n(RNAformer)'
            # continue
        elif model == 'rnafoldN100':
            # g.loc[:, 'algorithm'] = 'AF3shs\n(RNAfold)'
            continue
        else:
            raise UserWarning(f'Unknown algorithm {model} in rf vs af')
        plotting_data.append(g)
    
    plot_df = pd.concat(plotting_data)

    print(len(plot_df['algorithm'].unique()))

    n_bars = len(plot_df['algorithm'].unique())
    height=5.0 
    bar_width=0.3 
    bar_spacing=0.4
    width_add = 3.0
    width_add_inset = 4.0
    title_pad = 10  #  55
    annot_offset = 0.04
    bracket_arm = 0.0
    step_factor = 2.8

    fig_size = suggest_figsize(n_bars, height=height, bar_width=bar_width, bar_spacing=bar_spacing)
    fig_size_inset = suggest_figsize_with_inset(
        n_bars=n_bars,
        height=8.0, 
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        has_title=True,
        has_nsamples=False,
        inset=True,
        inset_outside=True,
        inset_rel_width=0.42,      # a bit wider than before
        inset_gap=1.0
    )

    fig_size_inset = (fig_size_inset[0]+width_add_inset, fig_size_inset[1])
    fig_size = (fig_size[0]+width_add, fig_size[1])

    print(fig_size)

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_RMSD',
        order=['SHS\n(RNAformer)', 'MSA'],  # , 'AF3shs\n(SPOT-RNA)', 'AF3shs\n(RNAfold)'],
        palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex RMSD',
        title=dataset.replace('; a2021', '').replace('; b2021', '').replace('w/ Homologs', 'w/\nHomologs').replace('w/o Homologs', 'w/o\nHomologs'),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        # yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 15), # (0, 15)           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(RNAformer)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        yticks=[0, 2, 4, 6],
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_af3_rnaformershs_bar_means_Complex_RMSD.svg"),
    )

    # plot_bar_means_fixed_panel_and_inset(
    #     df=plot_df,
    #     group_col='algorithm',
    #     metric_col='Complex_RMSD',
    #     order=['SHS\n(RNAformer)', 'MSA'],
    #     palette=['#9bd0ff', '#bfbfbf'],
    #     fig_size=fig_size_inset,
    #     bar_width=bar_width, bar_spacing=bar_spacing,
    #     ylabel='Complex RMSD',
    #     title=dataset.replace('; a2021','').replace('; b2021','')
    #                  .replace('w/ Homologs','w/\nHomologs')
    #                  .replace('w/o Homologs','w/o\nHomologs'),
    #     title_size=24,
    #     title_pad=10,
    #     ylabel_size=24, 
    #     xlabel_size=24,
    #     xtick_label_size=18, 
    #     ytick_label_size=18,
    #     ylim=(0, 8),
    #     annotate=True, 
    #     ref_group='SHS\n(RNAformer)',
    #     error='sem', 
    #     annot_offset=0.04, 
    #     bracket_arm=0.0, 
    #     step_factor=2.8,
    #     inset_enabled=True,
    #     inset_groups=['SHS\n(RNAformer)', 'MSA'],
    #     inset_auto=True,                 # let it pick tight y-lims
    #     inset_size=(0.34, 0.55),        # smaller & tidy
    #     inset_loc='upper right',
    #     inset_box_color='#79c34a',
    #     outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_af3_rnaformershs_bar_means_Complex_RMSD_inset_inside.svg"),
    # )

    # plot_bar_means_fixed_panel_and_inset(
    #     df=plot_df,
    #     group_col='algorithm',
    #     metric_col='Complex_RMSD',
    #     order=['SHS\n(RNAformer)', 'MSA'],
    #     palette=['#9bd0ff', '#bfbfbf'],
    #     fig_size=fig_size_inset,
    #     bar_width=bar_width, bar_spacing=bar_spacing,
    #     ylabel='Complex RMSD',
    #     title=dataset.replace('; a2021','').replace('; b2021','')
    #                  .replace('w/ Homologs','w/\nHomologs')
    #                  .replace('w/o Homologs','w/o\nHomologs'),
    #     title_size=24,
    #     title_pad=10,
    #     ylabel_size=24, 
    #     xlabel_size=24,
    #     xtick_label_size=18, 
    #     ytick_label_size=18,
    #     ylim=(0, 8),
    #     annotate=True, 
    #     ref_group='SHS\n(RNAformer)',
    #     error='sem', 
    #     annot_offset=0.04, 
    #     bracket_arm=0.0, 
    #     step_factor=2.8,
    #     inset_enabled=True,
    #     inset_groups=['SHS\n(RNAformer)', 'MSA'],
    #     inset_auto=True,                 # let it pick tight y-lims
    #     inset_size=(0.40, 0.70),        # smaller & tidy
    #     inset_loc='outside right',
    #     inset_box_color='#79c34a',
    #     outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_af3_rnaformershs_bar_means_Complex_RMSD_inset_outside.svg"),
    # )


    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_LDDT',
        order=['SHS\n(RNAformer)', 'MSA'],  # , 'AF3shs\n(SPOT-RNA)', 'AF3shs\n(RNAfold)'],
        palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex LDDT',
        title=dataset.replace('; a2021', '').replace('; b2021', '').replace('w/ Homologs', 'w/\nHomologs').replace('w/o Homologs', 'w/o\nHomologs'),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        # yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 1.0),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(RNAformer)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        # scale_in_0_1=True,
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_af3_rnaformershs_bar_means_Complex_LDDT.svg"),
    )

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_TM',
        order=['SHS\n(RNAformer)', 'MSA'],  # , 'AF3shs\n(SPOT-RNA)', 'AF3shs\n(RNAfold)'],
        palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex TM',
        title=dataset.replace('; a2021', '').replace('; b2021', '').replace('w/ Homologs', 'w/\nHomologs').replace('w/o Homologs', 'w/o\nHomologs'),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        # yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 1.0),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(RNAformer)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        # scale_in_0_1=True,
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_af3_rnaformershs_bar_means_Complex_TM.svg"),
    )
    

###########################################################################################################################################
# RNAformer vs AF3 number of RNAs analysis

    print('Plotting', dataset)

    # combined_df = pd.concat([af, rf, spot, rnafold])
    combined_df = pd.concat([af, n100, spotN100, rnafoldN100])
    print(len(combined_df), combined_df['algorithm'].unique())
    
    plotting_data = []
    
    for model, group in combined_df.groupby('algorithm'):
        g = group.copy()
        print(model)
        print(model, f"n={len(g)}", np.round(np.mean(g['Complex_RMSD']), 3), np.round(np.mean(g['Complex_LDDT']), 3), np.round(np.mean(g['Complex_TM']), 3))
        if model == 'alphafold':
            g.loc[:, 'algorithm'] = 'MSA'
        elif model == 'spotrnaN100':
            # g.loc[:, 'algorithm'] = 'AF3shs\n(SPOT-RNA)'
            continue
        elif model == 'rnaformerN100':
            g.loc[:, 'algorithm'] = 'SHS (RNAformer)'  # 'Synthetic MSA\n(RNAformer)'
            # continue
        elif model == 'rnafoldN100':
            # g.loc[:, 'algorithm'] = 'AF3shs\n(RNAfold)'
            continue
        else:
            raise UserWarning(f'Unknown algorithm {model} in rf vs af')
        plotting_data.append(g)
    
    plot_df = pd.concat(plotting_data)

    
    numRNA_classes = ['1 RNA Chain', '2 RNA Chains', '3 RNA Chains', '4 RNA Chains']
    numRNA_bins = [0, 1, 2, 3, np.inf]
        
    label_df = plot_df.copy()
    numRNA_labels = []
    label_df['numRNA_class'] = pd.cut(label_df['NumberRNAs'], 
                           bins=numRNA_bins,
                           labels=numRNA_classes)
    total = []
    for name, group in label_df.groupby('numRNA_class'):
        if len(group) == 0:
            total.append(0)
            numRNA_labels.append(f"{name}\n(n=0)")
        else:
            numRNA_labels.append(f"{name}\n(n={int(len(group) / len(group['algorithm'].unique()))})")
            total.append(int(len(group) / len(group['algorithm'].unique())))
        
    total = int(np.sum(total))

    for m in ['Complex_RMSD', 'Complex_LDDT', 'Complex_TM']:    

        plot_grouped_bar_with_error_and_stats(
            df=plot_df,
            x_col='NumberRNAs',
            hue_col='algorithm',
            value_col=m,
            x_bins=numRNA_bins,
            x_bin_labels=numRNA_labels,
            x_order=numRNA_labels,
            palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['#1f77b4', '#ff7f0e', '#2ca02c'],  # Custom colors for the bars
            fig_size=(16, 6),
            bar_width=0.5,
            bar_spacing=1.5,
            title=dataset.replace('; a2021', '').replace('; b2021', '').replace('w/ Homologs', 'w/\nHomologs').replace('w/o Homologs', 'w/o\nHomologs'),
            title_size=24,
            # xlabel='RNA Length (nt)',
            # xlabel_size=16,
            ylabel=m.replace('_', ' '),
            ylabel_size=24,
            xtick_label_size=18,
            ytick_label_size=18,
            # legend_title='Algorithm',
            # legend_title_size=16,
            legend_label_size=24,
            remove_spines=True,
            annotate=True,
            ref_hue='SHS\n(RNAformer)',
            test_method='wilcoxon',
            sig_levels=[(0.001, '***'), (0.01, '**'), (0.05, '*')],
            step_factor=5.5,
            # legend_position='below',
            legend_position='above',
            # legend_title_size=16,
            # legend_label_size=14,
            legend_offset = 0.95,
            yticks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] if m in ['Complex_LDDT', 'Complex_TM'] else None,
            limit_y_range=True if m in ['Complex_LDDT', 'Complex_TM'] else False,
            outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_af3_rnaformershs_bar_numRNAs_means_{m}.svg"),
        )


###########################################################################################################################################
# RNAformer vs AF3 length analysis

    print('Plotting', dataset)

    # combined_df = pd.concat([af, rf, spot, rnafold])
    combined_df = pd.concat([af, n100, spotN100, rnafoldN100])
    print(len(combined_df), combined_df['algorithm'].unique())
    
    plotting_data = []
    
    for model, group in combined_df.groupby('algorithm'):
        g = group.copy()
        print(model)
        print(model, f"n={len(g)}", np.round(np.mean(g['Complex_RMSD']), 3), np.round(np.mean(g['Complex_LDDT']), 3), np.round(np.mean(g['Complex_TM']), 3))
        if model == 'alphafold':
            g.loc[:, 'algorithm'] = 'MSA'
        elif model == 'spotrnaN100':
            # g.loc[:, 'algorithm'] = 'AF3shs\n(SPOT-RNA)'
            continue
        elif model == 'rnaformerN100':
            g.loc[:, 'algorithm'] = 'SHS (RNAformer)'  # 'Synthetic MSA\n(RNAformer)'
            # continue
        elif model == 'rnafoldN100':
            # g.loc[:, 'algorithm'] = 'AF3shs\n(RNAfold)'
            continue
        else:
            raise UserWarning(f'Unknown algorithm {model} in rf vs af')
        plotting_data.append(g)
    
    plot_df = pd.concat(plotting_data)


    length_classes = ['Length 0-20', 'Length 20-50', 'Length 50-100', 'Length 100-150', 'Length >150']
    length_bins = [0, 20, 50, 100, 150, np.inf]
        
    label_df = plot_df.copy()
    length_labels = []
    label_df['length_class'] = pd.cut(label_df['RNALength'], 
                           bins=length_bins,
                           labels=length_classes)
    total = []
    for name, group in label_df.groupby('length_class'):
        if len(group) == 0:
            total.append(0)
            length_labels.append(f"{name}\n(n=0)")
        else:
            length_labels.append(f"{name}\n(n={int(len(group) / len(group['algorithm'].unique()))})")
            total.append(int(len(group) / len(group['algorithm'].unique())))
        
    total = int(np.sum(total))

    for m in ['Complex_RMSD', 'Complex_LDDT', 'Complex_TM']:    

        plot_grouped_bar_with_error_and_stats(
            df=plot_df,
            x_col='RNALength',
            hue_col='algorithm',
            value_col=m,
            x_bins=length_bins,
            x_bin_labels=length_labels,
            x_order=length_labels,
            palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['#1f77b4', '#ff7f0e', '#2ca02c'],  # Custom colors for the bars
            fig_size=(16, 6),
            bar_width=0.5,
            bar_spacing=1.5,
            title=dataset.replace('; a2021', '').replace('; b2021', '').replace('w/ Homologs', 'w/\nHomologs').replace('w/o Homologs', 'w/o\nHomologs'),
            title_size=24,
            # xlabel='RNA Length (nt)',
            # xlabel_size=16,
            ylabel=m.replace('_', ' '),
            ylabel_size=24,
            xtick_label_size=18,
            ytick_label_size=18,
            # legend_title='Algorithm',
            # legend_title_size=16,
            legend_label_size=24,
            remove_spines=True,
            annotate=True,
            ref_hue='SHS\n(RNAformer)',
            test_method='wilcoxon',
            sig_levels=[(0.001, '***'), (0.01, '**'), (0.05, '*')],
            step_factor=5.5,
            # legend_position='below',
            legend_position='above',
            # legend_title_size=16,
            # legend_label_size=14,
            legend_offset = 0.95,
            yticks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] if m in ['Complex_LDDT', 'Complex_TM'] else None,
            limit_y_range=True if m in ['Complex_LDDT', 'Complex_TM'] else False,
            outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_af3_rnaformershs_bar_lengths_means_{m}.svg"),
        )

###########################################################################################################################################
# RNAformer vs AF3 msa size analysis

    print('Plotting', dataset)

    for tmp_df in (af, n100):
        tmp_df.drop(columns=['RNA_msa_size'], errors='ignore', inplace=True)

    n100 = n100.merge(
        msa_sizes[['exp_db_id', 'RNA_msa_size']],
        on='exp_db_id',
        how='left'
    )    # combined_df = pd.concat([af, rf, spot, rnafold])

    af = af.merge(
        msa_sizes[['exp_db_id', 'RNA_msa_size']],
        on='exp_db_id',
        how='left'
    )    # combined_df = pd.concat([af, rf, spot, rnafold])

    # n100 = n100.merge(
    #     af[['exp_db_id', 'RNA_msa_size']],
    #     on='exp_db_id',
    #     how='left'
    # )    # combined_df = pd.concat([af, rf, spot, rnafold])
    
    af['msa_size'] = af.apply(get_rna_msa_size, axis=1)
    n100['msa_size'] = n100.apply(get_rna_msa_size, axis=1)
    

    combined_df = pd.concat([af, n100, spotN100, rnafoldN100])
    print(len(combined_df), combined_df['algorithm'].unique())


    
    plotting_data = []
    
    for model, group in combined_df.groupby('algorithm'):
        g = group.copy()
        print(model)
        print(model, f"n={len(g)}", np.round(np.mean(g['Complex_RMSD']), 3), np.round(np.mean(g['Complex_LDDT']), 3), np.round(np.mean(g['Complex_TM']), 3))
        if model == 'alphafold':
            g.loc[:, 'algorithm'] = 'MSA'
        elif model == 'spotrnaN100':
            # g.loc[:, 'algorithm'] = 'AF3shs\n(SPOT-RNA)'
            continue
        elif model == 'rnaformerN100':
            g.loc[:, 'algorithm'] = 'SHS (RNAformer)'  # 'Synthetic MSA\n(RNAformer)'
            # continue
        elif model == 'rnafoldN100':
            # g.loc[:, 'algorithm'] = 'AF3shs\n(RNAfold)'
            continue
        else:
            raise UserWarning(f'Unknown algorithm {model} in rf vs af')
        plotting_data.append(g)
    
    plot_df = pd.concat(plotting_data)


    msa_classes = ['no_msa', 'Depth 1-20', 'Depth 21-100', 'Depth 101-500', 'Depth 501-1000', 'Depth >1000']
    msa_bins = [-np.inf, 0.9, 20, 100, 500, 1000, np.inf]
        
    label_df = plot_df.copy()
    msa_labels = []
    label_df['msa_class'] = pd.cut(label_df['msa_size'], 
                           bins=msa_bins,
                           labels=msa_classes)
    total = []
    for name, group in label_df.groupby('msa_class'):
        if len(group) == 0:
            total.append(0)
            msa_labels.append(f"{name}\n(n=0)")
        else:
            msa_labels.append(f"{name}\n(n={int(len(group) / len(group['algorithm'].unique()))})")
            total.append(int(len(group) / len(group['algorithm'].unique())))
        
    total = int(np.sum(total))

    for m in ['Complex_RMSD', 'Complex_LDDT', 'Complex_TM']:    

        plot_grouped_bar_with_error_and_stats(
            df=plot_df,
            x_col='msa_size',
            hue_col='algorithm',
            value_col=m,
            x_bins=msa_bins,
            x_bin_labels=msa_labels,
            x_order=msa_labels,
            palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['#1f77b4', '#ff7f0e', '#2ca02c'],  # Custom colors for the bars
            fig_size=(16, 6),
            bar_width=0.5,
            bar_spacing=1.5,
            title=dataset.replace('; a2021', '').replace('; b2021', '').replace('w/ Homologs', 'w/\nHomologs').replace('w/o Homologs', 'w/o\nHomologs'),
            title_size=24,
            # xlabel='Number of Sequences in MSA',
            # xlabel_size=16,
            ylabel=m.replace('_', ' '),
            ylabel_size=24,
            xtick_label_size=18,
            ytick_label_size=18,
            # legend_title='Algorithm',
            # legend_title_size=16,
            legend_label_size=24,
            remove_spines=True,
            annotate=True,
            ref_hue='SHS\n(RNAformer)',
            test_method='wilcoxon',
            sig_levels=[(0.001, '***'), (0.01, '**'), (0.05, '*')],
            step_factor=5.5,
            # legend_position='below',
            legend_position='above',
            # legend_title_size=16,
            # legend_label_size=14,
            legend_offset = 0.95,
            yticks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] if m in ['Complex_LDDT', 'Complex_TM'] else None,
            limit_y_range=True if m in ['Complex_LDDT', 'Complex_TM'] else False,
            outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_af3_rnaformershs_bar_msa_size_means_{m}.svg"),
        )


###########################################################################################################################################
# shs comparison

    print('Plotting', dataset)

    # combined_df = pd.concat([af, rf, spot, rnafold])
    combined_df = pd.concat([af, n100, spotN100, rnafoldN100])
    print(len(combined_df), combined_df['algorithm'].unique())
    
    plotting_data = []
    
    for model, group in combined_df.groupby('algorithm'):
        g = group.copy()
        print(model)
        print(model, f"n={len(g)}", np.round(np.mean(g['Complex_RMSD']), 3), np.round(np.mean(g['Complex_LDDT']), 3), np.round(np.std(g['Complex_LDDT']), 3), np.round(np.median(g['Complex_LDDT']), 3), np.round(np.mean(g['Complex_TM']), 3), np.round(np.std(g['Complex_TM']), 3), np.round(np.median(g['Complex_TM']), 3), g[g['Complex_TM'] >= 0.6].shape[0], np.round(np.mean(g['iLDDT']), 3), np.round(np.std(g['iLDDT']), 3), np.round(np.median(g['iLDDT']), 3))
        if model == 'alphafold':
            # g.loc[:, 'algorithm'] = 'AlphaFold 3'
            continue
        elif model == 'spotrnaN100':
            g.loc[:, 'algorithm'] = 'SHS\n(SPOT-RNA)'
            # continue
        elif model == 'rnaformerN100':
            g.loc[:, 'algorithm'] = 'SHS\n(RNAformer)'  # 'Synthetic MSA\n(RNAformer)'
            # continue
        elif model == 'rnafoldN100':
            g.loc[:, 'algorithm'] = 'SHS\n(RNAfold)'
            # continue
        else:
            raise UserWarning(f'Unknown algorithm {model} in shs comparison')
        plotting_data.append(g)
    
    plot_df = pd.concat(plotting_data)

    print(len(plot_df['algorithm'].unique()))

    n_bars = len(plot_df['algorithm'].unique())
    height=5.0 
    bar_width=0.3 
    bar_spacing=0.4
    width_add = 3.0
    title_pad = 30  # 55
    annot_offset = 0.04
    bracket_arm = 0.0
    step_factor = 2.8

    fig_size = suggest_figsize(n_bars, height=height, bar_width=bar_width, bar_spacing=bar_spacing)

    fig_size = (fig_size[0]+width_add, fig_size[1])

    print(fig_size)

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_RMSD',
        order=['SHS\n(RNAformer)', 'SHS\n(SPOT-RNA)', 'SHS\n(RNAfold)'],
        palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex RMSD',
        title=dataset.replace('; a2021', '').replace('; b2021', '').replace('w/ Homologs', 'w/\nHomologs').replace('w/o Homologs', 'w/o\nHomologs'),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 15),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(RNAformer)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_shscomparison_bar_means_Complex_RMSD.svg"),
    )

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_LDDT',
        order=['SHS\n(RNAformer)', 'SHS\n(SPOT-RNA)', 'SHS\n(RNAfold)'],
        palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex LDDT',
        title=dataset.replace('; a2021', '').replace('; b2021', '').replace('w/ Homologs', 'w/\nHomologs').replace('w/o Homologs', 'w/o\nHomologs'),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        # yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 1.0),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(RNAformer)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        # scale_in_0_1=True,
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_shscomparison_bar_means_Complex_LDDT.svg"),
    )

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_TM',
        order=['SHS\n(RNAformer)', 'SHS\n(SPOT-RNA)', 'SHS\n(RNAfold)'],
        palette=['#9bd0ff', '#bfbfbf', '#bfbfbf', '#bfbfbf'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex TM',
        title=dataset.replace('; a2021', '').replace('; b2021', '').replace('w/ Homologs', 'w/\nHomologs').replace('w/o Homologs', 'w/o\nHomologs'),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        # yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 1.0),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(RNAformer)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        # scale_in_0_1=True,
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_shscomparison_bar_means_Complex_TM.svg"),
    )


###########################################################################################################################################
# scaling comparison

    print('Plotting', dataset)

    combined_df = pd.concat([n100, n5000, n10000, n20000])
    print(len(combined_df), combined_df['algorithm'].unique())
    
    plotting_data = []
    
    for model, group in combined_df.groupby('algorithm'):
        g = group.copy()
        print(model)
        print(model, f"n={len(g)}", np.round(np.mean(g['Complex_RMSD']), 3), np.round(np.mean(g['Complex_LDDT']), 3), np.round(np.std(g['Complex_LDDT']), 3), np.round(np.median(g['Complex_LDDT']), 3), np.round(np.mean(g['Complex_TM']), 3), np.round(np.std(g['Complex_TM']), 3), np.round(np.median(g['Complex_TM']), 3), g[g['Complex_TM'] >= 0.6].shape[0], np.round(np.mean(g['iLDDT']), 3), np.round(np.std(g['iLDDT']), 3), np.round(np.median(g['iLDDT']), 3))
        if 'spotrna' in model:
            continue
        if 'rnafold' in model:
            continue
        if 'N10000' in model:
            g.loc[:, 'algorithm'] = 'SHS\n(N=10000)'
            # continue
        elif 'N5000' in model:
            g.loc[:, 'algorithm'] = 'SHS\n(N=5000)'
            # continue
        elif 'N100' in model:
            g.loc[:, 'algorithm'] = 'SHS\n(N=100)'  # 'Synthetic MSA\n(RNAformer)'
            # continue
        elif 'N20000' in model:
            g.loc[:, 'algorithm'] = 'SHS\n(N=20000)'
            # continue
        else:
            raise UserWarning(f'Unknown algorithm {model} in scaling comparison')
        plotting_data.append(g)
    
    plot_df = pd.concat(plotting_data)

    print(len(plot_df['algorithm'].unique()))

    n_bars = len(plot_df['algorithm'].unique())
    height=5.0 
    bar_width=2.5 
    bar_spacing=3.5
    width_add = 1.5
    title_pad = 55
    annot_offset = 0.04
    bracket_arm = 0.0
    step_factor = 2.8

    fig_size = suggest_figsize(n_bars, height=height, bar_width=bar_width, bar_spacing=bar_spacing)

    fig_size = (fig_size[0]+width_add, fig_size[1])

    print(fig_size)

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_RMSD',
        order=['SHS\n(N=100)', 'SHS\n(N=5000)', 'SHS\n(N=10000)', 'SHS\n(N=20000)'],
        palette=['#bfbfbf', '#999999', '#666666', '#3D3D3D'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex RMSD',
        title=dataset.replace('; a2021', '').replace('; b2021', ''),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 15),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(N=100)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_scaling_comparison_bar_means_Complex_RMSD.svg"),
    )

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_LDDT',
        order=['SHS\n(N=100)', 'SHS\n(N=5000)', 'SHS\n(N=10000)', 'SHS\n(N=20000)'],
        palette=['#bfbfbf', '#999999', '#666666', '#3D3D3D'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex LDDT',
        title=dataset.replace('; a2021', '').replace('; b2021', ''),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        # yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 1.0),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(N=100)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        # scale_in_0_1=True,
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_scaling_comparison_bar_means_Complex_LDDT.svg"),
    )

    plot_bar_means_fixed_panel(
        df=plot_df,
        group_col='algorithm',
        metric_col='Complex_TM',
        order=['SHS\n(N=100)', 'SHS\n(N=5000)', 'SHS\n(N=10000)', 'SHS\n(N=20000)'],
        palette=['#bfbfbf', '#999999', '#666666', '#3D3D3D'],  # ['skyblue', '#666666', '#3D3D3D', '#999999'],  # extend if 4 bars
        fig_size=fig_size,
        bar_width=bar_width,
        bar_spacing=bar_spacing,
        ylabel='Complex TM',
        title=dataset.replace('; a2021', '').replace('; b2021', ''),
        title_size=24,
        title_pad=title_pad,
        ylabel_size=24,
        xlabel_size=24,
        xtick_label_size=18,
        ytick_label_size=18,
        # yticks=None,              # e.g. [0, 5, 10, 15]
        xlim=None,                # e.g. (-0.5, 3.5)
        ylim=(0, 1.0),           # keep y ≥ 0; or set (0, 20)
        remove_spines=True,
        y_as_percent=False,         # keep % axis
        scale_in_0_1=False,        # set True if your metric is in [0,1]
        annotate=True,
        ref_group='SHS\n(N=100)',
        test_method='wilcoxon',   # auto-fallback if unpaired
        error='sem',              # or 'sem' if you want tighter bars
        annot_offset=annot_offset,
        bracket_arm=bracket_arm,
        step_factor=step_factor,
        floor_at_zero=True,
        show_nsamples=False,
        # scale_in_0_1=True,
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        outpath=Path(plotting_dir, f"{'_'.join(dataset.split()).replace('/', '_')}_scaling_comparison_bar_means_Complex_TM.svg"),
    )
