from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def _get_algo_color(algo, i, colors, color_cycle):
    """Helper to get a color for a given algorithm."""
    if isinstance(colors, dict):
        if algo in colors:
            return colors[algo]
    elif isinstance(colors, (list, tuple)):
        if len(colors) > 0:
            return colors[i % len(colors)]
    # fallback to matplotlib default cycle
    if color_cycle:
        return color_cycle[i % len(color_cycle)]
    return None


def plot_metric_scatter(
    results,
    algorithms_to_plot=('rnaformer', 'rnafold', 'spotrna'),
    x_metric='mcc',
    y_metrics=('RNA_TM', 'RNA_LDDT', 'RNA_RMSD'),
    figsize=(5, 4),
    marker_size=20,
    alpha=0.7,
    title_fontsize=16,
    label_fontsize=14,
    tick_fontsize=12,
    legend_fontsize=12,
    diagonal_line=True,
    grid=True,
    colors=None,
    save=False,
    save_dir='plots',
    dpi=300,
):
    """Joined scatterplots: multiple algorithms in one plot, separate figure per y_metric."""
    save_dir = Path(save_dir)
    if save:
        save_dir.mkdir(parents=True, exist_ok=True)

    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', [])

    for y_metric in y_metrics:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        all_x_vals = []
        all_y_vals = []

        for i, algo in enumerate(algorithms_to_plot):
            if algo not in results:
                continue

            df = results[algo]
            if x_metric not in df.columns or y_metric not in df.columns:
                continue

            x = df[x_metric]
            y = df[y_metric]

            all_x_vals.extend(x.values)
            all_y_vals.extend(y.values)

            color = _get_algo_color(algo, i, colors, color_cycle)

            ax.scatter(
                x,
                y,
                label=algo,
                s=marker_size,
                alpha=alpha,
                edgecolor='black',
                linewidth=0.5,
                marker='o',
                color=color,
            )

        if diagonal_line and all_x_vals and all_y_vals:
            xy_min = min(min(all_x_vals), min(all_y_vals))
            xy_max = max(max(all_x_vals), max(all_y_vals))
            ax.plot(
                [xy_min, xy_max],
                [xy_min, xy_max],
                linestyle='--',
                linewidth=1.0,
                color='gray',
                label='y = x',
            )

        ax.set_xlabel(x_metric, fontsize=label_fontsize)
        ax.set_ylabel(y_metric, fontsize=label_fontsize)
        ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)

        if grid:
            ax.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)

        ax.legend(
            fontsize=legend_fontsize,
            frameon=True,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
        )

        ax.set_title(f'{y_metric} vs {x_metric} (joined)', fontsize=title_fontsize)

        fig.tight_layout()

        if save:
            out_file = save_dir / f'scatter_joined_{x_metric}_vs_{y_metric}.svg'
            fig.savefig(out_file, dpi=dpi, bbox_inches='tight')

        plt.show()


def plot_metric_scatter_per_method(
    results,
    algorithms_to_plot=('rnaformer', 'rnafold', 'spotrna'),
    x_metric='mcc',
    y_metrics=('RNA_TM', 'RNA_LDDT', 'RNA_RMSD'),
    figsize=(5, 4),
    marker_size=20,
    alpha=0.7,
    title_fontsize=16,
    label_fontsize=14,
    tick_fontsize=12,
    diagonal_line=True,
    grid=True,
    colors=None,
    save=False,
    save_dir='plots',
    dpi=300,
):
    """Individual scatterplots: one algorithm per figure, separate figure per y_metric."""
    save_dir = Path(save_dir)
    if save:
        save_dir.mkdir(parents=True, exist_ok=True)

    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', [])

    for i, algo in enumerate(algorithms_to_plot):
        if algo not in results:
            continue

        df = results[algo]
        color = _get_algo_color(algo, i, colors, color_cycle)

        for y_metric in y_metrics:
            if x_metric not in df.columns or y_metric not in df.columns:
                continue

            x = df[x_metric]
            y = df[y_metric]

            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)

            ax.scatter(
                x,
                y,
                s=marker_size,
                alpha=alpha,
                edgecolor='black',
                linewidth=0.5,
                marker='o',
                color=color,
            )

            if diagonal_line and not x.empty and not y.empty:
                xy_min = min(x.min(), y.min())
                xy_max = max(x.max(), y.max())
                ax.plot(
                    [xy_min, xy_max],
                    [xy_min, xy_max],
                    linestyle='--',
                    linewidth=1.0,
                    color='gray',
                    label='y = x',
                )
                ax.legend(frameon=True, fontsize=tick_fontsize)

            ax.set_xlabel(x_metric, fontsize=label_fontsize)
            ax.set_ylabel(y_metric, fontsize=label_fontsize)
            ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)

            if grid:
                ax.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)

            ax.set_title(f'{algo}: {y_metric} vs {x_metric}', fontsize=title_fontsize)

            fig.tight_layout()

            if save:
                out_file = save_dir / f'scatter_{algo}_{x_metric}_vs_{y_metric}_individual.svg'
                fig.savefig(out_file, dpi=dpi, bbox_inches='tight')

            plt.show()


# NEW: per-method bin-distribution plot styled like plot_algorithm_distributions
def plot_method_interval_distribution(
    df: pd.DataFrame,
    method_name: str,
    sec_metric: str = 'mcc',
    threed_metric: str = 'RNA_TM',
    intervals=((-1.0, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, np.inf)),
    kind: str = 'box',             # 'box' or 'violin'
    palette: list = None,          # list of colors, one per bin
    fig_size: tuple = (6, 4),
    title: str = None,
    title_size: int = 14,
    xlabel: str = None,
    xlabel_size: int = 12,
    ylabel: str = None,
    ylabel_size: int = 12,
    xtick_label_size: int = 10,
    ytick_label_size: int = 10,
    xtick_label_rotation: float = 0,
    remove_spines: bool = True,
    box_width: float = 0.6,
    spacing: float = 1.0,
    ylim_clip_upper: float = 1.0,  # like your function: clip top to 1.0 by default
    outpath: Path = None,
    dpi: int = 300,
):
    """
    For a single method, plot distribution of a 3D metric across MCC bins.

    X labels: "MCC x–y\n(n = ...)" exactly as tick labels (not extra text).
    Visual style closely follows plot_algorithm_distributions.
    """
    data = df.copy()

    # Build data per interval
    bin_values = []
    bin_labels = []
    n_per_bin = []

    for low, high in intervals:
        sec_vals = data[sec_metric]
        th_vals = data[threed_metric]

        if np.isinf(low) and low < 0:
            mask = sec_vals < high
        elif np.isinf(high) and high > 0:
            mask = sec_vals >= low
        else:
            mask = (sec_vals >= low) & (sec_vals < high)

        vals = th_vals[mask].dropna().values
        if len(vals) == 0:
            continue

        bin_values.append(vals)
        n_per_bin.append(len(vals))

        # Label: MCC x–y, with ∞ handling
        if np.isinf(low) and low < 0:
            # label_top = f"{sec_metric.upper()} < {high:.2f}"
            label_top = f"< {high:.1f}"
        elif np.isinf(high) and high > 0:
            # label_top = f"{sec_metric.upper()} ≥ {low:.2f}"
            label_top = f"≥ {low:.1f}"
        else:
            # label_top = f"{sec_metric.upper()} {low:.2f}–{high:.2f}"
            label_top = f"{low:.1f} – {high:.1f}"

        bin_labels.append(f"{label_top}\n(n = {len(vals)})")

    if not bin_values:
        print(f"No data for method {method_name} with given intervals.")
        return

    n_bins = len(bin_values)
    positions = np.arange(n_bins) * spacing

    # Colors
    if palette is None:
        palette = ['skyblue' for _ in range(n_bins)]
    elif len(palette) < n_bins:
        # simple safeguard: repeat if too short
        palette = (palette * ((n_bins // len(palette)) + 1))[:n_bins]

    # Print means
    for label, vals in zip(bin_labels, bin_values):
        clean_label = label.replace("\n", " ")
        print(f"{method_name} | {clean_label}: mean={np.mean(vals):.3f}, n={len(vals)}")

    fig, ax = plt.subplots(figsize=fig_size)

    if kind == 'violin':
        parts = ax.violinplot(bin_values, positions=positions, widths=box_width,
                              showmedians=True)
        for body, color in zip(parts['bodies'], palette):
            body.set_facecolor(color)
            body.set_edgecolor('black')
            body.set_alpha(0.8)
        if 'cmedians' in parts:
            parts['cmedians'].set(color='firebrick', linewidth=2)
    else:
        parts = ax.boxplot(bin_values, positions=positions, widths=box_width,
                           patch_artist=True, showfliers=False)
        for patch, color in zip(parts['boxes'], palette):
            patch.set_facecolor(color)
            patch.set_edgecolor('black')
            patch.set_linewidth(1)
        for med in parts['medians']:
            med.set(color='firebrick', linewidth=2)

    if remove_spines:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # y-limits
    ymin, ymax = ax.get_ylim()
    if ylim_clip_upper is not None:
        ax.set_ylim(ymin, ylim_clip_upper)
        ymax = ylim_clip_upper
    y_min, y_max = ymin, ymax

    # labels
    ax.set_xticks(positions)
    ax.set_xticklabels(bin_labels, rotation=xtick_label_rotation,
                       ha='right' if xtick_label_rotation else 'center',
                       fontsize=xtick_label_size)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=xlabel_size)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=ylabel_size)
    ax.tick_params(axis='y', labelsize=ytick_label_size)

    if title:
        ax.set_title(title, loc='left', fontsize=title_size, pad=10)

    plt.tight_layout()

    if outpath is not None:
        outpath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, dpi=dpi)
    plt.show()
    plt.close(fig)


def main():
    algorithms = [
        'rnaformer',
        'rnafold',
        'spotrna',
        'rnaformerN100_dssr',
        'rnafold_dssr',
        'spotrna_dssr',
        'alphafold_dssr',
        'rnafoldN100_dssr',
        'spotrnaN100_dssr',
        'rnaformerN100',
        'rnafoldN100',
        'spotrnaN100',
    ]
    sec_results = pd.read_csv('RNA_Monomer_2D_and_3D_metrics.csv')
    sec_results = sec_results.dropna()
    mask = sec_results['method'] == 'rnaformer'
    dup = sec_results[mask].copy()
    dup['method'] = 'rnaformerN100'
    sec_results = pd.concat([sec_results, dup], ignore_index=True)

    mask = sec_results['method'] == 'rnafold'
    dup = sec_results[mask].copy()
    dup['method'] = 'rnafoldN100'
    sec_results = pd.concat([sec_results, dup], ignore_index=True)

    mask = sec_results['method'] == 'spotrna'
    dup = sec_results[mask].copy()
    dup['method'] = 'spotrnaN100'
    sec_results = pd.concat([sec_results, dup], ignore_index=True)

    print('Evaluating', len(sec_results), 'RNA monomer structures with both 2D and 3D metrics.')

    results = {}
    for algo in algorithms:
        sec_algo = sec_results[sec_results['method'] == algo]
        threed_algo = pd.read_csv(f"results/csvs/RNA_Monomers_{algo.split('_')[0]}.csv")
        threed_algo['pdbid'] = threed_algo['exp_db_id'].str.upper()
        results[algo] = pd.merge(sec_algo, threed_algo, on='pdbid', how='inner')

    for k, v in results.items():
        print(f"Algorithm: {k}")
        print('MCC:', v['mcc'].mean())
        print('RNA_LDDT:', v['RNA_LDDT'].mean())
        print('RNA_TM:', v['RNA_TM'].mean())
        print("\n")

    algo_colors = {
        'rnaformer': '#1f77b4',  # blue
        'rnafold':   '#ff7f0e',  # orange
        'spotrna':   '#2ca02c',  # green
    }

    # # Example scatter calls (unchanged)
    # plot_metric_scatter(
    #     results,
    #     algorithms_to_plot=('rnaformer', 'rnafold', 'spotrna'),
    #     x_metric='mcc',
    #     y_metrics=('RNA_TM', 'RNA_LDDT', 'RNA_RMSD'),
    #     figsize=(5, 4),
    #     marker_size=25,
    #     alpha=0.8,
    #     title_fontsize=16,
    #     label_fontsize=14,
    #     tick_fontsize=12,
    #     legend_fontsize=12,
    #     diagonal_line=True,
    #     grid=True,
    #     colors=algo_colors,
    #     save=False,
    #     save_dir='plots',
    #     dpi=300,
    # )

    # plot_metric_scatter_per_method(
    #     results,
    #     algorithms_to_plot=('rnaformer', 'rnafold', 'spotrna'),
    #     x_metric='mcc',
    #     y_metrics=('RNA_TM', 'RNA_LDDT', 'RNA_RMSD'),
    #     figsize=(5, 4),
    #     marker_size=25,
    #     alpha=0.8,
    #     title_fontsize=16,
    #     label_fontsize=14,
    #     tick_fontsize=12,
    #     diagonal_line=True,
    #     grid=True,
    #     colors=algo_colors,
    #     save=False,
    #     save_dir='plots',
    #     dpi=300,
    # )

    # === NEW: individual distribution plots for each method ===
    intervals = (
        (-np.inf, 0.5),
        (0.5, 0.7),
        (0.7, 0.9),
        (0.9, np.inf),
    )

    # rnaformer
    plot_method_interval_distribution(
        results['rnaformerN100'],
        method_name='rnaformerN100',
        sec_metric='mcc',
        threed_metric='RNA_TM',
        intervals=intervals,
        kind='box',  # or 'violin'
        # palette=['#c6e2ff', '#9ec5ff', '#6ea8fe', '#3d8bfd'],  # example per-bin colors
        palette=['skyblue', 'skyblue', 'skyblue', 'skyblue'],
        fig_size=(6, 4),
        # title='rnaformer: RNA_LDDT by MCC bin',
        title='SHS (RNAformer)',
        xtick_label_size=18,
        xlabel_size=24,
        ylabel_size=24,
        title_size=24,
        xlabel='MCC',
        ylabel='TM Score',
        outpath=Path('plots/rnaformerN100_mcc_bins_rnatm.svg'),
    )

    # rnafold
    plot_method_interval_distribution(
        results['rnafoldN100'],
        method_name='rnafoldN100',
        sec_metric='mcc',
        threed_metric='RNA_TM',
        intervals=intervals,
        kind='box',  # or 'violin'
        # palette=['#c6e2ff', '#9ec5ff', '#6ea8fe', '#3d8bfd'],  # example per-bin colors
        palette=['skyblue', 'skyblue', 'skyblue', 'skyblue'],
        fig_size=(6, 4),
        # title='rnaformer: RNA_LDDT by MCC bin',
        title='SHS (RNAfold)',
        xtick_label_size=18,
        xlabel_size=24,
        ylabel_size=24,
        title_size=24,
        xlabel='MCC',
        ylabel='TM Score',
        outpath=Path('plots/rnafoldN100_mcc_bins_rnatm.svg'),
    )

    # spotrna
    plot_method_interval_distribution(
        results['spotrnaN100'],
        method_name='spotrnaN100',
        sec_metric='mcc',
        threed_metric='RNA_TM',
        intervals=intervals,
        kind='box',  # or 'violin'
        # palette=['#c6e2ff', '#9ec5ff', '#6ea8fe', '#3d8bfd'],  # example per-bin colors
        palette=['skyblue', 'skyblue', 'skyblue', 'skyblue'],
        fig_size=(6, 4),
        # title='rnaformer: RNA_LDDT by MCC bin',
        title='SHS (SPOT-RNA)',
        xtick_label_size=18,
        xlabel_size=24,
        ylabel_size=24,
        title_size=24,
        xlabel='MCC',
        ylabel='TM Score',
        outpath=Path('plots/spotrnaN100_mcc_bins_rnatm.svg'),
    )


if __name__ == "__main__":
    main()
