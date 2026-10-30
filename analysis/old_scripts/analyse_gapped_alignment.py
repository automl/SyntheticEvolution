from pathlib import Path
from Bio import AlignIO
import logomaker as lm
import matplotlib.pyplot as plt
import itertools
import seaborn as sns

import numpy as np
import pandas as pd

#                                                                       
# INSTALL (once):
# pip install numpy pandas matplotlib scipy scikit-bio scikit-learn umap-learn
# (and you already have biopython & logomaker)
#                                                                        

import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
import subprocess

from scipy.stats import entropy
from sklearn.manifold import TSNE
import umap

from skbio.diversity import beta_diversity
from skbio import DistanceMatrix


# 1) Mutual-Information heatmaps
def compute_mi_matrix(seqs, pseudocount=1e-9):
    """
    seqs: list of equal-length strings
    returns: L×L numpy array of pairwise MI
    """
    arr = np.array([list(s) for s in seqs])
    L = arr.shape[1]
    mi = np.zeros((L, L))
    for i in range(L):
        xi = arr[:, i]
        for j in range(i+1, L):
            xj = arr[:, j]
            joint = pd.crosstab(xi, xj, normalize=True).values + pseudocount
            px = joint.sum(axis=1)
            py = joint.sum(axis=0)
            Hx = entropy(px, base=2)
            Hy = entropy(py, base=2)
            Hxy = entropy(joint.flatten(), base=2)
            mi[i, j] = mi[j, i] = Hx + Hy - Hxy
    return mi

def plot_mi_heatmaps(seqs1, seqs2, cmap='magma'):
    """
    Plots two side-by-side MI heatmaps for seqs1 vs. seqs2.
    """
    mi1 = compute_mi_matrix(seqs1)
    mi2 = compute_mi_matrix(seqs2)
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))
    im1 = ax1.imshow(mi1, origin='lower', cmap=cmap)
    ax1.set_title('MSA1 MI'); fig.colorbar(im1, ax=ax1)
    im2 = ax2.imshow(mi2, origin='lower', cmap=cmap)
    ax2.set_title('MSA2 MI'); fig.colorbar(im2, ax=ax2)
    plt.tight_layout()
    plt.show()


def plot_mi_with_upper_origin(mi, cmap='magma', tick_step=20, plotting_dir='plots/MSA', id=None):
    """
    Plot MI matrix so that (0,0) is in the top-left, and
    bottom-left shows (0, L-1), matching your mask orientation.
    """
    L = mi.shape[0]
    fig, ax = plt.subplots(figsize=(6,6))
    # origin='upper' places row 0 at the top
    im = ax.imshow(mi, cmap=cmap, origin='upper')
    
    # set ticks every tick_step
    ticks = np.arange(0, L, tick_step)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    
    # optionally label them
    ax.set_xticklabels(ticks)
    ax.set_yticklabels(ticks)
    
    # ax.set_xlabel('Position')
    # ax.set_ylabel('Position')
    # ax.set_title('MI heatmap (origin=upper)')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    print('Writing MI Matrix plot to:', f"{plotting_dir}/{id}_MI_matrix.png")
    plt.tight_layout()
    plt.savefig(f"{plotting_dir}/{id}_MI_matrix.png", dpi=300)
    plt.close()
    print('done.')
    # plt.show()

# 2) Covariation arc plot
def plot_covar_arcs(mi, threshold=0.1, figsize=(8,2)):
    """
    mi: L×L mutual information matrix
    threshold: draw arcs for mi>=threshold
    """
    L = mi.shape[0]
    fig, ax = plt.subplots(1,1, figsize=figsize)
    ax.set_xlim(0, L)
    ax.set_ylim(0, L/2)
    ax.axis('off')
    for i in range(L):
        for j in range(i+1, L):
            if mi[i,j] >= threshold:
                mid = (i+j)/2
                w = j - i
                arc = plt.Circle((mid, w/2), w/2,
                                 fill=False,
                                 linewidth=mi[i,j]*5,
                                 alpha=0.7)
                ax.add_patch(arc)
    plt.show()


# 3) R-scape wrapper
def run_rscape(stockholm_path, output_prefix=None):
    """
    Calls R-scape on a Stockholm file.
    Requires R-scape in your PATH.
    """
    out = output_prefix or stockholm_path.rsplit('.',1)[0]
    cmd = ['R-scape', stockholm_path, '-o', out]
    subprocess.run(cmd, check=True)
    print(f"R-scape results → {out}.rc, {out}.covariation.csv")


# 4) Per-column entropy & plot
def column_entropy(pwm):
    """
    pwm: pandas.DataFrame, rows=positions, cols=letters, values=probabilities
    returns: pandas.Series of length L
    """
    return pwm.apply(lambda row: entropy(row+1e-9, base=2), axis=1)

def plot_entropy(ent1, ent2, labels=('MSA1','MSA2')):
    """
    ent1, ent2: pandas.Series of per-position entropy
    """
    plt.figure(figsize=(10,3))
    plt.plot(ent1.values, label=labels[0])
    plt.plot(ent2.values, label=labels[1])
    plt.xlabel('Position')
    plt.ylabel('Entropy (bits)')
    plt.legend()
    plt.tight_layout()
    plt.show()


# 5) Sequence-distance histograms
def compute_hamming_dm(seqs):
    """
    Returns a skbio DistanceMatrix of pairwise Hamming distances.
    """
    return beta_diversity('hamming', seqs, ids=[str(i) for i in range(len(seqs))])

def plot_distance_histogram(dm1, dm2, bins=30):
    """
    dm1, dm2: skbio DistanceMatrix
    """
    d1 = dm1.condensed_form()
    d2 = dm2.condensed_form()
    plt.figure(figsize=(6,4))
    plt.hist(d1, bins=bins, alpha=0.6, label='MSA1')
    plt.hist(d2, bins=bins, alpha=0.6, label='MSA2')
    plt.xlabel('Hamming distance')
    plt.ylabel('Count')
    plt.legend()
    plt.tight_layout()
    plt.show()


# 6) t-SNE & UMAP embeddings
def one_hot_encode(seqs):
    """
    Simple A/C/G/U/- → 5-dim one-hot vectors, returns N×(L*5) array.
    """
    alphabet = ['A','C','G','U','-']
    mapping = {a:i for i,a in enumerate(alphabet)}
    arr = np.zeros((len(seqs), len(seqs[0]), 5), dtype=int)
    for i,s in enumerate(seqs):
        for j,ch in enumerate(s):
            k = mapping.get(ch,4)
            arr[i,j,k] = 1
    return arr.reshape(len(seqs), -1)

def plot_tsne_umap(seqs1, seqs2, perplexity=30, n_neighbors=15):
    """
    Renders t-SNE and UMAP plots overlaying points from both MSAs.
    """
    # prepare
    data = one_hot_encode(seqs1) 
    data2 = one_hot_encode(seqs2)
    X = np.vstack([data, data2])
    labels = (['MSA1']*len(data)) + (['MSA2']*len(data2))

    # t-SNE
    ts = TSNE(perplexity=perplexity)
    Z_ts = ts.fit_transform(X)
    # UMAP
    um = umap.UMAP(n_neighbors=n_neighbors)
    Z_um = um.fit_transform(X)

    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))
    for lab, col in zip(['MSA1','MSA2'], ['C0','C1']):
        mask = np.array(labels)==lab
        ax1.scatter(Z_ts[mask,0], Z_ts[mask,1], label=lab, alpha=0.6)
        ax2.scatter(Z_um[mask,0], Z_um[mask,1], label=lab, alpha=0.6)
    ax1.set_title('t-SNE'); ax2.set_title('UMAP')
    for ax in (ax1, ax2):
        ax.legend(); ax.axis('off')
    plt.tight_layout()
    plt.show()


def trim_to_reference(alignment, ref_id=None, gap_char='-'):
    """
    Given a Bio.AlignIO.MultipleSeqAlignment, return a list of
    sequences (strings) *only* keeping columns where the reference
    has a non-gap.  If ref_id is None, uses the first record.
    """
    # pick the reference
    if ref_id is None:
        ref = alignment[0]
    else:
        ref = next(rec for rec in alignment if rec.id == ref_id)
    ref_seq = str(ref.seq)
    
    # find columns where ref != gap
    keep = [i for i, c in enumerate(ref_seq) if c != gap_char]
    
    # extract those columns from every sequence
    trimmed = []
    for rec in alignment:
        seq = str(rec.seq)
        trimmed.append(''.join(seq[i] for i in keep))
    return trimmed

def pwm_from_seqs(sequences, ignore_gaps=True):
    """
    Build a position‐weight (probability) matrix from a list of equal-length
    sequences, summing each column to 1.
    """
    counts = lm.alignment_to_matrix(
        sequences,
        to_type='counts',
        characters_to_ignore='-' if ignore_gaps else None
    )
    # normalize rows → probabilities
    pwm = counts.div(counts.sum(axis=1), axis=0)
    return pwm

def plot_side_by_side_logos(
    msa1_path,
    msa2_path,
    fmt='fasta',
    ref_id=None,
    figsize=(20, 6),
    color_scheme='classic',
    show=False,
    # NEW: where & how to save the individual logos
    plotting_dir: Path = None,
    id1: str = None,
    id2: str = None,
    # OPTIONAL: also save the combined (side-by-side) figure
    save_side_by_side: bool = False,
    side_by_side_name: str = None,
):
    """
    Load two alignments, trim to the same reference coordinates,
    build PWMs, and plot them one above the other.

    NEW:
      - Saves *individual* logos for MSA1 and MSA2 to plotting_dir using id1/id2.
      - Optionally saves the combined (side-by-side) panel as well.
    """
    # -------- load --------
    aln1 = AlignIO.read(msa1_path, fmt)
    aln2 = AlignIO.read(msa2_path, fmt)

    # -------- trim to reference coords --------
    seqs1 = trim_to_reference(aln1, ref_id)
    seqs2 = trim_to_reference(aln2, ref_id)

    # -------- PWMs (rows sum to 1) --------
    pwm1 = pwm_from_seqs(seqs1)
    pwm2 = pwm_from_seqs(seqs2)

    # -------- ensure output dir & file names --------
    if plotting_dir is not None:
        plotting_dir = Path(plotting_dir)
        plotting_dir.mkdir(exist_ok=True, parents=True)

    stem1 = Path(msa1_path).stem
    stem2 = Path(msa2_path).stem
    file_id1 = id1 if id1 is not None else stem1
    file_id2 = id2 if id2 is not None else stem2

    # -------- individual logo: MSA1 (save only) --------
    if plotting_dir is not None:
        fig1, ax1 = plt.subplots(1, 1, figsize=(min(figsize[0], 12), 3))
        lm.Logo(pwm1, ax=ax1, color_scheme=color_scheme, stack_order='small_on_top')
        # ax1.set_title(stem1)
        ax1.set_ylim(0, 1)
        ax1.axis('off')
        out1 = plotting_dir / f"{file_id1}_logo.png"
        fig1.tight_layout()
        fig1.savefig(out1, dpi=300)
        plt.close(fig1)

    # -------- individual logo: MSA2 (save only) --------
    if plotting_dir is not None:
        fig2, ax2 = plt.subplots(1, 1, figsize=(min(figsize[0], 12), 3))
        lm.Logo(pwm2, ax=ax2, color_scheme=color_scheme, stack_order='small_on_top')
        # ax2.set_title(stem2)
        ax2.set_ylim(0, 1)
        ax2.axis('off')
        out2 = plotting_dir / f"{file_id2}_logo.png"
        fig2.tight_layout()
        fig2.savefig(out2, dpi=300)
        plt.close(fig2)

    # -------- side-by-side panel (as before) --------
    fig, (axA, axB) = plt.subplots(2, 1, sharex=True, figsize=figsize)
    lm.Logo(pwm1, ax=axA, color_scheme=color_scheme, stack_order='small_on_top')
    axA.set_title(stem1)
    axA.set_ylim(0, 1)
    axA.axis('off')

    lm.Logo(pwm2, ax=axB, color_scheme=color_scheme, stack_order='small_on_top')
    axB.set_title(stem2)
    axB.set_ylim(0, 1)
    axB.axis('off')

    fig.tight_layout()

    if save_side_by_side and plotting_dir is not None:
        panel_name = side_by_side_name or f"{file_id1}_vs_{file_id2}_logos.png"
        fig.savefig(plotting_dir / panel_name, dpi=300)

    if show:
        plt.show()
    plt.close(fig)


def plot_msa_nuc_gap_heatmap(
    msa_path,
    fmt="fasta",
    ref_id=None,
    color_map="viridis",
    figsize=(12, 3),
    plotting_dir=None,
    file_id=None,
    show=True
):
    """
    Plot a heatmap of nucleotide + gap occurrences per position in a gapped MSA,
    restricted to positions of the reference sequence (ref_id).
    
    Rows: A, C, G, U, -
    Cols: Positions in reference (ungapped coords)
    Color: Frequency of each base/gap at each position.
    """
    # Load MSA
    aln = AlignIO.read(msa_path, fmt)
    aln_df = pd.DataFrame([list(str(rec.seq).upper()) for rec in aln], index=[rec.id for rec in aln])

    # Determine reference sequence (query)
    if ref_id is None:
        ref_seq = aln[0].seq.upper()
    else:
        ref_seq = [rec.seq.upper() for rec in aln if rec.id == ref_id][0]

    # Identify columns where the reference has a non-gap character
    ref_mask = np.array([base != "-" for base in ref_seq])

    # Trim MSA to reference coordinates (like logos)
    trimmed_df = aln_df.loc[:, ref_mask]

    # Define nucleotide + gap alphabet
    alphabet = ["A", "C", "G", "U", "-"]

    # Count occurrences at each position
    counts = []
    for base in alphabet:
        counts.append((trimmed_df == base).sum(axis=0).values)
    counts = np.array(counts, dtype=float)

    # Normalize to frequencies
    freqs = counts / counts.sum(axis=0, keepdims=True)

    # Make DataFrame for plotting
    pos_labels = [f"{i+1}" for i in range(freqs.shape[1])]
    df_plot = pd.DataFrame(freqs, index=alphabet, columns=pos_labels)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(df_plot, cmap=color_map, ax=ax, cbar_kws={"label": "Frequency"}, vmin=0, vmax=1)
    ax.set_xlabel("Position in reference sequence")
    ax.set_ylabel("Nucleotide / Gap")
    ax.set_title(Path(msa_path).stem)

    plt.tight_layout()

    # Save if requested
    if plotting_dir is not None:
        plotting_dir = Path(plotting_dir)
        plotting_dir.mkdir(exist_ok=True, parents=True)
        file_id = file_id or Path(msa_path).stem
        fig.savefig(plotting_dir / f"{file_id}_nuc_gap_heatmap.png", dpi=300)

    if show:
        plt.show()
    plt.close(fig)


def plot_difference_logo(msa1_path, msa2_path, fmt='fasta',
                         ref_id=None, figsize=(20,3), show=True):
    """
    Compute pwm1 − pwm2 and plot the difference-logo in a single row.
    Positive letters mean enriched in MSA1, negatives enriched in MSA2.
    """
    # 1) load
    aln1 = AlignIO.read(msa1_path, fmt)
    aln2 = AlignIO.read(msa2_path, fmt)

    # 2) trim to reference coords
    seqs1 = trim_to_reference(aln1, ref_id)
    seqs2 = trim_to_reference(aln2, ref_id)

    # 3) build PWMs
    pwm1 = pwm_from_seqs(seqs1)
    pwm2 = pwm_from_seqs(seqs2)

    # 4) unify columns (so there are no NaNs)
    all_letters = sorted(set(pwm1.columns).union(pwm2.columns))
    pwm1 = pwm1.reindex(columns=all_letters, fill_value=0)
    pwm2 = pwm2.reindex(columns=all_letters, fill_value=0)

    # 5) subtract
    diff = pwm1 - pwm2

    # 6) plot
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    lm.Logo(
        diff,
        ax=ax,
        color_scheme='classic',
        stack_order='small_on_top',
        # if you *do* want to see gaps or ambig codes as blanks, you could:
        # allow_nan=True
    )
    ax.set_ylim(-1, 1)
    ax.axis('off')
    ax.set_title(f"{Path(msa1_path).stem} – {Path(msa2_path).stem}")
    plt.tight_layout()
    if show:
        plt.show()

from pathlib import Path
from Bio import AlignIO
import logomaker as lm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_logo_and_nuc_gap_heatmap(
    msa_path,
    fmt="fasta",
    ref_id=None,                   # if None, use the first sequence
    color_scheme="classic",        # logomaker color scheme
    heatmap_cmap="viridis",        # heatmap colormap
    combined_figsize=(20, 6),
    logo_figsize=(20, 3),
    heatmap_figsize=(20, 3),
    plotting_dir="plots/MSA",
    file_id=None,                  # defaults to msa_path.stem
    show=False,
    save_individual=True,
    xtick_step=20                   # tick every N positions
):
    """
    Make a combined (logo + nucleotide/gap heatmap) figure for a gapped MSA,
    trimming columns to positions where the reference (query) has non-gaps.

    Saves:
      <file_id>_logo+heatmap.png
      <file_id>_logo.png
      <file_id>_nuc_gap_heatmap.png
    """
    plotting_dir = Path(plotting_dir)
    plotting_dir.mkdir(exist_ok=True, parents=True)
    file_id = file_id or Path(msa_path).stem

    # ---------- load alignment ----------
    aln = AlignIO.read(msa_path, fmt)

    # ---------- choose reference (default = first) ----------
    if ref_id is None:
        ref_seq = str(aln[0].seq).upper()
    else:
        ref_seq = None
        for rec in aln:
            if rec.id == ref_id:
                ref_seq = str(rec.seq).upper()
                break
        if ref_seq is None:
            raise ValueError(f"ref_id '{ref_id}' not found in {msa_path}")

    # ---------- build trimmed sequences ----------
    keep_cols = [i for i, c in enumerate(ref_seq) if c != "-"]
    trimmed_seqs = [
        "".join(str(rec.seq).upper()[i] for i in keep_cols)
        for rec in aln
    ]

    # ---------- PWM for logo ----------
    alphabet = ['A', 'C', 'G', 'U']  # gaps ignored for logo
    counts = []
    for pos in range(len(trimmed_seqs[0])):
        col = [s[pos] for s in trimmed_seqs]
        counts.append({base: col.count(base) for base in alphabet})
    counts_df = pd.DataFrame(counts, columns=alphabet)
    pwm = counts_df.div(counts_df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    # ---------- nucleotide+gap frequency heatmap ----------
    aln_df = pd.DataFrame([list(s) for s in trimmed_seqs])
    alphabet5 = ['A', 'C', 'G', 'U', '-']
    mat_counts = np.vstack([(aln_df == base).sum(axis=0).values for base in alphabet5]).astype(float)
    col_sums = mat_counts.sum(axis=0, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        freqs = np.divide(mat_counts, col_sums, out=np.zeros_like(mat_counts), where=(col_sums != 0))
    pos_labels = [str(i+1) for i in range(freqs.shape[1])]
    heat_df = pd.DataFrame(freqs, index=alphabet5, columns=pos_labels)

    # ---------- combined figure ----------
    fig, (ax_logo, ax_heat) = plt.subplots(
        2, 1, figsize=combined_figsize, gridspec_kw={'height_ratios': [2, 1]}, sharex=True
    )

    # Logo
    lm.Logo(pwm, ax=ax_logo, color_scheme=color_scheme, stack_order='small_on_top')
    ax_logo.set_ylim(0, 1)
    ax_logo.set_ylabel("Freq.")
    for spine in ('top', 'right'):
        ax_logo.spines[spine].set_visible(False)
    ax_logo.tick_params(axis='y', length=0)
    ax_logo.set_xticks([])  # hide x ticks for top logo

    # Heatmap
    sns.heatmap(
        heat_df,
        cmap=heatmap_cmap,
        ax=ax_heat,
        cbar_kws={'label': 'Frequency'},
        vmin=0, vmax=1
    )
    ax_heat.set_xlabel("Position in reference sequence")
    ax_heat.set_ylabel("Nucleotide / Gap")

    # set x ticks every xtick_step
    total_len = heat_df.shape[1]
    xticks = list(range(0, total_len, xtick_step))
    ax_heat.set_xticks([x+0.5 for x in xticks])  # center ticks
    ax_heat.set_xticklabels([str(x+1) for x in xticks])

    plt.tight_layout()
    combined_out = plotting_dir / f"{file_id}_logo+heatmap.png"
    fig.savefig(combined_out, dpi=300)
    if show:
        plt.show()
    plt.close(fig)

    # ---------- individual saves (optional) ----------
    if save_individual:
        # Logo only
        fig_l, ax_l = plt.subplots(1, 1, figsize=logo_figsize)
        lm.Logo(pwm, ax=ax_l, color_scheme=color_scheme, stack_order='small_on_top')
        ax_l.set_ylim(0, 1)
        ax_l.set_ylabel("Freq.")
        for spine in ('top', 'right'):
            ax_l.spines[spine].set_visible(False)
        ax_l.tick_params(axis='y', length=0)
        ax_l.set_xticks([])
        plt.tight_layout()
        fig_l.savefig(plotting_dir / f"{file_id}_logo.png", dpi=300)
        plt.close(fig_l)

        # Heatmap only
        fig_h, ax_h = plt.subplots(1, 1, figsize=heatmap_figsize)
        sns.heatmap(
            heat_df,
            cmap=heatmap_cmap,
            ax=ax_h,
            cbar_kws={'label': 'Frequency'},
            vmin=0, vmax=1
        )
        ax_h.set_xlabel("Position in reference sequence")
        ax_h.set_ylabel("Nucleotide / Gap")
        total_len = heat_df.shape[1]
        xticks = list(range(0, total_len, xtick_step))
        ax_h.set_xticks([x+0.5 for x in xticks])
        ax_h.set_xticklabels([str(x+1) for x in xticks])
        plt.tight_layout()
        fig_h.savefig(plotting_dir / f"{file_id}_nuc_gap_heatmap.png", dpi=300)
        plt.close(fig_h)

    print(f"Saved combined:   {combined_out}")
    if save_individual:
        print(f"Saved logo:       {plotting_dir / f'{file_id}_logo.png'}")
        print(f"Saved heatmap:    {plotting_dir / f'{file_id}_nuc_gap_heatmap.png'}")


def compare_sequence_diversity(seqs1, seqs2, bins=30):
    """
    Compare the sequence‐space diversity of two MSAs by computing
    all pairwise Hamming distances.

    Parameters
    ----------
    seqs1, seqs2 : list of str
        Lists of equal‐length aligned sequences (e.g. your trimmed seqs1, seqs2).
    bins : int
        Number of bins for the histogram.

    Returns
    -------
    dict
        {
          'distances1': 1D numpy array of pairwise distances from seqs1,
          'distances2': 1D numpy array of pairwise distances from seqs2
        }
    """
    def pairwise_hamming(seqs):
        n = len(seqs)
        if n < 2:
            return np.array([])
        L = len(seqs[0])
        dists = []
        for i, j in itertools.combinations(range(n), 2):
            # fraction of mismatches
            mismatches = sum(c1 != c2 for c1, c2 in zip(seqs[i], seqs[j]))
            dists.append(mismatches / L)
        return np.array(dists)

    # compute
    d1 = pairwise_hamming(seqs1)
    d2 = pairwise_hamming(seqs2)

    # summary stats
    if d1.size:
        print(f"MSA1: mean distance = {d1.mean():.4f}, std = {d1.std():.4f}")
    else:
        print("MSA1: not enough sequences to compute distances")
    if d2.size:
        print(f"MSA2: mean distance = {d2.mean():.4f}, std = {d2.std():.4f}")
    else:
        print("MSA2: not enough sequences to compute distances")

    # plot
    plt.figure(figsize=(6, 4))
    if d1.size:
        plt.hist(d1, bins=bins, alpha=0.6, label='MSA1')
    if d2.size:
        plt.hist(d2, bins=bins, alpha=0.6, label='MSA2')
    plt.xlabel('Hamming distance')
    plt.ylabel('Count')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return {'distances1': d1, 'distances2': d2}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Analyse gapped RNA alignments.")
    parser.add_argument('--fasta1', type=Path, help='Path to first FASTA file')
    parser.add_argument('--fasta2', type=Path, help='Path to second FASTA file')
    parser.add_argument('--plotting_dir', type=Path, default='plots/MSA', help='Output directory for plots')
    parser.add_argument('--fmt', type=str, default='fasta',
                        help='Format of the input files (default: fasta)')
    parser.add_argument('--ref_id', type=str, default=None,
                        help='Reference sequence ID to trim to (default: first record)')
    parser.add_argument('--cmap', type=str, default='PRGn',
                        help='Colormap for MI heatmaps (default: PRGn)')
    parser.add_argument('--bins', type=int, default=40,
                        help='Number of bins for distance histograms (default: 30)')

    args = parser.parse_args()

    fasta_file1 = args.fasta1
    fasta_file2 = args.fasta2

    plotting_dir = args.plotting_dir
    plotting_dir.mkdir(exist_ok=True, parents=True)

    fmt = args.fmt
    ref_id = args.ref_id
    cmap = args.cmap
    bins = args.bins

    # # json_file = Path('1EUQ_RNAformerafc23_data.json')
    # fasta_file1 = Path('1SJ4_rnaprotein_rnaformerafc23_rnaprotein_gapped_rna_alignment.fasta')
    # # fasta_file2 = Path('1SJ4_rnaprotein_alphafold3_gapped_rna_alignment.fasta')
    # # fasta_file1 = Path('4p8z_gapped_rna_alignment_rnaformer.fasta')
    # # fasta_file2 = Path('4p8z_gapped_rna_alignment_spotrna.fasta')
    # fasta_file2 = Path('1sj4_gapped_rna_alignment_rnafold.fasta')
    # # data_dir = Path('RNAformerafc23_rnaprotein_nat_meth_data_files')
    # data_dir1 = Path('../synthetic_msa_gapped_rna_alignments/rnaprotein/rnaformerafc23_rnaprotein')
    # # data_dir1 = Path('gapped_alignments')
    # data_dir2 = Path('gapped_alignments')
    # plotting_dir = Path('msa_test_plots')
    # plotting_dir.mkdir(exist_ok=True, parents=True)
    # algorithm = 'sMSA'  # 'AlphaFold 3'
    # # algorithm = 'AlphaFold 3'  # 'AlphaFold 3'
    
    # data = load_json(data_dir / json_file)
    # msa = get_msa(data)

    plot_logo_and_nuc_gap_heatmap(
        Path(fasta_file1),
        fmt="fasta",
        ref_id=None,                   # if None, use the first sequence
        color_scheme="classic",        # logomaker color scheme
        heatmap_cmap="viridis",        # heatmap colormap
        combined_figsize=(20, 6),
        logo_figsize=(20, 3),
        heatmap_figsize=(20, 3),
        plotting_dir=args.plotting_dir,
        file_id=None,                  # defaults to msa_path.stem
        show=False,
        save_individual=True,
        xtick_step=20
    )

    plot_logo_and_nuc_gap_heatmap(
        Path(fasta_file2),
        fmt="fasta",
        ref_id=None,                   # if None, use the first sequence
        color_scheme="classic",        # logomaker color scheme
        heatmap_cmap="viridis",        # heatmap colormap
        combined_figsize=(20, 6),
        logo_figsize=(20, 3),
        heatmap_figsize=(20, 3),
        plotting_dir=args.plotting_dir,
        file_id=None,                  # defaults to msa_path.stem
        show=False,
        save_individual=True,
        xtick_step=20
    )



    # plot_side_by_side_logos(
    #     Path(fasta_file1),
    #     Path(fasta_file2),
    #     fmt=fmt,
    #     # ref_id='1EUY'  # or None to take first record
    #     plotting_dir=args.plotting_dir,
    #     id1=fasta_file1.stem,
    #     id2=fasta_file2.stem,
    # )

    # difference logo:
    plot_difference_logo(
        Path(fasta_file1),
        Path(fasta_file2),
        fmt=fmt,
        # ref_id='1EUY'  # or None to take first record
    )

    aln1 = AlignIO.read(Path(fasta_file1), fmt)
    aln2 = AlignIO.read(Path(fasta_file2), fmt)

    AlignIO.write(aln1, Path(plotting_dir, f"{Path(fasta_file1).stem}.stk"), "stockholm")
    AlignIO.write(aln2, Path(plotting_dir, f"{Path(fasta_file2).stem}.stk"), "stockholm")
    
    # 2) Trim every sequence to the non-gap columns of your reference
    #     (so that both seqs1 and seqs2 line up on the exact same query coords)
    seqs1 = trim_to_reference(aln1, ref_id=None, gap_char='-')
    seqs2 = trim_to_reference(aln2, ref_id=None, gap_char='-')

    print(len(seqs1), "sequences in MSA1")
    print(len(seqs2), "sequences in MSA2")

    # MI heatmaps
    # plot_mi_heatmaps(seqs1, seqs2)

    id1 = fasta_file1.stem
    id2 = fasta_file2.stem

    plot_mi_with_upper_origin(compute_mi_matrix(seqs1), cmap=cmap, plotting_dir=args.plotting_dir, id=id1)  #  cmap='PiYG')  # 'magma')
    plot_mi_with_upper_origin(compute_mi_matrix(seqs2), cmap=cmap, plotting_dir=args.plotting_dir, id=id2)   # ,   # cmap='PiYG')  # 'magma')
        
    # t-SNE & UMAP
    # plot_tsne_umap(seqs1, seqs2)

    # results = compare_sequence_diversity(seqs1, seqs2, bins=bins)
    # Access raw distances if needed:
    # d1 = results['distances1']
    # d2 = results['distances2']
