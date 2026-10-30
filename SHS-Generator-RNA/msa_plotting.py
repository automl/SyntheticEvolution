"""MSA feature visualisation for shs_generator.

Kept separate so the core generator does not depend on matplotlib/numpy: those
are only needed when --plot is used, and shs_generator imports this module
lazily at that point.
"""
import logging
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt


def pad_lowercase(seq: str, target_length: int) -> str:
    return ''.join([c for c in seq if not c.islower()]).ljust(target_length, '-')


def bool_left_deletion(seq: str) -> np.ndarray:
    return np.array([int(s == '-') for s in seq], dtype=bool)


def get_deletion_values(seq: np.ndarray) -> np.ndarray:
    return np.array(
        [(2 / np.pi) * np.arctan(np.sum(seq[:i]) / 3) if i > 0 else (2 / np.pi) * np.arctan(seq[-1] / 3)
         for i in range(len(seq))],
        dtype=np.float16
    )


def plot_final_features(msa: List[str], rna_seq: str, pairs: Dict[int, int],
                        pdb_id: str, show_plot: bool = False) -> Path:
    """Render the MSA feature panels and save them to
    msa_final_plots/msa_final_<pdb_id>.png (relative to the working directory).
    Returns the written path."""
    msa_sub = msa if len(msa) <= 100 else msa[:100]
    L = len(pad_lowercase(msa_sub[0], len(rna_seq)))
    cleaned_msa = [pad_lowercase(seq, len(rna_seq)) for seq in msa_sub]

    del_bool = np.array([bool_left_deletion(s) for s in cleaned_msa])

    del_smooth = np.array([get_deletion_values(row) for row in del_bool])

    query = cleaned_msa[0]
    nseq = len(cleaned_msa)
    diff_matrix = np.zeros((nseq, L), dtype=int)
    for i, seq in enumerate(cleaned_msa):
        for j in range(L):
            if seq[j] != query[j]:
                diff_matrix[i, j] = 1

    conservation = np.mean(diff_matrix == 0, axis=0)

    stem_mask = [1 if i in pairs else 0 for i in range(L)]
    stem_mask_row = np.array([stem_mask])

    nucleotides = ['A', 'U', 'G', 'C', '-']
    onehot = np.zeros((len(nucleotides), L))
    for j in range(L):
        col = [s[j] for s in cleaned_msa]
        total = len(col)
        for idx, nt in enumerate(nucleotides):
            count = sum(1 for c in col if c.upper() == nt)
            onehot[idx, j] = count / total


    fig, axs = plt.subplots(3, 2, figsize=(15, 15))

    im_a = axs[0, 0].imshow(del_bool, aspect='auto', interpolation='nearest', cmap='gray')
    axs[0, 0].set_title("Deletion Boolean Heatmap")
    axs[0, 0].set_ylabel("Sequence index")
    axs[0, 0].set_xlabel("Position")
    fig.colorbar(im_a, ax=axs[0, 0], orientation='vertical')

    im_b = axs[0, 1].imshow(del_smooth, aspect='auto', interpolation='nearest', cmap='viridis')
    axs[0, 1].set_title("Smoothed Deletion Heatmap")
    axs[0, 1].set_xlabel("Position")
    fig.colorbar(im_b, ax=axs[0, 1], orientation='vertical')

    im_c = axs[1, 0].imshow(diff_matrix, aspect='auto', interpolation='nearest', cmap='viridis')
    axs[1, 0].set_title("Substitution Differences Heatmap")
    axs[1, 0].set_ylabel("Sequence index")
    axs[1, 0].set_xlabel("Position")
    fig.colorbar(im_c, ax=axs[1, 0], orientation='vertical')

    axs[1, 1].plot(range(L), conservation, label="Conservation", color='blue')
    axs[1, 1].set_title("Conservation Score per Position")
    axs[1, 1].set_xlabel("Position")
    axs[1, 1].set_ylabel("Conservation")
    for j in range(L):
        if stem_mask[j] == 1:
            axs[1, 1].axvspan(j - 0.5, j + 0.5, color='orange', alpha=0.6)
    axs[1, 1].legend()

    im_e = axs[2, 0].imshow(stem_mask_row, aspect='auto', interpolation='nearest', cmap='Reds')
    axs[2, 0].set_title("Stem Mask Row (1 = Stem, 0 = Loop)")
    axs[2, 0].set_yticks([])
    axs[2, 0].set_xlabel("Position")
    fig.colorbar(im_e, ax=axs[2, 0], orientation='vertical')

    im_f = axs[2, 1].imshow(onehot, aspect='auto', interpolation='nearest', cmap='plasma')
    axs[2, 1].set_title("One-Hot Frequency Heatmap")
    axs[2, 1].set_yticks(range(len(nucleotides)))
    axs[2, 1].set_yticklabels(nucleotides)
    axs[2, 1].set_xlabel("Position")
    fig.colorbar(im_f, ax=axs[2, 1], orientation='vertical')

    plt.tight_layout()
    if show_plot:
        plt.show()
    final_dir = Path("msa_final_plots")
    final_dir.mkdir(parents=True, exist_ok=True)
    final_filename = final_dir / f"msa_final_{pdb_id}.png"
    plt.savefig(final_filename)
    plt.close()
    logging.info("Final MSA features plot saved to: %s", final_filename)
    return final_filename
