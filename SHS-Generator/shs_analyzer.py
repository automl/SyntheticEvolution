import argparse
import logging
import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional, Sequence, Tuple

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def parse_input(path: str) -> Tuple[str, np.ndarray]:
    with open(path, "r") as f:
        content = f.read()
    try:
        data = json.loads(content)
        if not isinstance(data, dict):
            raise json.decoder.JSONDecodeError("not an object", content, 0)
        all_fasta = [msa for c in data.get("sequences", []) if "rna" in c and (msa := c["rna"].get("unpairedMsa"))]
        if len(all_fasta) == 0:
            logging.error("No unpaired Msa found for any rna in the json file %s", input)
            raise ValueError("Invalid Input, aborting.")
        if len(all_fasta) > 1:
            logging.warning("Input json contains multiple rna chains with unpaired msa, using the first one")
        fasta = all_fasta[0]
    except json.decoder.JSONDecodeError:
        fasta = content
    seq, *raw_msa = fasta.splitlines()[1::2]
    if len(raw_msa) == 0:
        logging.error("Recovered MSA is empty")
        raise ValueError("Invalid Input, aborting.")
    clean_msa = [[c for c in seq if not c.islower()] for seq in raw_msa]
    if not all(len(s) == len(clean_msa[0]) for s in clean_msa):
        logging.error("Sequences in MSA have different lengths or the json was invalid")
        raise ValueError("Invalid Input, aborting.")
    msa_array = np.array(clean_msa)
    if not np.all(np.isin(msa_array, list("AGCU-"))):
        invalid_chars = np.unique(msa_array[~np.isin(msa_array, list("AGCU-"))])
        logging.error(f"Invalid characters in MSA: {invalid_chars}")
        raise ValueError("Invalid Input, aborting.")
    return (seq, msa_array)


def structure_matrix(pairs: Dict[int, int], seq_len: int,
                     paired_value, unpaired_value,
                     paired_diag, unpaired_diag) -> np.ndarray:
    mat = np.full((seq_len, seq_len), unpaired_value, dtype=float)
    for i, j in pairs.items():
        if 0 <= i < seq_len and 0 <= j < seq_len:
            mat[i, j] = paired_value
    for i in range(seq_len):
        mat[i, i] = paired_diag if i in pairs else unpaired_diag
    return mat


def compute_cov(input_path: str, use_corrcoef: bool = False) -> Tuple[str, np.ndarray, np.ndarray]:
    seq, msa = parse_input(input_path)
    seq_array = np.array(list(seq))
    similarities = (msa == seq_array).astype(int)
    if use_corrcoef:
        covariance = np.corrcoef(similarities, rowvar=False)
    else:
        covariance = np.cov(similarities, rowvar=False)
    return seq, msa, covariance

def plot_matrix(mat: np.ndarray, title: str = "Recovered (covariation)",
                show: bool = True, out: Optional[str] = None,
                vmin: Optional[float] = None, vmax: Optional[float] = None,
                show_values = False):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(mat.shape[0], mat.shape[1]))
    im = ax.imshow(mat, cmap="viridis", vmin=vmin, vmax=vmax)

    if show_values:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                text = ax.text(j, i, round(mat[i, j], 2), ha="center", va="center", color="w")
    
    ax.set_title(title)
    # fig.colorbar(im, ax=ax) 
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        logging.info("Figure written to: %s", out)
    if show:
        plt.show()
    return fig


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Recover an interaction matrix from a SHS.")
    p.add_argument("-i", "--input", help="Path to AF3 JSON input to a FASTA MSA File.", required=True)
    p.add_argument("-o", "--out", help="Output PNG path.")
    p.add_argument("-s", "--show", action="store_true", help="Display the figure.")
    p.add_argument("-t", "--title", default="", help="Figure title.")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    seq, msa, covariance = compute_cov(args.input)
    print(covariance.shape)
    title = args.title or "Recovered (covariation)"
    plot_matrix(covariance, title=title, show=args.show, out=args.out)


if __name__ == "__main__":
    main()
