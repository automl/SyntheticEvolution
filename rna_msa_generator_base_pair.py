#!/usr/bin/env python
import argparse
import random
import json
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import subprocess

import RnaBench

from pathlib import Path
from typing import Dict, List, Any

from RnaBench.lib.rna_folding_algorithms.rnafold import RNAFold
from RnaBench.lib.rna_folding_algorithms.DL.spotrna import SpotRna
from RnaBench.lib.utils import pairs2db

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

PAIR_MUTATIONS: Dict[str, List[str]] = {
    'CU': ['GU', 'AU', 'CG'],
    'CA': ['UA', 'CG'],
    'GA': ['UA', 'GC', 'GU'],
    'CC': ['CG', 'GC'],
    'AA': ['AU', 'UA'],
    'UU': ['AU', 'UA'],
    'GG': ['GC', 'CG', 'GU', 'UG'],
    'GC': ['AU', 'CG', 'GC'],
    'CG': ['GC', 'AU'],
    'AU': ['UA', 'GC'],
    'UA': ['AU', 'CG'],
    'GU': ['GC', 'AU', 'UG'],
    'UG': ['GC', 'AU', 'GU']
}
WC_PAIRS: List[str] = ['AU', 'UA', 'GC', 'CG']

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


def cluster_deletion_mean(deletion_values: np.ndarray) -> np.ndarray:
    return np.array(
        [(2 / np.pi) * np.arctan(np.mean(deletion_values[:, i]) / 3)
         for i in range(deletion_values.shape[1])],
        dtype=np.float16
    )


def load_json(json_path: str) -> Any:
    try:
        with open(json_path, "r") as f:
            return json.load(f)
    except Exception as e:
        logging.error("Failed to load JSON from %s: %s", json_path, e)
        raise


def convert_pred_pairs(pred_pairs: Any) -> Dict[int, int]:
    if isinstance(pred_pairs, dict):
        return pred_pairs
    elif isinstance(pred_pairs, list):
        new_dict = {}
        for item in pred_pairs:
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                i, j = int(item[0]), int(item[1])
                new_dict[i] = j
                new_dict[j] = i
        return new_dict
    else:
        raise ValueError("Unknown format for predicted pairs: {}".format(type(pred_pairs)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate RNA MSA and AF3-compatible JSON.")
    # I/O and structure parameters.
    io_group = parser.add_argument_group("I/O and Structure")
    io_group.add_argument('--structure_predictor', type=str, default=None,
                          help="Predictor for secondary structure (base pairs output). Options: rnafold, spotrna, rnaformer, dssr")
    io_group.add_argument('--structure', type=str, required=False, default=None,
                          help="Dot-bracket secondary structure (if provided, will be converted to base pairs)")
    io_group.add_argument('--rna-seq', type=str, required=False,
                          help="RNA sequence (used for MSA + JSON)")
    io_group.add_argument('--protein-seq', type=str, required=False,
                          help="Protein sequence (for JSON)")
    io_group.add_argument('--input_json_path', type=str, required=False,
                          help="A JSON file in Alphafold webservice format that will be adapted for custom RNA MSA")
    io_group.add_argument('--pdb_id', type=str, required=False,
                          help="PDB ID (used in JSON name)")
    io_group.add_argument('--output_json_dir', type=str, default="custom_msa_json_output")
    # Mutation parameters. TODO: Default to config 80
    mut_group = parser.add_argument_group("Mutation Parameters")
    mut_group.add_argument('-N', type=int, default=20, help="Number of sequences in the MSA")
    mut_group.add_argument('--insertion-prob-loop', type=float, default=0.2)
    mut_group.add_argument('--deletion-prob-loop', type=float, default=0.8)
    mut_group.add_argument('--insertion-prob-stem', type=float, default=0.01)
    mut_group.add_argument('--deletion-prob-stem', type=float, default=0.01)
    mut_group.add_argument('--wobble-prob', type=float, default=0.1)
    mut_group.add_argument('--long-insertion-prob', type=float, default=0.05)
    mut_group.add_argument('--long-deletion-prob', type=float, default=0.05)
    mut_group.add_argument('--max-insertion-fraction', type=float, default=0.1,
                           help="Max insertion length as fraction of RNA length")
    mut_group.add_argument('--max-deletion-fraction', type=float, default=0.1,
                           help="Max deletion length as fraction of RNA length")
    # Parameter to help preserve stem structure.
    mut_group.add_argument('--stem-keep-prob', type=float, default=0.99,
                           help="Probability of keeping paired (stem) residues unchanged to provide strong structural hints")
    # Additional options.
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--max_chains', type=int, default=None)
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--print_msa', action='store_true')
    parser.add_argument('--show_plot', action='store_true')
    return parser.parse_args()

class MsaGenerator:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        if self.args.seed is not None:
            random.seed(self.args.seed)
        Path(self.args.output_json_dir).mkdir(parents=True, exist_ok=True)

    def pair_indices(self, dotbracket: str) -> Dict[int, int]:
        stack = []
        pairs: Dict[int, int] = {}
        for i, char in enumerate(dotbracket):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if not stack:
                    logging.error("Unbalanced structure at position %d", i)
                    continue
                j = stack.pop()
                pairs[i] = j
                pairs[j] = i
        return pairs

    def get_structure(self, sequence: str) -> Dict[int, int]:
        seq = sequence.upper()
        # Use provided structure predictor if available.
        if self.args.structure_predictor:
            if self.args.structure:
                logging.warning("Both structure predictor and structure provided. Using provided structure.")
                logging.info("Provided dot-bracket structure: %s", self.args.structure)
                return self.pair_indices(self.args.structure)
            else:
                predictor = self.args.structure_predictor.lower()
                pdb_id = self.args.pdb_id.lower()[:4] if self.args.pdb_id else None
                if predictor == "rnafold":
                    model = RNAFold()
                    pred_pairs, _ = model(seq)
                elif predictor == "spotrna":
                    model = SpotRna()
                    pred_pairs = model(seq)
                elif predictor == "rnaformer":
                    df = pd.read_pickle('data/rnaformer_predictions.pkl')
                    d = df[df['pdb_id'].str.lower() == pdb_id]
                    if d.empty and df[df['sequence'] == seq]['pairs'].empty:
                        sample = {'Id': pdb_id.upper(), 'sequence': seq}
                        df = pd.DataFrame([sample])
                        
                        import pickle
                        
                        with open(Path(f'RNAformer/datasets/{pdb_id.upper()}.pkl'), 'wb') as f:
                            pickle.dump(df, f)
                        
                        subprocess.run(['./RNAformer/run_cpu.sh', 'infer_RNAformer.py', '-c', '1', '-m', 'models/af3_like_finetune', '-p', f'datasets/{pdb_id.upper()}.pkl'])
                        
                        pred = pd.read_pickle(f'RNAformer/datasets/{pdb_id.upper()}_processed.pkl')
                
                        pred_pairs = sorted(pred['pairs'].values[0], key=lambda x: x[0])
                    elif not d.empty:
                        pred_pairs = df[df['pdb_id'].str.lower() == pdb_id]['pairs'].values[0]
                    else:
                        pred_pairs = df[df['sequence'] == seq]['pairs'].values[0]
                    pred_pairs = [[p1, p2, 0] for p1, p2 in pred_pairs]
                else:
                    raise ValueError(f"Unknown structure predictor: {self.args.structure_predictor}")
                # Convert the predicted pairs to a dictionary, if necessary.
                pred_pairs = convert_pred_pairs(pred_pairs)
                logging.info("%s predicted %d unique base pairs.",
                             self.args.structure_predictor,
                             len({(min(i, j), max(i, j)) for i, j in pred_pairs.items() if i < j}))
                return pred_pairs
        elif self.args.structure:
            logging.info("Using provided dot-bracket structure: %s", self.args.structure)
            return self.pair_indices(self.args.structure)
        else:
            raise ValueError("Either a structure predictor or a structure must be provided.")

    def mutate_pair(self, nt1: str, nt2: str) -> str:
        original = nt1 + nt2
        candidates = PAIR_MUTATIONS.get(original, WC_PAIRS)
        chosen = random.choice(candidates)
        if random.random() < self.args.wobble_prob:
            return random.choice(['GU', 'UG'])
        return chosen

    def mutate_unpaired(self, nt: str) -> str:
        return random.choice(['A', 'U', 'G', 'C'])

    def mutate_sequence(self, seq: str, pairs: Dict[int, int]) -> str:
        mutated = list(seq)
        i = 0
        while i < len(seq):
            if i not in pairs:
                # Unpaired (loop) region: allow insertions and deletions.
                insertion = ''
                if random.random() < self.args.long_insertion_prob:
                    insertion_len = random.randint(2, self.max_insertion_length) if self.max_insertion_length > 2 else random.randint(2, 5)  # 5 chosen at random
                    insertion = ''.join(random.choice('augc') for _ in range(insertion_len))
                elif random.random() < self.args.insertion_prob_loop:
                    insertion = random.choice('augc')
                if random.random() < self.args.long_deletion_prob:
                    del_len = random.randint(2, self.max_deletion_length) if self.max_deletion_length > 2 else random.randint(2, 5)  # 5 chosen at random
                    for j in range(i, min(i + del_len, len(seq))):
                        mutated[j] = '-'
                    i += del_len
                    continue
                elif random.random() < self.args.deletion_prob_loop:
                    mutated[i] = '-'
                elif mutated[i] != '-':
                    mutated[i] = self.mutate_unpaired(seq[i])
                mutated[i] = insertion + mutated[i]
            elif i < pairs[i]:
                # Paired (stem) region: use stem-keep probability.
                j = pairs[i]
                if random.random() < self.args.stem_keep_prob:
                    mutated[i] = seq[i]
                    mutated[j] = seq[j]
                else:
                    new_pair = self.mutate_pair(seq[i], seq[j])
                    mutated[i] = new_pair[0]
                    mutated[j] = new_pair[1]
            i += 1
        return ''.join(mutated)

    def generate_msa(self, rna_seq: str, pairs: Dict[int, int]) -> List[str]:
        # Set maximum insertion/deletion lengths.
        self.max_insertion_length = int(len(rna_seq) * self.args.max_insertion_fraction)
        self.max_deletion_length = int(len(rna_seq) * self.args.max_deletion_fraction)
        msa = [rna_seq]
        for _ in range(self.args.N - 1):
            mutated = self.mutate_sequence(rna_seq, pairs)
            msa.append(mutated)
        return msa

    def plot_final_features(self, msa: List[str], rna_seq: str, pairs: Dict[int, int]) -> None:
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
        if self.args.show_plot:
            plt.show()
        final_dir = Path("msa_final_plots")
        final_dir.mkdir(parents=True, exist_ok=True)
        final_filename = final_dir / f"msa_final_{self.args.pdb_id}.png"
        plt.savefig(final_filename)
        plt.close()
        logging.info("Final MSA features plot saved to: %s", final_filename)

    def process(self) -> None:
        json_output: Dict[str, Any] = {}
        if self.args.input_json_path:
            data = load_json(self.args.input_json_path)
            if self.args.max_chains is not None and "sequences" in data:
                if len(data["sequences"]) > self.args.max_chains:
                    logging.warning("Input JSON contains %d chains, but --max_chains=%d is set. Skipping %s...",
                                    len(data["sequences"]), self.args.max_chains,
                                    Path(self.args.input_json_path).name)
                    return
            updated_sequences = []
            for key, value in data.items():
                if key == "sequences":
                    for chain in value:
                        if "rna" in chain:
                            updated_chain = chain.copy()
                            rna_seq = chain["rna"]["sequence"]
                            logging.info("Processing RNA sequence from JSON: %s", rna_seq)
                            pairs = self.get_structure(rna_seq)
                            logging.info("Predicted pairs: %s", pairs)
                            unique_pairs = {(min(i, j), max(i, j)) for i, j in pairs.items() if i < j}
                            logging.info("Secondary structure has %d unique base pairs.", len(unique_pairs))
                            msa = self.generate_msa(rna_seq, pairs)
                            updated_chain["rna"]["unpairedMsa"] = "\n".join(
                                [">query\n" + msa[0]] +
                                [f">sample_{i}\n{seq}" for i, seq in enumerate(msa[1:])]
                            )
                            if self.args.plot:
                                self.plot_final_features(msa, rna_seq, pairs)
                            updated_sequences.append(updated_chain)
                        if "protein" in chain:
                            updated_sequences.append(chain)
                else:
                    json_output[key] = value
            json_output["sequences"] = updated_sequences
            ins_len = getattr(self, "max_insertion_length", "NA")
            del_len = getattr(self, "max_deletion_length", "NA")
            name = f"{self.args.pdb_id}_custom_rnamsa_N{self.args.N}_seed{self.args.seed}_insl_{self.args.insertion_prob_loop}_dell_{self.args.deletion_prob_loop}_inss_{self.args.insertion_prob_stem}_dels_{self.args.deletion_prob_stem}_lins_{self.args.long_insertion_prob}_ldels_{self.args.long_deletion_prob}_maxinslen_{ins_len}_maxdellen_{del_len}_wobble_{self.args.wobble_prob}_stemkeep_{self.args.stem_keep_prob}_{self.args.structure_predictor}"
            json_output["name"] = name
        else:
            if not (self.args.rna_seq and self.args.protein_seq):
                raise ValueError("Either --rna-seq and --protein-seq or --input_json_path must be provided.")
            rna_seq = self.args.rna_seq.upper()
            logging.info("Input RNA sequence: %s", rna_seq)
            pairs = self.get_structure(rna_seq)
            unique_pairs = {(min(i, j), max(i, j)) for i, j in pairs.items() if i < j}
            logging.info("Secondary structure has %d unique base pairs.", len(unique_pairs))
            msa = self.generate_msa(rna_seq, pairs)
            if self.args.print_msa:
                logging.info("MSA:")
                for row in msa:
                    logging.info(row)
            if self.args.plot:
                self.plot_final_features(msa, rna_seq, pairs)
            msa_fasta = [">query\n" + msa[0]]
            msa_fasta += [f">sample_{i}\n{seq}" for i, seq in enumerate(msa[1:])]
            unpaired_msa = "\n".join(msa_fasta)
            name = f"{self.args.pdb_id}_custom_rnamsa_N{self.args.N}_seed{self.args.seed}_insl_{self.args.insertion_prob_loop}_dell_{self.args.deletion_prob_loop}_inss_{self.args.insertion_prob_stem}_dels_{self.args.deletion_prob_stem}_lins_{self.args.long_insertion_prob}_ldels_{self.args.long_deletion_prob}_maxinslen_{self.max_insertion_length}_maxdellen_{self.max_deletion_length}_wobble_{self.args.wobble_prob}_{self.args.structure_predictor}"
            json_output = {
                "name": name,
                "modelSeeds": [1],
                "sequences": [
                    {
                        "rna": {
                            "sequence": rna_seq,
                            "modifications": [],
                            "unpairedMsa": unpaired_msa,
                            "id": "A"
                        }
                    },
                    {
                        "protein": {
                            "sequence": self.args.protein_seq,
                            "modifications": [],
                            "id": "B"
                        }
                    }
                ],
                "dialect": "alphafold3",
                "version": 1
            }
        output_path = Path(self.args.output_json_dir, f"{name}.json")
        with open(output_path, "w") as f:
            json.dump(json_output, f, indent=2)
        logging.info("✅ JSON written to: %s", output_path)


def main() -> None:
    args = parse_args()
    generator = MsaGenerator(args)
    try:
        generator.process()
    except Exception as e:
        logging.error("Error during processing: %s", e)


if __name__ == '__main__':
    main()
