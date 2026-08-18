#!/usr/bin/env python
import argparse
import random
import json
import logging
import sys
from pathlib import Path
import numpy as np

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

HERE_DIR = Path(__file__).resolve().parent
if str(HERE_DIR) not in sys.path:
    sys.path.insert(0, str(HERE_DIR))

import json_generator
from pair_map import PairMap

from typing import Dict, List, Any, Optional

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

PAIR_MUTATION_PROBABILITIES = {
    "A": {"U": 1.0},
    "C": {"G": 1.0},
    "G": {"C": 0.75, "U": 0.25},
    "U": {"A": 0.75, "G": 0.25},
}

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

def load_json(json_path: str) -> Any:
    try:
        with open(json_path, "r") as f:
            return json.load(f)
    except Exception as e:
        logging.error("Failed to load JSON from %s: %s", json_path, e)
        raise

def parse_json(json_string: str) -> Any:
    try:
        json_string = json_string.replace("(", "[")
        json_string = json_string.replace(")", "]")
        return json.loads(json_string)
    except Exception as e:
        logging.error("Failed to parse JSON from string: %s", e)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate RNA MSA and AF3-compatible JSON.")
    io_group = parser.add_argument_group("I/O and Structure")
    io_group.add_argument('--structure_predictor', type=str, default=None,
                          help="Predictor for secondary structure (base pairs output). Options: rnafold, spotrna, rnaformer, dssr")
    io_group.add_argument('--structure', type=str, required=False, default=None,
                          help="Dot-bracket secondary structure (if provided, will be converted to base pairs) or base pairs")
    io_group.add_argument('--interactions', type=str, required=False, default=None,
                              help="List of tuples with interaction rates for every pair. The default is 0 for missing tuples")
    io_group.add_argument('--rna-seq', type=str, required=False,
                          help="RNA sequence (used for MSA + JSON)")
    io_group.add_argument('--pair-mutation', type=str, required=False)
    io_group.add_argument('--protein-seq', type=str, required=False,
                          help="Protein sequence (for JSON)")
    io_group.add_argument('--input_json_path', type=str, required=False,
                          help="A JSON file in Alphafold webservice format that will be adapted for custom RNA MSA")
    io_group.add_argument('--pdb_id', type=str, required=False,
                          help="PDB ID (used in JSON name)")
    io_group.add_argument('--output_json_dir', type=str, default="custom_msa_json_output")

    mut_group = parser.add_argument_group("Mutation Parameters")
    mut_group.add_argument('-N', type=int, default=20, help="Number of sequences in the MSA")
    mut_group.add_argument('--mutation-rate-unpaired', type=float, default=0.2)
    mut_group.add_argument('--mutation-rate-paired', type=float, default=0.2)
    mut_group.add_argument('--mutation-rates', type=str, required=False, 
                        help="Input a list with a mutation rate for every position in the sequence. Overrides mutation-rate-unpaired "
                        "and mutation-rate-paired if both are provided")
    mut_group.add_argument('--pair-mutation-approach', type=str, default="watson_crick_cov",
                           choices=["watson_crick", "covariance", "watson_crick_cov", "original", "none", "potts"],
                           help="Choose the method for the mutation of base pairs. 'original' corresponds "
                                "to the approach we used before the rework has problems. 'watson_crick' chooses"
                                " a random watson crick base pair with hardcoded probabilities. 'covariance' "
                                "chooses random base pairs but ensures both partners are always changing together. "
                                "'potts' samples the whole substitution layer jointly from a Potts model whose "
                                "fields come from the mutation rates and whose couplings come from the "
                                "interaction map; indels are unaffected")

    potts_group = parser.add_argument_group("Potts Parameters (--pair-mutation-approach potts)")
    potts_group.add_argument('--potts-coupling', type=float, default=6.0,
                             help="Coupling strength gamma. 0 removes all covariation, large values force "
                                  "every substituted pair to stay canonical")
    potts_group.add_argument('--potts-wobble', type=float, default=0.85,
                             help="Pairing score of GU/UG relative to a canonical pair. Controls the wobble "
                                  "fraction among substituted pairs")
    potts_group.add_argument('--potts-no-rate-matching', action='store_true',
                             help="Skip solving the fields for the requested mutation rates. Coupling suppresses "
                                  "substitution, so realised rates then fall below the requested ones")
    potts_group.add_argument('--potts-max-enumerate', type=int, default=7,
                             help="Largest coupled component tabulated exactly. Larger components fall back to "
                                  "Gibbs sampling")
    potts_group.add_argument('--potts-max-rate-match', type=int, default=4,
                             help="Largest coupled component for which the fields are solved to hit the requested "
                                  "mutation rates. Larger components use raw fields and undershoot")
    mut_group.add_argument('--stem_single_insertion_prob', type=float, default=0.05)
    mut_group.add_argument('--stem_long_insertion_prob', type=float, default=0.01)
    mut_group.add_argument('--stem_single_deletion_prob', type=float, default=0.01)
    mut_group.add_argument('--stem_pair_deletion_prob', type=float, default=0.01)
    mut_group.add_argument('--loop_single_insertion_prob', type=float, default=0.1)
    mut_group.add_argument('--loop_single_deletion_prob', type=float, default=0.1)
    mut_group.add_argument('--loop_long_insertion_prob', type=float, default=0.02)
    mut_group.add_argument('--loop_long_deletion_prob', type=float, default=0.02,
                           help="Note that long deletions override mutations and single deletions if they appear.")
    mut_group.add_argument('--wobble-prob', type=float, default=0.1)
    mut_group.add_argument('--max-insertion-fraction', type=float, default=0.1,
                           help="Max insertion length as fraction of RNA length")
    mut_group.add_argument('--max-deletion-fraction', type=float, default=0.1,
                           help="Max deletion length as fraction of RNA length")
    
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--max_chains', type=int, default=None)
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--print_msa', action='store_true')
    parser.add_argument('--show_plot', action='store_true')
    return parser.parse_args()


class MsaGenerator:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.pair_map: Optional[PairMap] = None
        self.potts = None
        if self.args.seed is not None:
            random.seed(self.args.seed)

    def get_structure(self, sequence: str) -> PairMap:
        seq = sequence.upper()

        if self.args.interactions:
            if self.args.structure:
                logging.warning("Both structure and interactions provided. Using provided interactions.")
            if self.args.structure_predictor:
                logging.warning("Both structure predictor and interactions provided. Using provided interactions.")
            logging.info("Using provided interactions: %s", self.args.interactions)
            interactions = parse_json(self.args.interactions)
            if self.args.mutation_rates:
                self.args.mutation_rate_paired = ""
                self.args.mutation_rate_unpaired = ""
                interactions = [(int(i), int(j), float(val)) for i, j, val in interactions]
                mutation_rates = parse_json(self.args.mutation_rates)
                self.args.mutation_rates = "custom"
                return PairMap.from_interactions(len(seq), interactions, mutation_rates)
            return PairMap.from_raw(interactions, len(seq), self.args.mutation_rate_paired, self.args.mutation_rate_unpaired)

        if self.args.structure:
            if self.args.structure_predictor:
                logging.warning("Both structure predictor and structure provided. Using provided structure.")
            logging.info("Using provided structure: %s", self.args.structure)
            return PairMap.from_raw(self.args.structure, len(seq), self.args.mutation_rate_paired, self.args.mutation_rate_unpaired)

        if self.args.structure_predictor:
            import structure_predictor  # Imported lazily so the no-predict path never pays the RnaBench cost.
            pdb_id = self.args.pdb_id.lower()[:4] if self.args.pdb_id else None
            pred_pairs = structure_predictor.predict(self.args.structure_predictor, seq, pdb_id)
            logging.info("%s predicted structure.", self.args.structure_predictor)
            return PairMap.from_raw(pred_pairs, len(seq), self.args.mutation_rate_paired, self.args.mutation_rate_unpaired)

        raise ValueError("Either a structure predictor a structure or interactions must be provided.")

    def loop_insertion(self) -> str:
        """Long insertions take priority, then single insertions. If all long insertions are disregarded
        the single insertion rate can be recovered accurately. Any insertion is randomly selected"""
        if random.random() < self.args.loop_long_insertion_prob:
            insertion_len = random.randint(2, self.max_insertion_length)
            return ''.join(random.choice('augc') for _ in range(insertion_len))
        if random.random() < self.args.loop_single_insertion_prob:
            return random.choice('acgu')
        return ""

    def stem_insertion(self) -> str:
        """Long insertions take priority, then single insertions. If all long insertions are disregarded
        the single insertion rate can be recovered accurately. Any insertion is randomly selected"""
        if random.random() < self.args.stem_long_insertion_prob:
            insertion_len = random.randint(2, self.max_insertion_length)
            return ''.join(random.choice('augc') for _ in range(insertion_len))
        if random.random() < self.args.stem_single_insertion_prob:
            return random.choice('acgu')
        return ""

    def paired_deletion(self, i, mutated_partner: str, partner_index: float) -> bool:
        single_del = self.args.stem_single_deletion_prob
        pair_del =  self.args.stem_pair_deletion_prob * self.pair_map.interaction(i, partner_index)
        if self.pair_map.is_multiplet_member(i):
            return random.random() < single_del
        if i < partner_index: # P(-)
            return random.random() < single_del + pair_del
        if mutated_partner == "-": # P(-|-)
            return random.random() * (single_del + pair_del) < pair_del
        return random.random() * (1 - single_del - pair_del) < single_del # P(-|!-)

    def mutate_random(self, nt: str, prob: float):
        if random.random() < prob:
            return random.choice([c for c in 'AUGC' if c != nt])
        return nt
    
    def mutate_unpaired(self, nt: str, loop_long_del_len: int, mutation_rate: float,
                        substitute: Optional[str] = None) -> tuple[str, int]:
        """Long deletions take priority, then single deletions, then mutations. Therefore
        the mutation rate is can be recovered accurately if deletions are disregarded.
        A substitute, if given, replaces the per-position draw but leaves deletions untouched."""
        if random.random() < self.args.loop_long_deletion_prob:
            loop_long_del_len = random.randint(2, self.max_deletion_length)
        if loop_long_del_len > 0:
            return "-", loop_long_del_len - 1
        if random.random() < self.args.loop_single_deletion_prob:
            return "-", 0
        if substitute is not None:
            return substitute, 0
        return self.mutate_random(nt, mutation_rate), 0

    def mutate_cov(self, nt: str, partners_original: np.ndarray, partners_mutated: np.ndarray,
                    interactions: np.ndarray, partner_mutation_rates: np.ndarray, mutation_rate: float) -> str:
        """Guarantees that all partners get mutated together but randomly and therefore only increases covariance."""
        opts = list(zip((partners_original != partners_mutated), interactions, partner_mutation_rates))
        m = mutation_rate
        prob = m
        for mut, inter, pm in opts:
            if mut:
                prob += ((m * pm) + inter * m * (1 - pm)) / max(pm, 0.0001)
            else:
                prob += ((m * (1 - pm)) - inter * m * (1 - pm)) / (1 - min(pm, 0.9999))
        return self.mutate_random(nt, prob / (len(opts)+1))

    def mutate_wc(self, nt: str, partners_original: np.ndarray, partners_mutated: np.ndarray, 
                  interactions: np.ndarray, partner_mutation_rates: np.ndarray, mutation_rate: float, increase_cov: bool) -> str:
        """Try to maximize the number of watson crick base pairs while also trying to leave no partner.
        This uses the probabilities in PAIR_MUTATION_PROBABILITIES."""
        # first in the multiplet is random with mutation rate
        if len(partners_original) == 0:
            return self.mutate_random(nt, mutation_rate)
        # use mutate_cov to decide if it should mutate for high covariance or decide randomly if nt should mutate for low covariance
        if increase_cov:
            mutate = self.mutate_cov(nt, partners_original, partners_mutated, interactions, partner_mutation_rates, mutation_rate) != nt
        else:
            mutate = self.mutate_random(nt, mutation_rate) != nt
        if not mutate:
            return nt
        # chose with interactions, always mutate
        options = {"A": 0.00001, "U": 0.00001, "G": 0.00001, "C": 0.00001}
        for j, mut in enumerate(partners_mutated):
            for opt, prob in PAIR_MUTATION_PROBABILITIES.get(mut, {}).items():
                options[opt] += prob * interactions[j]
        options.pop(nt)
        return random.choices(list(options.keys()), list(options.values()))[0]
        

    def mutate_pair_original(self, p_nt: str, nt: str, mutation_rate: float, partner_mutation_rate: float) -> str:
        if not random.random() < (mutation_rate + partner_mutation_rate) / 2:
            return p_nt + nt
        if random.random() < self.args.wobble_prob:
            return random.choice(['GU', 'UG'])
        candidates = PAIR_MUTATIONS.get(p_nt + nt, WC_PAIRS)
        return random.choice(candidates)

    def mutate_sequence(self, seq: str) -> str:
        approach = self.args.pair_mutation_approach
        if approach not in ["covariance", "watson_crick", "watson_crick_cov", "original", "none", "potts"]:
            logging.error("Unknown input for --pair-mutation-approach, '%s', please use one of the provided options", self.args.pair_mutation_approach)
            raise ValueError()
        seq = np.array(list(seq))
        # The Potts model couples positions, so its whole substitution layer is drawn up front;
        # the loop below then applies the unchanged insertion and deletion layers on top.
        potts_draw = self.potts.sample() if approach == "potts" else None
        new_seq = np.empty(len(seq), str)
        insertions = []
        loop_long_del_len = 0
        for i, nt in enumerate(seq):
            new_nt = nt
            new_insertion = ""
            mutation_rate = self.pair_map.mutation_rate(i)

            if self.pair_map.is_unpaired(i):
                new_nt, loop_long_del_len = self.mutate_unpaired(
                    nt, loop_long_del_len, mutation_rate,
                    substitute=potts_draw[i] if potts_draw is not None else None)

            if self.pair_map.is_paired(i):
                loop_long_del_len = 0
                partners =  np.array(self.pair_map.partners(i), int)
                prev = partners[partners < i]
                interactions = np.array([self.pair_map.interaction(i, j) for j in prev])
                partner_mutation_rates = np.array([self.pair_map.mutation_rate(j) for j in prev])
                if approach == "none":
                    new_nt = self.mutate_random(nt, self.pair_map.mutation_rate(i))
                if approach == "covariance":
                    new_nt = self.mutate_cov(nt, seq[prev], new_seq[prev], interactions, partner_mutation_rates,mutation_rate)
                if approach == "watson_crick":
                    new_nt = self.mutate_wc(nt, seq[prev], new_seq[prev], interactions, partner_mutation_rates, mutation_rate,  False)
                if approach == "watson_crick_cov":
                    new_nt = self.mutate_wc(nt, seq[prev], new_seq[prev], interactions, partner_mutation_rates, mutation_rate, True)
                if approach == "original" :
                    if self.pair_map.is_basic_pair(i) and i > partners[0]:
                        new_seq[partners[0]], new_nt = self.mutate_pair_original(seq[partners[0]], nt, mutation_rate, partner_mutation_rates[0])
                    if self.pair_map.is_multiplet_member(i):
                        new_nt = self.mutate_wc(nt, seq[prev], new_seq[prev], interactions, partner_mutation_rates, mutation_rate, True)
                if approach == "potts":
                    new_nt = potts_draw[i]
                # optionally override mutation if a deletion happens
                if self.paired_deletion(i, new_seq[partners[0]], partners[0]):
                    new_nt = "-"

            if self.pair_map.is_paired(i-1) and self.pair_map.is_paired(i):
                new_insertion = self.stem_insertion()
            else:
                new_insertion = self.loop_insertion()

            new_seq[i] = new_nt
            insertions.append(new_insertion)
        return ''.join(np.ravel(list(zip(insertions, new_seq)))) + self.loop_insertion()

    def generate_msa(self, rna_seq: str) -> List[str]:
        self.max_insertion_length = max(int(len(rna_seq) * self.args.max_insertion_fraction), 2)
        self.max_deletion_length = max(int(len(rna_seq) * self.args.max_deletion_fraction), 2)
        if self.args.pair_mutation_approach == "potts":
            # Imported lazily so the other approaches never pay for it.
            from potts_shs.potts import PottsModel
            self.potts = PottsModel(
                rna_seq, self.pair_map,
                coupling=self.args.potts_coupling,
                wobble=self.args.potts_wobble,
                match_rates=not self.args.potts_no_rate_matching,
                max_enumerate=self.args.potts_max_enumerate,
                max_rate_match=self.args.potts_max_rate_match,
            )
        msa = [rna_seq]
        for _ in range(self.args.N - 1):
            msa.append(self.mutate_sequence(rna_seq))
        return msa

    def build_output_name(self) -> str:
        ins_len = getattr(self, "max_insertion_length", "NA")
        del_len = getattr(self, "max_deletion_length", "NA")
        structure_predictor = self.args.structure_predictor or "none"
        name_parts = [
            f"{self.args.pdb_id}_custom_rnamsa",
            f"N{self.args.N}",
            f"seed{self.args.seed}",
            f"mru_{self.args.mutation_rate_unpaired}",
            f"mrp_{self.args.mutation_rate_paired}",
            f"mrs_{self.args.mutation_rates}",
            f"pma_{self.args.pair_mutation_approach}",
            f"ssi_{self.args.stem_single_insertion_prob}",
            f"sli_{self.args.stem_long_insertion_prob}",
            f"spd_{self.args.stem_pair_deletion_prob}",
            f"lsi_{self.args.loop_single_insertion_prob}",
            f"lsd_{self.args.loop_single_deletion_prob}",
            f"lli_{self.args.loop_long_insertion_prob}",
            f"lld_{self.args.loop_long_deletion_prob}",
            f"mif_{self.args.max_insertion_fraction}",
            f"mdf_{self.args.max_deletion_fraction}",
            f"maxinslen_{ins_len}",
            f"maxdellen_{del_len}",
            f"wp_{self.args.wobble_prob}",
        ]
        if self.args.pair_mutation_approach == "potts":
            name_parts += [
                f"pg_{self.args.potts_coupling}",
                f"pw_{self.args.potts_wobble}",
                f"prm_{not self.args.potts_no_rate_matching}",
            ]
        name_parts.append(structure_predictor)
        return "_".join(str(part) for part in name_parts if part not in {None, ""})

    def _build_rna_chain(self, chain: Dict[str, Any]) -> Dict[str, Any]:
        """Return a copy of an RNA chain with a freshly generated custom MSA
        written into its unpairedMsa field."""
        updated_chain = chain.copy()
        rna_seq = chain["rna"]["sequence"]
        logging.info("Processing RNA sequence: %s", rna_seq)

        self.pair_map = self.get_structure(rna_seq)

        logging.info("Pairs: %s", self.pair_map.pairs)
        logging.info("Secondary structure has %d unique base pairs.", len(self.pair_map.unique_pairs))
        if self.pair_map.basic_pairs:
            logging.info("Using %d stem pairs: %s", len(self.pair_map.basic_pairs), self.pair_map.basic_pairs)
        if self.pair_map.multiplets:
            logging.info("Using %d multiplets: %s", len(self.pair_map.multiplets), self.pair_map.multiplets)

        msa = self.generate_msa(rna_seq)
        updated_chain["rna"]["unpairedMsa"] = "\n".join(
            [">query\n" + msa[0]] +
            [f">sample_{i}\n{seq}" for i, seq in enumerate(msa[1:])]
        )

        if self.args.print_msa:
            logging.info("MSA:")
            for row in msa:
                logging.info(row)

        if self.args.plot:
            # Imported lazily so the no-plot path never pays the matplotlib cost.
            import msa_plotting
            msa_plotting.plot_final_features(msa, rna_seq, self.pair_map.pairs,
                                             self.args.pdb_id, self.args.show_plot)
        return updated_chain

    def process(self, write: bool = True) -> Dict[str, Any]:
        if self.args.input_json_path:
            data = load_json(self.args.input_json_path)
        elif self.args.rna_seq and self.args.protein_seq:
            data = json_generator.build_input_json(self.args.rna_seq, self.args.protein_seq)
        else:
            raise ValueError("Either --rna-seq and --protein-seq or --input_json_path must be provided.")

        chains = data.get("sequences", [])
        if self.args.max_chains is not None and len(chains) > self.args.max_chains:
            source = Path(self.args.input_json_path).name if self.args.input_json_path else "generated input"
            logging.warning("Input JSON contains %d chains, but --max_chains=%d is set. Skipping %s...",
                            len(chains), self.args.max_chains, source)
            return

        updated_sequences = []
        for chain in chains:
            if "rna" in chain:
                updated_sequences.append(self._build_rna_chain(chain))
            if "protein" in chain:
                updated_sequences.append(chain)

        # Preserve any extra top-level keys, then set the processed sequences/name.
        json_output = {k: v for k, v in data.items() if k != "sequences"}
        json_output["sequences"] = updated_sequences
        name = self.build_output_name()
        json_output["name"] = name

        if write:
            output_path = Path(self.args.output_json_dir, f"{name}.json")
            with open(output_path, "w") as f:
                json.dump(json_output, f, indent=2)
            logging.info("✅ JSON written to: %s", output_path)
        return json_output


def main() -> None:
    args = parse_args()
    generator = MsaGenerator(args)
    try:
        generator.process()
    except Exception as e:
        logging.error("Error during processing: %s", e)


if __name__ == '__main__':
    main()
