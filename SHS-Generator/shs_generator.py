#!/usr/bin/env python
import argparse
import random
import json
import logging
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

HERE_DIR = Path(__file__).resolve().parent
if str(HERE_DIR) not in sys.path:
    sys.path.insert(0, str(HERE_DIR))

import json_generator

from collections import defaultdict, deque
from typing import Dict, List, Any, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

PAIR_MUTATION_PROBABILITIES = {
    "AU" : 0.2,
    "UA" : 0.2,
    "GC" : 0.2,
    "CG" : 0.2,
    "GU" : 0.1,
    "UG" : 0.1
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


def _iter_pairs(raw: Any):
    """Yield (i, j) int tuples from any supported pair representation — a dict,
    or a list of [i, j, ...] — in input order. Shared by convert_pred_pairs and
    _normalize_pairs so the two can't diverge on format handling."""
    if isinstance(raw, dict):
        items = raw.items()
    elif isinstance(raw, (list, tuple)):
        items = [(item[0], item[1]) for item in raw if isinstance(item, (list, tuple)) and len(item) >= 2]
    else:
        raise ValueError(f"Unknown pair format: {type(raw)}")
    for a, b in items:
        yield int(a), int(b)


def convert_pred_pairs(pred_pairs: Any) -> Dict[int, int]:
    """Build a bidirectional {i: j, j: i} dict from predicted pairs, preserving
    input order (so the last partner wins where a position pairs more than once)."""
    result: Dict[int, int] = {}
    for i, j in _iter_pairs(pred_pairs):
        result[i] = j
        result[j] = i
    return result


def _normalize_pairs(raw: Any) -> List[Tuple[int, int]]:
    """Flatten predicted pairs into a sorted, de-duplicated list of (i, j) with
    i < j. Self-pairs are dropped."""
    out: set = set()
    for a, b in _iter_pairs(raw):
        if a == b:
            continue
        out.add((a, b) if a < b else (b, a))
    return sorted(out)


def _unique_pairs(pairs: Dict[int, int]) -> set:
    """The set of unordered base pairs (i, j) with i < j in a pairs dict."""
    return {(i, j) for i, j in pairs.items() if i < j}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate RNA MSA and AF3-compatible JSON.")
    io_group = parser.add_argument_group("I/O and Structure")
    io_group.add_argument('--structure_predictor', type=str, default=None,
                          help="Predictor for secondary structure (base pairs output). Options: rnafold, spotrna, rnaformer, dssr")
    io_group.add_argument('--structure', type=str, required=False, default=None,
                          help="Dot-bracket secondary structure (if provided, will be converted to base pairs) or base pairs")
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
    mut_group.add_argument('--pair-mutation-approach', type=str, nargs='*', default="watson_crick", 
                           choices=["watson_crick", "covariance", "original"]
                           help="Choose the method for the mutation of base pairs. 'original' corresponds "
                                "to the approach we used before the rework has problems. 'watson_crick' chooses"
                                " a random watson crick base pair with hardcoded probabilities. 'covariance' "
                                "chooses random base pairs but ensures both partners are always changing together")
    mut_group.add_argument('--insertion-prob-loop', type=float, default=0.2)
    mut_group.add_argument('--deletion-prob-loop', type=float, default=0.8)
    mut_group.add_argument('--insertion-prob-stem', type=float, default=0.01)
    mut_group.add_argument('--deletion-prob-stem', type=float, default=0.01)
    mut_group.add_argument('--long-insertion-prob', type=float, default=0.05)
    mut_group.add_argument('--long-deletion-prob', type=float, default=0.05)
    mut_group.add_argument('--wobble-prob', type=float, default=0.1)
    mut_group.add_argument('--max-insertion-fraction', type=float, default=0.1,
                           help="Max insertion length as fraction of RNA length")
    mut_group.add_argument('--max-deletion-fraction', type=float, default=0.1,
                           help="Max deletion length as fraction of RNA length")
    # Triplet co-mutation parameters (extension of the pairwise mutation
    # scheme to three correlated positions). Default --triplet-prob=0.0
    # disables triplets, preserving the original pair-only behavior.
    mut_group.add_argument('--triplet-prob', type=float, default=0.0,
                           help="Probability that a base triple/multiplet (a connected component of "
                                ">=3 positions in the predicted base-pair graph, e.g. pairs (1,15) and "
                                "(15,18) -> {1,15,18}) is promoted and co-mutated. Set to 0 to disable.")
    mut_group.add_argument('--triplet-keep-prob', type=float, default=0.99,
                           help="Probability of keeping a multiplet's residues unchanged (analog of --stem-keep-prob).")
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--max_chains', type=int, default=None)
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--print_msa', action='store_true')
    parser.add_argument('--show_plot', action='store_true')
    return parser.parse_args()

class MsaGenerator:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.pairs: dict[int, int] = {}
        if self.args.seed is not None:
            random.seed(self.args.seed)

    def pair_indices(self, dotbracket: str) -> Dict[int, int]:
        # Support nested pseudoknot bracket families ()[]{}<> via one stack
        # per family; any other character is treated as unpaired. Behavior
        # for plain ()-only structures is unchanged.
        openers = {'(': ')', '[': ']', '{': '}', '<': '>'}
        close_to_open = {c: o for o, c in openers.items()}
        stacks: Dict[str, list] = {o: [] for o in openers}
        pairs: Dict[int, int] = {}
        for i, char in enumerate(dotbracket):
            if char in openers:
                stacks[char].append(i)
            elif char in close_to_open:
                o = close_to_open[char]
                if not stacks[o]:
                    logging.error("Unbalanced structure at position %d", i)
                    continue
                j = stacks[o].pop()
                pairs[i] = j
                pairs[j] = i
        return pairs

    def get_structure(self, sequence: str) -> Tuple[Dict[int, int], List[Dict[str, Any]]]:
        seq = sequence.upper()
        # A provided dot-bracket structure takes precedence over a predictor.
        if self.args.structure:
            if self.args.structure_predictor:
                logging.warning("Both structure predictor and structure provided. Using provided structure.")
            try:
                data = json.loads(self.args.structure)
                if not isinstance(data, dict):
                    raise json.decoder.JSONDecodeError("not an object", data, 0)
                logging.info("Using provided pairs as structure: %s", self.args.structure)
                pairs = {int(k): int(v) for k, v in data.items()}
            except json.decoder.JSONDecodeError:
                logging.info("Using provided dot-bracket structure: %s", self.args.structure)
                pairs = self.pair_indices(self.args.structure)
            return pairs, self.derive_multiplets(_normalize_pairs(pairs))
        if self.args.structure_predictor:
            import structure_predictor  # Imported lazily so the no-predict path never pays the RnaBench cost.
            pdb_id = self.args.pdb_id.lower()[:4] if self.args.pdb_id else None
            pred_pairs = structure_predictor.predict(self.args.structure_predictor, seq, pdb_id)
            # Preserve the raw pair list before the dict conversion: a position
            # may participate in more than one pair (base triple/multiplet),
            # which a dict would silently collapse.
            pair_list = _normalize_pairs(pred_pairs)
            pred_pairs = convert_pred_pairs(pred_pairs)
            logging.info("%s predicted %d unique base pairs.",
                         self.args.structure_predictor,
                         len(_unique_pairs(pred_pairs)))
            return pred_pairs, self.derive_multiplets(pair_list)
        raise ValueError("Either a structure predictor or a structure must be provided.")

    def derive_multiplets(self, pair_list: List[Tuple[int, int]]) -> List[Dict[str, Any]]:
        """Derive base triples / higher multiplets from the predicted pairs.

        The pairs form a graph (positions = nodes, pairs = edges): any connected
        component of size >= 3 is a multiplet (e.g. (1,15) and (15,18) share
        position 15 -> {1,15,18}); size-2 components stay ordinary pairs handled
        by the existing pair machinery. Each multiplet is kept with probability
        --triplet-prob, once per run and deterministic given --seed."""
        if self.args.triplet_prob <= 0.0 or not pair_list:
            return []

        adj: Dict[int, set] = defaultdict(set)
        edges: set = set()
        for a, b in pair_list:
            adj[a].add(b)
            adj[b].add(a)
            edges.add((a, b))

        seen: set = set()
        multiplets: List[Dict[str, Any]] = []
        for start in sorted(adj):
            if start in seen:
                continue
            # BFS to collect the connected component.
            comp: set = set()
            queue = deque([start])
            seen.add(start)
            while queue:
                u = queue.popleft()
                comp.add(u)
                for v in adj[u]:
                    if v not in seen:
                        seen.add(v)
                        queue.append(v)
            if len(comp) < 3:
                continue  # ordinary base pair -> handled via the pairs dict
            if random.random() >= self.args.triplet_prob:
                continue  # not promoted this run; falls back to pair handling
            multiplets.append({
                "positions": tuple(sorted(comp)),
                "edges": sorted(e for e in edges if e[0] in comp and e[1] in comp)
            })

        sizes = [len(m["positions"]) for m in multiplets]
        logging.info(
            "Derived %d multiplets from %d base pairs (triplet_prob=%.3f); sizes=%s",
            len(multiplets), len(edges), self.args.triplet_prob, sizes)
        return multiplets

    def _partner_for(self, nt: str) -> str:
        """Pick a base that pairs with `nt`: Watson-Crick by default, with a
        --wobble-prob chance of a G-U wobble where chemically valid."""
        if nt == 'G':
            return 'U' if random.random() < self.args.wobble_prob else 'C'
        if nt == 'U':
            return 'G' if random.random() < self.args.wobble_prob else 'A'
        if nt == 'A':
            return 'U'
        if nt == 'C':
            return 'G'
        return random.choice(['A', 'U', 'G', 'C'])

    def mutate_multiplet(self, seq: str, positions: Tuple[int, ...],
                         edges: List[Tuple[int, int]]) -> Dict[int, str]:
        """Co-mutate the residues of a base triple / multiplet so that every
        structural contact (graph edge) stays a valid pair.

        Every edge is a real base pair, so we walk the component's edges in
        sorted order keeping a position->base assignment: the first edge of a
        sub-tree is mutated as a free pair; once one endpoint is fixed (e.g. a
        hub shared between two pairs), the other is chosen as a valid partner of
        it. This keeps every edge valid and degrades to `mutate_pair` for a lone
        pair."""
        assigned: Dict[int, str] = {}
        for a, b in sorted(edges):
            a_set, b_set = a in assigned, b in assigned
            if not a_set and not b_set:
                assigned[a], assigned[b] = self.mutate_pair()
            elif a_set and not b_set:
                assigned[b] = self._partner_for(assigned[a])
            elif b_set and not a_set:
                assigned[a] = self._partner_for(assigned[b])
            # else: both already fixed (a cycle's closing edge) -> leave as is.
        for p in positions:
            assigned.setdefault(p, seq[p])
        return assigned

    def _mutate_anchor(self, mutated: List[str], seq: str, m: Dict[str, Any]) -> None:
        # Multiplet co-mutation: parallel to the pair branch but acts on the
        # whole connected component via --triplet-keep-prob.
        if random.random() < self.args.triplet_keep_prob:
            for p in m["positions"]:
                mutated[p] = seq[p]
        else:
            for p, base in self.mutate_multiplet(seq, m["positions"], m["edges"]).items():
                mutated[p] = base

    def _mutate_unpaired_position(self, mutated: List[str], seq: str, i: int,
                                  multiplet_members: set) -> int:
        # Unpaired (loop) region: allow insertions and deletions. Returns how
        # far to advance i (a long deletion skips several positions at once).
        insertion = ''
        if random.random() < self.args.long_insertion_prob:
            insertion_len = random.randint(2, self.max_insertion_length) if self.max_insertion_length > 2 else random.randint(2, 5)  # 5 chosen at random
            insertion = ''.join(random.choice('augc') for _ in range(insertion_len))
        elif random.random() < self.args.insertion_prob_loop:
            insertion = random.choice('augc')
        elif random.random() < self.args.long_deletion_prob:
            del_len = random.randint(2, self.max_deletion_length) if self.max_deletion_length > 2 else random.randint(2, 5)  # 5 chosen at random
            for j in range(i, min(i + del_len, len(seq))):
                if j not in multiplet_members:  # keep multiplet geometry intact
                    mutated[j] = '-'
            return del_len
        elif random.random() < self.args.deletion_prob_loop:
            mutated[i] = '-'
        elif random.random() < self.args.mutation_rate_unpaired:
            mutated[i] = random.choice([c for c in 'AUGC' if c != seq[i]])
        mutated[i] = insertion + mutated[i]
        return 1

    def mutate_stem_pair(self, mutated: List[str], seq: str, i: int, j: int):
        if self.args.pair_mutation_approach == "covariance":
            self._mutate_stem_pair_cov(mutated, seq, i, j)
        elif self.args.pair_mutation_approach == "watson_crick":
            self._mutate_stem_pair_wc(mutated, seq, i, j)
        elif self.args.pair_mutation_approach == "original":
            self._mutate_pair_original(mutated, seq, i, j)
        else:
            logging.error("Unknown input for --pair-mutation-approach, '%s', please use on of the provided options",self.args.pair_mutation_approach)
            raise ValueError()

    def _mutate_stem_pair_wc(self, mutated: List[str], seq: str, i: int, j: int) -> None:
        """Mutate the base pair to a random watson crick base pair based on the probabilities in PAIR_MUTATIONS.
        This also increases covariance but slightly less then the pure covariance approach"""
        options = dict(PAIR_MUTATION_PROBABILITIES)
        options.pop(seq[i] + seq[j], 0)
        if random.random() < self.args.mutation_rate_paired:
            mutated[i], mutated[j] = random.choices(list(options.keys()), list(options.values()))[0]

    def _mutate_pair_original(self, mutated: List[str], seq: str, i: int, j: int) -> None:
        if random.random() < self.args.mutation_rate_paired:
            candidates = PAIR_MUTATIONS.get(seq[i] + seq[j], WC_PAIRS)
            mutated[i], mutated[j] = random.choice(candidates)
            if random.random() < self.args.wobble_prob:
                mutated[i], mutated[j] = random.choice(['GU', 'UG'])

    def _mutate_stem_pair_cov(self, mutated: List[str], seq: str, i: int, j: int) -> None:
        """Mutate the pair in a way that only guarantees, that both partners get mutated
        and therefore only increases covariance"""
        if random.random() < self.args.mutation_rate_paired:
            mutated[i] = random.choice([c for c in 'AUGC' if c != seq[i]])
            mutated[j] = random.choice([c for c in 'AUGC' if c != seq[j]])

    def mutate_sequence(self, seq: str, pairs: Dict[int, int],
                        multiplets: Optional[List[Dict[str, Any]]] = None) -> str:
        if multiplets is None:
            multiplets = []
        # Each multiplet is anchored at its smallest position; the remaining
        # member positions are co-mutated with the anchor and skipped by the
        # per-position branches so we never insert/delete/repair inside a
        # base triple or higher multiplet.
        anchor_to_multiplet: Dict[int, Dict[str, Any]] = {
            m["positions"][0]: m for m in multiplets
        }
        multiplet_members: set = {p for m in multiplets for p in m["positions"]}
        mutated = list(seq)
        i = 0
        while i < len(seq):
            if i in anchor_to_multiplet:
                self._mutate_anchor(mutated, seq, anchor_to_multiplet[i])
            elif i in multiplet_members:
                pass  # non-anchor member: already set by its anchor
            elif i not in pairs:
                i += self._mutate_unpaired_position(mutated, seq, i, multiplet_members)
                continue
            elif i < pairs[i] and pairs[i] not in multiplet_members:
                self.mutate_stem_pair(mutated, seq, i, pairs[i])
            i += 1
        return ''.join(mutated)

    def generate_msa(self, rna_seq: str, pairs: Dict[int, int],
                     multiplets: Optional[List[Dict[str, Any]]] = None) -> List[str]:
        self.max_insertion_length = int(len(rna_seq) * self.args.max_insertion_fraction)
        self.max_deletion_length = int(len(rna_seq) * self.args.max_deletion_fraction)
        msa = [rna_seq]
        for _ in range(self.args.N - 1):
            mutated = self.mutate_sequence(rna_seq, pairs, multiplets)
            msa.append(mutated)
        return msa

    def build_output_name(self) -> str:
        ins_len = getattr(self, "max_insertion_length", "NA")
        del_len = getattr(self, "max_deletion_length", "NA")
        triplet_suffix = (
            f"_triplet_{self.args.triplet_prob}_keep_{self.args.triplet_keep_prob}"
            if self.args.triplet_prob > 0 else ""
        )
        return (
            f"{self.args.pdb_id}_custom_rnamsa_N{self.args.N}_seed{self.args.seed}"
            f"_insl_{self.args.insertion_prob_loop}_dell_{self.args.deletion_prob_loop}"
            f"_inss_{self.args.insertion_prob_stem}_dels_{self.args.deletion_prob_stem}"
            f"_lins_{self.args.long_insertion_prob}_ldels_{self.args.long_deletion_prob}"
            f"_maxinslen_{ins_len}_maxdellen_{del_len}_wobble_{self.args.wobble_prob}"
            f"{triplet_suffix}_{self.args.structure_predictor}"
        )


    def _build_rna_chain(self, chain: Dict[str, Any]) -> Dict[str, Any]:
        """Return a copy of an RNA chain with a freshly generated custom MSA
        written into its unpairedMsa field."""
        updated_chain = chain.copy()
        rna_seq = chain["rna"]["sequence"]
        logging.info("Processing RNA sequence: %s", rna_seq)
        pairs, multiplets = self.get_structure(rna_seq)
        # expose pairs for the analyzer
        self.pairs = pairs
        logging.info("Predicted pairs: %s", pairs)
        logging.info("Secondary structure has %d unique base pairs.", len(_unique_pairs(pairs)))
        if multiplets:
            logging.info("Using %d multiplets: %s", len(multiplets), multiplets)
        msa = self.generate_msa(rna_seq, pairs, multiplets)
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
            msa_plotting.plot_final_features(msa, rna_seq, pairs,
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
