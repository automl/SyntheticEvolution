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

# Triplet co-mutations. Key = seq[i] + seq[j] + seq[k] where (i, j) is a
# Watson-Crick / wobble pair (i < j) and k is the third interacting residue.
# Values are alternative 3-nt motifs that preserve the triplet's chemistry.
# For motifs absent from this table, mutate_triplet falls back to mutating
# (nt1, nt2) via PAIR_MUTATIONS and drawing nt3 uniformly from {A, U, G, C}.
TRIPLET_MUTATIONS: Dict[str, List[str]] = {
    # A-minor motifs: pair + adenosine in minor groove
    'GCA': ['GCA', 'CGA', 'AUA', 'UAA', 'GCG'],
    'CGA': ['CGA', 'GCA', 'AUA', 'UAA'],
    'AUA': ['AUA', 'UAA', 'GCA', 'CGA'],
    'UAA': ['UAA', 'AUA', 'GCA', 'CGA'],
    # Hoogsteen-like triples (U in the third position)
    'AUU': ['AUU', 'UAU', 'GCU'],
    'UAU': ['UAU', 'AUU'],
    'GCU': ['GCU', 'AUU', 'CGU'],
    'CGU': ['CGU', 'GCU', 'UAU'],
    # G-rich / quadruplex-adjacent triples
    'GCG': ['GCG', 'CGG', 'GCA'],
    'CGG': ['CGG', 'GCG'],
    # Wobble-containing triples
    'GUA': ['GUA', 'GCA', 'AUA'],
    'UGA': ['UGA', 'GUA', 'AUA'],
    'GUU': ['GUU', 'GCU', 'AUU'],
    'UGU': ['UGU', 'GUU', 'UAU'],
}

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


def _normalize_pairs(raw: Any) -> List[Tuple[int, int]]:
    """Flatten any supported pair representation (dict, list of [i, j],
    list of [i, j, ...]) into a sorted, de-duplicated list of (i, j) with
    i < j. Self-pairs are dropped."""
    out: set = set()
    if isinstance(raw, dict):
        items = raw.items()
    elif isinstance(raw, (list, tuple)):
        items = [
            (item[0], item[1])
            for item in raw
            if isinstance(item, (list, tuple)) and len(item) >= 2
        ]
    else:
        raise ValueError(f"Unknown pair format: {type(raw)}")
    for a, b in items:
        a, b = int(a), int(b)
        if a == b:
            continue
        out.add((a, b) if a < b else (b, a))
    return sorted(out)


def _unique_pairs(pairs: Dict[int, int]) -> set:
    """The set of unordered base pairs (i, j) with i < j in a pairs dict."""
    return {(i, j) for i, j in pairs.items() if i < j}


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
    # Triplet co-mutation parameters (extension of the pairwise mutation
    # scheme to three correlated positions). Default --triplet-prob=0.0
    # disables triplets, preserving the original pair-only behavior.
    mut_group.add_argument('--triplet-prob', type=float, default=0.0,
                           help="Probability that a base triple/multiplet (a connected component of "
                                ">=3 positions in the predicted base-pair graph, e.g. pairs (1,15) and "
                                "(15,18) -> {1,15,18}) is promoted and co-mutated. Set to 0 to disable.")
    mut_group.add_argument('--triplet-keep-prob', type=float, default=0.99,
                           help="Probability of keeping a multiplet's residues unchanged (analog of --stem-keep-prob).")
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
        # Support nested pseudoknot bracket families ()[]{}<> via one stack
        # per family; any other character is treated as unpaired. Behaviour
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
            logging.info("Using provided dot-bracket structure: %s", self.args.structure)
            pairs = self.pair_indices(self.args.structure)
            return pairs, self.derive_multiplets(_normalize_pairs(pairs), len(seq))
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
            return pred_pairs, self.derive_multiplets(pair_list, len(seq))
        raise ValueError("Either a structure predictor or a structure must be provided.")

    def derive_multiplets(self, pair_list: List[Tuple[int, int]],
                          seq_len: int) -> List[Dict[str, Any]]:
        """Derive base triples / higher multiplets directly from the predicted
        pairs, rather than fabricating a third base from a search window.

        The base pairs form a graph (positions = nodes, pairs = edges). A
        normal secondary-structure pair is an isolated edge: two nodes, each of
        degree 1. A *base triple* is what you get when one position pairs with
        two partners, e.g. pairs (1, 15) and (15, 18) share position 15 -> the
        connected component {1, 15, 18}. The general "further extended" case is
        just a larger connected component: a chain (1,15),(15,18),(18,30) ->
        {1,15,18,30}, or a hub paired with several partners. So: any connected
        component of size >= 3 is a multiplet; components of size 2 stay
        ordinary pairs and are handled by the existing pair machinery.

        Each surviving multiplet is kept with probability --triplet-prob. Like
        the pair sourcing, this runs once and is deterministic given --seed, so
        every MSA row sees the same multiplet set."""
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
            positions = tuple(sorted(comp))
            comp_edges = sorted(e for e in edges if e[0] in comp and e[1] in comp)
            multiplets.append({"positions": positions, "edges": comp_edges})

        sizes = [len(m["positions"]) for m in multiplets]
        logging.info(
            "Derived %d multiplets from %d base pairs (triplet_prob=%.3f); sizes=%s",
            len(multiplets), len(edges), self.args.triplet_prob, sizes)
        return multiplets

    def mutate_pair(self, nt1: str, nt2: str) -> str:
        original = nt1 + nt2
        candidates = PAIR_MUTATIONS.get(original, WC_PAIRS)
        chosen = random.choice(candidates)
        if random.random() < self.args.wobble_prob:
            return random.choice(['GU', 'UG'])
        return chosen

    def mutate_triplet(self, nt1: str, nt2: str, nt3: str) -> str:
        """Co-mutate a triplet (i, j, k). Mirrors mutate_pair: look up the
        original 3-nt motif in TRIPLET_MUTATIONS, pick a candidate, and apply
        the same wobble override on the (i, j) pair. For unseen motifs, fall
        back to mutating (nt1, nt2) as a pair and drawing nt3 uniformly."""
        original = nt1 + nt2 + nt3
        candidates = TRIPLET_MUTATIONS.get(original)
        if candidates is None:
            new_pair = self.mutate_pair(nt1, nt2)
            new_third = random.choice(['A', 'U', 'G', 'C'])
            return new_pair + new_third
        chosen = random.choice(candidates)
        if random.random() < self.args.wobble_prob:
            return random.choice(['GU', 'UG']) + chosen[2]
        return chosen

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

        Multiplets here are built from *known pairs*: each edge is a real base
        pair, including both edges of a triple (e.g. position 15 pairs with
        both 1 and 18 in {1,15,18}). So we cannot use the old TRIPLET_MUTATIONS
        table — that was for an A-minor-style model where the third base only
        loosely contacts a WC pair and need not pair with it, which would break
        the (15,18) contact here. Instead, for triples *and* the extended
        chain/hub case alike, we walk the component's edges in sorted order
        keeping a position->base assignment: the first edge of a sub-tree is
        mutated as a free pair; once one endpoint is fixed (e.g. the hub shared
        between two pairs), the other endpoint is chosen to be a valid partner
        of it. This keeps every edge a valid pair, propagates a consistent
        assignment across the whole component, and degrades to `mutate_pair`
        for a lone pair."""
        assigned: Dict[int, str] = {}
        for a, b in sorted(edges):
            a_set, b_set = a in assigned, b in assigned
            if not a_set and not b_set:
                assigned[a], assigned[b] = self.mutate_pair(seq[a], seq[b])
            elif a_set and not b_set:
                assigned[b] = self._partner_for(assigned[a])
            elif b_set and not a_set:
                assigned[a] = self._partner_for(assigned[b])
            # else: both already fixed (a cycle's closing edge) -> leave as is.
        for p in positions:
            assigned.setdefault(p, seq[p])
        return assigned

    def mutate_unpaired(self, nt: str) -> str:
        return random.choice(['A', 'U', 'G', 'C'])

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
        if random.random() < self.args.long_deletion_prob:
            del_len = random.randint(2, self.max_deletion_length) if self.max_deletion_length > 2 else random.randint(2, 5)  # 5 chosen at random
            for j in range(i, min(i + del_len, len(seq))):
                if j not in multiplet_members:  # keep multiplet geometry intact
                    mutated[j] = '-'
            return del_len
        elif random.random() < self.args.deletion_prob_loop:
            mutated[i] = '-'
        elif mutated[i] != '-':
            mutated[i] = self.mutate_unpaired(seq[i])
        mutated[i] = insertion + mutated[i]
        return 1

    def _mutate_stem_pair(self, mutated: List[str], seq: str, i: int, j: int) -> None:
        # Paired (stem) region: use stem-keep probability. Any pair reaching
        # here is an ordinary (size-2) base pair; multiplet members are handled
        # separately by the anchor branch.
        if random.random() < self.args.stem_keep_prob:
            mutated[i] = seq[i]
            mutated[j] = seq[j]
        else:
            new_pair = self.mutate_pair(seq[i], seq[j])
            mutated[i] = new_pair[0]
            mutated[j] = new_pair[1]

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
                self._mutate_stem_pair(mutated, seq, i, pairs[i])
            i += 1
        return ''.join(mutated)

    def generate_msa(self, rna_seq: str, pairs: Dict[int, int],
                     multiplets: Optional[List[Dict[str, Any]]] = None) -> List[str]:
        # Set maximum insertion/deletion lengths.
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
            f"_stemkeep_{self.args.stem_keep_prob}{triplet_suffix}_{self.args.structure_predictor}"
        )


    def _build_rna_chain(self, chain: Dict[str, Any]) -> Dict[str, Any]:
        """Return a copy of an RNA chain with a freshly generated custom MSA
        written into its unpairedMsa field."""
        updated_chain = chain.copy()
        rna_seq = chain["rna"]["sequence"]
        logging.info("Processing RNA sequence: %s", rna_seq)
        pairs, multiplets = self.get_structure(rna_seq)
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

    def process(self) -> None:
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
