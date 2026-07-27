"""Normalized RNA base-pair mapping and multiplet derivation."""
import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List
import numpy as np
from collections import deque
from functools import cached_property

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


@dataclass
class PairMap:
    """Unified representation for secondary structure of rna. Including mutation rates
    
    Attributes:
        pairs_mat (N, N): matrix for every combination of possible pairs.
            Main diagonal: mutation rate for each position
            Off diagonal: Interaction strength for each combination of two positions. 
                Interpreted differently for differently by different mutation approaches.
            The matrix is mirrored along the main diagonal.

    """
    def __init__(self, pairs_mat: np.ndarray):
        self._pairs_mat = pairs_mat.copy()
        self._pairs_mat.flags.writeable = False
        self._direct_partners_cache: Dict[int, List[int]] = {
            i : [j for j, val in enumerate(self._pairs_mat[i]) if j != i and val > 0]
            for i in range(self._pairs_mat.shape[0])
        }
        self._partners_cache: Dict[int, List[int]] = {
            i: sorted(self._compute_partners(i))
            for i in range(self._pairs_mat.shape[0])
        }

    @classmethod
    def from_interactions(cls, size: int, pairs: list[tuple[int, int, float]], mutation_rates: list[int]):
        pairs_mat = np.zeros((size, size))
        for i, j, val in pairs:
            if i >= size or j >= size or i < 0 or j < 0:
                logging.warning("Skipping invalid pair %s for seq of length %s.", (i, j), size)
            if i == j:
                logging.warning("Self pairs found in pairs at position %s.", i)
                continue
            if pairs_mat[i][j] != 0:
                logging.warning("Tried to reassign interaction strength for pair %s. Using last provided value.", (i, j))
            if val < 0 or val > 1:
                logging.warning("Interaction strength for pair %s is %s, using %s instead.", (i, j), val, min(max(0, val), 1))
                val = min(max(0, val), 1)
            pairs_mat[i][j] = val
            pairs_mat[j][i] = val
        if len(mutation_rates) > size:
            logging.warning("Provided too many mutation rates (%s). Using the first (%s)", len(mutation_rates), size)
        if len(mutation_rates) < size:
            logging.error("Not enough mutation rates provided for sequence with size %s.", size)
            raise ValueError()
        for i in range(size):
            val = mutation_rates[i]
            if val < 0 or val > 1:
                logging.warning("Mutation rate at position %s is %s, using %s instead.", i, val, min(max(0, val), 1))
                val = min(max(0, val), 1)
            pairs_mat[i][i] = val
        return cls(pairs_mat = pairs_mat)
        
    @classmethod
    def from_raw(cls, raw: Any, size: int, mutation_rate_paired: int, 
                mutation_rate_unpaired: int, interaction_strength: float = 1) -> "PairMap":
        """
        Input formats:
            Datatype | [Example]
            str      | "((((...))<..))>"
            dict     | "{0: 1, 1: 0}"
            list     | "[[0, 1], [5, 2]]" or "[(0, 1), (5, 2)]"
            tuple    | "([0, 1], [5, 2])" or "((0, 1), (5, 2))"
        """
        if isinstance(raw, str):
            pairs, mutation_rates = PairMap._interactions_from_dotbracket(raw, mutation_rate_paired, mutation_rate_unpaired, interaction_strength)
        elif isinstance(raw, (dict, list, tuple)):
            pairs, mutation_rates = PairMap._interactions_from_iterable(raw, size, mutation_rate_paired, mutation_rate_unpaired, interaction_strength)
        else:
            logging.error("Unknown pair format: %s", type(raw))
            # try the iterable method anyway, maybe it works idk. This is fine as all safety checks are performed later in cls.from_interactions
            pairs, mutation_rates = PairMap._interactions_from_iterable(raw, size, mutation_rate_paired, mutation_rate_unpaired, interaction_strength)
        return cls.from_interactions(size, pairs, mutation_rates)


    @staticmethod
    def _interactions_from_dotbracket(dotbracket: str, 
                                       mutation_rate_paired: float, 
                                       mutation_rate_unpaired: float, 
                                       interaction_strength: float) -> tuple[list[tuple[int, int, float]], list[int]]:
        openers = {'(': ')', '[': ']', '{': '}', '<': '>'}
        close_to_open = {c: o for o, c in openers.items()}
        stacks: Dict[str, list] = {o: [] for o in openers}

        pairs: List[tuple[int, int, float]] = []
        mutation_rates: list[int] = []

        for i, char in enumerate(dotbracket):
            if char in openers:
                stacks[char].append(i)
                mutation_rates.append(mutation_rate_paired)
            elif char in close_to_open:
                o = close_to_open[char]
                if not stacks[o]:
                    logging.error("Unbalanced structure at position %d", i)
                    continue
                j = stacks[o].pop()
                pairs.append((j, i, interaction_strength))
                mutation_rates.append(mutation_rate_paired)
            else:
                mutation_rates.append(mutation_rate_unpaired)
        return pairs, mutation_rates

    @staticmethod
    def _interactions_from_iterable(size: int, pairs: dict, 
                                    mutation_rate_paired: float, 
                                    mutation_rate_unpaired: float, 
                                    interaction_strength: float) -> tuple[list[tuple[int, int, float]], list[int]]:
        unique_pairs: set[tuple[int, int]] = set()
        is_paired: set[int] = set()
        mutation_rates: list[int] = []

        if isinstance(pairs, dict):
            items = pairs.items()
        items = [(item[0], item[1]) for item in pairs]

        for a, b in items:
            a, b = int(a), int(b)
            if a == b:
                continue
            unique_pairs.add((a, b) if a < b else (b, a))
            is_paired.add(a)
            is_paired.add(b)
        for i in range(size):
            mutation_rates.append(mutation_rate_paired if i in is_paired else mutation_rate_unpaired)
        return [(a, b, interaction_strength) for a, b in unique_pairs], mutation_rates

    def _compute_partners(self, i: int) -> set[int]:
        partners = set()
        queue = deque([i])
        while queue:
            j = queue.popleft()
            if j not in partners:
                partners.add(j)
                queue.extend(self.direct_partners(j))
        partners.remove(i)
        return partners

    def direct_partners(self, i: int) -> List[int]:
        return self._direct_partners_cache[i]

    def partners(self, i: int) -> List[int]:
        return self._partners_cache[i]

    def is_multiplet_member(self, i: int):
        return len(self.partners(i)) > 1

    def is_multiplet_anchor(self, i: int) -> bool:
        return self.is_multiplet_member(i) and i < min(self.partners(i))

    def is_unpaired(self, i: int) -> bool:
        return len(self.partners(i)) == 0

    def is_basic_pair(self, i: int) -> bool:
        return len(self.partners(i)) == 1

    def is_paired(self, i: int) -> bool:
        return len(self.partners(i)) >= 1

    def is_pair(self, i: int, j: int) -> bool:
        return self._pairs_mat[i][j] > 0

    @cached_property
    def multiplets(self) -> dict[int, list[int]]:
        multiplets:  dict[int, list[int]] = dict()
        for i in range(self._pairs_mat.shape[0]):
            partners = self.partners(i)
            if len(partners) <= 1:
                continue
            if not (i in multiplets or partners[0] in multiplets):
                multiplets[i] = partners
        return multiplets

    @cached_property
    def pairs(self) -> List[tuple[int, int]]:
        i_indices, j_indices = np.where((self._pairs_mat > 0) & (np.eye(self._pairs_mat.shape[0], dtype=bool) == False))
        return list(zip(i_indices, j_indices))

    @cached_property
    def unique_pairs(self) -> set:
        return np.sum(self._pairs_mat[np.triu_indices(self._pairs_mat.shape[0], k=1)] > 0)


