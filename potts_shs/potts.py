"""Potts (direct-coupling) substitution model for SHS generation: per-position fields taken
from the mutation rates and pairwise couplings taken from the predicted interaction map."""

import itertools
import logging
import math
import random
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

# The pair representation still lives with the generator it belongs to.
GENERATOR_DIR = Path(__file__).resolve().parents[1] / "SHS-Generator"
if str(GENERATOR_DIR) not in sys.path:
    sys.path.insert(0, str(GENERATOR_DIR))

from pair_map import PairMap

BASES: Tuple[str, ...] = ("A", "C", "G", "U")

_CANONICAL = frozenset({("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")})
_WOBBLE = frozenset({("G", "U"), ("U", "G")})

_EPS = 1e-9

# Edge of a coupled component: (local index u, local index v, interaction strength).
Edge = Tuple[int, int, float]


def pair_score(a: str, b: str, wobble: float) -> float:
    """Pairing score of an ordered base combination. Symmetric in (a, b)."""
    if (a, b) in _CANONICAL:
        return 1.0
    if (a, b) in _WOBBLE:
        return wobble
    return 0.0


def _field(t: float, is_query: bool) -> float:
    """Unnormalised field weight: t on the query base, the remainder split over the other three."""
    return t if is_query else (1.0 - t) / 3.0


class _Table:
    """Exactly tabulated joint distribution over one coupled component."""

    def __init__(self, assignments: List[Tuple[str, ...]], cum_weights: List[float]) -> None:
        self.assignments = assignments
        self.cum_weights = cum_weights

    def draw(self) -> Tuple[str, ...]:
        return random.choices(self.assignments, cum_weights=self.cum_weights, k=1)[0]


class _GibbsChain:
    """Fallback sampler for components too large to tabulate. Restarted from a random state on
    every draw so that generated rows stay independent."""

    def __init__(self, size: int, query: Sequence[str], mus: Sequence[float],
                 edges: Sequence[Edge], coupling: float, wobble: float, sweeps: int) -> None:
        self.size = size
        self.query = list(query)
        self.t = [_clamp01(1.0 - mu) for mu in mus]
        self.coupling = coupling
        self.wobble = wobble
        self.sweeps = sweeps
        self.neighbours: List[List[Tuple[int, float]]] = [[] for _ in range(size)]
        for u, v, strength in edges:
            self.neighbours[u].append((v, strength))
            self.neighbours[v].append((u, strength))

    def draw(self) -> Tuple[str, ...]:
        state = [random.choice(BASES) for _ in range(self.size)]
        for _ in range(self.sweeps):
            for p in range(self.size):
                weights = []
                for base in BASES:
                    log_coupling = sum(
                        self.coupling * strength * pair_score(base, state[q], self.wobble)
                        for q, strength in self.neighbours[p]
                    )
                    weights.append(_field(self.t[p], base == self.query[p]) * math.exp(log_coupling))
                if sum(weights) <= 0.0:
                    weights = [1.0] * len(BASES)
                state[p] = random.choices(BASES, weights=weights, k=1)[0]
        return tuple(state)


def _clamp01(value: float) -> float:
    return min(max(value, 0.0), 1.0)


class PottsModel:
    """Samples the substitution layer of an SHS row from a Potts model.

    The distribution over a coupled component of positions is

        P(x) proportional to  prod_i f_i(x_i) * exp( gamma * sum_ij p_ij * s(x_i, x_j) )

    where f_i places weight t_i on the query base, s scores canonical pairs 1 and wobbles
    ``wobble``, and p_ij is the interaction strength held on the PairMap off-diagonal.
    Indels are not modelled here; the caller keeps its own insertion and deletion layer.
    """

    def __init__(self, sequence: str, pair_map: PairMap, coupling: float = 6.0,
                 wobble: float = 0.85, match_rates: bool = True, max_enumerate: int = 7,
                 max_rate_match: int = 4, gibbs_sweeps: int = 20) -> None:
        self.sequence = sequence.upper()
        self.pair_map = pair_map
        self.coupling = coupling
        self.wobble = wobble
        self.match_rates = match_rates
        self.max_enumerate = max_enumerate
        self.max_rate_match = max_rate_match
        self.gibbs_sweeps = gibbs_sweeps

        self.components: List[List[int]] = _components(pair_map, len(self.sequence))
        self._cache: Dict[Tuple, _Table] = {}
        self._samplers: List[Tuple[List[int], object]] = []
        for positions in self.components:
            self._samplers.append((positions, self._build(positions)))

        sizes = [len(c) for c in self.components]
        logging.info("Potts model: %d components, largest %d, gamma=%s wobble=%s rate matching=%s",
                     len(sizes), max(sizes) if sizes else 0, coupling, wobble, match_rates)

    # -- construction ------------------------------------------------------

    def _build(self, positions: Sequence[int]):
        size = len(positions)
        query = tuple(self.sequence[p] for p in positions)
        mus = tuple(_clamp01(self.pair_map.mutation_rate(p)) for p in positions)
        edges = tuple(
            (u, v, self.pair_map.interaction(positions[u], positions[v]))
            for u in range(size) for v in range(u + 1, size)
            if self.pair_map.interaction(positions[u], positions[v]) > 0
        )

        if size > self.max_enumerate:
            logging.warning("Component of size %d exceeds --potts-max-enumerate=%d; falling back to "
                            "Gibbs sampling without rate matching.", size, self.max_enumerate)
            return _GibbsChain(size, query, mus, edges, self.coupling, self.wobble, self.gibbs_sweeps)

        key = (query, mus, edges)
        table = self._cache.get(key)
        if table is None:
            table = self._tabulate(size, query, mus, edges)
            self._cache[key] = table
        return table

    def _tabulate(self, size: int, query: Tuple[str, ...], mus: Tuple[float, ...],
                  edges: Tuple[Edge, ...]) -> _Table:
        assignments = list(itertools.product(BASES, repeat=size))
        is_query = [tuple(a[p] == query[p] for p in range(size)) for a in assignments]

        log_coupling = [
            sum(self.coupling * strength * pair_score(a[u], a[v], self.wobble) for u, v, strength in edges)
            for a in assignments
        ]
        shift = max(log_coupling) if log_coupling else 0.0
        coupling_weight = [math.exp(lc - shift) for lc in log_coupling]

        if self.match_rates and edges and size <= self.max_rate_match:
            t = _solve_fields(size, mus, is_query, coupling_weight)
        else:
            if self.match_rates and edges and size > self.max_rate_match:
                logging.warning("Component of size %d exceeds --potts-max-rate-match=%d; using raw "
                                "fields, so realised mutation rates will fall below the requested ones.",
                                size, self.max_rate_match)
            t = [_clamp01(1.0 - mu) for mu in mus]

        weights = _weights(t, is_query, coupling_weight)
        total = sum(weights)
        if total <= 0.0:
            logging.error("Degenerate Potts component (query=%s); falling back to a uniform draw.", query)
            weights = [1.0] * len(assignments)

        cum, running = [], 0.0
        for w in weights:
            running += w
            cum.append(running)
        return _Table(assignments, cum)

    # -- sampling ----------------------------------------------------------

    def sample(self) -> List[str]:
        """Draw one substituted sequence. Returns a list of bases, gaps excluded."""
        out = list(self.sequence)
        for positions, sampler in self._samplers:
            draw = sampler.draw()
            for local, p in enumerate(positions):
                out[p] = draw[local]
        return out


def _components(pair_map: PairMap, length: int) -> List[List[int]]:
    """Connected components of the interaction graph; unpaired positions are singletons."""
    seen = set()
    components: List[List[int]] = []
    for i in range(length):
        if i in seen:
            continue
        component = sorted({i, *pair_map.partners(i)})
        seen.update(component)
        components.append(component)
    return components


def _weights(t: Sequence[float], is_query: Sequence[Tuple[bool, ...]],
             coupling_weight: Sequence[float]) -> List[float]:
    out = []
    for flags, cw in zip(is_query, coupling_weight):
        w = cw
        for p, flag in enumerate(flags):
            w *= _field(t[p], flag)
        out.append(w)
    return out


def _marginal(position: int, t: Sequence[float], is_query: Sequence[Tuple[bool, ...]],
              coupling_weight: Sequence[float]) -> float:
    """P(x_position != query base) under the current fields."""
    total = 0.0
    off = 0.0
    for flags, cw in zip(is_query, coupling_weight):
        w = cw
        for p, flag in enumerate(flags):
            w *= _field(t[p], flag)
        total += w
        if not flags[position]:
            off += w
    return off / total if total > 0.0 else 0.0


def _solve_fields(size: int, mus: Sequence[float], is_query: Sequence[Tuple[bool, ...]],
                  coupling_weight: Sequence[float], rounds: int = 12,
                  steps: int = 40) -> List[float]:
    """Choose fields so that each position's realised mutation rate matches the requested one.

    Coupling suppresses substitution, so the raw field 1-mu undershoots. The marginal is monotone
    decreasing in t, so a coordinate-wise bisection converges quickly.
    """
    t = [_clamp01(1.0 - mu) for mu in mus]
    free = []
    for p, mu in enumerate(mus):
        if mu <= 0.0:
            t[p] = 1.0
        elif mu >= 1.0:
            t[p] = 0.0
        else:
            free.append(p)

    for _ in range(rounds):
        shift = 0.0
        for p in free:
            previous = t[p]
            lo, hi = _EPS, 1.0 - _EPS
            for _ in range(steps):
                mid = 0.5 * (lo + hi)
                t[p] = mid
                if _marginal(p, t, is_query, coupling_weight) > mus[p]:
                    lo = mid
                else:
                    hi = mid
            t[p] = 0.5 * (lo + hi)
            shift = max(shift, abs(t[p] - previous))
        if shift < 1e-8:
            break
    return t
