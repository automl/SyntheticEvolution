"""Normalized RNA base-pair mapping and multiplet derivation."""
import json
import logging
import random
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import numpy as np


def _iter_pairs(raw: Any):
    """Yield (i, j) int tuples from dicts or lists of pairs."""
    if isinstance(raw, dict):
        items = raw.items()
    elif isinstance(raw, (list, tuple)):
        items = [(item[0], item[1]) for item in raw
                 if isinstance(item, (list, tuple)) and len(item) >= 2]
    else:
        raise ValueError(f"Unknown pair format: {type(raw)}")

    for a, b in items:
        yield int(a), int(b)


def _normalize_pairs(raw: Any) -> List[Tuple[int, int]]:
    out: set = set()
    for a, b in _iter_pairs(raw):
        if a == b:
            continue
        out.add((a, b) if a < b else (b, a))
    return sorted(out)


def _convert_pairs(raw: Any) -> Dict[int, int]:
    result: Dict[int, int] = {}
    for i, j in _iter_pairs(raw):
        result[i] = j
        result[j] = i
    return result


def _derive_multiplets(pair_list: List[Tuple[int, int]], triplet_prob: float) -> List[Dict[str, Any]]:
    if triplet_prob <= 0.0 or not pair_list:
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
            continue
        if random.random() >= triplet_prob:
            continue

        multiplets.append({
            "positions": tuple(sorted(comp)),
            "edges": sorted(e for e in edges if e[0] in comp and e[1] in comp),
        })

    logging.info("Derived %d multiplets from %d base pairs (triplet_prob=%.3f)",
                 len(multiplets), len(edges), triplet_prob)
    return multiplets


@dataclass
class PairMap:
    pairs: Dict[int, int] = field(default_factory=dict)
    multiplets: List[Dict[str, Any]] = field(default_factory=list)

    @classmethod
    def from_dotbracket(cls, dotbracket: str, triplet_prob: float = 0.0) -> "PairMap":
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

        return cls(pairs=pairs, multiplets=_derive_multiplets(_normalize_pairs(pairs), triplet_prob))

    @classmethod
    def from_raw(cls, raw: Any, triplet_prob: float = 0.0) -> "PairMap":
        if isinstance(raw, PairMap):
            return raw

        if isinstance(raw, str):
            try:
                parsed = json.loads(raw)
            except json.JSONDecodeError:
                return cls.from_dotbracket(raw, triplet_prob=triplet_prob)
            return cls.from_raw(parsed, triplet_prob=triplet_prob)

        if isinstance(raw, dict) and "pairs" in raw and isinstance(raw["pairs"], (dict, list, tuple)):
            raw = raw["pairs"]

        pairs = _convert_pairs(raw)
        return cls(pairs=pairs, multiplets=_derive_multiplets(_normalize_pairs(pairs), triplet_prob))

    @classmethod
    def from_json_pairs(cls, data: Dict[str, Any], triplet_prob: float = 0.0) -> "PairMap":
        return cls.from_raw(data, triplet_prob=triplet_prob)

    @classmethod
    def from_predicted_pairs(cls, raw: Any, triplet_prob: float = 0.0) -> "PairMap":
        return cls.from_raw(raw, triplet_prob=triplet_prob)

    def partner(self, i: int) -> Optional[int]:
        return self.pairs.get(i)

    def is_paired(self, i: int) -> bool:
        return i in self.pairs

    @property
    def unique_pairs(self) -> set:
        return {(i, j) for i, j in self.pairs.items() if i < j}

    @property
    def anchor_to_multiplet(self) -> Dict[int, Dict[str, Any]]:
        return {m["positions"][0]: m for m in self.multiplets}

    @property
    def multiplet_members(self) -> set:
        return {p for m in self.multiplets for p in m["positions"]}