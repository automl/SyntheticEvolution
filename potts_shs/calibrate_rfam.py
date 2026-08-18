"""Calibrate the Potts coupling and wobble weight against Rfam seed alignments, by matching
covariation statistics at the depth the generator is actually used at.

On Helix the seed archive is already unpacked at evaluation/data/rfam/Rfam.seed.gz:

    python SHS-Generator/calibrate_rfam.py \
        --seed-file evaluation/data/rfam/Rfam.seed.gz \
        --out results/potts_calibration.json

Calibrating at a depth needs families with at least that many seed sequences, which is the
binding constraint: of the 4227 families in the archive, 706 qualify at depth 20, 314 at 50 and
164 at 100 within the 15-200 nt window. `--max-families 300` is therefore not reached at depth
100 -- the run uses every family that qualifies.

Timing at depth 100 (measured): 5s to scan and prepare the archive, then ~15s per grid point over
all 164 families at --repeats 3. The full default 9x6 grid is about 14 minutes on one core.
"""

import argparse
import gzip
import json
import logging
import math
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

HERE_DIR = Path(__file__).resolve().parent
GENERATOR_DIR = HERE_DIR.parent / "SHS-Generator"
for _path in (HERE_DIR, GENERATOR_DIR):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from pair_map import PairMap
from potts import PottsModel

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

RFAM_SEED_URL = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz"

BASES = ("A", "C", "G", "U")
CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE = {("G", "U"), ("U", "G")}

# WUSS: matched bracket pairs, plus upper/lower letter pairs for pseudoknots.
_WUSS_OPEN = {"<": ">", "(": ")", "[": "]", "{": "}"}
_WUSS_CLOSE = {v: k for k, v in _WUSS_OPEN.items()}
_WUSS_UNPAIRED = set(".,_-:~")


# -- Rfam input ------------------------------------------------------------

@dataclass
class Family:
    accession: str
    name: str
    sequences: List[str] = field(default_factory=list)
    ss_cons: str = ""


def read_stockholm(path: Path) -> Iterator[Family]:
    """Stream families out of a (possibly gzipped, possibly interleaved) Stockholm file."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        rows: Dict[str, List[str]] = {}
        order: List[str] = []
        ss_parts: List[str] = []
        accession = name = ""
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("//"):
                if rows:
                    yield Family(accession or "NA", name or "NA",
                                 ["".join(rows[k]) for k in order], "".join(ss_parts))
                rows, order, ss_parts = {}, [], []
                accession = name = ""
            elif line.startswith("#=GF AC"):
                accession = line.split(maxsplit=2)[2].strip()
            elif line.startswith("#=GF ID"):
                name = line.split(maxsplit=2)[2].strip()
            elif line.startswith("#=GC SS_cons"):
                ss_parts.append(line.split(maxsplit=2)[2])
            elif line.startswith("#") or not line.strip():
                continue
            else:
                parts = line.split(maxsplit=1)
                if len(parts) != 2:
                    continue
                seq_id, chunk = parts
                if seq_id not in rows:
                    rows[seq_id] = []
                    order.append(seq_id)
                rows[seq_id].append(chunk.strip())
        if rows:
            yield Family(accession or "NA", name or "NA",
                         ["".join(rows[k]) for k in order], "".join(ss_parts))


def parse_wuss(ss: str) -> List[Tuple[int, int]]:
    """Base pairs from a WUSS consensus string, including its pseudoknot letter pairs."""
    stacks: Dict[str, List[int]] = {}
    pairs: List[Tuple[int, int]] = []
    for i, char in enumerate(ss):
        if char in _WUSS_OPEN or (char.isalpha() and char.isupper()):
            stacks.setdefault(char, []).append(i)
        elif char in _WUSS_CLOSE:
            opener = _WUSS_CLOSE[char]
            if stacks.get(opener):
                pairs.append((stacks[opener].pop(), i))
        elif char.isalpha() and char.islower():
            opener = char.upper()
            if stacks.get(opener):
                pairs.append((stacks[opener].pop(), i))
        elif char in _WUSS_UNPAIRED:
            continue
    return sorted((min(a, b), max(a, b)) for a, b in pairs)


def normalise(row: str) -> str:
    out = []
    for char in row.upper():
        if char in ("-", ".", "~", "_"):
            out.append("-")
        elif char == "T":
            out.append("U")
        elif char in BASES:
            out.append(char)
        else:
            out.append("N")
    return "".join(out)


@dataclass
class Case:
    accession: str
    name: str
    query: str
    pairs: List[Tuple[int, int]]
    rows: List[str]


def prepare(family: Family, depth: int, min_len: int, max_len: int,
            rng: random.Random) -> Optional[Case]:
    """Trim to a reference row, remap the consensus pairs onto it, and subsample to `depth`."""
    if not family.sequences or not family.ss_cons:
        return None
    width = len(family.sequences[0])
    if len(family.ss_cons) != width or any(len(s) != width for s in family.sequences):
        return None

    rows = [normalise(s) for s in family.sequences]
    # Reference = the row with the most resolved (ACGU) columns; it defines the query coordinates.
    reference = max(rows, key=lambda r: sum(c in BASES for c in r))
    keep = [i for i, c in enumerate(reference) if c in BASES]
    if not (min_len <= len(keep) <= max_len):
        return None

    remap = {col: k for k, col in enumerate(keep)}
    pairs = [(remap[a], remap[b]) for a, b in parse_wuss(family.ss_cons)
             if a in remap and b in remap]
    if len(pairs) < 3:
        return None

    others = [r for r in rows if r is not reference]
    if len(others) < depth - 1:
        return None
    sampled = rng.sample(others, depth - 1)
    trimmed = ["".join(r[i] for i in keep) for r in [reference, *sampled]]
    return Case(family.accession, family.name, trimmed[0], pairs, trimmed)


# -- statistics ------------------------------------------------------------

def mi_at(rows: Sequence[str], columns: Sequence[Tuple[int, int]]) -> np.ndarray:
    """Mutual information in bits for the given column pairs only (A/C/G/U/-/N alphabet)."""
    alphabet = BASES + ("-", "N")
    lookup = {b: k for k, b in enumerate(alphabet)}
    codes = np.array([[lookup.get(c, 5) for c in row] for row in rows])
    n = codes.shape[0]
    out = np.zeros(len(columns))
    for k, (i, j) in enumerate(columns):
        joint = np.zeros((len(alphabet), len(alphabet)))
        np.add.at(joint, (codes[:, i], codes[:, j]), 1.0)
        joint /= n
        expected = np.outer(joint.sum(axis=1), joint.sum(axis=0))
        nz = joint > 0
        out[k] = float(np.sum(joint[nz] * np.log2(joint[nz] / expected[nz])))
    return out


def background_columns(length: int, pairs: Sequence[Tuple[int, int]], count: int,
                       rng: random.Random) -> List[Tuple[int, int]]:
    truth = {(min(a, b), max(a, b)) for a, b in pairs}
    out, attempts = [], 0
    while len(out) < count and attempts < count * 40:
        attempts += 1
        i, j = sorted(rng.sample(range(length), 2))
        if (i, j) not in truth:
            out.append((i, j))
    return out


def measure(rows: Sequence[str], query: str, pairs: Sequence[Tuple[int, int]],
            rng: random.Random) -> Dict[str, float]:
    """Covariation and composition statistics of one alignment."""
    arr = np.array([list(r) for r in rows[1:]])
    query_arr = np.array(list(query))
    resolved = np.isin(arr, list(BASES))
    differs = (arr != query_arr) & resolved

    paired_cols = sorted({p for pair in pairs for p in pair})
    unpaired_cols = [i for i in range(len(query)) if i not in set(paired_cols)]

    def rate(cols: Sequence[int]) -> float:
        if not cols:
            return float("nan")
        denom = resolved[:, cols].sum()
        return float(differs[:, cols].sum() / denom) if denom else float("nan")

    canonical = wobble = other = 0
    for row in rows[1:]:
        for i, j in pairs:
            combo = (row[i], row[j])
            if combo[0] not in BASES or combo[1] not in BASES:
                continue
            if combo in CANONICAL:
                canonical += 1
            elif combo in WOBBLE:
                wobble += 1
            else:
                other += 1
    total = canonical + wobble + other

    mi_pairs = mi_at(rows[1:], list(pairs))
    background = background_columns(len(query), pairs, min(len(pairs) * 5, 400), rng)
    mi_background = mi_at(rows[1:], background) if background else np.array([0.0])

    return {
        "mut_paired": rate(paired_cols),
        "mut_unpaired": rate(unpaired_cols),
        "canonical": canonical / total if total else float("nan"),
        "wobble": wobble / total if total else float("nan"),
        "non_pair": other / total if total else float("nan"),
        "mi_pairs": float(np.mean(mi_pairs)) if len(mi_pairs) else float("nan"),
        "mi_background": float(np.mean(mi_background)),
        "gap_fraction": float(np.mean(arr == "-")),
    }


def potts_measure(case: Case, mu_paired: float, mu_unpaired: float, gamma: float,
                  wobble: float, depth: int, repeats: int, rng: random.Random) -> Dict[str, float]:
    """Same statistics, for alignments the Potts sampler produces on this case's structure."""
    paired = {p for pair in case.pairs for p in pair}
    rates = [mu_paired if i in paired else mu_unpaired for i in range(len(case.query))]
    pair_map = PairMap.from_interactions(
        len(case.query), [(a, b, 1.0) for a, b in case.pairs], rates)
    model = PottsModel(case.query, pair_map, coupling=gamma, wobble=wobble)

    runs = []
    for _ in range(repeats):
        rows = [case.query] + ["".join(model.sample()) for _ in range(depth - 1)]
        runs.append(measure(rows, case.query, case.pairs, rng))
    return {k: float(np.nanmean([r[k] for r in runs])) for k in runs[0]}


# -- driver ----------------------------------------------------------------

def discrepancy(target: Dict[str, float], candidate: Dict[str, float]) -> float:
    """Match the covariation signal first, the wobble share second."""
    mi = abs(target["mi_pairs"] - candidate["mi_pairs"])
    wob = abs(target["wobble"] - candidate["wobble"])
    return mi + 2.0 * wob


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-file", type=Path, required=True,
                        help=f"Rfam seed Stockholm file (.stk or .stk.gz). Download: {RFAM_SEED_URL}")
    parser.add_argument("--exclude", type=Path,
                        help="File of Rfam accessions or family names to skip, one per line. Use this "
                             "to drop every family present in the 919-target benchmark")
    parser.add_argument("--depth", type=int, default=100,
                        help="Alignment depth to calibrate at; must match the depth you will generate at")
    parser.add_argument("--min-len", type=int, default=15)
    parser.add_argument("--max-len", type=int, default=200)
    parser.add_argument("--max-families", type=int, default=300,
                        help="Upper bound only; at depth 100 just 164 families qualify and all "
                             "of them are used")
    parser.add_argument("--repeats", type=int, default=3, help="Potts draws per family per grid point")
    parser.add_argument("--gamma-grid", type=float, nargs="+",
                        default=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0])
    parser.add_argument("--wobble-grid", type=float, nargs="+",
                        default=[0.70, 0.80, 0.85, 0.90, 0.95, 1.00])
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=Path, default=Path("results/potts_calibration.json"))
    args = parser.parse_args()

    excluded = set()
    if args.exclude and args.exclude.exists():
        excluded = {l.strip().upper() for l in args.exclude.read_text().splitlines() if l.strip()}
        logging.info("Excluding %d families from calibration", len(excluded))

    rng = random.Random(args.seed)
    cases: List[Case] = []
    scanned = 0
    for family in read_stockholm(args.seed_file):
        scanned += 1
        if family.accession.upper() in excluded or family.name.upper() in excluded:
            continue
        case = prepare(family, args.depth, args.min_len, args.max_len, rng)
        if case is not None:
            cases.append(case)
        if len(cases) >= args.max_families:
            break
    logging.info("Scanned %d families, kept %d at depth %d in [%d, %d] nt",
                 scanned, len(cases), args.depth, args.min_len, args.max_len)
    if not cases:
        logging.error("No usable families. Check --seed-file, --depth and the length window.")
        return 1

    targets = [measure(c.rows, c.query, c.pairs, rng) for c in cases]
    mu_paired = float(np.nanmedian([t["mut_paired"] for t in targets]))
    mu_unpaired = float(np.nanmedian([t["mut_unpaired"] for t in targets]))
    pooled_target = {k: float(np.nanmedian([t[k] for t in targets])) for k in targets[0]}

    print("\nRfam targets (median over %d families, depth %d)" % (len(cases), args.depth))
    for key in ("mut_paired", "mut_unpaired", "canonical", "wobble", "non_pair",
                "mi_pairs", "mi_background", "gap_fraction"):
        print(f"  {key:<16}{pooled_target[key]:.4f}")
    print(f"\n=> recommended --mutation-rate-paired {mu_paired:.3f} "
          f"--mutation-rate-unpaired {mu_unpaired:.3f}")

    print(f"\nSweeping {len(args.gamma_grid)}x{len(args.wobble_grid)} grid at those rates...")
    header = f"{'gamma':>7}{'wobble':>8}{'MI@pairs':>10}{'wobble_f':>10}{'canon':>8}{'non_pair':>10}{'score':>9}"
    print(header)
    print("-" * len(header))

    results = []
    for gamma in args.gamma_grid:
        for wob in args.wobble_grid:
            per_case = [potts_measure(c, mu_paired, mu_unpaired, gamma, wob,
                                      args.depth, args.repeats, rng) for c in cases]
            pooled = {k: float(np.nanmedian([p[k] for p in per_case])) for k in per_case[0]}
            score = discrepancy(pooled_target, pooled)
            # Grid parameters are prefixed so they cannot collide with the measured statistics,
            # which carry their own 'wobble' key (the realised wobble fraction).
            results.append({"grid_gamma": gamma, "grid_wobble": wob, "score": score, **pooled})
            print(f"{gamma:>7.1f}{wob:>8.2f}{pooled['mi_pairs']:>10.4f}{pooled['wobble']:>10.4f}"
                  f"{pooled['canonical']:>8.4f}{pooled['non_pair']:>10.4f}{score:>9.4f}")

    best = min(results, key=lambda r: r["score"])
    print(f"\nBest fit: --potts-coupling {best['grid_gamma']} --potts-wobble {best['grid_wobble']} "
          f"(score {best['score']:.4f})")
    print(f"  Rfam  MI@pairs={pooled_target['mi_pairs']:.4f} wobble_fraction={pooled_target['wobble']:.4f}")
    print(f"  Potts MI@pairs={best['mi_pairs']:.4f} wobble_fraction={best['wobble']:.4f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        "depth": args.depth,
        "n_families": len(cases),
        "families": [c.accession for c in cases],
        "rfam_target": pooled_target,
        "recommended": {
            "mutation_rate_paired": mu_paired,
            "mutation_rate_unpaired": mu_unpaired,
            "potts_coupling": best["grid_gamma"],
            "potts_wobble": best["grid_wobble"],
        },
        "grid": results,
    }, indent=2))
    print(f"\nWritten to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
