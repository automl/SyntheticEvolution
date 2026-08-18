"""Builds a synthetic Stockholm file that mimics Rfam seed formatting (WUSS SS_cons, interleaved
blocks, pseudoknot letters) so calibrate_rfam.py can be exercised without the real download."""

import argparse
import random
from pathlib import Path

BASES = "ACGU"
COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}
WOBBLE_OF = {"G": "U", "U": "G"}

FAMILIES = [
    # (accession, name, WUSS consensus). Deliberately varied: nested, multi-stem, pseudoknotted.
    ("RFX0001", "SYN_hairpin", "<<<<<<<<____________>>>>>>>>"),
    ("RFX0002", "SYN_twostem", "<<<<<<<____>>>>>>>,,<<<<<<____>>>>>>"),
    ("RFX0003", "SYN_pknot", "<<<<<<<___AAAAA___>>>>>>>____aaaaa__"),
    ("RFX0004", "SYN_long", "<<<<<<<<<<<<____________>>>>>>>>>>>>...."),
]

OPEN = {"<": ">", "(": ")", "[": "]", "{": "}"}
CLOSE = {v: k for k, v in OPEN.items()}


def parse(ss):
    stacks, pairs = {}, []
    for i, c in enumerate(ss):
        if c in OPEN or (c.isalpha() and c.isupper()):
            stacks.setdefault(c, []).append(i)
        elif c in CLOSE and stacks.get(CLOSE[c]):
            pairs.append((stacks[CLOSE[c]].pop(), i))
        elif c.isalpha() and c.islower() and stacks.get(c.upper()):
            pairs.append((stacks[c.upper()].pop(), i))
    return sorted((min(a, b), max(a, b)) for a, b in pairs)


def make_row(ss, pairs, rng, mut_paired, mut_unpaired, wobble_rate, gap_rate):
    """One 'homolog': paired columns covary, unpaired columns drift, some columns gap out."""
    row = [None] * len(ss)
    partner = {}
    for a, b in pairs:
        partner[a] = b
        partner[b] = a
    for a, b in pairs:
        if rng.random() < mut_paired:
            left = rng.choice(BASES)
            right = WOBBLE_OF[left] if (left in WOBBLE_OF and rng.random() < wobble_rate) \
                else COMPLEMENT[left]
        else:
            left, right = "G", "C"
        row[a], row[b] = left, right
    for i in range(len(ss)):
        if row[i] is None:
            row[i] = rng.choice(BASES) if rng.random() < mut_unpaired else "A"
        if rng.random() < gap_rate:
            row[i] = "-"
    return "".join(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("SHS-Generator/tests/synthetic_rfam.stk"))
    parser.add_argument("--depth", type=int, default=60)
    parser.add_argument("--seed", type=int, default=3)
    parser.add_argument("--mut-paired", type=float, default=0.35)
    parser.add_argument("--mut-unpaired", type=float, default=0.55)
    parser.add_argument("--wobble-rate", type=float, default=0.18)
    parser.add_argument("--gap-rate", type=float, default=0.02)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        for accession, name, ss in FAMILIES:
            pairs = parse(ss)
            fh.write("# STOCKHOLM 1.0\n")
            fh.write(f"#=GF AC   {accession}\n#=GF ID   {name}\n\n")
            width = max(20, len(f"{name}_seq{args.depth}") + 2)
            for k in range(args.depth):
                row = make_row(ss, pairs, rng, args.mut_paired, args.mut_unpaired,
                               args.wobble_rate, args.gap_rate)
                fh.write(f"{f'{name}_seq{k}':<{width}}{row}\n")
            fh.write(f"{'#=GC SS_cons':<{width}}{ss}\n")
            fh.write("//\n")
    print(f"Wrote {len(FAMILIES)} synthetic families x {args.depth} rows to {args.out}")
    print("NOTE: synthetic test fixture only. Calibrate against the real Rfam seed archive:")
    print("  https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz")


if __name__ == "__main__":
    main()
