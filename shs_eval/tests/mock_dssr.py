"""Stand-in for x3dna-dssr that emits a JSON in DSSR's shape (chains/nts/pairs) from the mock
mmCIF, so the pipeline's DSSR parsing and scoring can be tested without DSSR installed."""

import argparse
import json
import random
from pathlib import Path


def read_mock_cif(cif: Path):
    residues = []
    for line in cif.read_text().splitlines():
        if line.startswith("ATOM "):
            parts = line.split()
            residues.append({"comp": parts[3], "chain": parts[4], "index": int(parts[5])})
    return residues


def nested_pairs(sequence, rng, drop_rate):
    """Greedy nested pairing of complementary bases, then drop a fraction. Not a folding model —
    it exists only to produce a plausible, non-empty pair list."""
    complement = {"A": "U", "U": "A", "G": "C", "C": "G"}
    pairs, used = [], set()
    left, right = 0, len(sequence) - 1
    while left < right:
        if left in used or right in used:
            break
        if complement.get(sequence[left]) == sequence[right] and right - left > 3:
            pairs.append((left, right))
            used.add(left)
            used.add(right)
            left += 1
            right -= 1
        elif rng.random() < 0.5:
            left += 1
        else:
            right -= 1
    return [p for p in pairs if rng.random() > drop_rate]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="cif", required=True)
    parser.add_argument("-o", "--output", dest="out", required=True)
    args, _ = parser.parse_known_args()

    cif = Path(args.cif.split("=", 1)[-1] if args.cif.startswith("=") else args.cif)
    out = Path(args.out.split("=", 1)[-1] if args.out.startswith("=") else args.out)

    residues = read_mock_cif(cif)
    if not residues:
        raise SystemExit(f"mock_dssr: no ATOM records in {cif}")
    chain = residues[0]["chain"]
    sequence = "".join(r["comp"] for r in residues)

    truth_path = cif.parent / "_mock_truth.json"
    truth = json.loads(truth_path.read_text()) if truth_path.exists() else {}
    rng = random.Random(truth.get("seed", 0) ^ len(sequence))
    pairs = nested_pairs(sequence, rng, truth.get("drop_rate", 0.15))

    nts = [{"nt_id": f"{chain}.{r['comp']}{r['index']}", "nt_name": r["comp"],
            "nt_code": r["comp"], "chain_name": chain, "index_chain": r["index"]}
           for r in residues]
    by_index = {r["index"]: nts[k]["nt_id"] for k, r in enumerate(residues)}

    dot = ["."] * len(sequence)
    for i, j in pairs:
        dot[i], dot[j] = "(", ")"

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "chains": {chain: {"bseq": sequence, "sstr": "".join(dot)}},
        "nts": nts,
        "pairs": [{"nt1": by_index[i + 1], "nt2": by_index[j + 1], "name": "WC"}
                  for i, j in pairs],
    }, indent=2))
    print(f"mock_dssr: wrote {out} ({len(pairs)} pairs)")


if __name__ == "__main__":
    main()
