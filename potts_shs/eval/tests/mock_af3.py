"""Stand-in for AlphaFold 3 that writes a minimal mmCIF and confidence file, so the pipeline's
paths, resumability and parsing can be tested where AF3 itself is not installed."""

import argparse
import json
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json_path", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--drop-rate", type=float, default=0.15,
                        help="Fraction of seeded pairs the mock 'fails' to reproduce")
    args, _ = parser.parse_known_args()

    payload = json.loads(Path(args.json_path).read_text())
    name = payload["name"]
    seed = (payload.get("modelSeeds") or [1])[0]

    # Every RNA chain, not just the first: a complex can hold two different RNAs, and which one
    # the pipeline scores is exactly what the chain-selection logic has to get right.
    chains = []
    for entity in payload["sequences"]:
        if "rna" not in entity:
            continue
        ids = entity["rna"].get("id", "A")
        for chain_id in (ids if isinstance(ids, list) else [ids]):
            chains.append((chain_id, entity["rna"]["sequence"].upper()))

    out_dir = Path(args.output_dir) / name.lower()
    out_dir.mkdir(parents=True, exist_ok=True)

    # One C1' atom per residue is enough to carry chain, residue index and residue name.
    lines = [f"data_{name}", "#", "loop_",
             "_atom_site.group_PDB", "_atom_site.id", "_atom_site.label_atom_id",
             "_atom_site.label_comp_id", "_atom_site.label_asym_id", "_atom_site.label_seq_id",
             "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
    rng = random.Random(hash((name, seed)) & 0xFFFF)
    atom = 0
    for chain_id, sequence in chains:
        for index, base in enumerate(sequence, start=1):
            atom += 1
            x, y, z = (round(rng.uniform(-30, 30), 3) for _ in range(3))
            lines.append(f"ATOM {atom} C1' {base} {chain_id} {index} {x} {y} {z}")
    (out_dir / f"{name}_model.cif").write_text("\n".join(lines) + "\n")

    (out_dir / f"{name}_summary_confidences.json").write_text(json.dumps({
        "ptm": round(rng.uniform(0.4, 0.9), 3),
        "iptm": round(rng.uniform(0.3, 0.8), 3),
        "ranking_score": round(rng.uniform(0.3, 0.95), 3),
        "has_clash": 0.0,
    }, indent=2))

    # The mock cannot infer structure, so it carries the seeded pairs alongside the model and
    # lets mock_dssr drop a fraction of them. This exercises the scoring path with non-trivial
    # numbers; it simulates nothing about folding.
    (out_dir / "_mock_truth.json").write_text(json.dumps({
        "chains": [c for c, _ in chains], "drop_rate": args.drop_rate, "seed": seed,
    }))
    print(f"mock_af3: wrote {out_dir}")


if __name__ == "__main__":
    main()
