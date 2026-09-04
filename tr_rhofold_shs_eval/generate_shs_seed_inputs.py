#!/usr/bin/env python3
"""Regenerate one SHS variant's RhoFold inputs (FASTA + A3M) under a chosen SHS seed.

Joins the variant's base-pair set from RNA_Monomers_base_pairs.csv with the target
sequences, runs the SyntheticEvolution SHS generator per target, and writes the same
fasta/ + a3m/ + metadata/ layout that prepare_shs_inputs_from_af3_json.py produces.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

HERE = Path(__file__).resolve().parent
# This folder lives inside the SyntheticEvolution checkout, so the SHS generator
# is one level up, at the repository root.
REPO_ROOT = HERE.parent
DEFAULT_GENERATOR = (REPO_ROOT / "rna_msa_generator_base_pair.py").resolve()
# The generator imports RnaBench, which only loads cleanly in `synEvo`
# (the environment created by install.sh / env/environment.yml).
DEFAULT_GEN_PYTHON = str(Path.home() / ".conda" / "envs" / "synEvo" / "bin" / "python")

# Which RNA_Monomers_base_pairs.csv method drives each SHS variant, and where the
# variant's target sequences come from (the existing seed-42 prepared folder).
# `*_denoise` variants are the second-round arms: their structure source is the DSSR
# annotation of a structure folded from the first-round N100 SHS MSA.
VARIANTS: Dict[str, Dict[str, str]] = {
    "rnafold": {
        "method": "rnafold_pred",
        "source": "prepared_shs_inputs_from_rnafold_folder",
    },
    "spotrna": {
        "method": "spotrna_pred",
        "source": "prepared_shs_inputs_from_spotrna_folder",
    },
    "rnaformerN100": {
        "method": "rnaformer_pred",
        "source": "prepared_shs_inputs_from_rnaformerN100_folder",
    },
    "rna_monomers_dssrN100": {
        "method": "dssrN100_dssr",
        "source": "prepared_shs_inputs_from_rna_monomers_dssrN100_folder",
    },
    "rna_monomers_rnafold_denoise": {
        "method": "rnafoldN100_dssr",
        "source": "prepared_shs_inputs_from_rna_monomers_rnafold_denoise_folder",
    },
    "rna_monomers_rnaformerN100_denoise": {
        "method": "rnaformerN100_dssr",
        "source": "prepared_shs_inputs_from_rna_monomers_rnaformerN100_denoise_folder",
    },
    "rna_monomers_spotrna_denoise": {
        "method": "spotrnaN100_dssr",
        "source": "prepared_shs_inputs_from_rna_monomers_spotrna_denoise_folder",
    },
}


def prep_dir_for(variant: str, seed: int) -> str:
    return f"prepared_shs_inputs_from_{variant}_seed{seed}_folder"


def pred_dir_for(variant: str, seed: int) -> str:
    return f"predictions_rhofold_on_shs_from_{variant}_seed{seed}_folder"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate seed-specific SHS FASTA/A3M inputs for one variant.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--variant", choices=sorted(VARIANTS), help="SHS variant to regenerate.")
    parser.add_argument("--seed", type=int, help="SHS-generator seed.")
    parser.add_argument("--n-msa", type=int, default=100, help="Synthetic sequences per MSA.")
    parser.add_argument(
        "--pairs-csv",
        default=str(HERE / "RNA_Monomers_base_pairs.csv"),
        help="CSV with pdb_id, method and pair-list columns.",
    )
    parser.add_argument(
        "--pairs-column",
        default="pairs",
        choices=["pairs", "wc_pairs", "wobble_pairs", "nc_pairs", "lone_pairs", "multiplets"],
        help="Which pair column feeds the generator.",
    )
    parser.add_argument("--method", default=None, help="Override the variant's CSV method.")
    parser.add_argument("--fasta-dir", default=None, help="Override the sequence source dir (a fasta/ folder).")
    parser.add_argument("--out-dir", default=None, help="Override the output dir.")
    parser.add_argument("--generator", default=str(DEFAULT_GENERATOR), help="SHS generator script.")
    parser.add_argument("--gen-python", default=DEFAULT_GEN_PYTHON, help="Interpreter for the generator (needs RnaBench).")
    parser.add_argument(
        "--shs-params",
        default=None,
        help='JSON dict of extra generator knobs, e.g. \'{"wobble-prob": 0.2}\'. '
        "Empty => generator defaults, which is what the seed-42 folders used.",
    )
    parser.add_argument("--ids-file", default=None, help="Optional file with one PDB ID per line.")
    parser.add_argument("--skip-existing", action="store_true", help="Leave targets that already have an A3M alone.")
    parser.add_argument("--keep-raw-json", action="store_true", help="Keep the generator's raw output dirs.")
    parser.add_argument(
        "--print-paths",
        action="store_true",
        help="Print this variant/seed's default dirs as shell assignments and exit (no generation).",
    )
    parser.add_argument("--list-variants", action="store_true", help="Print the variant -> CSV method table and exit.")
    args = parser.parse_args()

    if args.list_variants:
        return args
    if not args.variant and not (args.method and args.fasta_dir and args.out_dir):
        parser.error("--variant is required (or give --method, --fasta-dir and --out-dir explicitly)")
    if args.seed is None:
        parser.error("--seed is required")

    args.shs_params = json.loads(args.shs_params) if args.shs_params else {}
    if not isinstance(args.shs_params, dict):
        parser.error("--shs-params must be a JSON object (dict)")
    return args


def shs_params_to_flags(params: Dict[str, object]) -> List[str]:
    """Expand {knob: value} into generator CLI flags (mirrors run_shs_af3_dssr_eval.py).

    MSA depth is owned by --n-msa, so any n/N key here is dropped to avoid a duplicate -N."""
    flags: List[str] = []
    for key, val in params.items():
        flag = "--" + str(key).strip().lstrip("-").replace("_", "-")
        if flag.lower() in ("--n", "--n-msa"):
            continue
        if val is None or val is False:
            continue
        if val is True:
            flags.append(flag)
        else:
            flags += [flag, str(val)]
    return flags


def parse_pairs(cell: str) -> List[Tuple[int, int]]:
    """Parse a CSV pair cell. The column holds a Python literal list of tuples
    (`[(1, 25), ...]`), 0-based, so json.loads would choke on it."""
    text = (cell or "").strip()
    if not text or text.lower() in ("nan", "none", "[]"):
        return []
    parsed = ast.literal_eval(text)
    pairs = {(int(min(a, b)), int(max(a, b))) for a, b in parsed if int(a) != int(b)}
    return sorted(pairs)


def read_fasta_sequences(fasta_dir: Path) -> Dict[str, str]:
    """Map PDB ID -> sequence from a prepared fasta/ folder (one record per file)."""
    sequences: Dict[str, str] = {}
    for path in sorted(list(fasta_dir.glob("*.fasta")) + list(fasta_dir.glob("*.fa"))):
        seq = "".join(ln.strip() for ln in path.read_text().splitlines() if not ln.startswith(">"))
        seq = seq.upper().replace("T", "U")
        if seq:
            sequences[path.stem.upper()] = seq
    return sequences


def read_id_filter(ids_file: str | None) -> set[str] | None:
    if not ids_file:
        return None
    return {ln.strip().upper() for ln in Path(ids_file).read_text().splitlines() if ln.strip()}


def read_pairs_by_id(pairs_csv: Path, method: str, column: str) -> Dict[str, str]:
    """Raw pair cells for one method, keyed by PDB ID."""
    with open(pairs_csv, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or column not in rows[0]:
        raise SystemExit(f"{pairs_csv} has no '{column}' column")
    methods = {row["method"] for row in rows}
    if method not in methods:
        raise SystemExit(f"method '{method}' not in {pairs_csv} (have: {sorted(methods)})")
    return {row["pdb_id"].upper(): row[column] for row in rows if row["method"] == method}


def run_generator(args: argparse.Namespace, pdb_id: str, sequence: str,
                  pairs: List[Tuple[int, int]], work_dir: Path) -> Path:
    """Build the minimal AF3 input JSON and run the SHS generator on it.
    Returns the AF3 job JSON the generator wrote (unpairedMsa = synthetic homologs)."""
    af3_input = {
        "name": pdb_id,
        "modelSeeds": [1],
        "sequences": [{"rna": {"sequence": sequence, "modifications": [], "id": "A"}}],
        "dialect": "alphafold3",
        "version": 1,
    }
    input_path = work_dir / "af3_input.json"
    input_path.write_text(json.dumps(af3_input, indent=2))

    out_dir = work_dir / "shs_json"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        args.gen_python, str(args.generator),
        "--input_json_path", str(input_path),
        "--pairs", json.dumps([[a, b] for a, b in pairs]),
        "--output_json_dir", str(out_dir),
        "--pdb_id", pdb_id,
        "-N", str(args.n_msa),
        "--seed", str(args.seed),
    ] + shs_params_to_flags(args.shs_params)
    subprocess.run(cmd, check=True, cwd=str(Path(args.generator).parent))

    produced = sorted(out_dir.glob("*.json"), key=os.path.getmtime)
    if not produced:
        raise RuntimeError(f"[{pdb_id}] SHS generator produced no JSON in {out_dir}")
    return produced[-1]


def extract_fasta_a3m(job_json: Path, pdb_id: str) -> Tuple[str, str]:
    """Read sequence + synthetic MSA back out of the generator's AF3 job JSON
    (same read-out as prepare_shs_inputs_from_af3_json.py)."""
    data = json.loads(job_json.read_text())
    rna = next((item["rna"] for item in data.get("sequences", []) if "rna" in item), None)
    if rna is None:
        raise RuntimeError(f"{job_json} has no RNA sequence entry")
    seq = (rna.get("sequence") or "").strip().upper().replace("T", "U")
    if not seq:
        raise RuntimeError(f"{job_json} has an empty RNA sequence")
    msa = (rna.get("unpairedMsa") or "").strip()
    if not msa:
        msa = f">{pdb_id}\n{seq}"
    return f">{pdb_id}\n{seq}\n", msa + "\n"


def main() -> None:
    args = parse_args()

    if args.list_variants:
        print(f"{'variant':38s} {'csv method':22s} sequence source")
        for name, cfg in sorted(VARIANTS.items()):
            print(f"{name:38s} {cfg['method']:22s} {cfg['source']}")
        return

    cfg = VARIANTS.get(args.variant, {})
    method = args.method or cfg["method"]
    fasta_dir = Path(args.fasta_dir or (HERE / cfg["source"] / "fasta")).resolve()
    out_dir = Path(args.out_dir or (HERE / prep_dir_for(args.variant, args.seed))).resolve()

    if args.print_paths:
        print(f"DEFAULT_PREP_DIR={out_dir}")
        print(f"DEFAULT_OUT_DIR={HERE / pred_dir_for(args.variant, args.seed)}")
        print(f"METHOD={method}")
        return

    generator = Path(args.generator).resolve()
    if not generator.is_file():
        raise SystemExit(f"SHS generator not found: {generator} (pass --generator)")
    if not fasta_dir.is_dir():
        raise SystemExit(f"Sequence source not found: {fasta_dir} (pass --fasta-dir)")
    args.generator = str(generator)

    sequences = read_fasta_sequences(fasta_dir)
    if not sequences:
        raise SystemExit(f"No FASTA records under {fasta_dir}")
    pairs_by_id = read_pairs_by_id(Path(args.pairs_csv).resolve(), method, args.pairs_column)
    id_filter = read_id_filter(args.ids_file)

    fasta_out = out_dir / "fasta"
    a3m_out = out_dir / "a3m"
    json_out = out_dir / "af3_json"
    meta_out = out_dir / "metadata"
    work_root = out_dir / "work"
    for d in (fasta_out, a3m_out, json_out, meta_out, work_root):
        d.mkdir(parents=True, exist_ok=True)

    print(f"Variant     : {args.variant}  (method={method}, column={args.pairs_column})")
    print(f"Seed / N    : {args.seed} / {args.n_msa}")
    print(f"Sequences   : {fasta_dir} ({len(sequences)} targets)")
    print(f"Output      : {out_dir}", flush=True)

    rows: List[Dict[str, str]] = []
    exported_ids: List[str] = []
    exported = skipped = failed = 0

    for pdb_id in sorted(sequences):
        if id_filter is not None and pdb_id not in id_filter:
            continue
        sequence = sequences[pdb_id]
        status, reason = "skipped", ""
        seq_len, msa_nseq, n_pairs, source_json = str(len(sequence)), "0", "0", ""

        try:
            if args.skip_existing and (a3m_out / f"{pdb_id}.a3m").is_file():
                # Resumed run: the A3M is already there, so keep it and record it as
                # present, otherwise the rewritten manifest would under-report the folder.
                existing = (a3m_out / f"{pdb_id}.a3m").read_text()
                status, reason = "exported", "reused_existing"
                msa_nseq = str(sum(1 for ln in existing.splitlines() if ln.startswith(">")))
                exported_ids.append(pdb_id)
                skipped += 1
            elif pdb_id not in pairs_by_id:
                # e.g. SPOT-RNA has no prediction for 8TKK, DSSR none for 7EFG/7VFT.
                reason = f"no_{method}_row_in_pairs_csv"
                skipped += 1
            else:
                pairs = parse_pairs(pairs_by_id[pdb_id])
                out_of_range = [p for p in pairs if p[1] >= len(sequence)]
                if out_of_range:
                    raise ValueError(
                        f"pair index beyond sequence length {len(sequence)}: {out_of_range[:3]}")

                work_dir = work_root / pdb_id
                work_dir.mkdir(parents=True, exist_ok=True)
                produced = run_generator(args, pdb_id, sequence, pairs, work_dir)
                fasta_text, a3m_text = extract_fasta_a3m(produced, pdb_id)

                shutil.copy(produced, json_out / f"{pdb_id}.json")
                (fasta_out / f"{pdb_id}.fasta").write_text(fasta_text)
                (a3m_out / f"{pdb_id}.a3m").write_text(a3m_text)
                if not args.keep_raw_json:
                    shutil.rmtree(work_dir, ignore_errors=True)

                status, reason = "exported", f"n_pairs={len(pairs)}"
                n_pairs = str(len(pairs))
                msa_nseq = str(sum(1 for ln in a3m_text.splitlines() if ln.startswith(">")))
                source_json = produced.name
                exported += 1
                exported_ids.append(pdb_id)
        except Exception as exc:  # noqa: BLE001
            status, reason = "failed", f"error:{exc}"
            failed += 1
            print(f"[{pdb_id}] FAILED: {exc}", file=sys.stderr, flush=True)

        print(f"[{pdb_id}] {status} ({reason})", flush=True)
        rows.append({
            "pdb_id": pdb_id,
            "status": status,
            "reason": reason,
            "seq_len": seq_len,
            "msa_nseq": msa_nseq,
            "n_pairs": n_pairs,
            "generator_json": source_json,
        })

    with open(meta_out / "manifest.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["pdb_id", "status", "reason", "seq_len", "msa_nseq", "n_pairs", "generator_json"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    (meta_out / "exported_ids.txt").write_text("".join(f"{i}\n" for i in exported_ids))
    (meta_out / "run_info.json").write_text(json.dumps({
        "variant": args.variant,
        "method": method,
        "pairs_column": args.pairs_column,
        "pairs_csv": str(Path(args.pairs_csv).resolve()),
        "seed": args.seed,
        "n_msa": args.n_msa,
        "shs_params": args.shs_params,
        "generator": args.generator,
        "sequence_source": str(fasta_dir),
    }, indent=2))
    if not args.keep_raw_json:
        shutil.rmtree(work_root, ignore_errors=True)

    print(f"Exported {exported}, skipped {skipped}, failed {failed} -> {out_dir}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
