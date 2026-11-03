#!/usr/bin/env python3
"""Make RhoFold+ outputs scoreable by eval/run_eval_array.sh (eval env, Biopython).

  python eval/rhofold/rhofold_to_eval_layout.py --src <rhofold_dir> --dest <dir> --tag <tag>

RhoFold+ writes <src>/<PDBID>/{relaxed_1000_model.pdb, unrelaxed_model.pdb, ...}. The chosen model is written as
mmCIF plus an _entity_poly loop (chains typed 'polyribonucleotide'; required by the prediction parser) into
the AF3-like layout   <dest>/<pdbid>_rhofold/seed-1_sample-0/model.cif
and the shard line for eval/run_eval_array.sh is printed. RhoFold+ runs once per input (no seeds); the paper
uses the relaxed model.
"""
import argparse, io as _io
from pathlib import Path
from Bio.PDB import MMCIFIO, PDBParser

ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True, type=Path)
ap.add_argument("--dest", required=True, type=Path)
ap.add_argument("--tag", required=True)
ap.add_argument("--model", default="relaxed_1000_model.pdb")
a = ap.parse_args()

n = 0
for pdb in sorted(a.src.glob(f"*/{a.model}")):
    pid = pdb.parent.name[:4].lower()
    s = PDBParser(QUIET=True).get_structure(pid, pdb)
    buf = _io.StringIO(); w = MMCIFIO(); w.set_structure(s); w.save(buf)
    chains = [c.id for c in next(iter(s))]
    loop = ["#", "loop_", "_entity_poly.entity_id", "_entity_poly.type", "_entity_poly.pdbx_strand_id"]
    loop += [f"{k + 1} polyribonucleotide {c}" for k, c in enumerate(chains)]
    out = a.dest / f"{pid}_rhofold" / "seed-1_sample-0" / "model.cif"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(buf.getvalue().rstrip("\n") + "\n" + "\n".join(loop) + "\n#\n")
    n += 1
print(f"# converted {n} models")
print(f"{a.tag}\t{a.dest}\t1\t0\tsingle\trna_rna")
