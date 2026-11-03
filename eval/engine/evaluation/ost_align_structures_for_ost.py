#!/usr/bin/env python3
"""
Align predicted (model) structure to reference so chain IDs and residue indices conform
to OpenStructure expectations. Run for rows where OST_fail = 1 so that ost compare-structures
can compute iLDDT, QS, DockQ, ICS without "Could not map model chain to chem group".

Uses chain mapping from sequence alignment (ref_chain -> model_chain) then:
- Renames model chains to match reference chain IDs.
- Renumbers residue indices in the model to 1..N per chain (matching reference length).

Output: CIF files (or PDB if gemmi unavailable) to parent/af_res/ (model) and parent/pdb_res/ (ref).
Same filenames as originals; never overwrites af/ or pdb/.

Requires: Bio.PDB (MMCIFParser, PDBIO). Optional: gemmi for CIF output.
Usage:
  python align_structures_for_ost.py [--dry-run]
"""

import os
import sys
import argparse
import logging

from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import Array

from database.startConfig import StartConfig

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Standard RNA residue names
RNA_NAMES = {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}

# Standard amino acid 3-letter -> 1-letter (avoids Bio.PDB.Polypeptide.three_to_one which was removed in newer Biopython)
AA3_TO_1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H",
    "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q",
    "ARG": "R", "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def get_sequence_from_chain_bio(chain):
    """Return one-letter sequence and list of (residue, res_id) for a Bio.PDB chain."""
    seq = []
    res_list = []
    for res in chain:
        if is_aa(res):
            name = res.get_resname().strip().upper()
            seq.append(AA3_TO_1.get(name, "X"))
            res_list.append((res, res.id))
        else:
            name = res.get_resname().strip().upper()
            if name in RNA_NAMES or len(name) == 1:
                seq.append(name[0] if name in RNA_NAMES else "N")
            else:
                seq.append("N")
            res_list.append((res, res.id))
    return "".join(seq), res_list


def get_chain_sequences_bio(structure):
    """Return dict chain_id -> (sequence, list of (residue, res_id))."""
    out = {}
    for model in structure:
        for chain in model:
            seq, res_list = get_sequence_from_chain_bio(chain)
            if seq:
                out[chain.id] = (seq, res_list)
    return out


def align_sequences_and_match(ref_seqs, model_seqs):
    """
    Return dict ref_chain_id -> model_chain_id by best sequence alignment.
    ref_seqs, model_seqs: dict chain_id -> (sequence_string, res_list).
    Uses Bio.Align.PairwiseAligner (global, match=1, mismatch=0, no gap penalty).
    """
    # Unique alphabet: 20 standard aa + X (unknown) + U (RNA uracil)
    letters = "ACDEFGHIKLMNPQRSTVWYXU"
    identity_data = {(a, b): (1.0 if a == b else 0.0) for a in letters for b in letters}
    matrix = None
    try:
        mat = Array(letters, dims=2)
        if hasattr(mat, "update"):
            mat.update(identity_data)
        else:
            for (a, b), v in identity_data.items():
                mat[(a, b)] = v
        matrix = mat
    except Exception:
        try:
            matrix = Array(data=identity_data)
        except Exception:
            pass

    aligner = PairwiseAligner(mode="global")
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    if matrix is not None:
        aligner.substitution_matrix = matrix
    # else: use default matrix (still gives best chain match by score)

    mapping = {}
    for ref_cid, (ref_seq, _) in ref_seqs.items():
        best_mdl = None
        best_score = -1
        for mdl_cid, (mdl_seq, _) in model_seqs.items():
            if mdl_cid in mapping.values():
                continue
            try:
                alignments = aligner.align(ref_seq, mdl_seq)
                if alignments:
                    score = alignments[0].score
                    if score > best_score:
                        best_score = score
                        best_mdl = mdl_cid
            except Exception:
                pass
        if best_mdl is not None:
            mapping[ref_cid] = best_mdl
    return mapping


def conform_model_bio_simple(structure_ref, structure_model, chain_mapping):
    """
    Rename model chains to ref IDs and set residue IDs to match reference exactly.
    Removes from the model any chain not in the reference so OST can map all chains.
    Output chain order matches the reference (OST maps by position).
    Modifies structure_model in place. Handles overlapping ref/model chain IDs via temp ids.
    """
    ref_seqs = get_chain_sequences_bio(structure_ref)
    model_seqs = get_chain_sequences_bio(structure_model)
    model = structure_model[0]
    ref_model_obj = structure_ref[0]
    ref_cids = set(chain_mapping.keys())

    # Reference chain order (as in file) so we can write model in same order for OST
    ref_chain_order = []
    for ch in ref_model_obj:
        if ch.id in ref_cids:
            ref_chain_order.append(ch.id)

    def do_rename(mid, rid):
        if mid not in model.child_dict:
            return
        if mid == rid:
            return  # no-op
        chain = model.child_dict[mid]
        ref_res_list = ref_seqs.get(rid, (None, []))[1]  # list of (res, res.id)
        res_list = list(chain.get_residues())
        # Assign reference residue ids without triggering Entity.id setter (avoids KeyError when
        # ref reuses residue numbers or we overwrite keys). Set _id directly then rebuild child_dict.
        n = min(len(res_list), len(ref_res_list))
        for i in range(n):
            res_list[i]._id = ref_res_list[i][1]
        # Important: if the model chain is longer than the reference chain, truncating avoids
        # residue-id collisions (same (resseq, icode) written twice) which can later make
        # mmCIF output invalid and crash parsing elsewhere (e.g. "Atom P defined twice...").
        res_list = res_list[:n]
        chain.child_dict.clear()
        chain.child_list = res_list
        for res in res_list:
            chain.child_dict[res._id] = res
        # Update model dict and chain id
        model.child_dict[rid] = model.child_dict.pop(mid)
        chain._id = rid

    # Phase 1: renames where model_cid is not a ref_cid (no cycle)
    for ref_cid, model_cid in chain_mapping.items():
        if model_cid not in ref_cids:
            do_rename(model_cid, ref_cid)

    # Phase 2: cycles (model_cid is also a ref_cid): use temp id to avoid overwriting
    remaining = [(r, m) for r, m in chain_mapping.items() if m in ref_cids]
    for ref_cid, model_cid in remaining:
        tmp_id = "_tmp_" + ref_cid
        do_rename(model_cid, tmp_id)
    for ref_cid, model_cid in remaining:
        tmp_id = "_tmp_" + ref_cid
        do_rename(tmp_id, ref_cid)

    # Do not remove model chains: keep all so ref and pred with same number of chains (e.g. 2 RNA) are preserved.
    # Order: matched chains in ref order first, then any unmatched model chains (so they are still written).
    matched_in_order = [model[c] for c in ref_chain_order if c in model.child_dict]
    unmatched = [model[c] for c in list(model.child_dict.keys()) if c not in ref_cids]
    model.child_list = matched_in_order + unmatched
    return structure_model


def save_structure_bio(structure, path):
    """Save with Bio.PDB PDBIO (PDB format)."""
    io = PDBIO()
    io.set_structure(structure)
    io.save(path)


def save_structure_cif(structure_bio, path):
    """Save as CIF using gemmi (reads from a temp PDB written by Bio.PDB). Returns True if CIF written, False otherwise."""
    try:
        import gemmi
        import tempfile
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            tmp_pdb = f.name
        try:
            save_structure_bio(structure_bio, tmp_pdb)
            st = gemmi.read_structure(tmp_pdb)
            doc = st.make_mmcif_document()
            doc.write_file(path)
            normalize_cif_for_ost(path)
            return True
        finally:
            try:
                os.remove(tmp_pdb)
            except OSError:
                pass
    except Exception:
        return False


def _count_atoms(structure_bio) -> int:
    """Count atoms in a Bio.PDB structure object (used to avoid writing empty outputs)."""
    n = 0
    for model in structure_bio:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    n += 1
    return n


def normalize_cif_for_ost(cif_path):
    """
    Ensure CIF file is valid for OpenStructure: strip leading blank/comment lines
    so the file starts with a 'data_' block. OST and our pipeline (clean_cif_for_openstructure)
    can choke on a leading comment or BOM; writing a canonical mmCIF start avoids that.
    """
    try:
        with open(cif_path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
        start = 0
        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.startswith("data_"):
                start = i
                break
            start = i
            break
        if start == 0:
            return
        with open(cif_path, "w", encoding="utf-8") as f:
            f.writelines(lines[start:])
    except Exception as e:
        logger.debug("Could not normalize CIF for OST %s: %s", cif_path, e)


def copy_ref_to_pdb_res(ref_path, pdb_res_path):
    """Copy reference CIF to pdb_res/ (same filename)."""
    import shutil
    os.makedirs(os.path.dirname(pdb_res_path), exist_ok=True)
    shutil.copy2(ref_path, pdb_res_path)


def process_one_pair(ref_path, model_path, pdb_res_path, af_res_path, cif_parser, dry_run=False):
    """
    Load ref and model, align chains, conform model to ref, write to pdb_res_path and af_res_path.
    Returns True on success, False on skip/failure.
    """
    ref_filename = os.path.basename(ref_path)
    model_filename = os.path.basename(model_path)
    if not os.path.exists(ref_path) or not os.path.exists(model_path):
        logger.warning("Skip: ref or model missing (%s, %s)", ref_path, model_path)
        return False
    if dry_run:
        logger.info("Would process: %s vs %s", model_filename, ref_filename)
        return True
    try:
        ref_struct = cif_parser.get_structure("ref", ref_path)
        model_struct = cif_parser.get_structure("model", model_path)
    except Exception as e:
        logger.warning("Parse failed for %s / %s: %s", ref_filename, model_filename, e)
        return False
    ref_seqs = get_chain_sequences_bio(ref_struct)
    model_seqs = get_chain_sequences_bio(model_struct)
    if not ref_seqs or not model_seqs:
        logger.warning("No chains in ref or model")
        return False
    chain_mapping = align_sequences_and_match(ref_seqs, model_seqs)
    if not chain_mapping:
        logger.warning("No chain mapping")
        return False
    conform_model_bio_simple(ref_struct, model_struct, chain_mapping)

    # Safety: never overwrite inputs (this has caused partial/empty CIF corruption in your run).
    if os.path.abspath(pdb_res_path) == os.path.abspath(ref_path):
        logger.error("Refusing to overwrite reference CIF: %s", ref_path)
        return False
    if os.path.abspath(af_res_path) == os.path.abspath(model_path):
        logger.error("Refusing to overwrite model CIF: %s", model_path)
        return False

    # Safety: avoid writing empty conformations (can yield CIF without atom_site loops).
    if _count_atoms(model_struct) == 0:
        logger.error("Conformed model has 0 atoms; skipping write for model=%s", model_path)
        return False

    try:
        copy_ref_to_pdb_res(ref_path, pdb_res_path)
    except Exception as e:
        logger.error("Copy ref to pdb_res failed: %s", e)
        return False
    written_cif = save_structure_cif(model_struct, af_res_path)
    if not written_cif:
        af_res_out = os.path.join(os.path.dirname(af_res_path), os.path.splitext(model_filename)[0] + ".pdb")
        save_structure_bio(model_struct, af_res_out)
        logger.warning("Wrote PDB (CIF failed or gemmi unavailable): %s", af_res_out)
    else:
        af_res_out = af_res_path
        logger.info("Wrote CIF: %s", af_res_out)
    logger.info("Conformed %s -> %s (ref -> %s), chain_mapping=%s", model_filename, af_res_out, pdb_res_path, chain_mapping)
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Align pred structures to ref for OpenStructure (chain mapping / OST_fail)."
    )
    parser.add_argument("--dry-run", action="store_true", help="Only list files that would be processed, do not write")
    parser.add_argument(
        "--af",
        dest="af_folder",
        metavar="DIR",
        default=None,
        help="AF/prediction folder. If set, --pdb is required; process all pairs found in these folders (no DB). Output written to <af>/af_res/ and <pdb>/pdb_res/.",
    )
    parser.add_argument(
        "--pdb",
        dest="pdb_folder",
        metavar="DIR",
        default=None,
        help="Reference (PDB) folder. Required when using --af.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="When using DB mode (no --af/--pdb), process all pred_protein_rna rows with FileName not NULL, not just those with OST_fail != 0.",
    )
    parser.add_argument(
        "--exp_db_id",
        dest="exp_db_id",
        default=None,
        help="When using DB mode, only process the specified exp_db_id (e.g. 7UJ1).",
    )
    args = parser.parse_args()

    cif_parser = MMCIFParser(QUIET=True)

    if args.af_folder is not None or args.pdb_folder is not None:
        if not args.af_folder or not args.pdb_folder:
            logger.error("When using folder mode, both --af and --pdb are required.")
            sys.exit(1)
        af_folder = os.path.abspath(args.af_folder)
        pdb_folder = os.path.abspath(args.pdb_folder)
        if not os.path.isdir(af_folder):
            logger.error("AF folder not found: %s", af_folder)
            sys.exit(1)
        if not os.path.isdir(pdb_folder):
            logger.error("PDB folder not found: %s", pdb_folder)
            sys.exit(1)
        af_res_folder = os.path.join(af_folder, "af_res")
        pdb_res_folder = os.path.join(pdb_folder, "pdb_res")
        os.makedirs(af_res_folder, exist_ok=True)
        os.makedirs(pdb_res_folder, exist_ok=True)
        # For each .cif in af folder, derive pdb_id (e.g. fold_8cbw_s1_model_0.cif -> 8cbw) and pair with pdb_folder/8cbw.cif
        processed = 0
        for name in sorted(os.listdir(af_folder)):
            if not name.endswith(".cif"):
                continue
            model_path = os.path.join(af_folder, name)
            if not os.path.isfile(model_path):
                continue
            pdb_id = name.split("_")[1] if "_" in name else name[:4]
            ref_path = os.path.join(pdb_folder, f"{pdb_id}.cif")
            if not os.path.isfile(ref_path):
                logger.debug("No ref %s for %s", ref_path, name)
                continue
            pdb_res_path = os.path.join(pdb_res_folder, f"{pdb_id}.cif")
            af_res_path = os.path.join(af_res_folder, name)
            if process_one_pair(ref_path, model_path, pdb_res_path, af_res_path, cif_parser, dry_run=args.dry_run):
                processed += 1
        logger.info("Done. Processed %d pair(s). Output in %s and %s", processed, af_res_folder, pdb_res_folder)
        sys.exit(0)

    config = StartConfig()
    parent = config.get_parent_folder()
    if not parent:
        logger.error("parent_folder not set in config")
        sys.exit(1)
    pdb_folder = os.path.join(parent, "pdb")
    af_folder = os.path.join(parent, "af")
    pdb_res_folder = os.path.join(parent, "pdb_res")
    af_res_folder = os.path.join(parent, "af_res")
    os.makedirs(pdb_res_folder, exist_ok=True)
    os.makedirs(af_res_folder, exist_ok=True)

    import sqlite3
    db_path = config.get_database_path()
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    try:
        if args.all:
            cursor.execute(
                """SELECT FileName, exp_db_id FROM pred_protein_rna
                   WHERE FileName IS NOT NULL"""
            )
        else:
            cursor.execute(
                """SELECT FileName, exp_db_id FROM pred_protein_rna
                   WHERE FileName IS NOT NULL
                     AND (OST_fail = 1 OR (OST_fail IS NOT NULL AND OST_fail != '0'))"""
            )
        rows = cursor.fetchall()
    except sqlite3.OperationalError as e:
        conn.close()
        if "OST_fail" in str(e) or "no such column" in str(e).lower():
            logger.error("Column OST_fail not found. Run extractAllAtomLDDT pred_protein_rna first to add and fill OST_fail.")
        else:
            logger.error("Database error: %s", e)
        sys.exit(1)
    conn.close()

    if not rows:
        logger.info("No rows with OST_fail (failure). Ensure extractAllAtomLDDT has been run to add and fill OST_fail.")
        return

    logger.info("Found %d row(s) with OST_fail (chain mapping failure)", len(rows))

    for row in rows:
        file_name = row["FileName"]
        exp_db_id = row["exp_db_id"]
        if args.exp_db_id and str(exp_db_id).upper() != str(args.exp_db_id).upper():
            continue
        pdb_id = exp_db_id if exp_db_id else (file_name.split("_")[1] if "_" in file_name else file_name[:4])
        ref_path = os.path.join(pdb_folder, f"{pdb_id}.cif")
        model_path = os.path.join(af_folder, file_name)
        pdb_res_path = os.path.join(pdb_res_folder, os.path.basename(ref_path))
        af_res_path = os.path.join(af_res_folder, os.path.basename(model_path))
        process_one_pair(ref_path, model_path, pdb_res_path, af_res_path, cif_parser, dry_run=args.dry_run)

    logger.info("Done. Use af_res/ and pdb_res/ for compare-structures (same filenames).")


if __name__ == "__main__":
    main()
