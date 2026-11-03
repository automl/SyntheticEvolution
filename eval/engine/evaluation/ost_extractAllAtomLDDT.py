#!/usr/bin/env python3
"""
1. conda activate ost (if not present, check with conda info --envs)
2. Extract All-Atom LDDT Scores using OpenStructure

   python extractAllAtomLDDT.py pred_protein_rna [--residue-number-alignment] [--enforceMapping]

   Defaults align with SWISS-MODEL Structure Assessment (sequence chain mapping; no
   --residue-number-alignment).    Optional env (same as ost_rbp_interfaces.py): SWISSMODEL_OST_RESNUM_ALIGN=1 forces
   --residue-number-alignment on all compare-structures calls.

3. conda activate RNABench
"""

import argparse
import ast
import os
import sys
import sqlite3
import tempfile
import traceback
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.Polypeptide import is_aa
from database.startConfig import StartConfig
from database.databaseMethods import DatabaseMethods
import logging
import shlex
from collections import defaultdict
from typing import Dict, Optional, Set
from ost import mol, io
from ost.mol.alg import lddt

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

# Initialize database connection
config = StartConfig()
db = DatabaseMethods()

def ensure_rna_lddt_column(table_name):
    """
    Ensure the RNA all-atom LDDT column exists in the specified table
    
    Args:
        table_name (str): Name of the table to check/update
    """
    try:
        cursor = db.connection.cursor()
        
        # Check if RNA_LDDT_allAtom column exists and add it if it doesn't
        column_name = 'RNA_LDDT_allAtom'
        try:
            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} REAL")
            logger.info(f"Added column {column_name} to {table_name}")
        except sqlite3.OperationalError as e:
            if "duplicate column name" in str(e):
                logger.info(f"Column {column_name} already exists in {table_name}")
            else:
                raise e
        
        db.connection.commit()
        
    except Exception as e:
        logger.error(f"Error ensuring RNA LDDT column in {table_name}: {str(e)}")
        raise

def ensure_protein_complex_columns(table_name):
    """Ensure Protein_LDDT_allAtom, Complex_LDDT_allAtom, iLDDT, QS_global, Dockq and ICS columns exist.
    (OST_fail is owned by ost_calcRNA_AllChains.py and is not created/written here.)"""
    try:
        cursor = db.connection.cursor()
        for column_name in ['Protein_LDDT_allAtom', 'Complex_LDDT_allAtom', 'iLDDT', 'QS_global', 'Dockq', 'ICS']:
            try:
                cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} REAL")
                logger.info(f"Added column {column_name} to {table_name}")
            except sqlite3.OperationalError as e:
                if "duplicate column name" in str(e):
                    logger.debug(f"Column {column_name} already exists in {table_name}")
                else:
                    raise e
        db.connection.commit()
    except Exception as e:
        logger.error(f"Error ensuring Protein/Complex columns in {table_name}: {str(e)}")
        raise


def get_structure_files(file_name_cif_pdb, file_name_cif_af, pdb_folder, af_folder):
    """
    Get reference and predicted structure file paths

    Args:
        file_name_cif_pdb (str): Reference structure filename
        file_name_cif_af (str): Predicted structure filename
        pdb_folder (str): Reference structures folder
        af_folder (str): Predicted structures folder

    Returns:
        tuple: (ref_path, model_path) or (None, None) if files don't exist
    """
    ref_path = os.path.join(pdb_folder, file_name_cif_pdb)
    model_path = os.path.join(af_folder, file_name_cif_af)

    if not os.path.exists(ref_path):
        logger.warning(f"Reference file not found: {ref_path}")
        return None, None

    if not os.path.exists(model_path):
        logger.warning(f"Model file not found: {model_path}")
        return None, None

    return ref_path, model_path


def convert_cif_to_pdb(cif_path, output_path):
    """
    Convert CIF file to PDB format

    Args:
        cif_path (str): Path to CIF file
        output_path (str): Path for output PDB file
    """
    try:
        cif_parser = MMCIFParser(QUIET=True)
        structure = cif_parser.get_structure('structure', cif_path)

        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(output_path)

    except Exception as e:
        logger.error(f"Error converting CIF to PDB: {str(e)}")
        raise

def _parse_pdb_residues(path: str) -> Dict[str, Set[str]]:
    """
    Parse a PDB file and return a map chain_id -> set of residue identifiers.
    Residue identifier is the resseq with insertion code appended if present (e.g. "123A").
    """
    chains = defaultdict(set)
    with open(path, 'r') as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            # chain id is column 22 (index 21), resseq is cols 23-26 (22:26), iCode at col27 (index 26)
            chain = line[21].strip() or " "  # empty chain ID -> use space
            resseq_str = line[22:26].strip()
            icode = line[26].strip() if len(line) > 26 else ""
            # build residue id preserving insertion code if present
            res_id = f"{resseq_str}{icode}" if icode else resseq_str
            if res_id:
                chains[chain].add(res_id)
    return chains


def _parse_cif_residues(path: str) -> Dict[str, Set[str]]:
    """
    Basic mmCIF atom_site loop parser to extract chain and residue number.
    It looks for common tags: auth_asym_id/auth_seq_id then falls back to label_asym_id/label_seq_id.
    Returns a map chain_id -> set of residue identifiers (strings).
    Uses shlex.split to respect quoted fields.
    """
    chains = defaultdict(set)
    with open(path, 'r') as fh:
        lines = fh.readlines()

    in_atom_site = False
    atom_site_tags = []
    data_started = False
    # indices we'll try to find
    chain_tag_idx = None
    seq_tag_idx = None

    i = 0
    while i < len(lines):
        line = lines[i].rstrip('\n')
        stripped = line.strip()
        if stripped.startswith("loop_"):
            # start of a loop - check if following tags are atom_site tags
            j = i + 1
            tags = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                tags.append(lines[j].strip())
                j += 1
            # if any atom_site tags present, enter atom_site parsing
            if any(t.startswith("_atom_site.") for t in tags):
                # collect the atom_site tags order
                atom_site_tags = [t for t in tags if t.startswith("_atom_site.")]
                in_atom_site = True
                data_started = False
                # find indices for chain/seq
                # prefer auth_asym_id/auth_seq_id, else label_asym_id/label_seq_id
                def find_tag_index(pref):
                    for idx, t in enumerate(atom_site_tags):
                        if t.endswith(pref):
                            return idx
                    return None

                chain_tag_idx = find_tag_index("auth_asym_id") or find_tag_index("label_asym_id") \
                                or find_tag_index("auth_asym") or find_tag_index("label_asym")
                seq_tag_idx = find_tag_index("auth_seq_id") or find_tag_index("label_seq_id") \
                              or find_tag_index("auth_seq") or find_tag_index("label_seq")
                # move i to the first data line after tags
                i = j
                # now parse data rows until blank line or next loop_ or line starting with "_"
                while i < len(lines):
                    dat_line = lines[i].strip()
                    if dat_line == "" or dat_line.startswith("loop_") or dat_line.startswith("data_"):
                        break
                    if dat_line.startswith("_"):
                        # another tag found -> stop data section
                        break
                    # parse data tokens (respect quotes)
                    try:
                        tokens = shlex.split(lines[i])
                    except ValueError:
                        tokens = lines[i].split()
                    if chain_tag_idx is not None and seq_tag_idx is not None and len(tokens) > max(chain_tag_idx, seq_tag_idx):
                        chain = tokens[chain_tag_idx].strip()
                        seq = tokens[seq_tag_idx].strip()
                        if chain and seq and seq != "." and seq != "?":
                            chains[chain].add(seq)
                    i += 1
                # finished this atom_site loop, continue outer loop
                in_atom_site = False
                atom_site_tags = []
                continue  # skip i += 1 at bottom because we've already moved i
        i += 1

    # Fallback: if we couldn't find atom_site entries (rare), try a crude search for lines with _atom_site.label_seq_id values
    if not chains:
        # Attempt a very simple fallback: search for lines beginning with "_atom_site." to get tag positions
        # then read subsequent lines that don't start with "_" or "data_" or "loop_"
        # (this is less robust but better than nothing)
        # This fallback will be simplified due to mmCIF complexity and should rarely be needed.
        logger.debug("No atom_site-detected chains from mmCIF parsing; returning empty chain map.")
    return chains


def _parse_residues(path: str, is_pdb: bool) -> Dict[str, Set[str]]:
    if is_pdb:
        return _parse_pdb_residues(path)
    else:
        return _parse_cif_residues(path)


def should_normalize_resseq(ref_path: str, model_path: str, is_pdb: bool, threshold: float = 0.9) -> bool:
    """
    Compare per-chain residue-number overlap between reference and model files (PDB or mmCIF).

    Returns True if normalization SHOULD be performed (i.e. overlap < threshold).
    The overlap metric used is: (total matched residues across all chains) / (total residues in reference).
    - ref_path: path to reference structure (PDB or mmCIF)
    - model_path: path to model structure (PDB or mmCIF)
    - is_pdb: True if the files are in PDB format, False for mmCIF
    - threshold: fraction of reference residues that must be present in model to skip normalization

    Note: This function is conservative — missing chains in the model count as 0 matches.
    """
    ref_chains = _parse_residues(ref_path, is_pdb)
    model_chains = _parse_residues(model_path, is_pdb)

    if not ref_chains:
        logger.warning("No residues found in reference file; recommending normalization by default.")
        return True

    total_ref_res = sum(len(s) for s in ref_chains.values())
    if total_ref_res == 0:
        logger.warning("Reference has zero residues; recommending normalization by default.")
        return True

    total_matched = 0
    for chain_id, ref_set in ref_chains.items():
        model_set = model_chains.get(chain_id, set())
        matched = ref_set.intersection(model_set)
        total_matched += len(matched)

    coverage = total_matched / total_ref_res
    logger.info(f"Residue-number coverage (ref -> model): {coverage:.3f} ({total_matched}/{total_ref_res})")

    # If coverage below threshold, we should normalize
    return coverage < threshold

def calculate_allatom_lddt_openstructure(ref_path, model_path, selection=None, normalize_resseq=False):
    """
    Calculate all-atom LDDT using OpenStructure
    
    Args:
        ref_path (str): Path to reference structure
        model_path (str): Path to model structure
        selection (str, optional): Selection string for specific atoms
        
    Returns:
        float: LDDT score or None if calculation fails
    """

    try:
        from ost import mol, io
        from ost.mol.alg import lddt

        # Optionally pre-normalize residue numbers to avoid negative/resnum-offset mismatches
        if normalize_resseq:
            def _normalize_resseq_numbers(in_path):
                try:
                    import warnings as _warnings
                    from Bio import BiopythonWarning
                    from Bio.PDB import PDBParser, PDBIO
                    parser = PDBParser(QUIET=True)
                    structure = parser.get_structure('x', in_path)
                    # Sequential renumbering transiently reuses a sibling's id, which makes Bio.PDB
                    # warn per residue. The final numbering is correct, so silence that noise.
                    with _warnings.catch_warnings():
                        _warnings.simplefilter("ignore", BiopythonWarning)
                        for model_i in structure:
                            for chain in model_i:
                                new_idx = 1
                                for residue in chain:
                                    residue.id = (' ', new_idx, ' ')
                                    new_idx += 1
                    import tempfile, os as _os
                    tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
                    tmp_path = tmp.name
                    tmp.close()
                    io_writer = PDBIO()
                    io_writer.set_structure(structure)
                    io_writer.save(tmp_path)
                    return tmp_path
                except Exception:
                    return in_path

            ref_load_path = _normalize_resseq_numbers(ref_path)
            model_load_path = _normalize_resseq_numbers(model_path)
        else:
            ref_load_path = ref_path
            model_load_path = model_path

        if ref_load_path.endswith(".cif"):

            def clean_cif_for_openstructure(in_path, out_path):
                with open(in_path) as fin, open(out_path, "w") as fout:
                    next(fin)  # skip the first line
                    for line in fin:
                        fout.write(line)

            # Convert CIF to PDB for OpenStructure
            with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb", delete=False) as temp_ref_pdb, \
                    tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb", delete=False) as temp_model_pdb, \
                    tempfile.NamedTemporaryFile(mode="w+", suffix=".cif", delete=False) as cleaned_model_cif:

                # Clean model CIF (remove leading comments)
                clean_cif_for_openstructure(model_load_path, cleaned_model_cif.name)

                # Convert CIF -> PDB
                convert_cif_to_pdb(ref_load_path, temp_ref_pdb.name)
                convert_cif_to_pdb(cleaned_model_cif.name, temp_model_pdb.name)

                # Load into OpenStructure
                ref = io.LoadPDB(temp_ref_pdb.name)
                model = io.LoadPDB(temp_model_pdb.name)

                # with open(temp_model_pdb.name) as f:
                #     for i in range(5):
                #         print(repr(f.readline()))

        else:
            # Load reference and model structures
            ref = io.LoadPDB(ref_load_path)
            model = io.LoadPDB(model_load_path)

        # lDDT is superposition-free; skip superposition entirely to avoid failures from residue numbering
        ref_atoms = ref.GetAtomCount()
        model_atoms = model.GetAtomCount()
        if ref_atoms == 0 or model_atoms == 0:
            logger.warning(f"No atoms for LDDT: ref={ref_atoms}, model={model_atoms}")
            return None

        # # Create chain mapping if multiple chains exist
        # chain_mapping = None
        # ref_chains = [chain.GetName() for chain in ref.GetChainList()]
        # model_chains = [chain.GetName() for chain in model.GetChainList()]
        #
        # if len(ref_chains) > 1 or len(model_chains) > 1:
        #     # Create chain mapping by matching chains with similar residue sequences
        #     chain_mapping = {}
        #
        #     for model_chain_name in model_chains:
        #         model_chain = model.FindChain(model_chain_name)
        #         if model_chain is None:
        #             continue
        #
        #         # Get model chain sequence (first few residues for matching)
        #         model_residues = [res.GetName() for res in model_chain.GetResidueList()]
        #         model_seq_start = model_residues[:5] if len(model_residues) >= 5 else model_residues
        #
        #         best_match = None
        #         best_score = 0
        #
        #         for ref_chain_name in ref_chains:
        #             ref_chain = ref.FindChain(ref_chain_name)
        #             if ref_chain is None:
        #                 continue
        #
        #             # Get reference chain sequence
        #             ref_residues = [res.GetName() for res in ref_chain.GetResidueList()]
        #             ref_seq_start = ref_residues[:5] if len(ref_residues) >= 5 else ref_residues
        #
        #             # Calculate similarity score
        #             if len(model_seq_start) == len(ref_seq_start):
        #                 matches = sum(1 for m, r in zip(model_seq_start, ref_seq_start) if m == r)
        #                 score = matches / len(model_seq_start) if model_seq_start else 0
        #
        #                 if score > best_score and score > 0.6:  # Require at least 60% match
        #                     best_score = score
        #                     best_match = ref_chain_name
        #
        #         if best_match:
        #             chain_mapping[model_chain_name] = best_match
        #             logger.info(f"Mapped model chain {model_chain_name} to reference chain {best_match} (score: {best_score:.2f})")
        #
        #     if chain_mapping:
        #         logger.info(f"Created chain mapping: {chain_mapping}")
        #     else:
        #         logger.warning("No valid chain mappings found - using single chain mode")
        #         chain_mapping = None
        
        # Build a robust chain mapping to handle differing chain names (e.g., G vs B)
        chain_mapping = None
        try:
            ref_chains = [ch.GetName() for ch in ref.GetChainList()]
            model_chains = [ch.GetName() for ch in model.GetChainList()]
            if len(ref_chains) == 1 and len(model_chains) == 1:
                chain_mapping = {model_chains[0]: ref_chains[0]}
            elif len(ref_chains) and len(model_chains):
                # print("£%@£%££££££", len(ref_chains), len(model_chains))
                # Fallback: map by index up to min length
                chain_mapping = {m: r for m, r in zip(model_chains, ref_chains)}
        except Exception:
            chain_mapping = None

        # Try different OpenStructure LDDT APIs
        try:
            # Method 1: Using lDDTScorer class with correct API
            if hasattr(lddt, 'lDDTScorer'):
                # Create scorer with reference as target
                scorer = lddt.lDDTScorer(ref)
                if selection:
                    # For RNA, exclude hydrogens - use selection string directly
                    # Use the lDDT method, not Score method
                    global_score, _ = scorer.lDDT(
                        model,
                        thresholds=[0.5, 1.0, 2.0, 4.0],
                        chain_mapping=chain_mapping,
                        check_resnames=False,
                    )
                else:
                    global_score, _ = scorer.lDDT(
                        model,
                        thresholds=[0.5, 1.0, 2.0, 4.0],
                        chain_mapping=chain_mapping,
                        check_resnames=False,
                    )
                lddt_value = global_score
                # logger.info("Used Method 1: lDDTScorer.lDDT method")
            else:
                raise AttributeError("lDDTScorer class not found")
                
        except (AttributeError, TypeError) as e1:
            logger.info(f"Method 1 failed: {str(e1)}")
            try:
                # Method 2: Try lddt function if available
                if hasattr(lddt, 'lddt'):
                    if selection:
                        lddt_value = lddt.lddt(ref, model, selection)
                    else:
                        lddt_value = lddt.lddt(ref, model)
                    logger.info("Used Method 2: lddt.lddt function")
                else:
                    raise AttributeError("lddt.lddt function not found")
                    
            except (AttributeError, TypeError) as e2:
                logger.info(f"Method 2 failed: {str(e2)}")
                try:
                    # Method 3: Try direct module call
                    if selection:
                        lddt_value = lddt(ref, model, selection)
                    else:
                        lddt_value = lddt(ref, model)
                    logger.info("Used Method 3: Direct lddt module call")
                except (AttributeError, TypeError) as e3:
                    logger.error(f"All OpenStructure LDDT methods failed: {str(e3)}")
                    # Fallback to a simple implementation
                    # logger.info("Falling back to simple distance-based calculation")
                    # return calculate_simple_lddt(ref, model, selection)
        
        # Treat zero score as None per requirement
        if lddt_value is None:
            return None
        # lddt_value = round(lddt_value, 4)
        return None if lddt_value == 0.0 else lddt_value
        
    except ImportError:
        logger.error("OpenStructure not available. Please install OpenStructure to use this function.")
        return None
    except Exception as e:
        logger.error(f"Error calculating all-atom LDDT: {str(e)}")
        return None

import subprocess
import json
import os
import re

def _chain_seqs_from_aln(aln_list):
    """Build ref_chain_to_seq and pred_chain_to_seq from ost compare-structures JSON 'aln'."""
    ref_seqs = {}
    pred_seqs = {}
    if not aln_list:
        return ref_seqs, pred_seqs
    for block in aln_list:
        if not isinstance(block, str):
            continue
        # Each block can be ">reference:X\nSEQ\n>model:Y\nSEQ" or ">reference:X\nSEQ"
        part = block.strip()
        while part:
            if part.startswith(">reference:"):
                idx = part.find("\n")
                if idx == -1:
                    break
                header = part[:idx]
                rest = part[idx + 1:]
                next_ref = rest.find(">reference:")
                next_mdl = rest.find(">model:")
                if next_ref >= 0 and (next_mdl < 0 or next_ref < next_mdl):
                    seq_part = rest[:next_ref].strip()
                    part = rest[next_ref:]
                elif next_mdl >= 0:
                    seq_part = rest[:next_mdl].strip()
                    part = rest[next_mdl:]
                else:
                    seq_part = rest.replace("\n", "").replace(" ", "")
                    part = ""
                chain = header.replace(">reference:", "").strip()
                ref_seqs[chain] = seq_part
            elif part.startswith(">model:"):
                idx = part.find("\n")
                if idx == -1:
                    break
                header = part[:idx]
                rest = part[idx + 1:]
                next_ref = rest.find(">reference:")
                next_mdl = rest.find(">model:")
                if next_ref >= 0 and (next_mdl < 0 or next_ref < next_mdl):
                    seq_part = rest[:next_ref].strip()
                    part = rest[next_ref:]
                elif next_mdl >= 0:
                    seq_part = rest[:next_mdl].strip()
                    part = rest[next_mdl:]
                else:
                    seq_part = rest.replace("\n", "").replace(" ", "")
                    part = ""
                chain = header.replace(">model:", "").strip()
                pred_seqs[chain] = seq_part
            else:
                break
    return ref_seqs, pred_seqs

# OST_fail triggers (what sets 0 vs 1 vs list of chain IDs):
# - Success path (subprocess exit 0, output_json exists):
#   - unmapped = _unmapped_from_run(stderr, stdout, log). If unmapped non-empty -> os_fail = unmapped; else 0.
#   - We only parse "Could not map model chain X to any chem group" and "Chain X of the model was not mapped".
#   - If OST uses different wording or only puts the message in JSON log with different key, unmapped stays [] -> 0.
#   - Fallback: if os_fail is still 0 but lddt > 0 and ilddt/qs/dockq/ics all 0 -> treat as mapping failure -> os_fail = 1.
# - Retry path: after retry with chain mapping, if retry_unmapped empty -> os_fail = 0; else os_fail = retry_unmapped.
# - No output_json path (exit 0 but no file): we return os_fail = 0 (could be wrong if OST wrote only to stderr).
# - CalledProcessError path: we parse e.stderr (and e.stdout if we add it) for unmapped chain IDs; else os_fail = 1.
#   We do not have the JSON log on exception, so messages only in log are missed.

# Patterns for OST unmapped-chain messages (OST may use either wording)
_RE_UNMAPPED_MODEL_CHAIN = re.compile(r"Could not map model chain (\S+) to any chem group", re.I)
_RE_CHAIN_NOT_MAPPED = re.compile(r"Chain (\S+) of the model was not mapped", re.I)

def _unmapped_model_chain_ids(stderr_text, log_entries):
    """Return list of model chain IDs from OST unmapped-chain messages in stderr or JSON log."""
    ids = []
    def scan(text):
        if not text:
            return
        for m in _RE_UNMAPPED_MODEL_CHAIN.finditer(text):
            ids.append(m.group(1))
        for m in _RE_CHAIN_NOT_MAPPED.finditer(text):
            ids.append(m.group(1))
    if stderr_text:
        scan(stderr_text)
    if log_entries:
        for e in log_entries:
            if isinstance(e, dict):
                scan(str(e.get("message", "")))
    return list(dict.fromkeys(ids))  # preserve order, dedupe


def _unmapped_from_run(stderr_text, stdout_text, log_entries):
    """Collect unmapped model chain IDs from stderr, stdout, and JSON log (OST may print to either stream)."""
    ids = _unmapped_model_chain_ids(stderr_text or "", log_entries)
    if stdout_text:
        ids.extend(_unmapped_model_chain_ids(stdout_text, None))
    return list(dict.fromkeys(ids))

def _serialize_os_fail_for_db(os_fail):
    """Serialize os_fail for DB: 0 or [] -> '0'; one chain -> 'X'; multiple -> JSON list string."""
    if not os_fail:
        return "0"
    if isinstance(os_fail, list):
        return os_fail[0] if len(os_fail) == 1 else json.dumps(os_fail)
    return str(os_fail)

def _has_could_not_map(stderr_text, log_entries):
    """True if stderr or JSON log contains unmapped-chain / chem group warning."""
    for text in (stderr_text or "",):
        if text and ("Could not map" in text or "was not mapped" in text or "chem group" in text):
            return True
    if log_entries:
        for e in log_entries:
            if isinstance(e, dict):
                msg = str(e.get("message", ""))
                if "Could not map" in msg or "was not mapped" in msg or "chem group" in msg:
                    return True
    return False

def _chain_mapping_from_aln(aln_list):
    """
    Build ref_chain -> model_chain mapping from ost compare-structures JSON 'aln'.
    Each block has ">reference:X" then ">model:Y"; return dict {ref_chain: model_chain}.
    """
    mapping = {}
    if not aln_list:
        return mapping
    for block in aln_list:
        if not isinstance(block, str):
            continue
        part = block.strip()
        ref_chain = None
        for line in part.split("\n"):
            line = line.strip()
            if line.startswith(">reference:"):
                ref_chain = line.replace(">reference:", "").strip()
            elif line.startswith(">model:") and ref_chain is not None:
                model_chain = line.replace(">model:", "").strip()
                mapping[ref_chain] = model_chain
                ref_chain = None
    return mapping


def _entity_chain_names(entity):
    """Polymer chain names from an OST entity (same loader compare-structures uses)."""
    names = []
    try:
        for ch in entity.chains:
            names.append(ch.GetName())
        return names
    except Exception:
        pass
    try:
        for ch in entity.GetChainList():
            names.append(ch.GetName())
        return names
    except Exception:
        pass
    return []


def _chain_mapping_is_pure_identity(mapping: Optional[Dict[str, str]]) -> bool:
    """True if every pair is ref_chain == model_chain (e.g. A:A, B:B)."""
    if not mapping:
        return False
    return all(str(ref_k).strip() == str(mdl_v).strip() for ref_k, mdl_v in mapping.items())


def _parse_db_chain_field(field) -> list:
    """Normalize ProteinChainIDs / RNAChainIDs from DB to a list of chain id strings."""
    if field is None or field == "":
        return []
    if isinstance(field, (list, tuple)):
        return [str(x).strip() for x in field if str(x).strip()]
    s = str(field).strip()
    if not s:
        return []
    try:
        v = ast.literal_eval(s)
        if isinstance(v, (list, tuple)):
            return [str(x).strip() for x in v if str(x).strip()]
        return [str(v).strip()]
    except (SyntaxError, ValueError, TypeError):
        return [s]


_MIN_PARTITION_PROT_RES = 8
_MIN_PARTITION_NUC_RES = 4


def _bio_chain_sequences_partitioned(structure):
    """
    Split Bio.PDB structure chains into protein vs nucleic sequence dicts
    (same values as get_chain_sequences_bio). Skips short/ligand-like chains so they do not
    steal alignments from real protein/RNA (e.g. 7xw2 ligand vs polymer).
    """
    _eval_dir = os.path.dirname(os.path.abspath(__file__))
    if _eval_dir not in sys.path:
        sys.path.insert(0, _eval_dir)
    from ost_align_structures import RNA_NAMES, get_sequence_from_chain_bio

    prot = {}
    nuc = {}
    for model in structure:
        for chain in model:
            n_aa = 0
            n_na = 0
            for res in chain:
                name = res.get_resname().strip().upper()
                if is_aa(res):
                    n_aa += 1
                elif name in RNA_NAMES:
                    n_na += 1
                elif len(name) == 1 and name in ("A", "C", "G", "U", "T"):
                    n_na += 1
            seq, res_list = get_sequence_from_chain_bio(chain)
            if not seq:
                continue
            if n_aa >= max(n_na, 1) and n_aa >= _MIN_PARTITION_PROT_RES:
                prot[chain.id] = (seq, res_list)
            elif n_na >= _MIN_PARTITION_NUC_RES and n_na > n_aa:
                nuc[chain.id] = (seq, res_list)
            elif n_aa >= _MIN_PARTITION_PROT_RES:
                prot[chain.id] = (seq, res_list)
            elif n_na >= _MIN_PARTITION_NUC_RES:
                nuc[chain.id] = (seq, res_list)
    return prot, nuc


def _expected_ref_db_chain_ids_in_cif(pdb_id: str, ref_path: str) -> Set[str]:
    """DB-listed protein + RNA reference chains that actually appear in the reference mmCIF."""
    from structures.protein import Protein
    from structures.rna import RNA

    _eval_dir = os.path.dirname(os.path.abspath(__file__))
    if _eval_dir not in sys.path:
        sys.path.insert(0, _eval_dir)
    from ost_align_structures import get_chain_sequences_bio

    prot_o = Protein.get_protein_from_db(pdb_id.upper())
    rna_o = RNA.get_rna_from_db(pdb_id.upper())
    listed = _parse_db_chain_field(prot_o.get_protein_chain_IDs()) + _parse_db_chain_field(
        rna_o.get_rna_chain_IDs()
    )
    if not listed:
        return set()
    parser = MMCIFParser(QUIET=True)
    ref_s = parser.get_structure("ref", ref_path)
    ref_keys = set(get_chain_sequences_bio(ref_s).keys())
    return set(listed) & ref_keys


def _db_mapping_covers_ref_chains(mapping: Dict[str, str], pdb_id: str, ref_path: str) -> bool:
    want = _expected_ref_db_chain_ids_in_cif(pdb_id, ref_path)
    if not want:
        return False
    return want <= set(mapping.keys())


def _chain_mapping_from_db_reference(ref_path: str, model_path: str, pdb_id: str) -> Dict[str, str]:
    """
    Map reference protein/RNA chains (from experimental DB ProteinChainIDs & RNAChainIDs) to
    model chains by aligning separately against model protein-only and nucleic-only pools.
    """
    from structures.protein import Protein
    from structures.rna import RNA

    _eval_dir = os.path.dirname(os.path.abspath(__file__))
    if _eval_dir not in sys.path:
        sys.path.insert(0, _eval_dir)
    from ost_align_structures import align_sequences_and_match, get_chain_sequences_bio

    prot_o = Protein.get_protein_from_db(pdb_id.upper())
    rna_o = RNA.get_rna_from_db(pdb_id.upper())
    if not getattr(prot_o, "is_protein", False) or not getattr(rna_o, "is_rna", False):
        return {}

    ref_p_ids = _parse_db_chain_field(prot_o.get_protein_chain_IDs())
    ref_r_ids = _parse_db_chain_field(rna_o.get_rna_chain_IDs())
    parser = MMCIFParser(QUIET=True)
    ref_s = parser.get_structure("ref", ref_path)
    mdl_s = parser.get_structure("mdl", model_path)
    ref_all = get_chain_sequences_bio(ref_s)
    mdl_prot, mdl_nuc = _bio_chain_sequences_partitioned(mdl_s)

    out = {}
    rprot = {c: ref_all[c] for c in ref_p_ids if c in ref_all}
    rrna = {c: ref_all[c] for c in ref_r_ids if c in ref_all}
    if rprot and mdl_prot:
        out.update(align_sequences_and_match(rprot, mdl_prot))
    if rrna and mdl_nuc:
        out.update(align_sequences_and_match(rrna, mdl_nuc))
    return out


def _chain_mapping_sequence_align_partitioned(ref_path: str, model_path: str) -> Dict[str, str]:
    """ref -> model alignment with protein and nucleic chains matched in separate pools."""
    _eval_dir = os.path.dirname(os.path.abspath(__file__))
    if _eval_dir not in sys.path:
        sys.path.insert(0, _eval_dir)
    from ost_align_structures import align_sequences_and_match

    parser = MMCIFParser(QUIET=True)
    ref_s = parser.get_structure("ref", ref_path)
    mdl_s = parser.get_structure("mdl", model_path)
    ref_p, ref_n = _bio_chain_sequences_partitioned(ref_s)
    mdl_p, mdl_n = _bio_chain_sequences_partitioned(mdl_s)
    out = {}
    if ref_p and mdl_p:
        out.update(align_sequences_and_match(ref_p, mdl_p))
    if ref_n and mdl_n:
        out.update(align_sequences_and_match(ref_n, mdl_n))
    return out


def _chain_mapping_sequence_align(ref_path: str, model_path: str) -> Dict[str, str]:
    """
    ref_chain -> model_chain via sequence alignment; prefers protein/nucleic partitioning so
    extra chains (ligands, ions) do not consume model chains needed for A/B polymers.
    """
    part = _chain_mapping_sequence_align_partitioned(ref_path, model_path)
    if part:
        return part
    _eval_dir = os.path.dirname(os.path.abspath(__file__))
    if _eval_dir not in sys.path:
        sys.path.insert(0, _eval_dir)
    from ost_align_structures import align_sequences_and_match, get_chain_sequences_bio

    parser = MMCIFParser(QUIET=True)
    ref_s = parser.get_structure("ref", ref_path)
    mdl_s = parser.get_structure("mdl", model_path)
    ref_seqs = get_chain_sequences_bio(ref_s)
    model_seqs = get_chain_sequences_bio(mdl_s)
    if not ref_seqs or not model_seqs:
        return {}
    return align_sequences_and_match(ref_seqs, model_seqs)


def _enforce_chain_mapping_for_compare(
    ref_path: str, model_path: str, pdb_id: Optional[str] = None
) -> Dict[str, str]:
    """
    Build ref -> model mapping for compare-structures -c.

    Order: identity chain names; experimental DB protein/RNA chain IDs + partitioned model
    alignment; partitioned alignment without DB; legacy all-chain alignment; partial identity.
    """
    ref_ost: Set[str] = set()
    mdl_ost: Set[str] = set()
    try:
        ref_ent = io.LoadEntity(ref_path)
        mdl_ent = io.LoadEntity(model_path)
        ref_ost = set(_entity_chain_names(ref_ent))
        mdl_ost = set(_entity_chain_names(mdl_ent))
    except Exception as e:
        logger.debug("OST LoadEntity for chain names failed (%s); will try sequence alignment.", e)

    common = ref_ost & mdl_ost
    identity_map = {c: c for c in sorted(common)}

    if ref_ost and mdl_ost and ref_ost == mdl_ost and len(identity_map) == len(ref_ost):
        logger.info("enforceMapping: identity -c (matching chain names): %s", identity_map)
        return identity_map

    if pdb_id:
        db_map = _chain_mapping_from_db_reference(ref_path, model_path, pdb_id)
        if db_map and _db_mapping_covers_ref_chains(db_map, pdb_id, ref_path):
            logger.info("enforceMapping: DB ProteinChainIDs/RNAChainIDs ref->model: %s", db_map)
            return db_map
        if db_map:
            logger.debug("DB chain mapping incomplete vs mmCIF (%s); trying partitioned alignment", db_map)

    try:
        seq_map = _chain_mapping_sequence_align(ref_path, model_path)
        if seq_map:
            logger.info("enforceMapping: sequence_alignment ref->model for -c: %s", seq_map)
            return seq_map
    except Exception as e:
        logger.warning("enforceMapping: sequence_alignment chain mapping failed: %s", e)

    if not identity_map:
        ref_chains = _parse_cif_residues(ref_path)
        mdl_chains = _parse_cif_residues(model_path)
        common_txt = set(ref_chains.keys()) & set(mdl_chains.keys())
        identity_map = {c: c for c in sorted(common_txt)}
    if identity_map:
        logger.info("enforceMapping: partial identity -c (fallback): %s", identity_map)
    return identity_map


def _cli_chain_mapping_from_ost_json(data):
    """
    Build ref_chain -> model_chain for CLI -c flags from OST compare-structures JSON.
    Prefer top-level chain_mapping (same keys as get_interface_metrics); fall back to aln blocks.
    """
    if not data or not isinstance(data, dict):
        return {}
    cm = data.get("chain_mapping")
    if isinstance(cm, dict) and cm:
        out = {}
        for k, v in cm.items():
            if k is None or v is None:
                continue
            rk, rv = str(k).strip(), str(v).strip()
            if rk and rv:
                out[rk] = rv
        if out:
            return out
    return _chain_mapping_from_aln(data.get("aln"))


def _stderr_important_only(stderr_text):
    """Return only lines that contain unmapped-chain or chem group / chem mapping."""
    if not stderr_text:
        return ""
    lines = [ln for ln in stderr_text.strip().split("\n") if ln.strip()]
    keep = [
        ln for ln in lines
        if "Could not map" in ln or "was not mapped" in ln or "chem group" in ln or "chem mapping" in ln
    ]
    return "\n".join(keep) if keep else ""


def get_interface_metrics_for_chain_pair(data, ref_chain, pred_chain):
    """
    From full OST compare-structures JSON data, return (ilddt, qs, dockq, ics) for the
    interface that involves ref_chain (reference) and pred_chain (model). Uses
    qs_interfaces / per_interface_*, dockq_interfaces / dockq, contact_reference_interfaces / per_interface_ics.
    data: dict from OST JSON (or None).
    Returns (ilddt, qs, dockq, ics); each can be None if not found or missing.
    """
    if not data:
        return None, None, None, None
    ilddt_global = data.get("ilddt")
    # qs_interfaces: list of [ref_ch1, ref_ch2, mdl_ch1, mdl_ch2]
    qs_if = data.get("qs_interfaces") or []
    per_qs = data.get("per_interface_qs_global") or []
    dockq_if = data.get("dockq_interfaces") or []
    dockq_scores = data.get("dockq") or []
    contact_ref = data.get("contact_reference_interfaces") or []
    per_ics = data.get("per_interface_ics") or []

    def _idx_for_pair(interfaces, ref_c, pred_c):
        for k, iface in enumerate(interfaces):
            if len(iface) < 4:
                continue
            ref_ch1, ref_ch2, mdl_ch1, mdl_ch2 = iface[0], iface[1], iface[2], iface[3]
            if (ref_ch1 == ref_c and mdl_ch1 == pred_c) or (ref_ch2 == ref_c and mdl_ch2 == pred_c):
                return k
        return -1

    k_qs = _idx_for_pair(qs_if, ref_chain, pred_chain)
    k_dockq = _idx_for_pair(dockq_if, ref_chain, pred_chain)
    # contact_reference_interfaces: list of [ref_ch1, ref_ch2] (reference side); model side from chain_mapping
    chain_mapping = data.get("chain_mapping") or {}
    k_ics = -1
    for k, iface in enumerate(contact_ref):
        if len(iface) >= 2 and ref_chain in (iface[0], iface[1]):
            mdl_ch1 = chain_mapping.get(iface[0])
            mdl_ch2 = chain_mapping.get(iface[1])
            if pred_chain in (mdl_ch1, mdl_ch2):
                k_ics = k
                break

    qs_val = per_qs[k_qs] if k_qs >= 0 and k_qs < len(per_qs) else None
    dockq_val = dockq_scores[k_dockq] if k_dockq >= 0 and k_dockq < len(dockq_scores) else None
    ics_val = per_ics[k_ics] if k_ics >= 0 and k_ics < len(per_ics) else None
    return ilddt_global, qs_val, dockq_val, ics_val


def write_interface_pdbs_ost(ref_cif_path, pred_cif_path, ref_chain_ids, pred_chain_ids, ref_out_path, pred_out_path):
    """
    Write interface PDBs using OST (load CIF, select chains, save PDB) so format matches
    what ost compare-structures expects. Returns True on success, False otherwise.
    """
    try:
        ref_ent = io.LoadEntity(ref_cif_path)
        pred_ent = io.LoadEntity(pred_cif_path)
        ref_sel = " or ".join(f"cname={c}" for c in ref_chain_ids)
        pred_sel = " or ".join(f"cname={c}" for c in pred_chain_ids)
        ref_view = ref_ent.Select(ref_sel)
        pred_view = pred_ent.Select(pred_sel)
        if ref_view.GetChainCount() == 0 or pred_view.GetChainCount() == 0:
            logger.warning("OST interface export: empty chain selection (ref=%s pred=%s)", ref_chain_ids, pred_chain_ids)
            return False
        io.SavePDB(ref_view, ref_out_path)
        io.SavePDB(pred_view, pred_out_path)
        return True
    except Exception as e:
        logger.warning("OST write_interface_pdbs failed: %s", e)
        return False


def compute_lddt_with_cli(
    ref_path,
    model_path,
    lddt=None,
    chain_mapping=None,
    enforce_interface_mapping=False,
    pdb_id=None,
    log_prefix="",
    residue_number_alignment=False,
):
    """
    Compute LDDT (and optionally ilddt, qs, dockq, ics) using ost compare-structures CLI.
    chain_mapping: optional dict ref_chain -> model_chain (e.g. {"A": "A", "C": "B"}).

    When lddt is None (interface bundle: lddt, iLDDT, QS, DockQ, ICS), chain_mapping is normally
    not passed to OST. If enforce_interface_mapping is True, the first run uses -c ref:model:
    caller-supplied chain_mapping if set, otherwise ref->model mapping from matching chain names,
    experimental DB ProteinChainIDs/RNAChainIDs (when pdb_id is passed) with protein/RNA pools on
    the model, else partitioned sequence alignment.
    Explicit -c ref:model is passed only when the mapping is not pure identity (ref != model for some pair).
    Pure identity (e.g. {'A':'A','B':'B'}) omits -c so OST uses the same path as without --enforceMapping;
    passing identity -c zeros ICS for some complexes (e.g. 4qg3).
    residue_number_alignment: if False (default), omit --residue-number-alignment so OST uses
    sequence-based mapping like SWISS-MODEL Assessment / ost_rbp_interfaces.py (including chain-ID
    permutations when sequences imply them). If True, add --residue-number-alignment.

    When lddt is None, returns (lddt_comp, ilddt, qs, dockq, ics, os_fail, data).
    os_fail lists unmapped model chains when mapping fails.

    log_prefix: prepended to compare-structures log lines so callers can tell runs apart (e.g. protein
    PDB vs full complex CIF).
    """
    base_dir = os.path.dirname(os.path.dirname(ref_path))
    lddt_dir = os.path.join(base_dir, "lddt")
    os.makedirs(lddt_dir, exist_ok=True)

    model_name = os.path.splitext(os.path.basename(model_path))[0]
    output_json = os.path.join(lddt_dir, f"{model_name}_lddt.json")

    def build_cmd(custom_chain_mapping=None):
        # SWISS-MODEL OST API default: no residue-number-alignment → NW / sequence chain mapping
        # (see ost_rbp_interfaces.py docstring). Enables Model A ↔ Reference B style swaps.
        if lddt:
            c = [
                "ost", "compare-structures",
                "-r", ref_path,
                "-m", model_path,
                "--lddt",
            ]
            if residue_number_alignment:
                c.append("--residue-number-alignment")
            c.extend(["-o", output_json])
        else:
            c = [
                "ost", "compare-structures",
                "-m", model_path,
                "-r", ref_path,
                "-o", output_json,
                "--lddt",
                "--ilddt",
                "--qs-score",
                "--dockq",
            ]
            if residue_number_alignment:
                c.append("--residue-number-alignment")
            c.append("--ics")
            # Interface mode: include short chains in mapping (OST defaults exclude e.g. peptides < 6 residues)
            c.extend(["--min-pep-length", "1", "--min-nuc-length", "1"])
        if custom_chain_mapping:
            for ref_c, mdl_c in custom_chain_mapping.items():
                c.extend(["-c", f"{ref_c}:{mdl_c}"])
        return c

    try:
        lp = log_prefix or ""
        identity_enforce_omit_c = False
        # Interface bundle (lddt is None): with enforce_mapping, pass -c on first run when it helps.
        if lddt is None:
            if enforce_interface_mapping:
                cmd_chain_mapping = (
                    chain_mapping
                    if chain_mapping
                    else _enforce_chain_mapping_for_compare(ref_path, model_path, pdb_id=pdb_id)
                )
                if not cmd_chain_mapping:
                    cmd_chain_mapping = None
                elif _chain_mapping_is_pure_identity(cmd_chain_mapping):
                    # Explicit -c with ref==model breaks ICS / interface metrics for some complexes
                    # (e.g. 4qg3: ICS ~0.94 without -c, 0 with -c A:A B:B). Default ChainMapper matches.
                    logger.info(
                        "enforceMapping: mapping is identity (ref==model); omitting -c "
                        "(non-identity enforced maps still pass -c)."
                    )
                    identity_enforce_omit_c = True
                    cmd_chain_mapping = None
            else:
                cmd_chain_mapping = None
        else:
            cmd_chain_mapping = chain_mapping
            if cmd_chain_mapping and _chain_mapping_is_pure_identity(cmd_chain_mapping):
                logger.info(
                    "%scompare-structures: identity chain mapping -> omitting -c for lddt run",
                    lp,
                )
                cmd_chain_mapping = None
        pass1_had_explicit_c = bool(cmd_chain_mapping)
        cmd = build_cmd(cmd_chain_mapping)
        if lddt is None:
            logger.info("%sost compare-structures: %s", lp, cmd)
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        stderr_short = _stderr_important_only(result.stderr)
        stdout_short = _stderr_important_only(result.stdout or "")  # OST may print warnings to stdout
        if stdout_short and stderr_short:
            stderr_short = stderr_short + "\n" + stdout_short
        elif stdout_short:
            stderr_short = stdout_short
        had_could_not_map = _has_could_not_map(result.stderr, None) or _has_could_not_map(result.stdout, None)

        if os.path.exists(output_json) and lddt == True:
            with open(output_json) as f:
                data = json.load(f)
            if stderr_short:
                logger.info("%sost compare-structures: %s", lp, stderr_short.strip())
            return data.get("lddt", None)

        if os.path.exists(output_json) and lddt == None:
            with open(output_json) as f:
                data = json.load(f)
            if not had_could_not_map:
                had_could_not_map = _has_could_not_map(None, data.get("log")) or bool(result.stdout and "Could not map" in result.stdout)
            lddt_comp = data.get("lddt", None)
            ilddt = data.get("ilddt", None)
            qs = data.get("qs_global", None)
            dockq = data.get("dockq_wave", None)
            ics = data.get("ics", None)
            # Collect unmapped from stderr, stdout and log (OST may print to either when exit 0)
            unmapped = _unmapped_from_run(result.stderr, result.stdout, data.get("log"))
            os_fail = unmapped if unmapped else 0

            # Fallback: LDDT present but all interface metrics zero => mapping failure even if no message parsed
            def _num(x):
                if x is None:
                    return 0.0
                try:
                    return float(x)
                except (TypeError, ValueError):
                    return 0.0
            if not os_fail and _num(lddt_comp) > 0:
                if _num(ilddt) == 0 and _num(qs) == 0 and _num(dockq) == 0 and _num(ics) == 0:
                    os_fail = 1

            # Second pass only if first run had no -c and we did not deliberately omit identity -c.
            if (
                lddt is None
                and enforce_interface_mapping
                and not pass1_had_explicit_c
                and not identity_enforce_omit_c
            ):
                retry_mapping = _cli_chain_mapping_from_ost_json(data)
                if not retry_mapping:
                    logger.warning(
                        "enforceMapping: no chain_mapping or parsable aln in OST JSON; "
                        "second pass with -c skipped (%s)",
                        output_json,
                    )
                if retry_mapping and _chain_mapping_is_pure_identity(retry_mapping):
                    logger.info(
                        "%senforceMapping: JSON chain_mapping is identity; skip second pass -c (same as pass 1)",
                        lp,
                    )
                    retry_mapping = None
                if retry_mapping:
                    try:
                        cmd_retry = build_cmd(retry_mapping)
                        logger.info("%sost compare-structures (retry with -c from JSON): %s", lp, cmd_retry)
                        retry = subprocess.run(
                            cmd_retry, check=True, capture_output=True, text=True
                        )
                        if os.path.exists(output_json):
                            with open(output_json) as f:
                                data2 = json.load(f)
                            retry_unmapped = _unmapped_from_run(
                                retry.stderr, retry.stdout, data2.get("log")
                            )
                            lddt_comp = data2.get("lddt", lddt_comp)
                            ilddt = data2.get("ilddt", ilddt)
                            qs = data2.get("qs_global", qs)
                            dockq = data2.get("dockq_wave", dockq)
                            ics = data2.get("ics", ics)
                            data = data2
                            os_fail = retry_unmapped if retry_unmapped else 0
                            if not os_fail and _num(lddt_comp) > 0:
                                if (
                                    _num(ilddt) == 0
                                    and _num(qs) == 0
                                    and _num(dockq) == 0
                                    and _num(ics) == 0
                                ):
                                    os_fail = 1
                    except subprocess.CalledProcessError as e_retry:
                        logger.warning(
                            "compare-structures retry with -c chain mapping failed: %s",
                            _stderr_important_only(e_retry.stderr or ""),
                        )

            if stderr_short:
                logger.info("%sost compare-structures: %s", lp, stderr_short.strip())
            return lddt_comp, ilddt, qs, dockq, ics, os_fail, data

        if stderr_short:
            logger.info("%sost compare-structures: %s", lp, stderr_short.strip())
        if lddt is None:
            return None, None, None, None, None, 0, None
        return None

    except subprocess.CalledProcessError as e:
        err_short = _stderr_important_only(e.stderr or "")
        stdout_short = _stderr_important_only(e.stdout or "")
        if err_short or stdout_short:
            logger.warning("%sost compare-structures failed: %s", lp, (err_short + "\n" + stdout_short).strip())
        else:
            logger.warning("%sost compare-structures failed (no detail)", lp)
        if lddt is None:
            chains = _unmapped_from_run(e.stderr or "", e.stdout or "", None)
            os_fail = chains if chains else 1
            return None, None, None, None, None, os_fail, None
        return None



def extract_rna_lddt_scores(file_name_cif_af, aligned_folder, residue_number_alignment=False):
    """
    Extract RNA all-atom LDDT scores for a given structure
    
    Args:
        file_name_cif_af (str): Predicted structure filename (e.g., fold_6fql_s1_model_0.cif)
        aligned_folder (str): Path to aligned folder containing RNA PDB files
        residue_number_alignment (bool): Passed to ost compare-structures (default False).
        
    Returns:
        dict: Dictionary containing RNA LDDT score
    """
    try:
        # Extract PDB ID from filename (e.g., fold_6fql_s1_model_0.cif -> 6fql)
        pdb_id = file_name_cif_af.split('_')[1] if '_' in file_name_cif_af else file_name_cif_af[:4]
        
        # Use explicitly aligned single-chain RNA PDBs
        ref_rna_file = os.path.join(aligned_folder, f"{pdb_id}_pdb_aligned_rna.pdb")
        model_rna_file = os.path.join(aligned_folder, f"{pdb_id}_af_aligned_rna.pdb")
        # print("model_rna_file", model_rna_file, " !!!!", ref_rna_file)
        
        if ref_rna_file is None or model_rna_file is None:
            logger.warning(f"Could not find RNA PDB files for {pdb_id} in {aligned_folder}")
            logger.warning(f"Expected files: {pdb_id}_pdb_aligned_rna.pdb and {pdb_id}_af_aligned_rna.pdb")
            return None
        
        logger.info(f"Found RNA files: ref={os.path.basename(ref_rna_file)}, model={os.path.basename(model_rna_file)}")

        # Initialize scores - only RNA
        scores = {
            'rna_lddt_allatom': None
        }

        # Calculate RNA-specific LDDT (exclude hydrogens). Files are pre-aligned single-RNA-chain per instruction
        rna_lddt = compute_lddt_with_cli(
            ref_rna_file,
            model_rna_file,
            lddt=True,
            log_prefix="[RNA LDDT PDB] ",
            residue_number_alignment=residue_number_alignment,
        )

        if rna_lddt is not None:
            if rna_lddt >= 0.8:
                chosen_score = rna_lddt
            else:
                rna_lddt_all = calculate_allatom_lddt_openstructure(
                    ref_rna_file, model_rna_file, "aname!=H", normalize_resseq=True
                )
                chosen_score = max(rna_lddt, rna_lddt_all)
        else:
            rna_lddt_all = calculate_allatom_lddt_openstructure(
                ref_rna_file, model_rna_file, "aname!=H", normalize_resseq=True
            )
            chosen_score = rna_lddt_all

        # Final cleanup
        if chosen_score == 0.0:
            chosen_score = None

        scores['rna_lddt_allatom'] = chosen_score
        logger.info(f"RNA all-atom LDDT: {scores['rna_lddt_allatom']}")

        return scores
        
    except Exception as e:
        logger.error(f"Error extracting RNA LDDT scores for {file_name_cif_af}: {str(e)}")
        return None

def update_database_with_rna_lddt(table_name, file_name_cif_af, scores):
    """
    Update database with RNA all-atom LDDT scores
    
    Args:
        table_name (str): Database table name
        file_name_cif_af (str): Predicted structure filename
        scores (dict): Dictionary containing RNA LDDT score
    """
    try:
        cursor = db.connection.cursor()
        
        if table_name == 'pred_protein_rna':
            cursor.execute('''UPDATE pred_protein_rna 
                             SET RNA_LDDT_allAtom = ?
                             WHERE FileName = ?''',
                         (scores['rna_lddt_allatom'], file_name_cif_af))
        
        elif table_name == 'pred_rna_rna':
            # Also populate Complex_LDDT_allAtom for the rna_rna table. For single-chain RNA the whole
            # complex IS that RNA chain, so complex-level == chain-level LDDT; writing it now keeps the
            # column populated (the rna-only plotting falls back to RNA_LDDT_allAtom otherwise) and
            # future-proofs multichain RNA-RNA complexes where complex-level LDDT becomes distinct.
            try:
                cursor.execute("ALTER TABLE pred_rna_rna ADD COLUMN Complex_LDDT_allAtom REAL")
            except sqlite3.OperationalError as e:
                if "duplicate column name" not in str(e):
                    raise
            cursor.execute('''UPDATE pred_rna_rna
                             SET RNA_LDDT_allAtom = ?, Complex_LDDT_allAtom = ?
                             WHERE FileName = ?''',
                         (scores['rna_lddt_allatom'], scores['rna_lddt_allatom'], file_name_cif_af))
        
        db.connection.commit()
        logger.info(f"Updated {table_name} with RNA all-atom LDDT scores for {file_name_cif_af}")
        
    except Exception as e:
        logger.error(f"Error updating database: {str(e)}")
        raise

def update_protein_and_complex_lddt(
    file_name_cif_af,
    aligned_folder,
    enforce_interface_mapping=False,
    residue_number_alignment=False,
):
    """Compute Protein_LDDT_allAtom (largest protein chain) and Complex_LDDT_allAtom (all atoms) for protein_rna table.
    Assumes aligned_folder holds *_pdb_aligned_protein.pdb and *_af_aligned_protein.pdb and combined *_aligned_complex.pdb if available,
    otherwise falls back to original sources if present.

    residue_number_alignment: forwarded to ost compare-structures (default False = SWISS-MODEL Assessment default).
    """
    pdb_id = file_name_cif_af.split('_')[1] if '_' in file_name_cif_af else file_name_cif_af[:4]
    # Use explicitly aligned, single-chain protein PDBs
    ref_prot = os.path.join(aligned_folder, f"{pdb_id}_pdb_aligned_protein.pdb")
    model_prot = os.path.join(aligned_folder, f"{pdb_id}_af_aligned_protein.pdb")
    # Complex PDBs may not exist; try aligned complex names
    # ref_complex = os.path.join(aligned_folder, f"{pdb_id}_pdb_aligned_complex.pdb")
    # model_complex = os.path.join(aligned_folder, f"{pdb_id}_af_aligned_complex.pdb")

    cursor = db.connection.cursor()
    # Protein LDDT (heavy atoms only, prefer aligned files). Decide normalization based on residue number match
    prot_val = None
    if ref_prot and model_prot:
        prot_val = compute_lddt_with_cli(
            ref_prot,
            model_prot,
            lddt=True,
            log_prefix="[Protein LDDT PDB] ",
            residue_number_alignment=residue_number_alignment,
        )
    cursor.execute('''UPDATE pred_protein_rna SET Protein_LDDT_allAtom = ? WHERE FileName = ?''',
                   (prot_val, file_name_cif_af))

    # Complex LDDT from CIFs
    comp_val = ilddt = qs = dockq = ics = None
    os_fail = 0
    parent_folder = config.get_parent_folder()
    pdb_folder = os.path.join(parent_folder, 'pdb')
    af_folder = os.path.join(parent_folder, 'af')
    ref_cif, model_cif = get_structure_files(f"{pdb_id}.cif", file_name_cif_af, pdb_folder, af_folder)
    if ref_cif and model_cif:
        comp_val, ilddt, qs, dockq, ics, os_fail, _ = compute_lddt_with_cli(
            ref_cif,
            model_cif,
            enforce_interface_mapping=enforce_interface_mapping,
            pdb_id=pdb_id.lower(),
            log_prefix="[Complex mmCIF] ",
            residue_number_alignment=residue_number_alignment,
        )
    else:
        comp_val = ilddt = qs = dockq = ics = None
        os_fail = 0

    # OST_fail is owned by ost_calcRNA_AllChains.py; not written here.
    cursor.execute('''UPDATE pred_protein_rna
                      SET Complex_LDDT_allAtom = ?,
                          iLDDT = ?,
                          QS_global = ?,
                          Dockq = ?,
                          ICS = ?
                      WHERE FileName = ?''',
                   (comp_val, ilddt, qs, dockq, ics, file_name_cif_af))

    db.connection.commit()

def get_files_to_process(table_name, pdb_folder, af_folder, aligned_folder):
    """
    Get list of files to process from the database that actually exist in both af and pdb folders
    
    Args:
        table_name (str): Database table name
        pdb_folder (str): Path to PDB folder
        af_folder (str): Path to AF folder
        aligned_folder (str): Path to aligned folder
        
    Returns:
        list: List of (exp_db_id, file_name) tuples for files that exist
    """
    try:
        cursor = db.connection.cursor()
        
        if table_name == 'pred_protein_rna':
            cursor.execute("SELECT exp_db_id, FileName FROM pred_protein_rna WHERE FileName IS NOT NULL")
        elif table_name == 'pred_rna_rna':
            cursor.execute("SELECT exp_db_id, FileName FROM pred_rna_rna WHERE FileName IS NOT NULL")
        else:
            raise ValueError(f"Unsupported table: {table_name}")
        
        all_files = cursor.fetchall()
        
        # Filter to only include files that exist in both folders
        existing_files = []
        for exp_db_id, file_name_cif_af in all_files:
            # Extract PDB ID from filename
            pdb_id = file_name_cif_af.split('_')[1] if '_' in file_name_cif_af else file_name_cif_af[:4]
            file_name_cif_pdb = f"{pdb_id}.cif"
            
            # Check if CIF files exist in both folders
            ref_cif_path = os.path.join(pdb_folder, file_name_cif_pdb)
            model_cif_path = os.path.join(af_folder, file_name_cif_af)
            
            # Check if aligned RNA PDB files exist (required for RNA LDDT)
            ref_rna_pdb = os.path.join(aligned_folder, f"{pdb_id}_pdb_aligned_rna.pdb")
            model_rna_pdb = os.path.join(aligned_folder, f"{pdb_id}_af_aligned_rna.pdb")
            
            # Only include if both CIF files exist AND aligned RNA PDB files exist
            if (os.path.exists(ref_cif_path) and os.path.exists(model_cif_path) and
                os.path.exists(ref_rna_pdb) and os.path.exists(model_rna_pdb)):
                existing_files.append((exp_db_id, file_name_cif_af))
            else:
                logger.debug(f"Skipping {file_name_cif_af}: missing files - "
                           f"ref_cif={os.path.exists(ref_cif_path)}, "
                           f"model_cif={os.path.exists(model_cif_path)}, "
                           f"ref_rna_pdb={os.path.exists(ref_rna_pdb)}, "
                           f"model_rna_pdb={os.path.exists(model_rna_pdb)}")
        
        logger.info(f"Filtered {len(all_files)} database entries to {len(existing_files)} files with existing structure files")
        return existing_files
        
    except Exception as e:
        logger.error(f"Error getting files to process: {str(e)}")
        raise

def process_rna_lddt_scores(
    table_name,
    aligned_folder,
    pdb_folder,
    af_folder,
    enforce_interface_mapping=False,
    residue_number_alignment=False,
):
    """
    Process RNA all-atom LDDT scores for all structures in the specified table
    
    Args:
        table_name (str): Database table name ('pred_protein_rna' or 'pred_rna_rna')
        aligned_folder (str): Path to aligned folder containing RNA PDB files
        pdb_folder (str): Path to PDB folder
        af_folder (str): Path to AF folder
        enforce_interface_mapping (bool): Optional explicit -c chain mapping (see --enforceMapping).
        residue_number_alignment (bool): Pass --residue-number-alignment to compare-structures (default
            False = SWISS-MODEL Assessment style sequence mapping).
    """
    logger.info(f"Processing RNA all-atom LDDT scores for {table_name}")
    if table_name == "pred_protein_rna":
        logger.info(
            "compare-structures: residue_number_alignment=%s (False matches SWISS-MODEL Assessment default)",
            residue_number_alignment,
        )
    if enforce_interface_mapping:
        logger.info(
            "Interface bundle: --enforceMapping enabled (explicit -c when non-identity after DB/alignment)"
        )
    
    # Ensure RNA LDDT column exists
    ensure_rna_lddt_column(table_name)
    # Ensure protein/complex columns for protein_rna table
    if table_name == 'pred_protein_rna':
        ensure_protein_complex_columns(table_name)
    
    # Get files to process (only those that exist in both folders)
    files = get_files_to_process(table_name, pdb_folder, af_folder, aligned_folder)
    logger.info(f"Found {len(files)} files to process (with existing structure files)")
    
    processed_count = 0
    error_count = 0

    from tqdm import tqdm
    for exp_db_id, file_name_cif_af in tqdm(files, total=len(files), desc="structures", unit="struct"):
        try:
            logger.info(f"Processing {exp_db_id}: {file_name_cif_af}")

            # print("file_name_cif_af!!!", file_name_cif_af)
            # Extract RNA LDDT scores
            scores = extract_rna_lddt_scores(
                file_name_cif_af, aligned_folder, residue_number_alignment=residue_number_alignment
            )
            
            if scores:
                # Update database
                update_database_with_rna_lddt(table_name, file_name_cif_af, scores)
                # Optionally compute protein and complex LDDT for protein_rna
                if table_name == 'pred_protein_rna':
                    try:
                        update_protein_and_complex_lddt(
                            file_name_cif_af,
                            aligned_folder,
                            enforce_interface_mapping=enforce_interface_mapping,
                            residue_number_alignment=residue_number_alignment,
                        )
                    except Exception as e:
                        logger.warning(f"Protein/Complex LDDT update failed for {file_name_cif_af}: {e}")
                processed_count += 1
            else:
                error_count += 1
                logger.warning(f"Could not extract RNA LDDT scores for {file_name_cif_af}")
                
        except Exception as e:
            error_count += 1
            logger.error(f"Error processing {file_name_cif_af}: {str(e)}")
            logger.error(traceback.format_exc())
    
    logger.info(f"Processing complete. Processed: {processed_count}, Errors: {error_count}")

def main():
    """
    Main function to process RNA all-atom LDDT scores
    """
    parser = argparse.ArgumentParser(
        description=(
            "Extract RNA/protein/complex LDDT and interface scores via local ost compare-structures. "
            "Defaults match SWISS-MODEL Assessment (no --residue-number-alignment; automatic chain mapping)."
        ),
    )
    parser.add_argument(
        "table_name",
        choices=["pred_protein_rna", "pred_rna_rna"],
        help="Prediction table to process",
    )
    parser.add_argument(
        "--residue-number-alignment",
        "--ost-residue-number-alignment",
        "--resnum-align",
        action="store_true",
        dest="residue_number_alignment",
        help=(
            "Pass --residue-number-alignment to every ost compare-structures call. Default OFF matches "
            "SWISS-MODEL Structure Assessment (sequence-based chain mapping; same default as "
            "evaluation/ost_rbp_interfaces.py). Use ON only when reference/model numbering already match."
        ),
    )
    parser.add_argument(
        "--enforceMapping",
        action="store_true",
        help=(
            "Optional fallback for pred_protein_rna: force explicit -c ref:model when computed mapping "
            "is non-identity (DB ProteinChainIDs/RNAChainIDs + partitioned alignment). Default pipeline "
            "lets OST auto-map like SWISS-MODEL (no -c)."
        ),
    )
    args = parser.parse_args()
    table_name = args.table_name
    enforce_interface_mapping = bool(args.enforceMapping)
    residue_number_alignment = bool(getattr(args, "residue_number_alignment", False))
    if os.environ.get("SWISSMODEL_OST_RESNUM_ALIGN", "").strip().lower() in ("1", "true", "yes"):
        residue_number_alignment = True

    # Get folders from config
    parent_folder = config.get_parent_folder()
    aligned_folder = os.path.join(parent_folder, 'aligned')
    pdb_folder = os.path.join(parent_folder, 'pdb')
    af_folder = os.path.join(parent_folder, 'af')
    
    if not os.path.exists(aligned_folder):
        print(f"Error: Aligned directory not found: {aligned_folder}")
        sys.exit(1)
    
    if not os.path.exists(pdb_folder):
        print(f"Error: PDB directory not found: {pdb_folder}")
        sys.exit(1)
    
    if not os.path.exists(af_folder):
        print(f"Error: AF directory not found: {af_folder}")
        sys.exit(1)
    
    try:
        process_rna_lddt_scores(
            table_name,
            aligned_folder,
            pdb_folder,
            af_folder,
            enforce_interface_mapping=enforce_interface_mapping,
            residue_number_alignment=residue_number_alignment,
        )
        logger.info("RNA all-atom LDDT processing completed successfully")
        
    except Exception as e:
        logger.error(f"Error in main processing: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    # Quiet per-sample INFO chatter so the tqdm progress bar renders cleanly (warnings/errors stay).
    logging.getLogger().setLevel(logging.WARNING)
    main()
