"""
RNA AllChains runner. Writes three JSON keyed-dict columns to the predicted table (+ OST_fail):

  RNA_LDDT_AllChains   -> robust all-vs-all RNA-chain LDDT, best per predicted RNA chain, keyed by
                          the single predicted RNA chain: [{"D": 0.53}, {"E": 0.59}]
  RNA_iLDDT_AllChains  -> per protein-RNA interface, keyed by the predicted "protein:RNA" pair
  RNA_ICS_AllChains       (RBP_-style): [{"A:E": 0.73}, {"B:D": 0.58}]  (filter at evaluation)

RNA_LDDT does NOT use OST chain mapping (all-vs-all, then best per predicted chain). Interfaces use
one full-complex ``ost compare-structures`` (sequence mapping) for the ref->pred chain mapping +
per-interface ICS; iLDDT is computed by extracting each protein-RNA interface and reading its
(global) iLDDT. QS/DockQ dropped -> use RBP_QS / RBP_DockQ from ost_rbp_interfaces.py.

Modes: default (all three), ``--lddt`` (LDDT only), ``--interfaces`` (iLDDT/ICS + OST_fail).
"""
import os
import re
import sys
import json
import ast
import sqlite3
import tempfile
import traceback
import logging

from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

from database.startConfig import StartConfig
from structures.rna import RNA
from evaluation.chain_assignment import assigned_by_pred_chain, assigned_from_map
from evaluation.ost_extractAllAtomLDDT import (
    calculate_allatom_lddt_openstructure,
    compute_lddt_with_cli,
)

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def _safe_filename_segment(s: str) -> str:
    """Sanitize a chain id for use in a filename."""
    return re.sub(r"[^\w.\-]", "_", str(s))


# Basic classification of RNA residues (standard bases + common modifications)
RNA_BASES = {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}
RNA_MODIFIED = {
    "6MZ", "PSU", "5MC", "OMC", "4OC", "5MU", "OMU", "UR3", "A2M", "MA6",
    "2MG", "OMG", "7MG", "RSQ", "5CM", "C34", "5HC", "6OG", "6MA", "1CC",
    "8OG", "5FC", "3DR",
}


def is_rna_residue(residue):
    """Return True if residue looks like an RNA (or modified RNA) nucleotide."""
    resname = residue.get_resname().strip().upper()
    return resname in RNA_BASES or resname in RNA_MODIFIED


def get_rna_chain_lengths(structure):
    """
    Get RNA chain lengths from a Bio.PDB structure.
    Returns: dict(chain_id -> residue_count) for RNA-only chains.
    """
    chain_lengths = {}
    for model in structure:
        for chain in model:
            count = 0
            for residue in chain:
                if is_rna_residue(residue):
                    count += 1
            if count > 0:
                # If multiple models are present, accumulate counts per chain ID
                chain_lengths[chain.id] = chain_lengths.get(chain.id, 0) + count
    # print("chain_lengths", chain_lengths)
    return chain_lengths


def get_protein_chain_ids(structure):
    """
    Get chain IDs that contain protein (amino acid) residues.
    Returns: set of chain_id for chains that have at least one standard amino acid.
    """
    protein_chains = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue):
                    protein_chains.add(chain.id)
                    break
    return protein_chains


def _parse_sequence_chain_mapping(rna_sequence_val, rna_chain_ids_val):
    """
    Parse RNASequence and RNAChainIDs from DB into dict chain_id -> sequence.
    RNASequence can be '["seq1", "seq2"]' or 'seq1'; RNAChainIDs can be '["D", "E"]' or 'D'.
    """
    chain_to_seq = {}
    if not rna_sequence_val or not rna_chain_ids_val:
        return chain_to_seq
    try:
        seq_raw = rna_sequence_val.strip() if isinstance(rna_sequence_val, str) else rna_sequence_val
        ids_raw = rna_chain_ids_val.strip() if isinstance(rna_chain_ids_val, str) else rna_chain_ids_val
        if not seq_raw or not ids_raw:
            return chain_to_seq
        try:
            seqs = ast.literal_eval(seq_raw) if isinstance(seq_raw, str) else seq_raw
        except (ValueError, SyntaxError):
            seqs = [seq_raw] if isinstance(seq_raw, str) else []
        try:
            ids = ast.literal_eval(ids_raw) if isinstance(ids_raw, str) else ids_raw
        except (ValueError, SyntaxError):
            ids = [ids_raw] if ids_raw else []
        if not isinstance(seqs, (list, tuple)):
            seqs = [seqs]
        if not isinstance(ids, (list, tuple)):
            ids = [ids]
        for k, cid in enumerate(ids):
            if k < len(seqs) and cid is not None:
                chain_to_seq[str(cid).strip()] = str(seqs[k]).strip() if seqs[k] else ""
    except Exception:
        pass
    return chain_to_seq


def _parse_length_chain_mapping(rna_length_val, rna_chain_ids_val):
    """
    Parse RNALength and RNAChainIDs from DB into dict chain_id -> length (int).
    RNALength can be [8, 5, 8] or 8; RNAChainIDs can be '["F", "G", "H"]' or 'F'.
    """
    chain_to_len = {}
    if rna_length_val is None or not rna_chain_ids_val:
        return chain_to_len
    try:
        len_raw = rna_length_val
        ids_raw = rna_chain_ids_val.strip() if isinstance(rna_chain_ids_val, str) else rna_chain_ids_val
        if ids_raw is None:
            return chain_to_len
        try:
            lengths = ast.literal_eval(str(len_raw).strip()) if isinstance(len_raw, str) else len_raw
        except (ValueError, SyntaxError):
            lengths = [int(len_raw)] if len_raw is not None else []
        try:
            ids = ast.literal_eval(ids_raw) if isinstance(ids_raw, str) else ids_raw
        except (ValueError, SyntaxError):
            ids = [ids_raw] if ids_raw else []
        if not isinstance(lengths, (list, tuple)):
            lengths = [lengths]
        if not isinstance(ids, (list, tuple)):
            ids = [ids]
        for k, cid in enumerate(ids):
            if k < len(lengths) and cid is not None:
                try:
                    chain_to_len[str(cid).strip()] = int(lengths[k])
                except (TypeError, ValueError):
                    pass
    except Exception:
        pass
    return chain_to_len


def _parse_ordered_chain_data(rna_sequence_val, rna_chain_ids_val, rna_length_val):
    """
    Parse RNASequence, RNAChainIDs, RNALength from DB into ordered lists (same index = same chain).
    Returns (chain_ids_list, lengths_list, seqs_list) in DB order; missing data yields empty list for that slot.
    """
    chain_ids = []
    lengths = []
    seqs = []
    if not rna_chain_ids_val:
        return chain_ids, lengths, seqs
    try:
        ids_raw = rna_chain_ids_val.strip() if isinstance(rna_chain_ids_val, str) else rna_chain_ids_val
        ids = ast.literal_eval(ids_raw) if isinstance(ids_raw, str) else ids_raw
        if not isinstance(ids, (list, tuple)):
            ids = [ids]
        chain_ids = [str(cid).strip() for cid in ids if cid is not None]
    except Exception:
        return [], [], []

    if rna_length_val is not None:
        try:
            len_raw = rna_length_val
            lengths_parsed = ast.literal_eval(str(len_raw).strip()) if isinstance(len_raw, str) else len_raw
            if not isinstance(lengths_parsed, (list, tuple)):
                lengths_parsed = [lengths_parsed]
            for k in range(len(chain_ids)):
                lengths.append(int(lengths_parsed[k]) if k < len(lengths_parsed) else 0)
        except (ValueError, TypeError, SyntaxError):
            pass
    if len(lengths) != len(chain_ids):
        lengths = [0] * len(chain_ids)

    if rna_sequence_val:
        try:
            seq_raw = rna_sequence_val.strip() if isinstance(rna_sequence_val, str) else rna_sequence_val
            seqs_parsed = ast.literal_eval(seq_raw) if isinstance(seq_raw, str) else seq_raw
            if not isinstance(seqs_parsed, (list, tuple)):
                seqs_parsed = [seqs_parsed]
            for k in range(len(chain_ids)):
                seqs.append(str(seqs_parsed[k]).strip() if k < len(seqs_parsed) and seqs_parsed[k] else "")
        except (ValueError, SyntaxError):
            pass
    if len(seqs) != len(chain_ids):
        seqs = [""] * len(chain_ids)

    return chain_ids, lengths, seqs


def get_chain_to_sequence_from_db(cursor, ref_table: str, pred_table: str, pdb_id: str, file_name_cif_af: str):
    """
    Load ref and pred RNASequence / RNAChainIDs / RNALength and return
    (ref_chain_to_sequence, pred_chain_to_sequence, ref_chain_to_length, pred_chain_to_length,
     ref_ordered, pred_ordered).
    ref_ordered = (ref_chain_ids_list, ref_lengths_list, ref_seqs_list) in DB order (for matrix layout).
    pred_ordered = (pred_chain_ids_list, pred_lengths_list, pred_seqs_list) in DB order.
    Pairing is by RNALength and sequence: only chains with same length and same sequence are compared.
    """
    ref_chain_to_seq = {}
    pred_chain_to_seq = {}
    ref_chain_to_len = {}
    pred_chain_to_len = {}
    ref_ordered = ([], [], [])  # chain_ids, lengths, seqs
    pred_ordered = ([], [], [])
    pdb_id = pdb_id.upper()
    try:
        cursor.execute(
            f"SELECT RNASequence, RNAChainIDs, RNALength FROM {ref_table} WHERE PDBId = ? OR UniprotId = ? LIMIT 1",
            (pdb_id, pdb_id),
        )
        row = cursor.fetchone()
        if row:
            ref_chain_to_seq = _parse_sequence_chain_mapping(row[0], row[1])
            if len(row) > 2 and row[2] is not None:
                ref_chain_to_len = _parse_length_chain_mapping(row[2], row[1])
            ref_ordered = _parse_ordered_chain_data(row[0], row[1], row[2] if len(row) > 2 else None)
    except sqlite3.OperationalError:
        try:
            cursor.execute(
                f"SELECT RNASequence, RNAChainIDs FROM {ref_table} WHERE PDBId = ? OR UniprotId = ? LIMIT 1",
                (pdb_id, pdb_id),
            )
            row = cursor.fetchone()
            if row:
                ref_chain_to_seq = _parse_sequence_chain_mapping(row[0], row[1])
                ref_ordered = _parse_ordered_chain_data(row[0], row[1], None)
        except Exception as e:
            logger.debug("DB ref lookup failed for %s: %s", file_name_cif_af, e)
    except Exception as e:
        logger.debug("DB ref lookup failed for %s: %s", file_name_cif_af, e)

    try:
        cursor.execute(
            f"SELECT RNASequence, RNAChainIDs, RNALength FROM {pred_table} WHERE exp_db_id = ? AND FileName = ? LIMIT 1",
            (pdb_id, file_name_cif_af),
        )
        row = cursor.fetchone()
        if row:
            pred_chain_to_seq = _parse_sequence_chain_mapping(row[0], row[1])
            if len(row) > 2 and row[2] is not None:
                pred_chain_to_len = _parse_length_chain_mapping(row[2], row[1])
            pred_ordered = _parse_ordered_chain_data(row[0], row[1], row[2] if len(row) > 2 else None)
    except sqlite3.OperationalError:
        try:
            cursor.execute(
                f"SELECT RNASequence, RNAChainIDs FROM {pred_table} WHERE exp_db_id = ? AND FileName = ? LIMIT 1",
                (pdb_id, file_name_cif_af),
            )
            row = cursor.fetchone()
            if row:
                pred_chain_to_seq = _parse_sequence_chain_mapping(row[0], row[1])
                pred_ordered = _parse_ordered_chain_data(row[0], row[1], None)
        except Exception as e:
            logger.debug("DB pred lookup failed for %s: %s", file_name_cif_af, e)
    except Exception as e:
        logger.debug("DB pred lookup failed for %s: %s", file_name_cif_af, e)

    return ref_chain_to_seq, pred_chain_to_seq, ref_chain_to_len, pred_chain_to_len, ref_ordered, pred_ordered


class RNAChainSelect(Select):
    """Select a single RNA chain (and only RNA residues) for PDBIO.save."""

    def __init__(self, target_chain_id):
        self.target_chain_id = target_chain_id

    def accept_chain(self, chain):
        return chain.id == self.target_chain_id

    def accept_residue(self, residue):
        return is_rna_residue(residue)


class ContactInterfaceSelect(Select):
    """
    Select the given RNA chain and all chains in contact with it (protein and/or RNA within distance threshold).
    """

    def __init__(self, rna_chain_id, contact_chain_ids):
        self.rna_chain_id = rna_chain_id
        self.contact_chain_ids = set(contact_chain_ids) if contact_chain_ids else {rna_chain_id}

    def accept_chain(self, chain):
        return chain.id in self.contact_chain_ids

    def accept_residue(self, residue):
        return is_rna_residue(residue) or is_aa(residue)


# Column names for per-chain metric matrices (same logic as RNA_LDDT_AllChains)
# Kept columns (QS/DockQ dropped as obsolete — use RBP_QS/RBP_DockQ). Each is stored as a JSON
# keyed dict of the best (max) pairing per reference RNA chain, e.g. [{"C:A": 0.59}, {"D:B": 0.61}].
RNA_METRIC_ALLCHAINS_COLUMNS = [
    "RNA_LDDT_AllChains",
    "RNA_iLDDT_AllChains",
    "RNA_ICS_AllChains",
]
# Interface-only metrics (--interfaces): OST_fail is only set when these are evaluated
INTERFACE_METRIC_COLUMNS = [
    "RNA_iLDDT_AllChains",
    "RNA_ICS_AllChains",
]


def ensure_rna_metric_allchains_column(cursor, table_name: str, column_name: str):
    """
    Ensure a per-chain metric column exists (same pattern as RNA_TM_AllChains).
    """
    try:
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} REAL")
        logger.info(f"Added column {column_name} to {table_name}")
    except sqlite3.OperationalError as e:
        if "duplicate column name" in str(e):
            logger.debug(f"Column {column_name} already exists in {table_name}")
        else:
            raise


def ensure_rna_lddt_allchains_column(cursor, table_name: str, column_name: str = "RNA_LDDT_AllChains"):
    """Backward-compat: ensure RNA_LDDT_AllChains (REAL for legacy)."""
    try:
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} REAL")
        logger.info(f"Added column {column_name} to {table_name}")
    except sqlite3.OperationalError as e:
        if "duplicate column name" in str(e):
            logger.debug(f"Column {column_name} already exists in {table_name}")
        else:
            raise


def ensure_os_fail_column(cursor, table_name: str):
    """Ensure OST_fail column exists (TEXT: '0', 'X', or JSON list of chain IDs). Only for pred_protein_rna."""
    if table_name != "pred_protein_rna":
        return
    try:
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN OST_fail TEXT")
        logger.info(f"Added column OST_fail to {table_name}")
    except sqlite3.OperationalError as e:
        if "duplicate column name" in str(e):
            logger.debug(f"Column OST_fail already exists in {table_name}")
        else:
            raise


def _float_or_zero(val, decimals=None):
    """Return the RAW float value, or 0.0 if None/invalid.

    No rounding: raw full-precision values are stored so downstream statistics are not biased by
    truncation. (`decimals` is accepted for backward compatibility and ignored.)
    """
    if val is None:
        return 0.0
    try:
        return float(val)
    except (TypeError, ValueError):
        return 0.0


def _read_rna_chain_assignment(cursor, table_name, file_name_cif_af):
    """Return the {pred_chain: ref_chain} RNA assignment the TM step stored for this row, else None.

    Gracefully returns None when the RNA_chain_assignment column does not exist (older DBs) or is
    empty, so the LDDT step falls back to its own one-to-one assignment.
    """
    try:
        cursor.execute(
            f"SELECT RNA_chain_assignment FROM {table_name} WHERE FileName = ?",
            (file_name_cif_af,),
        )
        row = cursor.fetchone()
    except sqlite3.OperationalError:
        return None
    if not row or row[0] in (None, "", "null"):
        return None
    try:
        m = json.loads(row[0])
    except (ValueError, TypeError):
        return None
    return m if isinstance(m, dict) and m else None


def _serialize_os_fail_for_db(os_fail):
    """Serialize OST_fail for DB: 0 -> '0'; one chain -> 'X'; list -> JSON string."""
    if not os_fail:
        return "0"
    if isinstance(os_fail, list):
        return os_fail[0] if len(os_fail) == 1 else json.dumps(os_fail)
    return str(os_fail)


def _rna_lddt_from_aligned(aligned_folder: str, pdb_id: str):
    """
    Compute RNA LDDT using aligned RNA-only PDBs (same as extractAllAtomLDDT).
    Returns float or None if files missing / computation fails.
    So RNA_LDDT_AllChains can match RNA_LDDT_allAtom when aligned files exist.
    """
    ref_rna = os.path.join(aligned_folder, f"{pdb_id}_pdb_aligned_rna.pdb")
    model_rna = os.path.join(aligned_folder, f"{pdb_id}_af_aligned_rna.pdb")
    if not os.path.isfile(ref_rna) or not os.path.isfile(model_rna):
        return None
    try:
        rna_lddt = compute_lddt_with_cli(ref_rna, model_rna, lddt=True)
        if rna_lddt is not None:
            if rna_lddt >= 0.8:
                return _float_or_zero(rna_lddt)
            lddt_all = calculate_allatom_lddt_openstructure(
                ref_rna, model_rna, "aname!=H", normalize_resseq=True
            )
            return _float_or_zero(max(rna_lddt, lddt_all) if lddt_all is not None else rna_lddt)
        lddt_all = calculate_allatom_lddt_openstructure(
            ref_rna, model_rna, "aname!=H", normalize_resseq=True
        )
        return _float_or_zero(lddt_all) if lddt_all is not None else None
    except Exception as e:
        logger.debug("Aligned RNA LDDT failed for %s: %s", pdb_id, e)
        return None


def _rna_lddt_from_aligned_pair(
    aligned_folder: str, pdb_id: str, ref_chain: str, pred_chain: str
):
    """
    Compute RNA LDDT for one (ref_chain, pred_chain) pair using per-chain aligned PDBs.
    Looks for {pdb_id}_pdb_aligned_rna_{ref_chain}_{pred_chain}.pdb and
    {pdb_id}_af_aligned_rna_{ref_chain}_{pred_chain}.pdb.
    Returns float or None if files missing / computation fails.
    """
    safe_ref = ref_chain.replace(" ", "_")
    safe_pred = pred_chain.replace(" ", "_")
    ref_rna = os.path.join(
        aligned_folder, f"{pdb_id}_pdb_aligned_rna_{safe_ref}_{safe_pred}.pdb"
    )
    model_rna = os.path.join(
        aligned_folder, f"{pdb_id}_af_aligned_rna_{safe_ref}_{safe_pred}.pdb"
    )
    if not os.path.isfile(ref_rna) or not os.path.isfile(model_rna):
        return None
    try:
        rna_lddt = compute_lddt_with_cli(ref_rna, model_rna, lddt=True)
        if rna_lddt is not None:
            if rna_lddt >= 0.8:
                return _float_or_zero(rna_lddt)
            lddt_all = calculate_allatom_lddt_openstructure(
                ref_rna, model_rna, "aname!=H", normalize_resseq=True
            )
            return _float_or_zero(max(rna_lddt, lddt_all) if lddt_all is not None else rna_lddt)
        lddt_all = calculate_allatom_lddt_openstructure(
            ref_rna, model_rna, "aname!=H", normalize_resseq=True
        )
        return _float_or_zero(lddt_all) if lddt_all is not None else None
    except Exception as e:
        logger.debug(
            "Aligned RNA LDDT failed for %s chain pair %s/%s: %s",
            pdb_id, ref_chain, pred_chain, e,
        )
        return None


def _compute_rna_only_lddt(
    structure_reference,
    structure_predicted,
    ref_chain: str,
    pred_chain: str,
    pdb_id: str,
    parent_folder: str,
    file_name_cif_af: str,
):
    """
    Extract ref_chain and pred_chain (RNA only) to PDBs under parent_folder/lddt/ so that
    OST compare-structures can write its JSON there (base_dir = parent_folder). Return LDDT float or None.
    """
    lddt_dir = os.path.join(parent_folder, "lddt") if parent_folder else None
    if not lddt_dir:
        try:
            lddt_dir = tempfile.mkdtemp(prefix="lddt_")
        except Exception:
            return None
    os.makedirs(lddt_dir, exist_ok=True)
    safe_ref = ref_chain.replace(" ", "_").replace("/", "_")
    safe_pred = pred_chain.replace(" ", "_").replace("/", "_")
    base_name = f"{pdb_id}_ref_{safe_ref}_pred_{safe_pred}"
    ref_path = os.path.join(lddt_dir, f"{base_name}_ref.pdb")
    model_path = os.path.join(lddt_dir, f"{base_name}_model.pdb")
    try:
        io_ref = PDBIO()
        io_ref.set_structure(structure_reference)
        io_ref.save(ref_path, RNAChainSelect(ref_chain))
        io_pred = PDBIO()
        io_pred.set_structure(structure_predicted)
        io_pred.save(model_path, RNAChainSelect(pred_chain))
        logger.info(
            "[RNA LDDT AllChains] OST input: ref_chain=%s pred_chain=%s ref_path=%s model_path=%s chain_mapping=%s",
            ref_chain, pred_chain, ref_path, model_path, {ref_chain: pred_chain},
        )
        # OST matches chains by ID by default; ref has ref_chain and model has pred_chain, so pass explicit mapping
        lddt_val = compute_lddt_with_cli(ref_path, model_path, lddt=True, chain_mapping={ref_chain: pred_chain})
        if lddt_val is None:
            lddt_val = calculate_allatom_lddt_openstructure(
                ref_path, model_path, "aname!=H", normalize_resseq=True
            )
        logger.info(
            "[RNA LDDT AllChains] OST result: ref_chain=%s pred_chain=%s lddt=%s",
            ref_chain, pred_chain, lddt_val,
        )
        return float(lddt_val) if lddt_val is not None else None
    except Exception as e:
        logger.warning(
            "RNA-only LDDT failed for %s chain %s/%s: %s",
            file_name_cif_af, ref_chain, pred_chain, e,
        )
        return None
    finally:
        for p in (ref_path, model_path):
            if p and os.path.isfile(p):
                try:
                    os.remove(p)
                except OSError:
                    pass


def _lddt_by_pred_chain(matrix, ref_chains, pred_chains, chain_map=None):
    """
    Reduce the all-vs-all RNA LDDT matrix (rows = reference chains, cols = predicted chains) to one
    RAW entry per PREDICTED RNA chain via a one-to-one assignment, keyed by the predicted chain id ->
    [{"D": 0.53}, {"E": 0.59}, {"F": 0.52}]. Predicted chains left unmatched or with a non-positive
    score are skipped. Matches RNAChainIDs (predicted).

    If `chain_map` ({pred_chain: ref_chain}) is given -- the assignment the TM step already computed --
    LDDT is read for the SAME reference partner so TM and LDDT never refer to different GT chains.
    Otherwise a fresh one-to-one maximum-weight assignment is run on this LDDT matrix (still avoids the
    greedy per-column max that let several predicted chains claim the same reference chain).
    """
    if chain_map:
        return assigned_from_map(matrix, ref_chains, pred_chains, chain_map, positive_only=True)
    values, _ = assigned_by_pred_chain(matrix, ref_chains, pred_chains, positive_only=True)
    return values


def _classify_ref_chains_from_aln(aln):
    """
    Classify each REFERENCE chain as 'rna' or 'protein' from the OST compare-structures 'aln' field,
    which carries per-chain sequences in OST's own chain namespace (>reference:<chain>\\n<seq>).
    This avoids the Biopython (auth_asym_id) vs OST (label_asym_id) chain-id mismatch that otherwise
    breaks protein-RNA interface detection.
    """
    types = {}
    for entry in (aln or []):
        m = re.match(r">reference:(\S+)\s*\n([^\n>]+)", str(entry))
        if not m:
            continue
        ch = m.group(1)
        seq = re.sub(r"[^A-Za-z]", "", m.group(2)).upper()
        if not seq:
            continue
        types[ch] = "rna" if set(seq) <= set("ACGUN") else "protein"
    return types


def _ilddt_for_protein_rna_pair(ref_path, pred_path, ref_pair, pred_pair, pdb_id,
                                parent_folder, file_name_cif_af):
    """
    Extract the two-chain protein-RNA interface (ref_pair from the reference, pred_pair from the
    model) using OST's own loader so the chain ids match OST's namespace, run OST, and return the
    global iLDDT of that 2-chain sub-structure (= the interface iLDDT for this protein-RNA pair).
    ref_pair / pred_pair are (protein_chain, rna_chain) in OST's chain namespace. None on failure.
    """
    from ost import io as _ost_io

    lddt_dir = os.path.join(parent_folder, "lddt") if parent_folder else None
    if not lddt_dir:
        try:
            lddt_dir = tempfile.mkdtemp(prefix="ilddt_")
        except Exception:
            return None
    os.makedirs(lddt_dir, exist_ok=True)
    tag = f"{pdb_id}_iface_{_safe_filename_segment(ref_pair[0])}_{_safe_filename_segment(ref_pair[1])}"
    ref_out = os.path.join(lddt_dir, f"{tag}_ref.pdb")
    mdl_out = os.path.join(lddt_dir, f"{tag}_mdl.pdb")
    try:
        ref_ent = _ost_io.LoadMMCIF(ref_path)
        mdl_ent = _ost_io.LoadMMCIF(pred_path)
        ref_sel = ref_ent.Select(f"cname='{ref_pair[0]}' or cname='{ref_pair[1]}'")
        mdl_sel = mdl_ent.Select(f"cname='{pred_pair[0]}' or cname='{pred_pair[1]}'")
        _ost_io.SavePDB(ref_sel, ref_out)
        _ost_io.SavePDB(mdl_sel, mdl_out)
        result = compute_lddt_with_cli(
            ref_out, mdl_out, lddt=None,
            chain_mapping={ref_pair[0]: pred_pair[0], ref_pair[1]: pred_pair[1]},
        )
        if isinstance(result, (list, tuple)) and len(result) >= 2:
            return result[1]  # global iLDDT of the 2-chain interface
        return None
    except Exception as e:
        logger.warning("iLDDT pair failed %s %s/%s: %s", file_name_cif_af, ref_pair, pred_pair, e)
        return None
    finally:
        for p in (ref_out, mdl_out):
            if p and os.path.isfile(p):
                try:
                    os.remove(p)
                except OSError:
                    pass


def compute_protein_rna_interfaces(file_name_cif_pdb, file_name_cif_af, pdb_folder, af_folder,
                                   parent_folder=None, residue_number_alignment=False):
    """
    Per protein-RNA interface metrics, keyed by the PREDICTED "protein:RNA" chain pair (RBP_-style),
    e.g. [{"A:E": 0.73}, {"B:D": 0.58}].

    Everything is done in OST's chain namespace: chains are classified protein vs RNA from the
    'aln' sequences in the OST JSON (NOT from Biopython, whose author chain ids can disagree with
    OST's and would misclassify interfaces). ICS comes from per_interface_ics; iLDDT is obtained by
    extracting each protein-RNA interface with OST's loader and reading the 2-chain global iLDDT.

    Returns {"RNA_ICS_AllChains": [...], "RNA_iLDDT_AllChains": [...], "OST_fail": ...}.
    """
    pdb_id = file_name_cif_af.split("_")[1] if "_" in file_name_cif_af else file_name_cif_af[:4]
    ref_path = os.path.join(pdb_folder, file_name_cif_pdb)
    pred_path = os.path.join(af_folder, file_name_cif_af)
    _parent = parent_folder or (os.path.dirname(pdb_folder) if pdb_folder else None)

    out = {"RNA_ICS_AllChains": [], "RNA_iLDDT_AllChains": [], "OST_fail": 0}

    # One full-complex OST call: chain_mapping (ref->model) + per-interface ICS + per-chain sequences.
    result = compute_lddt_with_cli(
        ref_path, pred_path, lddt=None, residue_number_alignment=residue_number_alignment
    )
    ost_json = result[6] if isinstance(result, (list, tuple)) and len(result) >= 7 and isinstance(result[6], dict) else None
    os_fail = result[5] if isinstance(result, (list, tuple)) and len(result) >= 6 else 0
    if os_fail:
        out["OST_fail"] = os_fail if isinstance(os_fail, list) else 1
    if not ost_json:
        logger.warning("[RNA interfaces] %s: no OST JSON returned (OST could not map chains).", pdb_id)
        return out

    cmap = ost_json.get("chain_mapping") or {}
    ifaces = ost_json.get("contact_reference_interfaces") or []
    ics_vals = ost_json.get("per_interface_ics") or []
    chain_types = _classify_ref_chains_from_aln(ost_json.get("aln"))  # OST-namespace ref chain -> rna/protein
    logger.info("[RNA interfaces] %s chain_types=%s (%d interfaces)", pdb_id, chain_types, len(ifaces))

    ics_list, ilddt_list, seen = [], [], set()
    for k, iface in enumerate(ifaces):
        if not iface or len(iface) < 2:
            continue
        a, b = str(iface[0]), str(iface[1])
        ta, tb = chain_types.get(a), chain_types.get(b)
        # Keep only protein-RNA interfaces (one protein chain, one RNA chain on the reference side).
        if ta == "protein" and tb == "rna":
            ref_prot_c, ref_rna_c = a, b
        elif tb == "protein" and ta == "rna":
            ref_prot_c, ref_rna_c = b, a
        else:
            continue
        pred_prot, pred_rna = cmap.get(ref_prot_c), cmap.get(ref_rna_c)
        if not pred_prot or not pred_rna:
            continue  # unmapped chain -> cannot key by pred pair or extract
        key = f"{pred_prot}:{pred_rna}"
        if key in seen:
            continue
        seen.add(key)
        ics = ics_vals[k] if k < len(ics_vals) else None
        ics_list.append({key: _float_or_zero(ics)})
        ilddt = _ilddt_for_protein_rna_pair(
            ref_path, pred_path, (ref_prot_c, ref_rna_c), (pred_prot, pred_rna),
            pdb_id, _parent, file_name_cif_af,
        )
        ilddt_list.append({key: _float_or_zero(ilddt)})

    out["RNA_ICS_AllChains"] = ics_list
    out["RNA_iLDDT_AllChains"] = ilddt_list
    return out


def compute_rna_lddt_by_pred_chain(
    file_name_cif_pdb, file_name_cif_af, pdb_folder, af_folder,
    ref_chain_to_sequence, pred_chain_to_sequence,
    ref_chain_to_len=None, pred_chain_to_len=None,
    ref_ordered=None, pred_ordered=None,
    aligned_folder=None, parent_folder=None,
    chain_map=None,
):
    """
    All-vs-all RNA-chain LDDT (robust; does NOT use OST chain mapping). Only ref/pred chains with the
    same RNALength (and same sequence when both are known) are compared, via single-chain RNA-only OST
    LDDT. The matrix is reduced to one RAW value per PREDICTED RNA chain by a one-to-one assignment;
    if `chain_map` ({pred_chain: ref_chain}, from the TM step) is given it is reused so LDDT scores the
    same reference partner as TM.
    Returns {"RNA_LDDT_AllChains": [{"<pred_rna_chain>": lddt}, ...]}.
    """
    pdb_id = file_name_cif_af.split("_")[1] if "_" in file_name_cif_af else file_name_cif_af[:4]
    ref_chain_to_len = ref_chain_to_len or {}
    pred_chain_to_len = pred_chain_to_len or {}
    ref_ids_ord, ref_lens_ord, _ = ref_ordered or ([], [], [])
    pred_ids_ord, pred_lens_ord, _ = pred_ordered or ([], [], [])

    aligned_single = _rna_lddt_from_aligned(aligned_folder, pdb_id) if aligned_folder else None

    parser = MMCIFParser(QUIET=True)
    ref_path = os.path.join(pdb_folder, file_name_cif_pdb)
    pred_path = os.path.join(af_folder, file_name_cif_af)

    def _load(structure_id, primary_path, folder_path):
        # mmCIF may be truncated (missing _atom_site.id); retry ../done/<sub>/<name>.
        try:
            return parser.get_structure(structure_id, primary_path)
        except (FileNotFoundError, OSError, KeyError) as e:
            if isinstance(e, KeyError) and "_atom_site.id" not in str(e):
                raise
            done_dir = os.path.join(os.path.dirname(folder_path), "done", os.path.basename(folder_path))
            fb = os.path.join(done_dir, os.path.basename(primary_path))
            if os.path.isfile(fb):
                logger.warning("MMCIF load failed (%s); using fallback %s", os.path.basename(primary_path), fb)
                return parser.get_structure(structure_id, fb)
            raise

    structure_reference = _load(file_name_cif_pdb, ref_path, pdb_folder)
    structure_predicted = _load(file_name_cif_af, pred_path, af_folder)
    ref_counts = get_rna_chain_lengths(structure_reference)
    pred_counts = get_rna_chain_lengths(structure_predicted)

    if ref_ids_ord and pred_ids_ord:
        ref_chains = [(cid, ref_lens_ord[i] if i < len(ref_lens_ord) else ref_counts.get(cid, 0))
                      for i, cid in enumerate(ref_ids_ord)]
        pred_chains = [(cid, pred_lens_ord[i] if i < len(pred_lens_ord) else pred_counts.get(cid, 0))
                       for i, cid in enumerate(pred_ids_ord)]
    else:
        ref_chains = sorted(ref_counts.items())
        pred_chains = sorted(pred_counts.items())

    if not ref_chains or not pred_chains:
        return {"RNA_LDDT_AllChains": []}

    single_chain_mode = len(ref_chains) == 1 and len(pred_chains) == 1
    n_ref, n_pred = len(ref_chains), len(pred_chains)
    lddt_matrix = [[0.0] * n_pred for _ in range(n_ref)]
    _parent = parent_folder or (os.path.dirname(pdb_folder) if pdb_folder else None)

    logger.info("[RNA LDDT AllChains] %s ref_chains=%s pred_chains=%s", pdb_id, ref_chains, pred_chains)
    for i, (ref_chain, ref_len) in enumerate(ref_chains):
        for j, (pred_chain, pred_len) in enumerate(pred_chains):
            if not single_chain_mode:
                if ref_chain_to_len.get(ref_chain, ref_len) != pred_chain_to_len.get(pred_chain, pred_len):
                    continue  # length mismatch
                ref_seq = ref_chain_to_sequence.get(ref_chain) if ref_chain_to_sequence else None
                pred_seq = pred_chain_to_sequence.get(pred_chain) if pred_chain_to_sequence else None
                if ref_seq is not None and pred_seq is not None and ref_seq != pred_seq:
                    continue  # sequence mismatch
            if ref_chain not in ref_counts or pred_chain not in pred_counts:
                continue
            # Prefer aligned per-pair / single-chain LDDT files, else compute RNA-only LDDT via OST.
            val = _rna_lddt_from_aligned_pair(aligned_folder, pdb_id, ref_chain, pred_chain) if aligned_folder else None
            if val is None and single_chain_mode and aligned_single is not None:
                val = aligned_single
            if val is None:
                val = _compute_rna_only_lddt(structure_reference, structure_predicted,
                                             ref_chain, pred_chain, pdb_id, _parent, file_name_cif_af)
            lddt_matrix[i][j] = _float_or_zero(val)

    return {"RNA_LDDT_AllChains": _lddt_by_pred_chain(lddt_matrix, ref_chains, pred_chains, chain_map)}


def update_database_with_lddt_matrix(
    file_name_cif_pdb: str,
    file_name_cif_af: str,
    pdb_folder: str,
    af_folder: str,
    only_lddt: bool = False,
    only_interfaces: bool = False,
    residue_number_alignment: bool = False,
    **_legacy_ignored,
):
    """
    Compute and store the RNA AllChains metrics as JSON keyed dicts:
      RNA_LDDT_AllChains  -> robust all-vs-all RNA LDDT, best per predicted RNA chain: [{"D":0.53}, ...]
      RNA_iLDDT_AllChains / RNA_ICS_AllChains -> per protein-RNA interface, keyed by predicted
                                                 "protein:RNA" pair: [{"A:E":0.73}, ...]
    only_lddt -> only RNA_LDDT_AllChains (OST_fail not set). only_interfaces -> only iLDDT/ICS + OST_fail.
    Legacy interface-export kwargs are accepted and ignored (mapping now comes from OST directly).
    """
    file_pdb = os.path.splitext(file_name_cif_pdb)[0]
    file_af = os.path.splitext(file_name_cif_af)[0]
    pdb_id = file_name_cif_af.split("_")[1] if "_" in file_name_cif_af else file_name_cif_af[:4]

    rna = None
    if len(file_af) > 6:
        parts = file_name_cif_af.split("_")
        rna = RNA.get_rna_from_db(id=parts[1].upper(), file_name=file_name_cif_af)
    elif len(file_pdb) in [4, 6]:
        rna = RNA.get_rna_from_db(id=file_pdb.upper())

    if not rna or not rna.is_rna:
        logger.info("%s: no RNA present according to database; skipping.", file_name_cif_af)
        return

    config = StartConfig()
    database_path = config.get_database_path()
    table_name = config.get_predicted_table_name()
    ref_table = config.get_reference_table_name()
    parent_folder = config.get_parent_folder()
    aligned_folder = os.path.join(parent_folder, "aligned") if parent_folder else None
    if not aligned_folder or not os.path.isdir(aligned_folder):
        aligned_folder = None

    matrices = {"RNA_LDDT_AllChains": [], "RNA_iLDDT_AllChains": [], "RNA_ICS_AllChains": [], "OST_fail": 0}

    # RNA_LDDT: robust all-vs-all RNA comparison, reduced to the best score per predicted RNA chain.
    if not only_interfaces:
        conn = sqlite3.connect(database_path)
        cursor = conn.cursor()
        try:
            ref_chain_to_seq, pred_chain_to_seq, ref_chain_to_len, pred_chain_to_len, ref_ordered, pred_ordered = get_chain_to_sequence_from_db(
                cursor, ref_table, table_name, pdb_id, file_name_cif_af
            )
            # Reuse the RNA chain assignment the TM step (extractEvalMetrics) already wrote, so LDDT is
            # scored against the same reference partner as TM. Absent/old rows -> None -> LDDT does its
            # own one-to-one assignment on the LDDT matrix (still not the greedy per-column max).
            tm_chain_map = _read_rna_chain_assignment(cursor, table_name, file_name_cif_af)
        finally:
            conn.close()
        lddt_res = compute_rna_lddt_by_pred_chain(
            file_name_cif_pdb, file_name_cif_af, pdb_folder, af_folder,
            ref_chain_to_seq, pred_chain_to_seq,
            ref_chain_to_len=ref_chain_to_len, pred_chain_to_len=pred_chain_to_len,
            ref_ordered=ref_ordered, pred_ordered=pred_ordered,
            aligned_folder=aligned_folder, parent_folder=parent_folder,
            chain_map=tm_chain_map,
        )
        matrices["RNA_LDDT_AllChains"] = lddt_res.get("RNA_LDDT_AllChains", [])

    # RNA_iLDDT + RNA_ICS: per protein-RNA interface, keyed by predicted "protein:RNA" pair.
    if not only_lddt:
        iface = compute_protein_rna_interfaces(
            file_name_cif_pdb, file_name_cif_af, pdb_folder, af_folder,
            parent_folder=parent_folder, residue_number_alignment=residue_number_alignment,
        )
        matrices["RNA_iLDDT_AllChains"] = iface.get("RNA_iLDDT_AllChains", [])
        matrices["RNA_ICS_AllChains"] = iface.get("RNA_ICS_AllChains", [])
        matrices["OST_fail"] = iface.get("OST_fail", 0)

    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    try:
        ensure_rna_lddt_allchains_column(cursor, table_name, "RNA_LDDT_AllChains")
        if only_lddt:
            # Only LDDT; OST_fail is not set
            cursor.execute(
                f"UPDATE {table_name} SET RNA_LDDT_AllChains = ? WHERE FileName = ?",
                (json.dumps(matrices["RNA_LDDT_AllChains"]), file_name_cif_af),
            )
            logger.info(
                "Updated RNA_LDDT_AllChains only for FileName=%s (%d RNA chains)",
                file_name_cif_af, len(matrices["RNA_LDDT_AllChains"]),
            )
        elif only_interfaces:
            # Only interface metrics; OST_fail is set when interface metrics are evaluated
            for col in INTERFACE_METRIC_COLUMNS:
                ensure_rna_metric_allchains_column(cursor, table_name, col)
            ensure_os_fail_column(cursor, table_name)
            columns_to_update = list(INTERFACE_METRIC_COLUMNS) + (
                ["OST_fail"] if table_name == "pred_protein_rna" else []
            )
            placeholders = ", ".join(f"{col} = ?" for col in columns_to_update)
            values = [json.dumps(matrices[col]) for col in INTERFACE_METRIC_COLUMNS]
            if table_name == "pred_protein_rna":
                values.append(_serialize_os_fail_for_db(matrices.get("OST_fail", 0)))
            values.append(file_name_cif_af)
            cursor.execute(
                f"UPDATE {table_name} SET {placeholders} WHERE FileName = ?",
                values,
            )
            logger.info(
                "Updated RNA interface AllChains (iLDDT/ICS) and OST_fail for FileName=%s (%d iLDDT, %d ICS entries)",
                file_name_cif_af, len(matrices["RNA_iLDDT_AllChains"]), len(matrices["RNA_ICS_AllChains"]),
            )
        else:
            # Full run: all metrics and OST_fail
            for col in RNA_METRIC_ALLCHAINS_COLUMNS[1:]:  # iLDDT, ICS
                ensure_rna_metric_allchains_column(cursor, table_name, col)
            ensure_os_fail_column(cursor, table_name)

            columns_to_update = list(RNA_METRIC_ALLCHAINS_COLUMNS) + (
                ["OST_fail"] if table_name == "pred_protein_rna" else []
            )
            placeholders = ", ".join(f"{col} = ?" for col in columns_to_update)
            values = [json.dumps(matrices[col]) for col in RNA_METRIC_ALLCHAINS_COLUMNS]
            if table_name == "pred_protein_rna":
                values.append(_serialize_os_fail_for_db(matrices.get("OST_fail", 0)))
            values.append(file_name_cif_af)

            cursor.execute(
                f"UPDATE {table_name} SET {placeholders} WHERE FileName = ?",
                values,
            )
            logger.info(
                "Updated %s RNA_*_AllChains for FileName=%s (%d LDDT, %d iLDDT, %d ICS entries)",
                table_name, file_name_cif_af,
                len(matrices["RNA_LDDT_AllChains"]), len(matrices["RNA_iLDDT_AllChains"]),
                len(matrices["RNA_ICS_AllChains"]),
            )

        conn.commit()
    finally:
        conn.close()


def process_files(
    pdb_path: str,
    af_path: str,
    only_lddt: bool = False,
    only_interfaces: bool = False,
    residue_number_alignment: bool = False,
):
    """
    Iterate over predicted CIF files, find matching reference CIFs, and
    fill RNA_LDDT_AllChains / RNA_iLDDT_AllChains / RNA_ICS_AllChains for each (RNA-only).
    only_lddt -> only RNA_LDDT_AllChains (OST_fail not set); only_interfaces -> only iLDDT/ICS + OST_fail.
    Prints a progress bar (per-sample INFO logging is suppressed by the CLI).
    """
    from tqdm import tqdm

    cif_files = sorted(f for f in os.listdir(af_path) if f.endswith(".cif"))
    for file_name_cif_af in tqdm(cif_files, total=len(cif_files), desc="structures", unit="struct"):
        try:
            pdb_id = file_name_cif_af.split("_")[1]
            file_name_cif_pdb = f"{pdb_id}.cif"
            pdb_file_path = os.path.join(pdb_path, file_name_cif_pdb)

            if not os.path.exists(pdb_file_path):
                logger.warning("No corresponding PDB CIF found for %s", pdb_id)
                continue

            try:
                update_database_with_lddt_matrix(
                    file_name_cif_pdb,
                    file_name_cif_af,
                    pdb_path,
                    af_path,
                    only_lddt=only_lddt,
                    only_interfaces=only_interfaces,
                    residue_number_alignment=residue_number_alignment,
                )
            except Exception as e:
                tb = traceback.format_exc()
                logger.error("Error calculating RNA *_AllChains for %s: %s", pdb_id, e)
                logger.error("Traceback:\n%s", tb)
                continue

        except Exception as e:
            logger.error("Error processing %s: %s", file_name_cif_af, e)
            continue


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compute RNA LDDT/iLDDT/ICS AllChains (all-vs-all, reduced to per-chain max keyed dicts).")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--lddt",
        action="store_true",
        help="Only compute and update RNA_LDDT_AllChains (skip interface metrics; OST_fail not set).",
    )
    group.add_argument(
        "--interfaces",
        action="store_true",
        help="Only compute and update interface metrics (RNA_iLDDT_AllChains, RNA_ICS_AllChains) and OST_fail.",
    )
    parser.add_argument(
        "--residue-number-alignment",
        action="store_true",
        help=(
            "Pass --residue-number-alignment to OST's interface call. Default OFF: sequence-based "
            "mapping maps far more chains here (matches RBP_DockQ); resnum-align often maps only RNA."
        ),
    )
    args = parser.parse_args()

    # Quiet the per-sample chatter: only the progress bar (+ hard errors) should show. Root at
    # WARNING silences INFO; the OST helper is pushed to ERROR to hide its noisy per-call warnings
    # (e.g. "ost compare-structures failed (no detail)").
    logging.getLogger().setLevel(logging.WARNING)
    logging.getLogger("evaluation.ost_extractAllAtomLDDT").setLevel(logging.ERROR)

    config = StartConfig()
    parent_folder = config.get_parent_folder()
    if parent_folder is None:
        print("Error: parent_folder is not set in StartConfig/config.json")
        sys.exit(1)

    pdb_folder = os.path.join(parent_folder, "pdb")
    af_folder = os.path.join(parent_folder, "af")

    process_files(
        pdb_folder,
        af_folder,
        only_lddt=args.lddt,
        only_interfaces=args.interfaces,
        residue_number_alignment=args.residue_number_alignment,
    )

