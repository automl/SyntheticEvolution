import os
import sys
from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from pymol import cmd

# pymol 2.5 (e.g. the synEvo env) lacks the 'ions' selection keyword that pymol 3.x
# provides; every cmd.remove("ions") here is immediately followed by an explicit
# resn MG+K+NA+CA+ZN+CL+SO4 removal, so make the generic call a safe no-op when the
# keyword is absent (re-raising any other error) rather than crashing the scoring.
_pymol_orig_remove = cmd.remove
def _safe_remove(selection, *a, **k):
    try:
        return _pymol_orig_remove(selection, *a, **k)
    except Exception:
        if isinstance(selection, str) and "ion" in selection.lower():
            return None
        raise
cmd.remove = _safe_remove
import tempfile
from structures.rna import RNA
from structures.dna import DNA
from structures.protein import Protein
import re
from scipy.spatial.distance import pdist, squareform
import sqlite3
import json
import traceback
import contextlib
import numpy as np

from database.startConfig import StartConfig
from evaluation.chain_assignment import assigned_by_pred_chain


@contextlib.contextmanager
def _suppress_stdout_fd():
    """
    Redirect the process's stdout file descriptor (fd 1) to /dev/null so that C-level output --
    notably PyMOL's alignment feedback ('Match:', 'MatchAlign:', 'ExecutiveAlign:', 'ExecutiveRMS:',
    'Executive: RMSD ...') -- is silenced together with ordinary Python prints. stderr (fd 2, where
    the tqdm bar lives) and anything printed outside this block are unaffected.
    """
    sys.stdout.flush()
    saved_fd = os.dup(1)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull_fd, 1)
        yield
    finally:
        sys.stdout.flush()
        os.dup2(saved_fd, 1)
        os.close(saved_fd)
        os.close(devnull_fd)

_ALLOWED_PRED_TABLES = frozenset(
    {"pred_protein_rna", "pred_protein_rna_dna", "pred_protein_dna", "pred_rna_rna"}
)
_DEFAULT_EVAL_DB_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "database", "rbpDatabase.db"
)
_eval_metrics_db_path = _DEFAULT_EVAL_DB_PATH


def set_eval_metrics_database_path(path):
    """SQLite path used when writing scores in get_rmsd_tm_lddt_scores (default: bundled rbpDatabase.db)."""
    global _eval_metrics_db_path
    _eval_metrics_db_path = path


def load_predicted_entries_from_table(table_name, db_path, pdb_folder, af_folder):
    """
    Return sorted (exp_db_id, FileName) entries from table that have both reference and predicted CIF on disk.

    Mirrors the idea behind extractAllAtomLDDT: only rows in the chosen prediction table are processed,
    not every file in af/.
    """
    if table_name not in _ALLOWED_PRED_TABLES:
        raise ValueError(
            f"Unsupported prediction table '{table_name}'. "
            f"Use one of: {', '.join(sorted(_ALLOWED_PRED_TABLES))}"
        )

    conn = sqlite3.connect(db_path)
    try:
        cur = conn.cursor()
        cur.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
            (table_name,),
        )
        if not cur.fetchone():
            raise RuntimeError(f"Table not found in database: {table_name}")

        cur.execute(
            f"SELECT exp_db_id, FileName FROM {table_name} "
            "WHERE FileName IS NOT NULL AND TRIM(FileName) != ''"
        )
        rows = cur.fetchall()
    finally:
        conn.close()

    found = []
    missing = []
    seen = set()
    for exp_db_id, file_name_cif_af in rows:
        if not str(file_name_cif_af).endswith(".cif"):
            continue
        uniq_key = (exp_db_id, file_name_cif_af)
        if uniq_key in seen:
            continue
        seen.add(uniq_key)
        try:
            pdb_id = file_name_cif_af.split("_")[1] if "_" in file_name_cif_af else None
            if pdb_id is None:
                missing.append(file_name_cif_af)
                continue
        except (IndexError, AttributeError):
            missing.append(file_name_cif_af)
            continue

        file_name_cif_pdb = f"{pdb_id}.cif"
        pdb_fp = os.path.join(pdb_folder, file_name_cif_pdb)
        af_fp = os.path.join(af_folder, file_name_cif_af)
        if os.path.isfile(pdb_fp) and os.path.isfile(af_fp):
            found.append((exp_db_id, file_name_cif_af))
        else:
            print(
                f"Skipping DB row exp_db_id={exp_db_id}, FileName={file_name_cif_af}: "
                f"ref exists={os.path.isfile(pdb_fp)}, af exists={os.path.isfile(af_fp)}"
            )

    print(
        f"Table {table_name}: {len(found)} entries with paired CIF files "
        f"({len(rows)} DB rows)."
    )
    return sorted(found, key=lambda x: (str(x[0]), str(x[1])))


def count_residues_per_chain(selection):
    chain_residue_counts = {}
    model = cmd.get_model(selection)

    for atom in model.atom:
        chain = atom.chain
        resi = atom.resi
        if chain not in chain_residue_counts:
            chain_residue_counts[chain] = set()
        chain_residue_counts[chain].add(resi)

    for chain in chain_residue_counts:
        chain_residue_counts[chain] = len(chain_residue_counts[chain])
    return chain_residue_counts

def truncate_chain_ids(structure):
    # If chain_id = AAA or XXX
    for model in structure:
        for chain in model:
            if len(chain.id) > 1:
                chain.id = chain.id[0]

def _extract_tm_score_from_text(content):
    """TM-score normalized by Structure_1, rounded to 5dp (USalign's max output precision)."""
    match = re.search(r"TM-score=\s*([0-9.]+)\s*\(normalized by length of Structure_1", content)
    if match:
        return float(match.group(1))
    return None


def _add_column_if_missing(cursor, table_name, column_name, col_type="REAL"):
    """ALTER TABLE ... ADD COLUMN, silently ignoring the 'already exists' case."""
    try:
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} {col_type}")
    except sqlite3.OperationalError as e:
        if "duplicate column name" not in str(e):
            raise


# Per-composition mapping: (target table, non-LDDT columns, LDDT columns), each column with its
# SQLite type. LDDT columns are kept separate so they can be skipped entirely when LDDT computation
# is turned off (the -lddt flag). TEXT for values that may be stored as JSON strings (RMSD lists via
# handle_list_or_value, and the RNA_TM_AllChains matrix); REAL for plain scalar scores. Keep these in
# sync with the UPDATE statements at the end of get_rmsd_tm_lddt_scores.
def get_target_table_and_columns(is_protein, is_rna, is_dna, rna_chain_number=None):
    if is_protein and is_rna and not is_dna:
        return "pred_protein_rna", [
            ("Complex_RMSD", "TEXT"), ("Protein_RMSD", "TEXT"), ("RNA_RMSD", "TEXT"),
            ("Protein_TM", "REAL"), ("RNA_TM", "REAL"), ("Complex_TM", "REAL"),
            ("RNA_TM_AllChains", "TEXT"), ("RNA_chain_assignment", "TEXT"),
        ], [
            ("Protein_LDDT", "REAL"), ("RNA_LDDT", "REAL"), ("Complex_LDDT", "REAL"),
        ]
    if is_protein and is_rna and is_dna:
        return "pred_protein_rna_dna", [
            ("Complex_RMSD", "TEXT"), ("Protein_RMSD", "TEXT"), ("RNA_RMSD", "TEXT"), ("DNA_RMSD", "TEXT"),
            ("Protein_TM", "REAL"), ("RNA_TM", "REAL"), ("DNA_TM", "REAL"),
            ("Complex_TM", "REAL"), ("RNA_TM_AllChains", "TEXT"), ("RNA_chain_assignment", "TEXT"),
        ], [
            ("Protein_LDDT", "REAL"), ("RNA_LDDT", "REAL"), ("DNA_LDDT", "REAL"), ("Complex_LDDT", "REAL"),
        ]
    if is_protein and not is_rna and is_dna:
        return "pred_protein_dna", [
            ("Complex_RMSD", "TEXT"), ("Protein_RMSD", "TEXT"), ("DNA_RMSD", "TEXT"),
            ("Protein_TM", "REAL"), ("DNA_TM", "REAL"), ("Complex_TM", "REAL"),
        ], [
            ("Protein_LDDT", "REAL"), ("DNA_LDDT", "REAL"), ("Complex_LDDT", "REAL"),
        ]
    if not is_protein and is_rna and not is_dna:
        return "pred_rna_rna", [
            ("Complex_RMSD", "TEXT"), ("RNA_RMSD", "TEXT"),
            ("RNA_TM", "REAL"), ("Complex_TM", "REAL"), ("RNA_TM_AllChains", "TEXT"),
            ("RNA_chain_assignment", "TEXT"),
        ], [
            ("RNA_LDDT", "REAL"), ("Complex_LDDT", "REAL"),
        ]
    return None, [], []


def ensure_metric_columns(cursor, is_protein, is_rna, is_dna, rna_chain_number=None,
                          compute_lddt=True, compute_rnacomplex_tm=True):
    """
    Make sure every column the upcoming UPDATE will populate exists, creating any that are missing.
    The target table and column set depend on the complex composition (e.g.
    not is_protein and is_rna and not is_dna and rna_chain_number == 1 -> pred_rna_rna).
    LDDT columns are only ensured when compute_lddt is True; RNAComplex_TM (RNA tables only) only
    when compute_rnacomplex_tm is True.
    Returns the target table name (or None if the composition is unsupported).
    """
    table_name, columns, lddt_columns = get_target_table_and_columns(
        is_protein, is_rna, is_dna, rna_chain_number
    )
    if table_name is None:
        return None
    if compute_lddt:
        columns = columns + lddt_columns
    if compute_rnacomplex_tm and is_rna:
        columns = columns + [("RNAComplex_TM", "REAL")]
    for column_name, col_type in columns:
        _add_column_if_missing(cursor, table_name, column_name, col_type)
    return table_name


def compute_rna_tm_allchains(temp_ref_pdb_path, temp_tgt_pdb_path, af_folder, pdb_id):
    """
    Per-predicted-RNA-chain TM-score, reusing the already-written temp PDBs (no CIF re-parsing).
    Internally an all-vs-all TM matrix is built (TM-score run only for chain pairs of equal residue
    length; other cells stay 0.0), then reduced by a one-to-one maximum-weight assignment (Hungarian):
    every predicted chain is matched to a DISTINCT reference chain, and the RAW assigned TM is kept,
    keyed by the predicted chain id -- same format as RNA_LDDT_AllChains. This avoids the greedy
    per-column max that let several predicted chains claim the same reference chain and over-scored
    multichain / homomeric complexes.

    Returns (values, chain_map):
      values    -> list of single-key dicts, e.g. [{"F": 0.29}, {"G": 0.32}, {"H": 0.33}] (raw TM).
                   Predicted chains left unmatched or with a non-positive TM are skipped.
      chain_map -> {predicted_chain_id: reference_chain_id} for the matched pairs, so the OST LDDT/ICS
                   step can report against the SAME reference partner. Empty if no RNA on either side.
    """
    cmd.reinitialize()
    cmd.load(temp_ref_pdb_path, "reference")
    cmd.load(temp_tgt_pdb_path, "predicted")
    cmd.remove("solvent")
    cmd.remove("organic")
    cmd.remove("ions")
    cmd.remove("resn MG+K+NA+CA+ZN+CL+SO4")

    ref_chains = sorted(count_residues_per_chain("reference and polymer.nucleic").items())
    pred_chains = sorted(count_residues_per_chain("predicted and polymer.nucleic").items())
    if not ref_chains or not pred_chains:
        return [], {}

    # With a single RNA chain on each side there is no ambiguity, so compare them even if the
    # residue counts differ (e.g. missing residues in the experimental structure); USalign aligns
    # them regardless. For multi-chain, keep the equal-length filter to pair matching copies.
    single_chain_mode = len(ref_chains) == 1 and len(pred_chains) == 1

    tm_matrix = [[0.0 for _ in pred_chains] for _ in ref_chains]
    for i, (ref_chain, ref_len) in enumerate(ref_chains):
        for j, (pred_chain, pred_len) in enumerate(pred_chains):
            if not single_chain_mode and ref_len != pred_len:
                continue  # keep as 0.0
            with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb", delete=False) as ref_chain_pdb:
                ref_chain_path = ref_chain_pdb.name
                cmd.save(ref_chain_path, f"reference and chain {ref_chain} and polymer.nucleic")
            with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb", delete=False) as pred_chain_pdb:
                pred_chain_path = pred_chain_pdb.name
                cmd.save(pred_chain_path, f"predicted and chain {pred_chain} and polymer.nucleic")
            try:
                tm_output_path = os.path.join(af_folder, f"{pdb_id}_rna_tm_{ref_chain}_{pred_chain}.txt")
                os.system(
                    f"./USalign/USalign {ref_chain_path} {pred_chain_path} "
                    f"-mol RNA -mm 0 -ter 2 > {tm_output_path}"
                )
                with open(tm_output_path, "r") as fh:
                    tm_score = _extract_tm_score_from_text(fh.read())
                tm_matrix[i][j] = tm_score if tm_score is not None else 0.0
            except Exception as e:
                print(f"Error computing RNA_TM_AllChains for {pdb_id} {ref_chain}->{pred_chain}: {e}")
                tm_matrix[i][j] = 0.0
            finally:
                for p in (ref_chain_path, pred_chain_path):
                    try:
                        os.remove(p)
                    except OSError:
                        pass

    # Reduce the matrix by a one-to-one maximum-weight assignment (Hungarian): each predicted RNA
    # chain is matched to a distinct reference chain and the RAW assigned TM is kept, keyed by the
    # predicted chain id. Also return the {pred_chain: ref_chain} map so the OST step can score
    # LDDT/ICS against the same reference partner.
    out, chain_map = assigned_by_pred_chain(tm_matrix, ref_chains, pred_chains, positive_only=True)
    return out, chain_map


def compute_rna_complex_tm(temp_ref_pdb_path, temp_tgt_pdb_path, af_folder, pdb_id):
    """
    RNA-only complex TM-score: concatenate all RNA chains (dropping protein, DNA, ions, ligands and
    solvent) and run USalign in multimer mode (-mol RNA -mm 1 -ter 0). Reuses the already-written temp
    PDBs. Original script: RNAlddtextractEvalMetrics_complex.py.

    Only meaningful for multi-chain RNA; for a single RNA chain it is identical to RNA_TM.
    Returns a float or None.
    """
    def _save_rna_only(temp_pdb_path, obj_name, out_path):
        cmd.reinitialize()
        cmd.load(temp_pdb_path, obj_name)
        cmd.remove("solvent")
        cmd.remove("organic")
        cmd.remove("ions")
        cmd.remove("resn MG+K+NA+CA+ZN+CL+SO4")
        cmd.remove("polymer.protein")
        cmd.remove("resn DA+DC+DG+DT")  # drop canonical DNA, keep RNA
        cmd.save(out_path, f"{obj_name} and polymer.nucleic")

    with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb", delete=False) as ref_tmp:
        ref_complex_path = ref_tmp.name
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb", delete=False) as tgt_tmp:
        tgt_complex_path = tgt_tmp.name
    try:
        _save_rna_only(temp_ref_pdb_path, "ref_rna_complex", ref_complex_path)
        _save_rna_only(temp_tgt_pdb_path, "tgt_rna_complex", tgt_complex_path)

        tm_output_path = os.path.join(af_folder, f"{pdb_id}_rna_complex_tm_output.txt")
        os.system(
            f"./USalign/USalign {ref_complex_path} {tgt_complex_path} "
            f"-mol RNA -mm 1 -ter 0 > {tm_output_path}"
        )
        with open(tm_output_path, "r") as fh:
            return _extract_tm_score_from_text(fh.read())
    except Exception as e:
        print(f"Error computing RNAComplex_TM for {pdb_id}: {e}")
        return None
    finally:
        for p in (ref_complex_path, tgt_complex_path):
            try:
                os.remove(p)
            except OSError:
                pass
        cmd.reinitialize()


def get_rmsd_tm_lddt_scores(file_name_cif_pdb, file_name_cif_af, compute_lddt=True,
                            compute_rnacomplex_tm=True):
    cif_parser = MMCIFParser(QUIET=True)
    structure_reference = cif_parser.get_structure(file_name_cif_pdb, os.path.join(pdb_folder, file_name_cif_pdb))
    structure_predicted = cif_parser.get_structure(file_name_cif_af, os.path.join(af_folder, file_name_cif_af))
    # truncate_chain_ids(structure_predicted)  # We have a separate script to do this: renameChainsCIF in dataPrep

    # Convert structures to PDB format and save to temporary files
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file_ref:
        pdb_io = PDBIO()
        pdb_io.set_structure(structure_reference)
        pdb_io.save(temp_pdb_file_ref.name)
        temp_pdb_file_ref_path = temp_pdb_file_ref.name

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file_tgt:
        pdb_io.set_structure(structure_predicted)
        pdb_io.save(temp_pdb_file_tgt.name)
        temp_pdb_file_tgt_path = temp_pdb_file_tgt.name

    protein_rmsd = []
    rna_rmsd = []
    dna_rmsd = []
    complex_rmsd = []
    # LDDT scores are optional (skipped when compute_lddt is False); default to None so the
    # separate LDDT UPDATE can be conditionally skipped without NameErrors.
    protein_lddt_score = None
    rna_lddt_score = None
    dna_lddt_score = None
    complex_lddt = None
    # RNAComplex_TM is optional too (skipped when compute_rnacomplex_tm is False / -nRNA_TM).
    rnacomplex_tm = None
    cmd.load(temp_pdb_file_ref_path, 'reference')
    cmd.load(temp_pdb_file_tgt_path, 'predicted')
    cmd.remove('solvent')
    cmd.remove('organic') # Remove all ligands
    cmd.remove('ions') # Remove all ions - not sufficient
    cmd.remove('resn MG')
    cmd.remove('resn K')
    cmd.remove('resn NA')
    cmd.remove('resn CA')
    cmd.remove('resn ZN')
    cmd.remove('resn CL')
    cmd.remove('resn SO4')
    # # Rename 6MZ to A (adenosine)
    # cmd.alter("resn 6MZ", "resn='A'")
    # cmd.alter("resn 6MZ", "name='A'")
    # # Rename 5MC to C (cytidine)
    # cmd.alter("resn 5MC", "resn='C'")
    # cmd.alter("resn 5MC", "name='C'")
    # cmd.remove('hetatm') # cmd.remove('resn MG') to remove mg and other hetatm, problem: also removes 6MZ nucleotide modifications!!!!

    file_pdb = os.path.splitext(file_name_cif_pdb)[0]
    file_af = os.path.splitext(file_name_cif_af)[0]
    pdb_id = file_name_cif_pdb[:4]
    if len(file_pdb) in [4, 6]:
        rna = RNA.get_rna_from_db(id=file_pdb.upper())
        # print(rna.is_rna)
        if rna.is_rna:
            # pdb_rna_chain_id = rna.get_rna_chain_IDs()[0]
            pdb_rna_chain_id = rna.get_longest_chain_id()
            # print("pdb_rna_chain_id!!!", pdb_rna_chain_id)
        protein = Protein.get_protein_from_db(id=file_pdb.upper())
        if protein.is_protein:
            pdb_protein_chain_id = protein.get_longest_chain_id()
        dna = DNA.get_dna_from_db(id=file_pdb.upper())
        if dna.is_dna:
            pdb_dna_chain_id = dna.get_longest_chain_id()
    if len(file_af) > 6:
        parts = file_name_cif_af.split('_')
        rna = RNA.get_rna_from_db(id=parts[1].upper(), file_name=file_name_cif_af)
        if rna.is_rna:
            af_rna_chain_id = rna.get_longest_chain_id()
            # print("af_rna_chain_id", af_rna_chain_id)
        protein = Protein.get_protein_from_db(id=parts[1].upper(), file_name=file_name_cif_af)
        if protein.is_protein:
            af_protein_chain_id = protein.get_longest_chain_id()
        dna = DNA.get_dna_from_db(id=parts[1].upper(), file_name=file_name_cif_af)
        if dna.is_dna:
            af_dna_chain_id = dna.get_longest_chain_id()

    aligned_folder = os.path.join(os.path.dirname(af_folder), "aligned")
    if not os.path.exists(aligned_folder):
        os.makedirs(aligned_folder)

    # chains = cmd.get_chains("reference")
    # print("reference residues count", count_residues_per_chain("reference"))
    # print("predicted residues count", count_residues_per_chain("predicted"))

    def extract_score(content, score_type='TM-score'):
        if score_type == 'RMSD':
            # Extract the RMSD value specifically
            match = re.search(r'Aligned length=\s*\d+,\s*RMSD=\s*([0-9.]+)', content)
            if match:
                score = float(match.group(1))
                # score = round(float(match.group(1)), 2)
                print(f"Extracted RMSD: {score}")
                return score
        elif score_type == 'TM-score':
            # Extract the TM-score value specifically from the reference structure normalization
            match = re.search(r'TM-score=\s*([0-9.]+)\s*\(normalized by length of Structure_1', content)
            if match:
                # USalign prints TM-score to 5 decimals; keep full available precision.
                score = float(match.group(1))
                print(f"Extracted TM-score: {score}")
                return score
        return None

    # def calculate_distances(atoms):
    #     coordinates = [atom.get_coord() for atom in atoms]
    #     return squareform(pdist(coordinates))
    #
    # def lddt(ref_distances, pred_distances, cutoff=0.5):
    #     total_pairs = 0
    #     satisfied_pairs = 0
    #
    #     for i in range(len(ref_distances)):
    #         for j in range(i + 1, len(ref_distances)):
    #             total_pairs += 1
    #             try:
    #                 if abs(ref_distances[i][j] - pred_distances[i][j]) <= cutoff:
    #                     satisfied_pairs += 1
    #             except IndexError:
    #                 # to-do: Handle the case where the index is out of bounds
    #                 print(f"Warning: Index out of bounds at ref_distances[{i}][{j}] or pred_distances[{i}][{j}]")
    #                 continue  # Skip this pair and continue with the next one
    #     return satisfied_pairs / total_pairs if total_pairs > 0 else None

    parser = PDBParser(QUIET=True)

    # Define modified nucleotide and protein residue mappings to canonical equivalents
    modified_residues_nucleotide = {
        '6MZ': 'A', 'PSU': 'U', '5MC': 'C', 'OMC': 'C', '4OC': 'C', '5MU': 'U',
        'OMU': 'U', 'UR3': 'U', 'A2M': 'A', 'MA6': 'A', '2MG': 'G', 'OMG': 'G',
        '7MG': 'G', 'RSQ': 'G', '5CM': 'C', 'C34': 'C', '5HC': 'C', '6OG': 'G',
        '6MA': 'A', '1CC': 'C', '8OG': 'G', '5FC': 'C', '3DR': 'T'
    }
    
    modified_residues_protein = {
        'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'NEP': 'H', 'HIP': 'H', 'ALY': 'K',
        'MLY': 'K', 'M3L': 'K', 'MLZ': 'K', '2MR': 'R', 'AGM': 'R', 'MCS': 'C',
        'HYP': 'P', 'HY3': 'H', 'LYZ': 'K', 'AHB': 'A', 'P1L': 'P', 'SNN': 'S',
        'SNC': 'C', 'TRF': 'W', 'KCR': 'K', 'CIR': 'R', 'YHA': 'Y'
    }
    
    def normalize_sequence(sequence, is_protein=False):
        """Replace modified nucleotides/protein residues with their canonical equivalents and remove parentheses."""
        if not sequence:
            return sequence
        
        # Choose the appropriate mapping
        modified_residues = modified_residues_protein if is_protein else modified_residues_nucleotide
        
        # Replace known modifications with their canonical equivalents
        # Pattern: (MOD) where MOD is a known modification code
        normalized = sequence
        for mod, base in modified_residues.items():
            # Replace (MOD) with the canonical base
            normalized = re.sub(rf'\({re.escape(mod)}\)', base, normalized, flags=re.IGNORECASE)
        
        # Remove any remaining parentheses and their contents (unknown modifications)
        normalized = re.sub(r'\([^)]*\)', '', normalized)
        
        return normalized
    
    def contains_unhandled_modifications(sequence, is_protein=False):
        """Check if sequence contains parentheses (unhandled modifications) after normalization."""
        normalized = normalize_sequence(sequence, is_protein=is_protein)
        pattern = r"\(.*?\)"
        if re.search(pattern, normalized):
            return True  # Still has parentheses after normalization
        return False

    def rename_protein_chain(pdb_protein_chain_id, pdb_rna_chain_id):
        """Ensures protein chain ID is a single letter, avoiding conflicts with RNA chain IDs."""
        if len(pdb_protein_chain_id) == 1:
            return pdb_protein_chain_id

        available_letters = set(string.ascii_uppercase) - set(pdb_rna_chain_id)
        if not available_letters:
            raise ValueError("No available chain IDs left!")
        new_chain_id = sorted(available_letters)[0]
        print(f"Renaming {pdb_protein_chain_id} to {new_chain_id}")
        return new_chain_id

    if protein.is_protein:
        if rna.is_rna:
            # print("pdb_protein_chain_id, pdb_rna_chain_id", pdb_protein_chain_id, pdb_rna_chain_id)
            rename_protein_chain(pdb_protein_chain_id, pdb_rna_chain_id)
            # print("pdb_protein_chain_id, pdb_rna_chain_id", pdb_protein_chain_id, pdb_rna_chain_id)
    # Alter the chain identifiers in the predicted structure to match the reference structure
    # https://pymol.org/dokuwiki/doku.php?id=command:alter
    # https://pymol.org/pymol-command-ref.html#align
    #     cmd.alter(f'predicted and chain {af_protein_chain_id}', f'chain="{pdb_protein_chain_id}"')
        alignment_result = cmd.align(f'predicted and chain {af_protein_chain_id} and polymer.protein and name CA',
                                     f'reference and chain {pdb_protein_chain_id} and polymer.protein and name CA',
                                     quiet=0, reset=1)  # , object='alignment_rna')
        #rmsd_protein = round(alignment_result[0], 2)
        rmsd_protein = alignment_result[0]
        print(f"RMSD for protein structure from alignment = {rmsd_protein} Å")

        if rmsd_protein is not None:
            protein_rmsd.append(rmsd_protein)
        else:
            protein_rmsd.append("None")
        # num_protein_atoms = cmd.count_atoms(f"predicted and chain {af_protein_chain_id} and polymer.protein")

        protein_aligned_reference_pdb = os.path.join(aligned_folder, f'{pdb_id}_pdb_aligned_protein.pdb')
        protein_aligned_predicted_pdb = os.path.join(aligned_folder, f'{pdb_id}_af_aligned_protein.pdb')
        # print(f"Reference chain ID: {pdb_protein_chain_id}")
        # print(f"Predicted chain ID: {af_protein_chain_id}")
        cmd.save(protein_aligned_reference_pdb, f'reference and chain {pdb_protein_chain_id} and polymer.protein')
        cmd.save(protein_aligned_predicted_pdb, f'predicted and chain {af_protein_chain_id} and polymer.protein')

        if contains_unhandled_modifications(protein.protein_sequence, is_protein=True) == False:
            # Obtain TM-score from https://zhanggroup.org/TM-score/
            tm_output_path_protein = os.path.join(af_folder, f'{pdb_id}_protein_tm_output.txt')
            # -mol prot: Specifies that we're aligning protein molecules
            # -mm 0: For monomeric structures (single chain)
            # -ter 2: Only align the first chain (default for monomeric cases)
            os.system(
                f"./USalign/USalign {protein_aligned_reference_pdb} {protein_aligned_predicted_pdb} "
                f"-mol prot -mm 0 -ter 2 > {tm_output_path_protein}"
            )
            with open(tm_output_path_protein, 'r') as file:
                content = file.read()
                tm_score_protein = extract_score(content, 'TM-score')
                rmsd_score_protein = extract_score(content, 'RMSD')
                if rmsd_score_protein is not None:
                    protein_rmsd.append(rmsd_score_protein)
                else:
                    protein_rmsd.append("None")
        else:
            tm_score_protein = None
            rmsd_score_protein = None

        # Calculate protein LDDT using the complex LDDT function (optional)
        if compute_lddt:
            protein_lddt_score = calculate_complex_lddt(protein_aligned_reference_pdb, protein_aligned_predicted_pdb)
            print(f"Protein lDDT Score: {protein_lddt_score:.2f}")

        # alignment_result2 = cmd.align(f'predicted and chain {af_protein_chain_id} and polymer.protein',
        #                              f'reference and chain {pdb_protein_chain_id} and polymer.protein',
        #                              quiet=0)  # , object='alignment_rna', reset=1)
        # rmsd_protein3 = round(alignment_result2[0], 2)
        # print(f"All atom RMSD for protein structure from alignment = {rmsd_protein3} Å")
        #
        # if rmsd_protein3 is not None:
        #     protein_rmsd.append(rmsd_protein3)
        # else:
        #     protein_rmsd.append("None")

    if rna.is_rna:
        if protein.is_protein: # Problem: 7wl0 still contained protein chains despite selection
            cmd.remove('polymer.protein')
            rename_protein_chain(pdb_rna_chain_id, pdb_protein_chain_id)
        # cmd.alter(f'predicted and chain {af_rna_chain_id}', f'chain="{pdb_rna_chain_id}"')
        alignment_result = cmd.align(f'predicted and chain {af_rna_chain_id} and polymer.nucleic and name P',
                  f'reference and chain {pdb_rna_chain_id} and polymer.nucleic and name P',
                  quiet=0, reset=1) #, object='alignment_rna', reset=1)
        # rmsd_rna = round(alignment_result[0], 2)
        rmsd_rna = alignment_result[0]
        print(f"RMSD for RNA structure from alignment = {rmsd_rna} Å")
        if rmsd_rna is not None:
            rna_rmsd.append(rmsd_rna)
        else:
            rna_rmsd.append("None")
        # num_rna_atoms = cmd.count_atoms(f"predicted and chain {af_rna_chain_id} and polymer.nucleic")

        rna_aligned_reference_pdb = os.path.join(aligned_folder, f'{pdb_id}_pdb_aligned_rna.pdb')
        rna_aligned_predicted_pdb = os.path.join(aligned_folder, f'{pdb_id}_af_aligned_rna.pdb')

        # Ensure chain IDs match before saving
        # cmd.alter(f'reference and chain {pdb_rna_chain_id}', f'chain="{af_rna_chain_id}"')
        # cmd.save(rna_aligned_reference_pdb, f'reference and chain {af_rna_chain_id}')
        # cmd.save(rna_aligned_predicted_pdb, f'predicted and chain {af_rna_chain_id}')

        cmd.save(rna_aligned_reference_pdb, f'reference and chain {pdb_rna_chain_id}')
        cmd.save(rna_aligned_predicted_pdb, f'predicted and chain {af_rna_chain_id}')

        if contains_unhandled_modifications(rna.rna_sequence, is_protein=False) == False: # RNA with known modifications will be processed
            # Obtain TM-score from https://zhanggroup.org/TM-score/
            # https://aideepmed.com/US-align/help/
            # Try compiling: g + + -O3 - ffast - math - o RNAalign RNAalign.cpp
            # MAC: The file includes both stdlib.h and malloc.h. On macOS, we should remove the malloc.h include since stdlib.h already provides the necessary memory allocation functions.
            tm_output_path_rna = os.path.join(af_folder, f'{pdb_id}_rna_tm_output.txt')
            # os.system(f'./RNAalign/RNAalign {rna_aligned_reference_pdb} {rna_aligned_predicted_pdb} > {tm_output_path_rna}')
            # print("C1")
            # os.system(
            #     f"./USalign/USalign {rna_aligned_reference_pdb} {rna_aligned_predicted_pdb} "
            #     f"-mol RNA -mm 0 -atom \" C1'\" -ter 2 > {tm_output_path_rna}"
            # )
            #
            # with open(tm_output_path_rna, 'r') as file:
            #     content = file.read()
            #     tm_score_rna1 = extract_score(content, 'TM-score')
            #
            # print("C4")
            # os.system(
            #     f"./USalign/USalign {rna_aligned_reference_pdb} {rna_aligned_predicted_pdb} "
            #     f"-mol RNA -mm 0 -atom \" C4'\" -ter 2 > {tm_output_path_rna}"
            # )
            #
            # with open(tm_output_path_rna, 'r') as file:
            #     content = file.read()
            #     tm_score_rna2 = extract_score(content, 'TM-score')
            #
            # print("P")
            # os.system(
            #     f"./USalign/USalign {rna_aligned_reference_pdb} {rna_aligned_predicted_pdb} "
            #     f"-mol RNA -mm 0 -atom \" P  \" -ter 2 > {tm_output_path_rna}"
            # )
            #
            # with open(tm_output_path_rna, 'r') as file:
            #     content = file.read()
            #     tm_score_rna3 = extract_score(content, 'TM-score')

            # Default is " C3'" for RNA/DNA
            os.system(
                f"./USalign/USalign {rna_aligned_reference_pdb} {rna_aligned_predicted_pdb} "
                f"-mol RNA -mm 0 -ter 2 > {tm_output_path_rna}"
            )

            with open(tm_output_path_rna, 'r') as file:
                content = file.read()
                tm_score_rna = extract_score(content, 'TM-score')
                rmsd_score_rna = extract_score(content, 'RMSD')
                if rmsd_score_rna is not None:
                    rna_rmsd.append(rmsd_score_rna)
                else:
                    rna_rmsd.append("None")
        else:
            tm_score_rna = None
            rmsd_score_rna = None

        # Calculate RNA LDDT using the complex LDDT function (optional)
        if compute_lddt:
            rna_lddt_score = calculate_complex_lddt(rna_aligned_reference_pdb, rna_aligned_predicted_pdb)
            # rna_lddt_score = 0
            print(f"RNA lDDT Score: {rna_lddt_score}")

        # rmsd_rna2 = cmd.rms_cur(
        #     "predicted and polymer.nucleic",
        #     "reference and polymer.nucleic"
        # )
        #

        alignment_result2 = cmd.align(f'predicted and chain {af_rna_chain_id} and polymer.nucleic',
                                     f'reference and chain {pdb_rna_chain_id} and polymer.nucleic',
                                     quiet=0, reset=1)  # , object='alignment_rna', reset=1)
        # rmsd_rna3 = round(alignment_result2[0], 2)
        rmsd_rna3 = alignment_result2[0]
        print(f"All atom RMSD for RNA structure from alignment = {rmsd_rna3} Å")

        if rmsd_rna3 is not None:
            rna_rmsd.append(rmsd_rna3)
        else:
            rna_rmsd.append("None")

        # print("rms_cur rmsd_rna predicted and polymer.nucleic ", rmsd_rna2, "normal pymol", rmsd_rna, "TM",
        #       rmsd_score_rna, "alignment_result2", rmsd_rna3)

    if dna.is_dna:
        # cmd.alter(f'predicted and chain {af_dna_chain_id}', f'chain="{pdb_dna_chain_id}"') # Change chain label?
        # It seems that the other chain labelled as {pdb_dna_chain_id} is considered instead
        alignment_result = cmd.align(f'predicted and chain {af_dna_chain_id} and polymer.nucleic and name P',
                  f'reference and chain {pdb_dna_chain_id} and polymer.nucleic and name P',
                  quiet=0) #, object='alignment_rna', reset=1)
        # rmsd_dna = round(alignment_result[0], 2)
        rmsd_dna = alignment_result[0]
        print(f"RMSD for DNA structure from alignment = {rmsd_dna} Å")
        dna_rmsd.append(rmsd_dna)
        num_dna_atoms = cmd.count_atoms(f"predicted and chain {af_dna_chain_id} and polymer.nucleic")

        dna_aligned_reference_pdb = os.path.join(aligned_folder, f'{pdb_id}_pdb_aligned_dna.pdb')
        dna_aligned_predicted_pdb = os.path.join(aligned_folder, f'{pdb_id}_af_aligned_dna.pdb')
        cmd.save(dna_aligned_reference_pdb, f'reference and chain {pdb_dna_chain_id}')
        cmd.save(dna_aligned_predicted_pdb, f'predicted and chain {af_dna_chain_id}')

        if contains_unhandled_modifications(protein.protein_sequence, is_protein=True) == False:
            # Obtain TM-score from https://zhanggroup.org/TM-score/
            tm_output_path_dna = os.path.join(af_folder, f'{pdb_id}_dna_tm_output.txt')
            os.system(
                f"./USalign/USalign {dna_aligned_reference_pdb} {dna_aligned_predicted_pdb} "
                f"-mol DNA -mm 0 -ter 2 > {tm_output_path_dna}"
            )
            with open(tm_output_path_dna, 'r') as file:
                content = file.read()
                tm_score_dna = extract_score(content, 'TM-score')
                rmsd_score_dna = extract_score(content, 'RMSD')
                if rmsd_score_dna:
                    dna_rmsd.append(rmsd_score_dna)
        else:
            tm_score_dna = None
            rmsd_score_dna = None

        # Calculate DNA LDDT using the complex LDDT function (optional)
        if compute_lddt:
            dna_lddt_score = calculate_complex_lddt(dna_aligned_reference_pdb, dna_aligned_predicted_pdb)
            print(f"DNA lDDT Score: {dna_lddt_score:.2f}")

        alignment_result2 = cmd.align(f'predicted and chain {af_dna_chain_id} and polymer.nucleic',
                                     f'reference and chain {pdb_dna_chain_id} and polymer.nucleic',
                                     quiet=0, reset=1)  # , object='alignment_rna')
        rmsd_dna3 = alignment_result2[0]
        # rmsd_dna3 = round(alignment_result2[0], 2)
        print(f"RMSD for DNA structure from alignment = {rmsd_dna3} Å")
        dna_rmsd.append(rmsd_dna3)

    conn = sqlite3.connect(_eval_metrics_db_path)
    cursor = conn.cursor()

    # Make sure all columns the per-composition UPDATE below will write actually exist
    # (LDDT columns only when LDDT is computed; RNAComplex_TM only when it is computed).
    ensure_metric_columns(
        cursor,
        protein.is_protein, rna.is_rna, dna.is_dna,
        rna.get_rna_chain_number() if rna.is_rna else None,
        compute_lddt=compute_lddt,
        compute_rnacomplex_tm=compute_rnacomplex_tm,
    )

    def handle_list_or_value(item):
        if isinstance(item, list):
            if len(item) == 1:
                return item[0]
            return json.dumps(item)
        return item

    cmd.reinitialize()
    # Calculate complex RMSD and LDDT once for all cases
    cmd.load(temp_pdb_file_ref_path, 'referenceA')
    cmd.load(temp_pdb_file_tgt_path, 'predictedA')
    # Match interactive PyMOL: align polymer only, with solvents/ligands/ions stripped
    for sel in ('solvent', 'organic', 'ions'):
        cmd.remove(sel)
    for ion_resn in ('MG', 'K', 'NA', 'CA', 'ZN', 'CL', 'SO4'):
        cmd.remove(f'resn {ion_resn}')
    # chains = cmd.get_chains('referenceA')
    # print("Chains in referenceA:", chains)
    # chains = cmd.get_chains('predictedA')
    # print("Chains in referenceA:", chains)

    if not protein.is_protein and rna.is_rna and not dna.is_dna and rna.get_rna_chain_number() == 1:
        complex_lddt = 0
        complex_rmsd = 0
        tm_score_complex = 0
    else:
        # # Calculate LDDT for the entire complex (optional)
        if compute_lddt:
            complex_lddt = calculate_complex_lddt(temp_pdb_file_ref_path, temp_pdb_file_tgt_path)
            print(f"Complex LDDT = {complex_lddt}")

        alignment_result = cmd.align('predictedA and polymer', 'referenceA and polymer', quiet=0, reset=1)
        # pymol_rmsd = round(alignment_result[0], 2)
        pymol_rmsd = alignment_result[0]
        if pymol_rmsd is not None:
            complex_rmsd.append(pymol_rmsd)
        else:
            complex_rmsd.append("None")

        print(f"RMSD after alignment = {pymol_rmsd} Å")

        tm_output_path_complex = os.path.join(af_folder, f'{pdb_id}_complex_tm_output.txt')
        os.system(
            f"./USalign/USalign {temp_pdb_file_ref_path} {temp_pdb_file_tgt_path} "
            f"-mol auto -mm 1 -ter 0 > {tm_output_path_complex}"
        )

        with open(tm_output_path_complex, 'r') as file:
            content = file.read()
            tm_score_complex = extract_score(content, 'TM-score')
            rmsd_score_complex = extract_score(content, 'RMSD')

            if rmsd_score_complex is not None:
                complex_rmsd.append(rmsd_score_complex)
            else:
                complex_rmsd.append("None")

    # RNA_TM_AllChains: assigned TM per predicted RNA chain, keyed by predicted chain id, e.g.
    # [{"F": 0.29}, {"G": 0.32}] (raw values from a one-to-one assignment; same keyed-dict format as
    # RNA_LDDT_AllChains). RNA_chain_assignment stores the {pred_chain: ref_chain} map so the OST
    # LDDT/ICS step reports against the same reference partner.
    rna_tm_allchains_json = None
    rna_chain_assignment_json = None
    if rna.is_rna:
        rna_chain_map = {}
        if rna.get_rna_chain_number() == 1:
            # Single RNA chain: reuse RNA_TM (no recompute), keyed by the predicted RNA chain id.
            # Assignment is trivial (1x1); the OST step falls back to its own single-chain match.
            pred_chain = locals().get("af_rna_chain_id")
            if pred_chain and tm_score_rna is not None and tm_score_rna > 0:
                rna_tm_allchains = [{pred_chain: tm_score_rna}]
            else:
                rna_tm_allchains = []
        else:
            rna_tm_allchains, rna_chain_map = compute_rna_tm_allchains(
                temp_pdb_file_ref_path, temp_pdb_file_tgt_path, af_folder, pdb_id
            )
        rna_tm_allchains_json = json.dumps(rna_tm_allchains)
        rna_chain_assignment_json = json.dumps(rna_chain_map) if rna_chain_map else None
        print(f"RNA_TM_AllChains = {rna_tm_allchains_json}")
        print(f"RNA_chain_assignment = {rna_chain_assignment_json}")

    # RNAComplex_TM: RNA-only complex TM. Single RNA chain -> identical to RNA_TM;
    # >1 RNA chains -> dedicated multimer alignment of all RNA chains concatenated.
    if rna.is_rna and compute_rnacomplex_tm:
        if rna.get_rna_chain_number() == 1:
            rnacomplex_tm = float(tm_score_rna) if tm_score_rna is not None else None
        elif not contains_unhandled_modifications(rna.rna_sequence, is_protein=False):
            rnacomplex_tm = compute_rna_complex_tm(
                temp_pdb_file_ref_path, temp_pdb_file_tgt_path, af_folder, pdb_id
            )
        print(f"RNAComplex_TM = {rnacomplex_tm}")

    # Handle all RMSD values consistently
    complex_rmsd = handle_list_or_value(complex_rmsd)
    protein_rmsd = handle_list_or_value(protein_rmsd)
    rna_rmsd = handle_list_or_value(rna_rmsd)
    dna_rmsd = handle_list_or_value(dna_rmsd)

    if protein.is_protein and rna.is_rna and not dna.is_dna:
        cursor.execute('''UPDATE pred_protein_rna
                          SET Complex_RMSD = ?, Protein_RMSD = ?, RNA_RMSD = ?,
                              Protein_TM = ?, RNA_TM = ?, Complex_TM = ?,
                              RNA_TM_AllChains = ?, RNA_chain_assignment = ?
                          WHERE FileName = ?''',
                       (complex_rmsd, protein_rmsd, rna_rmsd,
                        tm_score_protein, tm_score_rna, tm_score_complex,
                        rna_tm_allchains_json, rna_chain_assignment_json, file_name_cif_af))

    if protein.is_protein and rna.is_rna and dna.is_dna:
        cursor.execute('''UPDATE pred_protein_rna_dna
                                  SET Complex_RMSD = ?, Protein_RMSD = ?, RNA_RMSD = ?, DNA_RMSD = ?,
                                  Protein_TM = ?, RNA_TM = ?, DNA_TM = ?,
                                  Complex_TM = ?, RNA_TM_AllChains = ?, RNA_chain_assignment = ?
                                  WHERE FileName = ?''',
                       (complex_rmsd, protein_rmsd, rna_rmsd, dna_rmsd,
                        tm_score_protein, tm_score_rna, tm_score_dna, tm_score_complex,
                        rna_tm_allchains_json, rna_chain_assignment_json, file_name_cif_af))

    if protein.is_protein and not rna.is_rna and dna.is_dna:
        cursor.execute('''UPDATE pred_protein_dna
                                  SET Complex_RMSD = ?, Protein_RMSD = ?, DNA_RMSD = ?,
                                  Protein_TM = ?, DNA_TM = ?, Complex_TM = ?
                                  WHERE FileName = ?''',
                       (complex_rmsd, protein_rmsd, dna_rmsd, tm_score_protein,
                        tm_score_dna, tm_score_complex, file_name_cif_af))

    if not protein.is_protein and rna.is_rna and not dna.is_dna:
        cursor.execute('''UPDATE pred_rna_rna
                                  SET Complex_RMSD = ?, RNA_RMSD = ?, RNA_TM = ?, Complex_TM = ?,
                                  RNA_TM_AllChains = ?, RNA_chain_assignment = ?
                                  WHERE FileName = ?''',
                       (complex_rmsd, rna_rmsd,
                        tm_score_rna, tm_score_complex, rna_tm_allchains_json,
                        rna_chain_assignment_json, file_name_cif_af))

    if not protein.is_protein and rna.is_rna and not dna.is_dna and rna.get_rna_chain_number() == 1:
        cursor.execute('''UPDATE pred_rna_rna
                                  SET RNA_RMSD = ?, RNA_TM = ?
                                  WHERE FileName = ?''',
                       (rna_rmsd, tm_score_rna, file_name_cif_af))

    # LDDT scores are written separately so the whole block can be skipped when LDDT is off.
    if compute_lddt:
        lddt_value_by_column = {
            "Protein_LDDT": protein_lddt_score,
            "RNA_LDDT": rna_lddt_score,
            "DNA_LDDT": dna_lddt_score,
            "Complex_LDDT": complex_lddt,
        }
        lddt_table, _non_lddt_cols, lddt_cols = get_target_table_and_columns(
            protein.is_protein, rna.is_rna, dna.is_dna,
            rna.get_rna_chain_number() if rna.is_rna else None,
        )
        if lddt_table is not None and lddt_cols:
            set_clause = ", ".join(f"{col} = ?" for col, _ in lddt_cols)
            values = [lddt_value_by_column[col] for col, _ in lddt_cols]
            values.append(file_name_cif_af)
            cursor.execute(
                f"UPDATE {lddt_table} SET {set_clause} WHERE FileName = ?",
                tuple(values),
            )

    # RNAComplex_TM written separately (RNA tables only) so it can be skipped with -nRNA_TM.
    if compute_rnacomplex_tm and rna.is_rna:
        rc_table, _rc_non_lddt, _rc_lddt = get_target_table_and_columns(
            protein.is_protein, rna.is_rna, dna.is_dna, rna.get_rna_chain_number(),
        )
        if rc_table is not None:
            cursor.execute(
                f"UPDATE {rc_table} SET RNAComplex_TM = ? WHERE FileName = ?",
                (rnacomplex_tm, file_name_cif_af),
            )

    conn.commit()
    conn.close()

    os.remove(temp_pdb_file_ref_path)
    os.remove(temp_pdb_file_tgt_path)
    cmd.reinitialize()

def backfill_rna_tm_allchains(file_name_cif_pdb, file_name_cif_af, table_name):
    """
    Recompute ONLY RNA_TM_AllChains for one (reference, predicted) pair and write just that column
    to `table_name` (matched by FileName). Used to backfill old databases; the full pipeline already
    computes RNA_TM_AllChains inline in get_rmsd_tm_lddt_scores, so this is not needed going forward.

    Uses the same compute_rna_tm_allchains() as the main path (all-vs-all matrix reduced to the best
    TM per predicted RNA chain); single- and multi-chain RNA are handled uniformly. Rows with no RNA
    yield an empty list.
    """
    pdb_id = file_name_cif_pdb[:4]
    cif_parser = MMCIFParser(QUIET=True)
    structure_reference = cif_parser.get_structure(
        file_name_cif_pdb, os.path.join(pdb_folder, file_name_cif_pdb)
    )
    structure_predicted = cif_parser.get_structure(
        file_name_cif_af, os.path.join(af_folder, file_name_cif_af)
    )
    pdb_io = PDBIO()
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_ref:
        pdb_io.set_structure(structure_reference)
        pdb_io.save(temp_ref.name)
        ref_tmp = temp_ref.name
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_tgt:
        pdb_io.set_structure(structure_predicted)
        pdb_io.save(temp_tgt.name)
        tgt_tmp = temp_tgt.name
    try:
        rna_tm_allchains, rna_chain_map = compute_rna_tm_allchains(ref_tmp, tgt_tmp, af_folder, pdb_id)
        rna_chain_assignment_json = json.dumps(rna_chain_map) if rna_chain_map else None
        conn = sqlite3.connect(_eval_metrics_db_path)
        cursor = conn.cursor()
        try:
            _add_column_if_missing(cursor, table_name, "RNA_TM_AllChains", "TEXT")
            _add_column_if_missing(cursor, table_name, "RNA_chain_assignment", "TEXT")
            cursor.execute(
                f"UPDATE {table_name} SET RNA_TM_AllChains = ?, RNA_chain_assignment = ? WHERE FileName = ?",
                (json.dumps(rna_tm_allchains), rna_chain_assignment_json, file_name_cif_af),
            )
            conn.commit()
        finally:
            conn.close()
        print(f"RNA_TM_AllChains = {json.dumps(rna_tm_allchains)}")
    finally:
        for p in (ref_tmp, tgt_tmp):
            try:
                os.remove(p)
            except OSError:
                pass
        cmd.reinitialize()


def process_files(pdb_path, af_path, entries, compute_lddt=True, compute_rnacomplex_tm=True,
                  only_rna_tm_allchains=False, table_name=None):
    """
    Extract RMSD/TM/LDDT and update the DB.

    entries is an iterable of (exp_db_id, FileName) table rows (from the config-selected table);
    only those rows are processed.
    If compute_lddt is False, LDDT scores are neither computed nor written.
    If compute_rnacomplex_tm is False, RNAComplex_TM is neither computed nor written.
    If only_rna_tm_allchains is True, ONLY RNA_TM_AllChains is (re)computed and written (backfill mode).
    """
    from tqdm import tqdm
    entries = list(entries)
    for exp_db_id, file_name_cif_af in tqdm(entries, total=len(entries), desc="structures", unit="struct"):
        if not str(file_name_cif_af).endswith(".cif"):
            continue
        try:
            pdb_id = file_name_cif_af.split("_")[1]
            file_name_cif_pdb = f"{pdb_id}.cif"
            pdb_file_path = os.path.join(pdb_path, file_name_cif_pdb)
            af_file_path = os.path.join(af_path, file_name_cif_af)

            if not os.path.exists(af_file_path):
                print(f"Warning: predicted CIF missing, skipping {file_name_cif_af}")
                continue

            # Check if PDB file exists
            if not os.path.exists(pdb_file_path):
                print(f"Warning: No corresponding file found in PDB folder for {pdb_id}")
                continue

            # Try to get scores, skip file if it fails
            try:
                # This function (and PyMOL) emit a lot of per-file output; suppress it at the fd level
                # so the tqdm bar stays clean. Errors below are printed outside the block and stay visible.
                with _suppress_stdout_fd():
                    if only_rna_tm_allchains:
                        backfill_rna_tm_allchains(file_name_cif_pdb, file_name_cif_af, table_name)
                    else:
                        get_rmsd_tm_lddt_scores(
                            file_name_cif_pdb, file_name_cif_af,
                            compute_lddt=compute_lddt, compute_rnacomplex_tm=compute_rnacomplex_tm,
                        )
            except Exception as e:
                tb = traceback.format_exc()
                print(f"Error calculating scores for {pdb_id}: {e}")
                print(f"Traceback:\n{tb}")
                continue

        except Exception as e:
            print(f"Error processing {file_name_cif_af}: {e}")
            continue

def calculate_complex_lddt(reference_path, predicted_path):
    """Calculate LDDT for the entire complex using all heavy atoms from polymers only."""
    from Bio.PDB import PDBParser
    import numpy as np
    from Bio.PDB.Polypeptide import is_aa

    def get_polymer_atom_coords(structure):
        """Get coordinates of 'CA' atoms for amino acids and 'C3\'' atoms for nucleotides."""
        coords = []
        mask = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        for atom in residue:
                            if atom.get_name() == 'CA':
                                coords.append(atom.get_coord())
                                mask.append(1.0)
                    else:
                        for atom in residue:
                            if atom.get_name() == "C3'":
                                coords.append(atom.get_coord())
                                mask.append(1.0)

        return np.array(coords), np.array(mask).reshape(-1, 1)

    def calculate_lddt_score(pred_points, true_points, true_mask, cutoff=15.0):
        """Calculate LDDT score following the simplified approach.
        https://github.com/google-deepmind/alphafold/blob/7c9114c8423ac9db981d8365168464bab09b3e54/alphafold/model/lddt.py#L19
        https://git.scicore.unibas.ch/schwede/openstructure/-/blob/master/modules/mol/alg/pymod/lddt.py
        https://github.com/MDAnalysis/mdanalysis/issues/4134
        https://openstructure.org/docs/2.4/mol/alg/lddt/
        https://swissmodel.expasy.org/lddt
        Mariani, V., Biasini, M., Barbato, A. & Schwede, T. lDDT: A local
        superposition-free score for comparing protein structures and models using
        distance difference tests. Bioinformatics 29, 2722–2728 (2013).
        """
        # Compute distance matrices
        dmat_true = np.sqrt(1e-10 + np.sum(
            (true_points[:, None] - true_points[None, :])**2, axis=-1))
        
        dmat_predicted = np.sqrt(1e-10 + np.sum(
            (pred_points[:, None] - pred_points[None, :])**2, axis=-1))
        
        # Create mask for valid distances to compare
        dists_to_score = (
            (dmat_true < cutoff).astype(np.float32) * true_mask *
            true_mask.T * (1. - np.eye(dmat_true.shape[0]))
        )
        
        # Calculate absolute differences
        dist_l1 = np.abs(dmat_true - dmat_predicted)
        
        # Score using the standard LDDT thresholds
        score = 0.25 * (
            (dist_l1 < 0.5).astype(np.float32) +
            (dist_l1 < 1.0).astype(np.float32) +
            (dist_l1 < 2.0).astype(np.float32) +
            (dist_l1 < 4.0).astype(np.float32)
        )
        
        # Normalize
        norm = 1.0 / (1e-10 + np.sum(dists_to_score))
        final_score = norm * (1e-10 + np.sum(dists_to_score * score))
        
        return final_score

    try:
        # # Parse structures
        parser = PDBParser(QUIET=True)
        ref_structure = parser.get_structure('reference', reference_path)
        pred_structure = parser.get_structure('predicted', predicted_path)

        # Get coordinates and masks for all heavy atoms
        ref_coords, ref_mask = get_polymer_atom_coords(ref_structure)
        pred_coords, _ = get_polymer_atom_coords(pred_structure)

        if len(ref_coords) == 0 or len(pred_coords) == 0:
            print("No backbone atoms found in one or both structures")
            return None

        # Ensure equal number of points
        min_length = min(len(ref_coords), len(pred_coords))
        ref_coords = ref_coords[:min_length]
        pred_coords = pred_coords[:min_length]
        ref_mask = ref_mask[:min_length]

        print(f"Calculating LDDT using {min_length} backbone atoms from polymers")

        # Calculate LDDT score
        lddt_score = calculate_lddt_score(pred_coords, ref_coords, ref_mask)
        
        if lddt_score is not None:
            # return round(lddt_score, 2)
            return lddt_score
        return None
    
    except Exception as e:
        print(f"Error in LDDT calculation: {str(e)}")
        return None

def _print_usage(script):
    print(f"Usage:")
    print(f"  python {script} [-lddt] [-nRNA_TM] [-RNA_TM_AllChains]")
    print(f"    Table and database/paths all come from config.json (classification, base_data_path, db_path).")
    print()
    print(f"Options:")
    print(f"  -h, --help         Show this help message and exit.")
    print(f"  -lddt              Skip LDDT computation (Protein/RNA/DNA/Complex LDDT not computed or written).")
    print(f"  -nRNA_TM           Skip RNAComplex_TM computation (not computed or written).")
    print(f"  -RNA_TM_AllChains  Backfill mode: (re)compute and write ONLY RNA_TM_AllChains (nothing else).")


if __name__ == "__main__":
    script = os.path.basename(sys.argv[0])
    args = sys.argv[1:]
    # -h/--help: show available options and exit successfully.
    if "-h" in args or "--help" in args:
        _print_usage(script)
        sys.exit(0)
    # Optional flags (may appear anywhere): when present the corresponding metric is NOT computed/written.
    compute_lddt = "-lddt" not in args
    compute_rnacomplex_tm = "-nRNA_TM" not in args
    only_rna_tm_allchains = "-RNA_TM_AllChains" in args
    args = [a for a in args if a not in ("-lddt", "-nRNA_TM", "-RNA_TM_AllChains")]
    if only_rna_tm_allchains:
        print("RNA_TM_AllChains-only backfill mode: only RNA_TM_AllChains is (re)computed and written.")
    elif not compute_lddt:
        print("LDDT computation disabled (-lddt): Protein_LDDT/RNA_LDDT/DNA_LDDT/Complex_LDDT will be skipped.")
    if not compute_rnacomplex_tm and not only_rna_tm_allchains:
        print("RNAComplex_TM computation disabled (-nRNA_TM): RNAComplex_TM will be skipped.")
    if args:
        print(f"Error: unrecognized argument(s): {' '.join(args)}\n")
        _print_usage(script)
        sys.exit(1)

    # Config-driven mode: table and database/paths all come from config.json.
    config = StartConfig()
    table_name = config.get_predicted_table_name()
    parent = config.get_parent_folder()
    if not parent or not os.path.isdir(parent):
        print(f"Error: parent_folder from config is missing or not a directory: {parent!r}")
        sys.exit(1)
    pdb_folder = os.path.join(parent, "pdb")
    af_folder = os.path.join(parent, "af")
    if not os.path.isdir(pdb_folder):
        print(f"Error: pdb folder not found: {pdb_folder}")
        sys.exit(1)
    if not os.path.isdir(af_folder):
        print(f"Error: af folder not found: {af_folder}")
        sys.exit(1)
    set_eval_metrics_database_path(config.get_database_path())
    print(f"Using database: {_eval_metrics_db_path}")
    print(f"Table: {table_name} | pdb: {pdb_folder} | af: {af_folder}")
    entries = load_predicted_entries_from_table(
        table_name, _eval_metrics_db_path, pdb_folder, af_folder
    )
    if not entries:
        print("No matching structures to process (check DB FileName rows and CIF paths).")
        sys.exit(0)
    process_files(
        pdb_folder, af_folder, entries=entries,
        compute_lddt=compute_lddt, compute_rnacomplex_tm=compute_rnacomplex_tm,
        only_rna_tm_allchains=only_rna_tm_allchains, table_name=table_name,
    )