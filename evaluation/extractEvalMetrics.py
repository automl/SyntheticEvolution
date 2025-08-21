import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from pymol import cmd
import tempfile
from structures.rna import RNA
from structures.dna import DNA
from structures.protein import Protein
import re
from scipy.spatial.distance import pdist, squareform
import sqlite3
import json
import traceback
import numpy as np

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

def get_rmsd_tm_lddt_scores(file_name_cif_pdb, file_name_cif_af):
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
    cmd.load(temp_pdb_file_ref_path, 'reference')
    cmd.load(temp_pdb_file_tgt_path, 'predicted')
    cmd.remove('solvent')
    cmd.remove('organic') # Remove all ligands
    # cmd.remove('ions') # Remove all ions - not sufficient
    cmd.remove('inorganic') # Remove all ions - not sufficient
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
    print("reference residues count", count_residues_per_chain("reference"))
    print("predicted residues count", count_residues_per_chain("predicted"))

    def extract_score(content, score_type='TM-score'):
        if score_type == 'RMSD':
            # Extract the RMSD value specifically
            match = re.search(r'Aligned length=\s*\d+,\s*RMSD=\s*([0-9.]+)', content)
            if match:
                score = round(float(match.group(1)), 2)
                print(f"Extracted RMSD: {score}")
                return score
        elif score_type == 'TM-score':
            # Extract the TM-score value specifically from the reference structure normalization
            match = re.search(r'TM-score=\s*([0-9.]+)\s*\(normalized by length of Structure_1', content)
            if match:
                score = round(float(match.group(1)), 2)
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

    def contains_parentheses(sequence):
        pattern = r"\(.*?\)"
        if re.search(pattern, sequence):
            return True  # Pattern found
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
            print("pdb_protein_chain_id, pdb_rna_chain_id", pdb_protein_chain_id, pdb_rna_chain_id)
            rename_protein_chain(pdb_protein_chain_id, pdb_rna_chain_id)
            print("pdb_protein_chain_id, pdb_rna_chain_id", pdb_protein_chain_id, pdb_rna_chain_id)
    # Alter the chain identifiers in the predicted structure to match the reference structure https://pymol.org/dokuwiki/doku.php?id=command:alter
    #     cmd.alter(f'predicted and chain {af_protein_chain_id}', f'chain="{pdb_protein_chain_id}"')
        alignment_result = cmd.align(f'predicted and chain {af_protein_chain_id} and polymer.protein and name CA',
                                     f'reference and chain {pdb_protein_chain_id} and polymer.protein and name CA',
                                     quiet=0)  # , object='alignment_rna', reset=1)
        rmsd_protein = round(alignment_result[0], 2)
        print(f"RMSD for protein structure from alignment = {rmsd_protein} Å")
        protein_rmsd.append(rmsd_protein)
        # num_protein_atoms = cmd.count_atoms(f"predicted and chain {af_protein_chain_id} and polymer.protein")

        protein_aligned_reference_pdb = os.path.join(aligned_folder, f'{pdb_id}_pdb_aligned_protein.pdb')
        protein_aligned_predicted_pdb = os.path.join(aligned_folder, f'{pdb_id}_af_aligned_protein.pdb')
        # print(f"Reference chain ID: {pdb_protein_chain_id}")
        # print(f"Predicted chain ID: {af_protein_chain_id}")
        cmd.save(protein_aligned_reference_pdb, f'reference and chain {pdb_protein_chain_id} and polymer.protein')
        cmd.save(protein_aligned_predicted_pdb, f'predicted and chain {af_protein_chain_id} and polymer.protein')

        if contains_parentheses(protein.protein_sequence) == False:
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
                if rmsd_score_protein:
                    protein_rmsd.append(rmsd_score_protein)
        else:
            tm_score_protein = None
            rmsd_score_protein = None

        # Calculate protein LDDT using the complex LDDT function
        protein_lddt_score = calculate_complex_lddt(protein_aligned_reference_pdb, protein_aligned_predicted_pdb)
        print(f"Protein lDDT Score: {protein_lddt_score:.2f}")


    if rna.is_rna:
        if protein.is_protein: # Problem: 7wl0 still contained protein chains despite selection
            cmd.remove('polymer.protein')
            rename_protein_chain(pdb_rna_chain_id, pdb_protein_chain_id)
        # cmd.alter(f'predicted and chain {af_rna_chain_id}', f'chain="{pdb_rna_chain_id}"')
        alignment_result = cmd.align(f'predicted and chain {af_rna_chain_id} and polymer.nucleic and name P',
                  f'reference and chain {pdb_rna_chain_id} and polymer.nucleic and name P',
                  quiet=0, reset=1) #, object='alignment_rna', reset=1)
        rmsd_rna = round(alignment_result[0], 2)
        print(f"RMSD for RNA structure from alignment = {rmsd_rna} Å")
        rna_rmsd.append(rmsd_rna)
        num_rna_atoms = cmd.count_atoms(f"predicted and chain {af_rna_chain_id} and polymer.nucleic")

        rna_aligned_reference_pdb = os.path.join(aligned_folder, f'{pdb_id}_pdb_aligned_rna.pdb')
        rna_aligned_predicted_pdb = os.path.join(aligned_folder, f'{pdb_id}_af_aligned_rna.pdb')

        # Ensure chain IDs match before saving
        # cmd.alter(f'reference and chain {pdb_rna_chain_id}', f'chain="{af_rna_chain_id}"')
        # cmd.save(rna_aligned_reference_pdb, f'reference and chain {af_rna_chain_id}')
        # cmd.save(rna_aligned_predicted_pdb, f'predicted and chain {af_rna_chain_id}')

        cmd.save(rna_aligned_reference_pdb, f'reference and chain {pdb_rna_chain_id}')
        cmd.save(rna_aligned_predicted_pdb, f'predicted and chain {af_rna_chain_id}')

        if contains_parentheses(rna.rna_sequence) == False: # RNA does not contain any modified nucleotides
            # Obtain TM-score from https://zhanggroup.org/TM-score/
            # Try compiling: g + + -O3 - ffast - math - o RNAalign RNAalign.cpp
            # MAC: The file includes both stdlib.h and malloc.h. On macOS, we should remove the malloc.h include since stdlib.h already provides the necessary memory allocation functions.
            tm_output_path_rna = os.path.join(af_folder, f'{pdb_id}_rna_tm_output.txt')
            # os.system(f'./RNAalign/RNAalign {rna_aligned_reference_pdb} {rna_aligned_predicted_pdb} > {tm_output_path_rna}')
            os.system(
                f"./USalign/USalign {rna_aligned_reference_pdb} {rna_aligned_predicted_pdb} "
                f"-mol RNA -mm 0 -ter 2 > {tm_output_path_rna}"
            )

            with open(tm_output_path_rna, 'r') as file:
                content = file.read()
                tm_score_rna = extract_score(content, 'TM-score')
                rmsd_score_rna = extract_score(content, 'RMSD')
                if rmsd_score_rna:
                    rna_rmsd.append(rmsd_score_rna)
        else:
            tm_score_rna = None
            rmsd_score_rna = None

        # Calculate RNA LDDT using the complex LDDT function
        rna_lddt_score = calculate_complex_lddt(rna_aligned_reference_pdb, rna_aligned_predicted_pdb)
        # rna_lddt_score = 0
        print(f"RNA lDDT Score: {rna_lddt_score:.2f}")

    if dna.is_dna:
        # cmd.alter(f'predicted and chain {af_dna_chain_id}', f'chain="{pdb_dna_chain_id}"') # Change chain label?
        # It seems that the other chain labelled as {pdb_dna_chain_id} is considered instead
        alignment_result = cmd.align(f'predicted and chain {af_dna_chain_id} and polymer.nucleic and name P',
                  f'reference and chain {pdb_dna_chain_id} and polymer.nucleic and name P',
                  quiet=0) #, object='alignment_rna', reset=1)
        rmsd_dna = round(alignment_result[0], 2)
        print(f"RMSD for DNA structure from alignment = {rmsd_dna} Å")
        dna_rmsd.append(rmsd_dna)
        num_dna_atoms = cmd.count_atoms(f"predicted and chain {af_dna_chain_id} and polymer.nucleic")

        dna_aligned_reference_pdb = os.path.join(aligned_folder, f'{pdb_id}_pdb_aligned_dna.pdb')
        dna_aligned_predicted_pdb = os.path.join(aligned_folder, f'{pdb_id}_af_aligned_dna.pdb')
        cmd.save(dna_aligned_reference_pdb, f'reference and chain {pdb_dna_chain_id}')
        cmd.save(dna_aligned_predicted_pdb, f'predicted and chain {af_dna_chain_id}')

        if contains_parentheses(protein.protein_sequence) == False:
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

        # Calculate DNA LDDT using the complex LDDT function
        dna_lddt_score = calculate_complex_lddt(dna_aligned_reference_pdb, dna_aligned_predicted_pdb)
        print(f"DNA lDDT Score: {dna_lddt_score:.2f}")

    database_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'database', 'rbpDatabase.db')
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    # try:
    #     cursor.execute('''ALTER TABLE pred_protein_rna
    #                      ADD COLUMN Complex_LDDT REAL''')
    # except sqlite3.OperationalError:
    #     pass  # Column already exists

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
    # chains = cmd.get_chains('referenceA')
    # print("Chains in referenceA:", chains)
    # chains = cmd.get_chains('predictedA')
    # print("Chains in referenceA:", chains)

    # # Calculate LDDT for the entire complex
    complex_lddt = calculate_complex_lddt(temp_pdb_file_ref_path, temp_pdb_file_tgt_path)
    print(f"Complex LDDT = {complex_lddt}")

    if protein.is_protein:
        alignment_result = cmd.align('predictedA and polymer and name CA', 'referenceA and polymer and name CA',
                                     quiet=0)
    else:
        alignment_result = cmd.align('predictedA and polymer and name P', 'referenceA and polymer and name P',
                                     quiet=0)

    pymol_rmsd = round(alignment_result[0], 2)
    print(f"PyMOL RMSD after alignment = {pymol_rmsd} Å")
    complex_rmsd = [pymol_rmsd]  # Initialize list with PyMOL RMSD

    tm_output_path_complex = os.path.join(af_folder, f'{pdb_id}_complex_tm_output.txt')
    os.system(
        f"./USalign/USalign {temp_pdb_file_ref_path} {temp_pdb_file_tgt_path} "
        f"-mol auto -mm 1 -ter 0 > {tm_output_path_complex}"
    )

    with open(tm_output_path_complex, 'r') as file:
        content = file.read()
        tm_score_complex = extract_score(content, 'TM-score')
        rmsd_score_complex = extract_score(content, 'RMSD')
        if rmsd_score_complex:
            complex_rmsd.append(rmsd_score_complex)

    # Handle all RMSD values consistently
    complex_rmsd = handle_list_or_value(complex_rmsd)
    protein_rmsd = handle_list_or_value(protein_rmsd)
    rna_rmsd = handle_list_or_value(rna_rmsd)
    dna_rmsd = handle_list_or_value(dna_rmsd)

    if protein.is_protein and rna.is_rna and not dna.is_dna:
        cursor.execute('''UPDATE pred_protein_rna
                          SET Complex_RMSD = ?, Protein_RMSD = ?, RNA_RMSD = ?, 
                              Protein_LDDT = ?, RNA_LDDT = ?, 
                              Protein_TM = ?, RNA_TM = ?, Complex_TM = ?,
                              Complex_LDDT = ?
                          WHERE FileName = ?''',
                       (complex_rmsd, protein_rmsd, rna_rmsd,
                        protein_lddt_score, rna_lddt_score, 
                        tm_score_protein, tm_score_rna,
                        tm_score_complex, complex_lddt,
                        file_name_cif_af))

    if protein.is_protein and rna.is_rna and dna.is_dna:
        cursor.execute('''UPDATE pred_protein_rna_dna
                                  SET Complex_RMSD = ?, Protein_RMSD = ?, RNA_RMSD = ?, DNA_RMSD = ?, Protein_LDDT = ?, 
                                  RNA_LDDT = ?, DNA_LDDT = ?, Protein_TM = ?, RNA_TM = ?, DNA_TM = ?, 
                                  Complex_TM = ?, Complex_LDDT = ?
                                  WHERE FileName = ?''',
                       (complex_rmsd, protein_rmsd, rna_rmsd, dna_rmsd, protein_lddt_score, rna_lddt_score, dna_lddt_score,
                        tm_score_protein, tm_score_rna, tm_score_dna, tm_score_complex, complex_lddt, file_name_cif_af))

    if protein.is_protein and not rna.is_rna and dna.is_dna:
        cursor.execute('''UPDATE pred_protein_dna
                                  SET Complex_RMSD = ?, Protein_RMSD = ?, DNA_RMSD = ?, Protein_LDDT = ?, DNA_LDDT = ?, 
                                  Protein_TM = ?, DNA_TM = ?, Complex_TM = ?, Complex_LDDT = ?
                                  WHERE FileName = ?''',
                       (complex_rmsd, protein_rmsd, dna_rmsd, protein_lddt_score, dna_lddt_score, tm_score_protein,
                        tm_score_dna, tm_score_complex, complex_lddt, file_name_cif_af))

    if not protein.is_protein and rna.is_rna and not dna.is_dna:
        cursor.execute('''UPDATE pred_rna_rna
                                  SET Complex_RMSD = ?, RNA_RMSD = ?, RNA_LDDT = ?, RNA_TM = ?, Complex_TM = ?, Complex_LDDT = ?
                                  WHERE FileName = ?''',
                       (complex_rmsd, rna_rmsd, rna_lddt_score,
                        tm_score_rna, tm_score_complex, complex_lddt, file_name_cif_af))

    conn.commit()
    conn.close()

    os.remove(temp_pdb_file_ref_path)
    os.remove(temp_pdb_file_tgt_path)
    cmd.reinitialize()

def process_files(pdb_path, af_path):
    for file_name_cif_af in os.listdir(af_path):
        if file_name_cif_af.endswith('.cif'):
            try:
                pdb_id = file_name_cif_af.split('_')[1]
                file_name_cif_pdb = f"{pdb_id}.cif"
                pdb_file_path = os.path.join(pdb_path, file_name_cif_pdb)
                
                # Check if PDB file exists
                if not os.path.exists(pdb_file_path):
                    print(f"Warning: No corresponding file found in PDB folder for {pdb_id}")
                    continue
                
                print(f"\nProcessing {pdb_id}...")

                # Try to get scores, skip file if it fails
                try:
                    get_rmsd_tm_lddt_scores(file_name_cif_pdb, file_name_cif_af)
                except Exception as e:
                    tb = traceback.format_exc()
                    print(f"Error calculating scores for {pdb_id}: {e}")
                    print(f"Traceback:\n{tb}")
                    # print(f"Skipping {pdb_id}...")
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
            return round(lddt_score, 2)
        return None
    
    except Exception as e:
        print(f"Error in LDDT calculation: {str(e)}")
        return None

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <reference_directory_path> <predicted_directory_path>")
        sys.exit(1)
    pdb_folder = sys.argv[1]
    af_folder = sys.argv[2]

    process_files(pdb_folder, af_folder)