import os
import sys
import sqlite3
from Bio.PDB import MMCIFParser, PDBIO
import tempfile
from pymol import cmd
from structures.rna import RNA
from structures.dna import DNA
from structures.protein import Protein

def calculate_complex_lddt(reference_path, predicted_path):
    """Calculate LDDT for the entire complex using all heavy atoms from polymers only."""
    from Bio.PDB import PDBParser
    import numpy as np
    
    def get_polymer_atom_coords(structure):
        """Get coordinates of all heavy atoms from polymers."""
        coords = []
        mask = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        # Only include heavy atoms (not H)
                        if atom.element != 'H':
                            coords.append(atom.get_coord())
                            mask.append(1.0)
        return np.array(coords), np.array(mask).reshape(-1, 1)

    def calculate_lddt_score(pred_points, true_points, true_mask, cutoff=15.0):
        """Calculate LDDT score following the simplified approach."""
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
        # Load structures into PyMOL and select only polymers
        cmd.load(reference_path, 'reference')
        cmd.load(predicted_path, 'predicted')
        cmd.select("polymers_only", "polymer")
        
        # Save the polymer-only structures to temporary files
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_ref:
            cmd.save(temp_ref.name, 'reference and polymers_only')
            temp_ref_path = temp_ref.name
            
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pred:
            cmd.save(temp_pred.name, 'predicted and polymers_only')
            temp_pred_path = temp_pred.name

        # Parse structures
        parser = PDBParser(QUIET=True)
        ref_structure = parser.get_structure('reference', temp_ref_path)
        pred_structure = parser.get_structure('predicted', temp_pred_path)

        # Get coordinates and masks for all heavy atoms
        ref_coords, ref_mask = get_polymer_atom_coords(ref_structure)
        pred_coords, _ = get_polymer_atom_coords(pred_structure)

        if len(ref_coords) == 0 or len(pred_coords) == 0:
            print("No heavy atoms found in one or both structures")
            return None

        # Ensure equal number of points, max 30000 atoms
        min_length = min(len(ref_coords), len(pred_coords), 30000)
        ref_coords = ref_coords[:min_length]
        pred_coords = pred_coords[:min_length]
        ref_mask = ref_mask[:min_length]

        print(f"Calculating LDDT using {min_length} heavy atoms from polymers")

        # Calculate LDDT score
        lddt_score = calculate_lddt_score(pred_coords, ref_coords, ref_mask)
        
        # Clean up temporary files
        os.remove(temp_ref_path)
        os.remove(temp_pred_path)
        cmd.reinitialize()
        
        if lddt_score is not None:
            return round(lddt_score, 2)
        return None
    
    except Exception as e:
        print(f"Error in LDDT calculation: {str(e)}")
        return None

def update_lddt_scores(pdb_folder, af_folder, database_path):
    """Update LDDT scores for protein, RNA, and DNA components."""
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()
    
    # Get all files that need updating
    cursor.execute('''
        SELECT FileName, exp_db_id 
        FROM pred_rna_rna 
        WHERE FileName IS NOT NULL
    ''')
    files = cursor.fetchall()
    
    for file_name_cif_af, exp_db_id in files:
        try:
            print(f"\nProcessing {exp_db_id}...")
            
            # Get corresponding PDB file
            file_name_cif_pdb = f"{exp_db_id}.cif"
            pdb_file_path = os.path.join(pdb_folder, file_name_cif_pdb)
            
            if not os.path.exists(pdb_file_path):
                print(f"Warning: No corresponding file found in PDB folder for {exp_db_id}")
                continue
            
            # Load structures
            cif_parser = MMCIFParser(QUIET=True)
            structure_reference = cif_parser.get_structure(file_name_cif_pdb, pdb_file_path)
            structure_predicted = cif_parser.get_structure(file_name_cif_af, os.path.join(af_folder, file_name_cif_af))
            
            # Create aligned folder if it doesn't exist
            aligned_folder = os.path.join(os.path.dirname(af_folder), "aligned")
            os.makedirs(aligned_folder, exist_ok=True)
            
            # Process each component
            protein_lddt_score = None
            rna_lddt_score = None
            dna_lddt_score = None
            
            # Get component information
            rna = RNA.get_rna_from_db(id=exp_db_id.upper())
            protein = Protein.get_protein_from_db(id=exp_db_id.upper())
            dna = DNA.get_dna_from_db(id=exp_db_id.upper())
            
            # Process protein if present
            if protein.is_protein:
                pdb_protein_chain_id = protein.get_longest_chain_id()
                af_protein_chain_id = pdb_protein_chain_id  # Assuming chain IDs match
                
                # Save aligned structures
                protein_aligned_reference_pdb = os.path.join(aligned_folder, f'{exp_db_id}_pdb_aligned_protein.pdb')
                protein_aligned_predicted_pdb = os.path.join(aligned_folder, f'{exp_db_id}_af_aligned_protein.pdb')
                
                cmd.load(pdb_file_path, 'reference')
                cmd.load(os.path.join(af_folder, file_name_cif_af), 'predicted')
                cmd.save(protein_aligned_reference_pdb, f'reference and chain {pdb_protein_chain_id} and polymer.protein')
                cmd.save(protein_aligned_predicted_pdb, f'predicted and chain {af_protein_chain_id} and polymer.protein')
                
                protein_lddt_score = calculate_complex_lddt(protein_aligned_reference_pdb, protein_aligned_predicted_pdb)
                print(f"Protein LDDT Score: {protein_lddt_score}")
            
            # Process RNA if present
            if rna.is_rna:
                pdb_rna_chain_id = rna.get_longest_chain_id()
                af_rna_chain_id = pdb_rna_chain_id  # Assuming chain IDs match
                
                # Save aligned structures
                rna_aligned_reference_pdb = os.path.join(aligned_folder, f'{exp_db_id}_pdb_aligned_rna.pdb')
                rna_aligned_predicted_pdb = os.path.join(aligned_folder, f'{exp_db_id}_af_aligned_rna.pdb')
                
                cmd.load(pdb_file_path, 'reference')
                cmd.load(os.path.join(af_folder, file_name_cif_af), 'predicted')
                cmd.save(rna_aligned_reference_pdb, f'reference and chain {pdb_rna_chain_id} and polymer.nucleic')
                cmd.save(rna_aligned_predicted_pdb, f'predicted and chain {af_rna_chain_id} and polymer.nucleic')
                
                rna_lddt_score = calculate_complex_lddt(rna_aligned_reference_pdb, rna_aligned_predicted_pdb)
                print(f"RNA LDDT Score: {rna_lddt_score}")
            
            # Process DNA if present
            if dna.is_dna:
                pdb_dna_chain_id = dna.get_longest_chain_id()
                af_dna_chain_id = pdb_dna_chain_id  # Assuming chain IDs match
                
                # Save aligned structures
                dna_aligned_reference_pdb = os.path.join(aligned_folder, f'{exp_db_id}_pdb_aligned_dna.pdb')
                dna_aligned_predicted_pdb = os.path.join(aligned_folder, f'{exp_db_id}_af_aligned_dna.pdb')
                
                cmd.load(pdb_file_path, 'reference')
                cmd.load(os.path.join(af_folder, file_name_cif_af), 'predicted')
                cmd.save(dna_aligned_reference_pdb, f'reference and chain {pdb_dna_chain_id} and polymer.nucleic')
                cmd.save(dna_aligned_predicted_pdb, f'predicted and chain {af_dna_chain_id} and polymer.nucleic')
                
                dna_lddt_score = calculate_complex_lddt(dna_aligned_reference_pdb, dna_aligned_predicted_pdb)
                print(f"DNA LDDT Score: {dna_lddt_score}")
            
            # Update database with new LDDT scores
            # cursor.execute('''
            #     UPDATE pred_rna_rna
            #     SET Protein_LDDT = ?, RNA_LDDT = ?, DNA_LDDT = ?
            #     WHERE FileName = ?
            # ''', (protein_lddt_score, rna_lddt_score, dna_lddt_score, file_name_cif_af))

            cursor.execute('''
                            UPDATE pred_rna_rna 
                            SET RNA_LDDT = ?
                            WHERE FileName = ?
                        ''', (rna_lddt_score, file_name_cif_af))
            
            conn.commit()
            print(f"Updated LDDT scores for {exp_db_id}")
            
        except Exception as e:
            print(f"Error processing {exp_db_id}: {str(e)}")
            continue
    
    conn.close()

def main():
    if len(sys.argv) != 4:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <reference_directory_path> <predicted_directory_path> <database_path>")
        sys.exit(1)
    
    pdb_folder = sys.argv[1]
    af_folder = sys.argv[2]
    database_path = sys.argv[3]
    
    update_lddt_scores(pdb_folder, af_folder, database_path)

if __name__ == "__main__":
    main() 