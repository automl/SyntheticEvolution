import os
import subprocess
from Bio import PDB
import glob
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig
import tempfile
import shutil
import argparse

# Initialize database connection
config = StartConfig()
db = DatabaseMethods()

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate electrostatic potential for RNA structures')
    parser.add_argument('--save-files', action='store_true', 
                      help='Save PQR and DX files (default: only calculate potential)')
    return parser.parse_args()

"""'conda install -c conda-forge pdb2pqr apbs' if conda is not installed"""
def check_dependencies():
    """Check if required software is installed"""
    # Check PDB2PQR
    pdb2pqr_found = False
    for cmd in ["pdb2pqr30", "pdb2pqr", "pdb2pqr.exe"]:
        try:
            result = subprocess.run(f"{cmd} --help", shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                pdb2pqr_found = True
                break
        except:
            continue
    
    if not pdb2pqr_found:
        print("ERROR: PDB2PQR not found. Please install with: conda install -c conda-forge pdb2pqr")
        return False
    
    # Check APBS
    try:
        result = subprocess.run("apbs --version", shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print("ERROR: APBS not found. Please install with: conda install -c conda-forge apbs")
            return False
    except:
        print("ERROR: APBS not found. Please install with: conda install -c conda-forge apbs")
        return False
    
    return True


def clean_pdb_for_pdb2pqr(input_pdb, output_pdb):
    """
    Clean PDB file to make it compatible with PDB2PQR.
    - Removes HETATM records
    - Ensures proper chain IDs and residue connectivity
    - Removes alternative conformations
    - Only processes standard RNA residues (A, U, G, C)
    - Handles terminal residues properly
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('temp', input_pdb)
        
        # Standard RNA residue mapping (no modifications)
        residue_map = {
            # Standard RNA residues only
            'A': 'RA', 'U': 'RU', 'G': 'RG', 'C': 'RC',
            # DNA residues (convert to RNA)
            'DA': 'RA', 'DT': 'RU', 'DG': 'RG', 'DC': 'RC',
            # Alternative names
            'ADE': 'RA', 'URI': 'RU', 'GUA': 'RG', 'CYT': 'RC',
            'THY': 'RU'
        }

        # Define atom mappings (old name -> new name)
        atom_map = {
            'O1P': 'OP1',
            'O2P': 'OP2',
            "O5'": "O5'",
            "C5'": "C5'",
            "C4'": "C4'",
            "O4'": "O4'",
            "C3'": "C3'",
            "O3'": "O3'",
            "C2'": "C2'",
            "O2'": "O2'",
            "C1'": "C1'",
            'OP1': 'OP1',
            'OP2': 'OP2'
        }

        # Base-specific atoms
        purine_atoms = {
            'N9', 'C8', 'N7', 'C5', 'C6', 'N6',  # For A
            'O6', 'N1', 'C2', 'N2', 'N3', 'C4'   # Additional for G
        }
        
        pyrimidine_atoms = {
            'N1', 'C2', 'O2', 'N3', 'C4', 'O4',  # For U
            'N4', 'C5', 'C6'                      # Additional for C
        }

        # Write cleaned PDB with ordered atoms
        with open(output_pdb, 'w') as out_f:
            model = structure[0]
            prev_chain_id = None
            atom_serial = 1
            
            for chain in model:
                chain_id = chain.id
                if prev_chain_id is not None and chain_id != prev_chain_id:
                    out_f.write("TER\n")
                prev_chain_id = chain_id
                
                # Get valid residues in this chain
                valid_residues = [res for res in chain if res.get_id()[0] == ' ' and res.get_resname() in residue_map]
                
                for i, residue in enumerate(valid_residues):
                    is_first = (i == 0)
                    is_last = (i == len(valid_residues) - 1)
                    resname = residue_map[residue.get_resname()]
                    
                    # Write backbone atoms first
                    if not is_first:  # Skip phosphate group for first residue
                        for old_name in ['O1P', 'O2P', 'OP1', 'OP2', 'P']:
                            if old_name in residue:
                                atom = residue[old_name]
                                new_name = atom_map.get(old_name, old_name)
                                x, y, z = atom.get_coord()
                                b_factor = atom.get_bfactor()
                                occupancy = atom.get_occupancy()
                                element = atom.element.strip()
                                out_f.write(f"ATOM  {atom_serial:5d} {new_name:^4s} {resname:3s} {chain_id:1s}{residue.get_id()[1]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{b_factor:6.2f}          {element:>2s}  \n")
                                atom_serial += 1

                    # Write sugar atoms
                    for old_name in ["O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]:
                        if old_name in residue:
                            atom = residue[old_name]
                            new_name = atom_map.get(old_name, old_name)
                            x, y, z = atom.get_coord()
                            b_factor = atom.get_bfactor()
                            occupancy = atom.get_occupancy()
                            element = atom.element.strip()
                            out_f.write(f"ATOM  {atom_serial:5d} {new_name:^4s} {resname:3s} {chain_id:1s}{residue.get_id()[1]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{b_factor:6.2f}          {element:>2s}  \n")
                            atom_serial += 1

                    # Write base atoms based on residue type
                    base_atoms = purine_atoms if resname in ['RA', 'RG'] else pyrimidine_atoms
                    for atom_name in sorted(base_atoms):
                        if atom_name in residue:
                            atom = residue[atom_name]
                            x, y, z = atom.get_coord()
                            b_factor = atom.get_bfactor()
                            occupancy = atom.get_occupancy()
                            element = atom.element.strip()
                            out_f.write(f"ATOM  {atom_serial:5d} {atom_name:^4s} {resname:3s} {chain_id:1s}{residue.get_id()[1]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{b_factor:6.2f}          {element:>2s}  \n")
                            atom_serial += 1
                
                out_f.write("TER\n")
            
            out_f.write("END\n")

        # Only verify if file exists and is not empty
        if os.path.exists(output_pdb) and os.path.getsize(output_pdb) > 0:
            return True
        else:
            print(f"ERROR: Failed to create cleaned PDB file: {output_pdb}")
            return False
            
    except Exception as e:
        print(f"ERROR: Error cleaning PDB file {input_pdb}: {str(e)}")
        return False


def generate_pqr(input_pdb, output_pqr, exp_ph=7.0):
    """
    Converts PDB file to PQR format using PDB2PQR.
    PQR files contain atomic charges and radii necessary for electrostatic calculations.
    PDB2PQR automatically adds missing hydrogens and assigns proper protonation states.

    Args:
        input_pdb (str): Path to input PDB file
        output_pqr (str): Path to output PQR file
        exp_ph (float): Experimental pH value to use for protonation states

    References:
    1. Dolinsky et al. (2004), PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann electrostatics calculations
       Nucleic Acids Research, 32(W665–W667), DOI: 10.1093/nar/gkh381
    2. Jurrus et al. (2018), Improvements to the APBS biomolecular solvation software suite
       Protein Science, 27(1), 112-128, DOI: 10.1002/pro.3280
    """
    try:
        pdb2pqr_cmd = None
        for cmd in ["pdb2pqr30", "pdb2pqr", "pdb2pqr.exe"]:
            try:
                result = subprocess.run(f"{cmd} --help", shell=True, capture_output=True, text=True)
                if result.returncode == 0:
                    pdb2pqr_cmd = cmd
                    break
            except:
                continue
        
        if not pdb2pqr_cmd:
            print("ERROR: PDB2PQR not found")
            return False

        with tempfile.TemporaryDirectory() as temp_dir:
            cleaned_pdb = os.path.join(temp_dir, "cleaned.pdb")
            if not clean_pdb_for_pdb2pqr(input_pdb, cleaned_pdb):
                return False

            force_fields = ["AMBER", "PARSE", "CHARMM", "TYL06"]
            success = False
            
            for ff in force_fields:
                try:
                    command = f"{pdb2pqr_cmd} --ff={ff} --with-ph={exp_ph} --drop-water {cleaned_pdb} {output_pqr}"
                    result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
                    success = True
                    break
                except subprocess.CalledProcessError:
                    if os.path.exists(output_pqr):
                        os.remove(output_pqr)
                    continue

            if not success:
                print(f"ERROR: Failed to convert {input_pdb} with any force field")
                return False

            if not os.path.exists(output_pqr) or os.path.getsize(output_pqr) == 0:
                print(f"ERROR: PQR file {output_pqr} was not created or is empty")
                return False

            return True

    except Exception as e:
        print(f"ERROR: Error in generate_pqr for {input_pdb}: {str(e)}")
        if os.path.exists(output_pqr):
            os.remove(output_pqr)
        return False


def run_apbs(pqr_file, output_dir, save_files=False):
    """
    Runs APBS (Adaptive Poisson-Boltzmann Solver) to compute electrostatic potential.
    Uses the nonlinear Poisson-Boltzmann equation with physiological ionic conditions.
    Optionally saves visualization outputs that can be viewed in PyMOL.

    Args:
        pqr_file (str): Path to input PQR file
        output_dir (str): Directory to save output files
        save_files (bool): Whether to save output files for visualization

    Returns:
        float: Electrostatic potential energy in kJ/mol, or None if calculation fails
    """
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            base_name = os.path.splitext(os.path.basename(pqr_file))[0]
            base_name = base_name.replace('.dx', '').replace('_potential', '')
            potential_dx = os.path.join(temp_dir if not save_files else output_dir, f"{base_name}.dx")

            apbs_input = f"""read
    mol pqr {pqr_file}
end
elec name rna
    mg-auto
    dime 193 193 193
    cglen 120.0 120.0 120.0
    fglen 100.0 100.0 100.0
    cgcent mol 1
    fgcent mol 1
    mol 1
    npbe
    bcfl sdh
    pdie 4.0
    sdie 80.0
    ion charge 1 conc 0.150 radius 2.0
    ion charge -1 conc 0.150 radius 2.0
    srfm mol
    chgm spl2
    sdens 10.00
    srad 1.40
    swin 0.30
    temp 298.15
    calcenergy total
    calcforce no
    write pot dx {potential_dx}
end
print elecEnergy rna end
quit
"""
            with tempfile.NamedTemporaryFile(mode='w', suffix='.in') as f:
                f.write(apbs_input)
                f.flush()
                
                result = subprocess.run(f"apbs {f.name} 2>/dev/null", shell=True, check=True, 
                                    capture_output=True, text=True)
                
                energy = None
                for line in result.stdout.split('\n'):
                    if "Global net ELEC energy" in line:
                        energy = float(line.split('=')[1].strip().split()[0])
                
                if energy is None:
                    print(f"ERROR: Could not find energy value in APBS output for {pqr_file}")
                    return None

                if save_files:
                    readme_path = os.path.join(output_dir, "README_visualization.txt")
                    with open(readme_path, 'w') as f:
                        f.write("PyMOL Visualization Instructions\n")
                        f.write("============================\n\n")
                        f.write(f"Output files are saved in: {output_dir}\n\n")
                        f.write("To visualize the electrostatic potential:\n\n")
                        f.write(f"1. Load the PQR file:\n   load {pqr_file}\n\n")
                        f.write(f"2. Load the potential map:\n   load {potential_dx}\n\n")
                        f.write("3. Show surface:\n   show surface\n\n")
                        f.write("4. Color the surface:\n   color white\n\n")
                        f.write("5. Show potentials (copy and paste these commands):\n")
                        f.write(f"   isosurface pos, {base_name}, 1\n")
                        f.write(f"   isosurface neg, {base_name}, -1\n")
                        f.write("   set surface_color, blue, neg\n")
                        f.write("   set surface_color, red, pos\n")
                        f.write("   set transparency, 0.4\n")
                
                return energy

    except subprocess.CalledProcessError as e:
        print(f"ERROR: Error running APBS for {pqr_file}")
        return None
    except Exception as e:
        print(f"ERROR: Unexpected error running APBS for {pqr_file}: {str(e)}")
        return None


def extract_electrostatic_potential(pdb_file, exp_ph=7.0, save_files=False):
    """
    Extracts electrostatic potential from a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file
        exp_ph (float): Experimental pH value to use for protonation states
        save_files (bool): Whether to save PQR and DX files for visualization
        
    Returns:
        float: Electrostatic potential energy in kJ/mol, or None if calculation fails
    """
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            parent_dir = os.path.dirname(os.path.dirname(pdb_file))
            output_dir = os.path.join(parent_dir, "el_potential_output") if save_files else temp_dir
            if save_files:
                os.makedirs(output_dir, exist_ok=True)
            
            base_name = os.path.splitext(os.path.basename(pdb_file))[0]
            base_name = base_name.replace('.dx', '')
            pqr_file = os.path.join(output_dir if save_files else temp_dir, f"{base_name}.pqr")
            
            if not generate_pqr(pdb_file, pqr_file, exp_ph):
                return None
                
            potential = run_apbs(pqr_file, output_dir, save_files)
            return potential
            
    except Exception as e:
        print(f"ERROR: Error processing {pdb_file}: {str(e)}")
        return None


def process_aligned_folder():
    """
    Processes all RNA-containing PDB files in the aligned folder and updates the database.
    Only processes structures with unmodified RNA (RNAModified = '[]').
    Uses experimental pH values from the database for both experimental and predicted structures.
    """
    # Parse command line arguments
    args = parse_args()

    # Get aligned folder path
    aligned_folder = os.path.join(config.get_parent_folder(), "aligned")
    if not os.path.exists(aligned_folder):
        print(f"ERROR: Aligned folder not found at {aligned_folder}")
        return

    # Find all RNA PDB files
    pdb_files = glob.glob(os.path.join(aligned_folder, "*rna*.pdb"))

    # Get list of unmodified structures and experimental pH
    query = "SELECT PDBId, RNAModified, exp_pH FROM exp_protein_rna"
    db_data = db.execute_query(query)
    rna_modifications = {row[0].lower(): row[1] for row in db_data}
    exp_ph_values = {row[0].lower(): row[2] if row[2] is not None else 7.0 for row in db_data}
    unmodified_pdbs = {pdb_id for pdb_id, mods in rna_modifications.items() if mods == '[]'}

    # Group files by PDB ID
    pdb_groups = {}
    for pdb_file in pdb_files:
        base_name = os.path.basename(pdb_file)
        pdb_id = base_name.split('_')[0].lower()
        if pdb_id not in pdb_groups:
            pdb_groups[pdb_id] = []
        pdb_groups[pdb_id].append(pdb_file)

    # Process each PDB ID group
    for pdb_id, group_files in pdb_groups.items():
        exp_ph = exp_ph_values.get(pdb_id, 7.0)
        
        if pdb_id not in unmodified_pdbs:
            continue

        # Process each file in the group
        for pdb_file in group_files:
            base_name = os.path.basename(pdb_file)
            is_predicted = '_af_' in base_name
            print(f"Processing {base_name} (PDB ID: {pdb_id}, Predicted: {is_predicted}, pH: {exp_ph})")
            
            potential = extract_electrostatic_potential(pdb_file, exp_ph, args.save_files)
            
            if potential is not None:
                try:
                    table = "pred_protein_rna" if is_predicted else "exp_protein_rna"
                    id_column = "exp_db_id" if is_predicted else "PDBId"
                    pdb_id = pdb_id.upper()
                    query = f"UPDATE {table} SET RNA_ElectrostaticPotential = ? WHERE {id_column} = ?"
                    db.execute_query(query, (potential, pdb_id))
                except Exception as e:
                    print(f"ERROR: Failed to update database for {pdb_id}: {str(e)}")


if __name__ == "__main__":
    process_aligned_folder()

# Usage:
# Calculate potentials and update database (minimal output)
# python extractElectrostaticPotential.py

# Calculate potentials and save PQR/DX files with visualization instructions
# python extractElectrostaticPotential.py --save-files