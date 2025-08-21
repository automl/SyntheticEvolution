import MDAnalysis as mda
from simtk.openmm import app
from simtk.openmm import unit
from simtk.openmm.app import PDBFile
from simtk.openmm.app import ForceField
from simtk.openmm.app import Simulation
from simtk.openmm.app import PME
from simtk.openmm import LangevinIntegrator
from simtk.openmm import Platform
import sqlite3


# Step 1: Load CIF File
def load_cif(file_path):
    """Loads the CIF file into an MDAnalysis Universe."""
    try:
        u = mda.Universe(file_path)
        print("CIF file loaded successfully.")
        return u
    except Exception as e:
        print(f"Error loading CIF file: {e}")
        return None


# Step 2: Extract Protein-RNA Complex
def extract_complex(universe):
    """Extracts protein and RNA components from the molecular structure."""
    protein = universe.select_atoms("protein")
    rna = universe.select_atoms("nucleic")
    print(f"Extracted Protein: {protein.n_atoms} atoms, RNA: {rna.n_atoms} atoms.")
    return protein, rna


# Step 3: Calculate Gibbs Free Energy
def calculate_gibbs_free_energy(protein, rna, temperature=300):
    """Estimates Gibbs free energy change using forcefield-based scoring."""
    # Save protein and RNA structures to temporary files
    protein.write("protein.pdb")
    rna.write("rna.pdb")

    # Load structures with OpenMM
    pdb_protein = PDBFile("protein.pdb")
    pdb_rna = PDBFile("rna.pdb")

    # Define the forcefield
    forcefield = ForceField("amber99sb.xml", "tip3p.xml")  # Common biomolecular FF

    # Create a combined system
    modeller = app.Modeller(pdb_protein.topology, pdb_protein.positions)
    modeller.add(pdb_rna.topology, pdb_rna.positions)

    # Add solvent
    modeller.addSolvent(forcefield, model='tip3p', boxSize=app.Vec3(10.0, 10.0, 10.0) * unit.nanometers)

    # Build system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds
    )

    # Define integrator and simulation
    integrator = LangevinIntegrator(temperature * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
    platform = Platform.getPlatformByName('CPU')
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    # Minimize energy
    simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy=True)
    free_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    print(f"Estimated Gibbs Free Energy: {free_energy:.2f} kJ/mol")
    return free_energy


# Step 4: Update Database
def update_database(db_path, complex_id, gibbs_free_energy):
    """Updates the database with calculated Gibbs Free Energy."""
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("""
            UPDATE pred_protein_rna
            SET GibbsFreeEnergy = ?
            WHERE Id = ?
        """, (gibbs_free_energy, complex_id))
        conn.commit()
        print(f"Updated Gibbs Free Energy for Complex ID {complex_id}.")
        conn.close()
    except Exception as e:
        print(f"Error updating database: {e}")


# Step 5: Main Function
def main():
    db_path = "protein_rna_predictions.db"  # Path to SQLite database
    cif_dir = "path/to/cif/files"  # Directory containing CIF files

    # Connect to database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Query table for complex data
    cursor.execute("SELECT Id, FileName FROM pred_protein_rna WHERE GibbsFreeEnergy IS NULL")
    complexes = cursor.fetchall()
    conn.close()

    for complex_id, file_name in complexes:
        print(f"Processing Complex ID: {complex_id}, File: {file_name}")

        # Load CIF file and calculate Gibbs Free Energy
        cif_path = f"{cif_dir}/{file_name}"
        universe = load_cif(cif_path)
        if universe:
            protein, rna = extract_complex(universe)
            gibbs_free_energy = calculate_gibbs_free_energy(protein, rna)

            # Update database
            update_database(db_path, complex_id, gibbs_free_energy)


if __name__ == "__main__":
    main()
