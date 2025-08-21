import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from database.startConfig import StartConfig
from database.databaseMethods import DatabaseMethods
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import sqlite3

allowed_ligands = {'ADP', 'ATP', 'AMP', 'GTP', 'GDP', 'FAD', 'NAD', 'NAP', 'NDP', 'HEM', 'HEC', 'PLM', 'OLA', 'MYR',
                   'CIT', 'CLA', 'CHL', 'BCL', 'BCB'}
allowed_ions = {'MG', 'ZN', 'CL', 'CA', 'NA', 'MN', 'K', 'FE', 'CU', 'CO'}

config = StartConfig()
db = DatabaseMethods()
ref_folder = config.ref_folder
exp_table = config.get_reference_table_name()

# Add columns if not present
def add_column_if_not_exists(conn, table, column, coltype):
    cursor = conn.cursor()
    cursor.execute(f"PRAGMA table_info({table})")
    columns = [row[1] for row in cursor.fetchall()]
    if column not in columns:
        cursor.execute(f"ALTER TABLE {table} ADD COLUMN {column} {coltype}")
        conn.commit()

add_column_if_not_exists(db.connection, exp_table, 'Ions', 'TEXT')
add_column_if_not_exists(db.connection, exp_table, 'Ligands', 'TEXT')

# Get structures from database
query = f"SELECT PDBId FROM {exp_table}"
structures = db.execute_query(query)

for (pdb_id,) in structures:
    cif_file = os.path.join(ref_folder, f"{pdb_id}.cif")
    if not os.path.isfile(cif_file):
        print(f"CIF file not found for {pdb_id}")
        continue
    try:
        pdb_info = MMCIF2Dict(cif_file)
        non_poly_ids = pdb_info.get('_pdbx_entity_nonpoly.comp_id', [])
        if isinstance(non_poly_ids, str):
            non_poly_ids = [non_poly_ids]
    except Exception as e:
        print(f"Error parsing {cif_file}: {e}")
        continue

    ions = []
    ligands = []
    for comp_id in non_poly_ids:
        if comp_id == 'HOH':
            continue  # Disregard HOH completely
        if comp_id == 'SO4':
            ions.append([0, 'SO4'])  # Treat SO4 as a disallowed ion
            continue
            
        # Check against allowed lists FIRST
        if comp_id in allowed_ions:
            ions.append([1, comp_id])
        elif comp_id in allowed_ligands:
            ligands.append([1, comp_id])
        # Fallback for DISALLOWED items based on length
        elif len(comp_id) <= 2:
            ions.append([0, comp_id])
        elif len(comp_id) == 3:
            ligands.append([0, comp_id])

    # If no ions/ligands, store empty list
    ions_str = str(ions) if ions else '[]'
    ligands_str = str(ligands) if ligands else '[]'

    # Update the database
    try:
        db.connection.execute(f"UPDATE {exp_table} SET Ions = ?, Ligands = ? WHERE PDBId = ?", (ions_str, ligands_str, pdb_id))
        db.connection.commit()
    except Exception as e:
        print(f"Error updating {pdb_id}: {e}")

db.close_connection()
print("Done updating Ions and Ligands columns.") 