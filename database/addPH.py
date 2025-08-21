import os
import glob
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from database.startConfig import StartConfig
from database.databaseMethods import DatabaseMethods
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

# Initialize database connection
config = StartConfig()
db = DatabaseMethods()

# Buffer pH mapping
BUFFER_PH_MAP = {
    'formate': 4.0,
    'tris': 8.0,  # Average of 7.0-9.0
    'hepes': 7.5,  # Average of 6.8-8.2
    'citrate': 4.5,  # Average of 3.0-6.0
    'mes': 6.1,  # Average of 5.5-6.7
    'phosphate': 7.25,  # Average of 6.5-8.0
    'pbs': 7.25,  # Same as phosphate
    'acetate': 4.8,  # Average of 4.5-5.1
    'mops': 7.2,  # Average of 6.5-7.9
    'tricine': 8.0,  # Average of 7.4-8.8
    'bis-tris': 6.5,  # Average of 5.8-7.2
    'cacodylate': 6.5  # Average of 6.0-7.0
}

def extract_ph_from_details(details):
    """Extract pH value from experimental details text."""
    if not details or details == '?' or details.lower() == 'null':
        return None
        
    # First try to find explicit pH value
    details = details.lower()
    if 'ph' in details:
        try:
            # Find the position of 'ph' and look for numbers after it
            ph_pos = details.find('ph')
            ph_text = details[ph_pos:ph_pos+10]  # Take a substring after 'ph'
            # Extract the first number after 'ph'
            import re
            ph_match = re.search(r'ph\s*(\d+\.?\d*)', ph_text)
            if ph_match:
                return float(ph_match.group(1))
        except:
            pass
    
    # If no explicit pH, look for buffer names
    for buffer, ph in BUFFER_PH_MAP.items():
        if buffer in details:
            return ph
            
    return None

def add_columns():
    """Add exp_pH, exp_Temp, and RNA_ElectrostaticPotential columns if they don't exist."""
    try:
        exp_table = config.get_reference_table_name()
        pred_table = config.get_predicted_table_name()
        
        # Check exp_protein_rna columns
        columns_query = f"PRAGMA table_info({exp_table})"
        columns = db.execute_query(columns_query)
        column_names = [col[1] for col in columns]
        
        if 'exp_pH' not in column_names:
            db.execute_query(f"ALTER TABLE {exp_table} ADD COLUMN exp_pH REAL")
            logger.info(f"Added exp_pH column to {exp_table}")
            
        if 'exp_Temp' not in column_names:
            db.execute_query(f"ALTER TABLE {exp_table} ADD COLUMN exp_Temp REAL")
            logger.info(f"Added exp_Temp column to {exp_table}")
            
        if 'RNA_ElectrostaticPotential' not in column_names:
            db.execute_query(f"ALTER TABLE {exp_table} ADD COLUMN RNA_ElectrostaticPotential REAL")
            logger.info(f"Added RNA_ElectrostaticPotential column to {exp_table}")
            
        # Check pred_protein_rna columns
        columns_query = f"PRAGMA table_info({pred_table})"
        columns = db.execute_query(columns_query)
        column_names = [col[1] for col in columns]
        
        if 'RNA_ElectrostaticPotential' not in column_names:
            db.execute_query(f"ALTER TABLE {pred_table} ADD COLUMN RNA_ElectrostaticPotential REAL")
            logger.info(f"Added RNA_ElectrostaticPotential column to {pred_table}")
            
    except Exception as e:
        logger.error(f"Error adding columns: {str(e)}")
        return False
    return True

def process_structure(pdb_id, experiment_type, cif_file):
    """Process a single structure and extract pH and temperature values."""
    try:
        pdb_info = MMCIF2Dict(cif_file)
        exp_ph = None
        exp_temp = None
        
        if experiment_type == 'X-RAY DIFFRACTION':
            # Try to get pH from crystal growth conditions
            ph_value = pdb_info.get('_exptl_crystal_grow.pH', ['?'])[0]
            if ph_value != '?' and ph_value.lower() != 'null':
                try:
                    exp_ph = float(ph_value)
                except:
                    pass
                    
            # If no pH value, check details
            if not exp_ph:
                details = pdb_info.get('_exptl_crystal_grow.pdbx_details', ['?'])[0]
                exp_ph = extract_ph_from_details(details)
                
            # Get temperature
            temp_value = pdb_info.get('_exptl_crystal_grow.temp', ['?'])[0]
            if temp_value != '?' and temp_value.lower() != 'null':
                try:
                    exp_temp = float(temp_value)
                except:
                    pass
                    
        elif experiment_type == 'ELECTRON MICROSCOPY':
            # Try different fields for pH
            ph_fields = ['_em_buffer.pH', '_pdbx_buffer_pH', '_em_experiment.pH']
            for field in ph_fields:
                ph_value = pdb_info.get(field, ['?'])[0]
                if ph_value != '?' and ph_value.lower() != 'null':
                    try:
                        exp_ph = float(ph_value)
                        break
                    except:
                        continue
                        
        # Set default pH if none found
        if exp_ph is None:
            exp_ph = 7.0
            
        return exp_ph, exp_temp
        
    except Exception as e:
        logger.error(f"Error processing {pdb_id}: {str(e)}")
        return None, None

def update_database():
    """Update database with experimental conditions from CIF files."""
    if not add_columns():
        return
        
    # Get reference folder path
    ref_folder = config.ref_folder
    exp_table = config.get_reference_table_name()
    
    # Get structures and their experiment types from database
    query = f"SELECT PDBId, Experiment FROM {exp_table}"
    structures = db.execute_query(query)
    
    for pdb_id, experiment_type in structures:
        if not pdb_id:
            continue
            
        # Find CIF file
        cif_file = os.path.join(ref_folder, f"{pdb_id.lower()}.cif")
        if not os.path.exists(cif_file):
            logger.warning(f"CIF file not found for {pdb_id}")
            continue
            
        # Process structure
        exp_ph, exp_temp = process_structure(pdb_id, experiment_type, cif_file)
        
        # Update database
        try:
            if experiment_type == 'X-RAY DIFFRACTION':
                query = f"""
                UPDATE {exp_table} 
                SET exp_pH = ?, exp_Temp = ?
                WHERE PDBId = ?
                """
                db.execute_query(query, (exp_ph, exp_temp, pdb_id))
            else:  # ELECTRON MICROSCOPY
                query = f"""
                UPDATE {exp_table} 
                SET exp_pH = ?
                WHERE PDBId = ?
                """
                db.execute_query(query, (exp_ph, pdb_id))
                
            logger.info(f"Updated {pdb_id}: pH={exp_ph}, Temp={exp_temp if experiment_type == 'X-RAY DIFFRACTION' else 'N/A'}")
            
        except Exception as e:
            logger.error(f"Error updating database for {pdb_id}: {str(e)}")

if __name__ == "__main__":
    update_database() 