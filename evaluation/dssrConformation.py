import os
import sys
import json
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig

class DSSRConformation:
    def __init__(self):
        self.database = DatabaseMethods()
        self.config = StartConfig()
        # Load config from JSON file
        self.config.load_config()

    def classify_chi_angle(self, chi_angle):
        """Classify glycosidic bond angle into anti/syn"""
        try:
            chi = float(chi_angle)
            if (90 <= chi <= 180) or (-180 <= chi <= -90):
                return 'anti'
            elif -90 <= chi <= 90:
                return 'syn'
        except (ValueError, TypeError):
            pass
        return 'other'

    def classify_sugar_pucker(self, pucker):
        """Classify sugar pucker conformation"""
        if pucker == "C3'-endo":
            return "C3'-endo"
        elif pucker == "C2'-endo":
            return "C2'-endo"
        elif pucker == "C3'-exo":
            return "C3'-exo"
        elif pucker == "C4'-exo":
            return "C4'-exo"
        elif pucker == "C2'-exo":
            return "C2'-exo"
        elif pucker == "O4'-endo":
            return "O4'-endo"
        elif pucker == "O4'-exo":
            return "O4'-exo"
        else:
            return 'other'

    def classify_gamma_angle(self, gamma):
        """Classify gamma angle into g+, g-, or t"""
        try:
            gamma = float(gamma)
            if 30 <= gamma <= 90:
                return 'g+'
            elif -90 <= gamma <= -30:
                return 'g-'
            elif 150 <= abs(gamma) <= 180:
                return 't'
        except (ValueError, TypeError):
            pass
        return 'other'

    def parse_dssr_file(self, dssr_file):
        """Parse DSSR output file and extract conformational parameters"""
        try:
            with open(dssr_file) as f:
                data = json.load(f)
        except Exception as e:
            print(f"Error reading DSSR file {dssr_file}: {e}")
            return None

        # Extract nucleotides data
        nts = data.get('nts', [])
        
        # Classify glycosidic bond orientations
        glycosidic = Counter()
        for nt in nts:
            chi = nt.get('chi', 'unknown')
            if chi != 'unknown':
                classification = self.classify_chi_angle(chi)
                glycosidic[classification] += 1
        
        # Classify sugar pucker conformations
        sugar_pucker = Counter()
        for nucleotide in nts:
            sugar_pucker_value = nucleotide.get('puckering', 'unknown')
            if sugar_pucker_value != 'unknown':
                classification = self.classify_sugar_pucker(sugar_pucker_value)
                sugar_pucker[classification] += 1
        
        # Classify gamma angles
        gamma_conf = Counter()
        for nt in nts:
            gamma = nt.get('gamma', 'unknown')
            if gamma != 'unknown':
                classification = self.classify_gamma_angle(gamma)
                gamma_conf[classification] += 1

        # Create dictionary representations (only counts)
        glycosidic_dict = {k: v for k, v in glycosidic.items() if v > 0}
        sugar_pucker_dict = {k: v for k, v in sugar_pucker.items() if v > 0}
        gamma_dict = {k: v for k, v in gamma_conf.items() if v > 0}

        parameters = {
            'glycosidic_bond': glycosidic_dict,
            'sugar_pucker': sugar_pucker_dict,
            'gamma_angles': gamma_dict
        }

        # Print detailed parameters with PDB ID
        pdb_id = os.path.basename(dssr_file).split('_')[0].upper()  # Extract and uppercase PDB ID
        print(f"\n{pdb_id}:")
        print(f"Glycosidic Bond Orientations: {glycosidic_dict}")
        print(f"Sugar Pucker Conformations: {sugar_pucker_dict}")
        print(f"Gamma Angle Conformations: {gamma_dict}")

        return parameters

    def store_parameters(self, identifier, parameters, is_predicted=True):
        """Store conformational parameters in database"""
        table_name = self.config.pred_table if is_predicted else self.config.ref_table
        id_column = 'exp_db_id' if is_predicted else 'PDBId'  # Correct case for column names
        
        # First check if the identifier exists
        check_query = f"SELECT {id_column} FROM {table_name} WHERE {id_column} = ?"
        self.database.cursor.execute(check_query, (identifier.upper(),))  # Convert identifier to uppercase
        if not self.database.cursor.fetchone():
            print(f"Warning: {id_column} {identifier.upper()} not found in {table_name}")
            return

        # print(f"\nStoring parameters for {identifier.upper()} in {table_name}")
        
        try:
            self.database.update_or_insert(
                table_name,
                ['RNA_GlycosidicBond', 'RNA_SugarPucker', 'RNA_GammaAngle'],
                (
                    str(parameters['glycosidic_bond']),
                    str(parameters['sugar_pucker']),
                    str(parameters['gamma_angles'])
                ),
                condition=f"{id_column} = '{identifier.upper()}'"  # Convert identifier to uppercase
            )
            # print("Parameters stored successfully")
        except Exception as e:
            print(f"Error updating {table_name} for {identifier.upper()}: {e}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python dssrConformation.py <dssr_folder>")
        print("Example: python dssrConformation.py ./dssr_files")
        sys.exit(1)

    dssr_folder = sys.argv[1]

    if not os.path.exists(dssr_folder):
        print(f"Error: Folder {dssr_folder} not found")
        sys.exit(1)

    if not os.path.isdir(dssr_folder):
        print(f"Error: {dssr_folder} is not a directory")
        sys.exit(1)

    processor = DSSRConformation()
    parent_folder = os.path.dirname(os.path.dirname(os.path.abspath(dssr_folder)))
    processor.config.set_config(parent_folder)
    processor.config.load_config()

    # Get all JSON files in the directory
    dssr_files = [f for f in os.listdir(dssr_folder) if f.endswith('_dssr.json')]
    
    if not dssr_files:
        print(f"Error: No DSSR JSON files found in {dssr_folder}")
        sys.exit(1)

    print(f"\nProcessing {len(dssr_files)} DSSR files...")
    
    for dssr_file in dssr_files:
        try:
            # Determine if it's predicted or experimental based on filename
            is_predicted = '_af_' in dssr_file
            
            # Extract identifier (xxxx) from xxxx_af_dssr.json or xxxx_pdb_dssr.json
            identifier = dssr_file.split('_')[0]
            
            # Process DSSR file
            parameters = processor.parse_dssr_file(os.path.join(dssr_folder, dssr_file))
            if parameters:
                processor.store_parameters(identifier, parameters, is_predicted=is_predicted)
                print(f"Processed: {dssr_file} ({'predicted' if is_predicted else 'experimental'})")
            
        except Exception as e:
            print(f"Error processing {dssr_file}: {e}")
            continue

    print("\nProcessing complete.")

if __name__ == "__main__":
    main() 