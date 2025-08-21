import os
from parseAFcif2DB import parse_cif_file as parse_af_cif_file, insert_pred_data
from parsePDBcif2DB import parse_cif_file as parse_pdb_cif_file, insert_exp_data

class CIFParser:
    def __init__(self, cif_file_path: str, database_path: str):
        self.cif_file_path = cif_file_path
        self.database_path = database_path

    def process_cif_folder(self, cif_file: str): # TO-DO: merge Alphafold and PDB
        folder_name = os.path.basename(os.path.dirname(cif_file)).lower()

        if folder_name == 'af':
            # Use parseAFcif2DB methods
            data = parse_af_cif_file(cif_file)
            insert_pred_data(self.database_path, data)
        elif folder_name == 'pdb':
            # Use parsePDBcif2DB methods
            data = parse_pdb_cif_file(cif_file)
            insert_exp_data(self.database_path, data)
        else:
            print(f"Unknown folder type for file: {cif_file}")

    def process_cif_files(self):
        cif_files = [os.path.join(self.cif_file_path, f) for f in os.listdir(self.cif_file_path) if f.endswith('.cif')]
        for cif_file in cif_files:
            self.process_cif_folder(cif_file)

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <directory_path>")
    else:
        directory_path = sys.argv[1]
        db_name = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'database', 'rbpDatabase.db')
        parser = CIFParser(cif_file_path=directory_path, database_path=db_name)
        parser.process_cif_files()