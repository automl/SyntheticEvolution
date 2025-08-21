# from Bio.PDB import MMCIFParser, PDBIO
# from pymol import cmd
#
# cif_parser = MMCIFParser(QUIET=True)
# structure_reference = cif_parser.get_structure(file_name_cif_pdb, os.path.join(pdb_folder, file_name_cif_pdb))
# structure_predicted = cif_parser.get_structure(file_name_cif_af, os.path.join(af_folder, file_name_cif_af))
#
# # Convert structures to PDB format and save to temporary files
# with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file_ref:
#     pdb_io = PDBIO()
#     pdb_io.set_structure(structure_reference)
#     pdb_io.save(temp_pdb_file_ref.name)
#     temp_pdb_file_ref_path = temp_pdb_file_ref.name
#
# with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file_tgt:
#     pdb_io.set_structure(structure_predicted)
#     pdb_io.save(temp_pdb_file_tgt.name)
#     temp_pdb_file_tgt_path = temp_pdb_file_tgt.name
#
# cmd.load(temp_pdb_file_ref_path, 'referenceA')
# cmd.load(temp_pdb_file_tgt_path, 'predictedA')
# alignment_result = cmd.align('predictedA and polymer and name CA', 'referenceA and polymer and name CA',
#                              quiet=0)
# complex_rmsd = round(alignment_result[0], 2)

# import os
#
# folder_path = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/af/rnaformer_config23_all_data"  # Change this to your actual path
#
# # Get list of filenames in the folder
# filenames = os.listdir(folder_path)
#
# # Extract the part before '_rnaformer2025_config23'
# pdb_ids = [name.split('_rnaformer2025_config23')[0] for name in filenames if '_rnaformer2025_config23' in name]
#
# print(pdb_ids)

# import os
#
# # List of PDB IDs to match (case-insensitive, expected in lowercase in filenames)
# pdb_ids = [
#     "1S03", "4OOG", "3CUL", "6CAE", "2ZJQ", "3HHN", "4V8D", "5AH5", "4V84",
#     "3EGZ", "2CZJ", "2ZZM", "2OZB", "3AKZ", "5DDP", "5M73", "5ZTM", "6SY4",
#     "1H4Q", "5O7H", "6SY6", "6DCB", "1F7U", "2DEU", "3ZJV", "1ZHO", "4YBB",
#     "4PKD", "3EPH", "6LAS", "4V9K", "5WLH", "6IV8", "4V8N", "6AAY", "6MWN",
#     "6JE9", "4TVX", "4C4W", "4X4V", "5D6G", "4WJ4", "4R8I", "5B63", "5WTK", "4ATO",
# ]
# #     "6R47", "4RUM", "1f1t", "6dvk", "6dlq", "3npq", "1u9s", "5t83"!, "4jf2", "3vrs", "1kpy"
#
# # Convert to lowercase and format to match _xxxx_ in filenames
# match_substrings = [f"_{pdb_id.lower()}_" for pdb_id in pdb_ids]
#
# # Folder to clean (replace with your target folder)
# target_folder = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/rna_rna/af"
#
# # Iterate through files and remove matches
# for filename in os.listdir(target_folder):
#     file_path = os.path.join(target_folder, filename)
#     if os.path.isfile(file_path):
#         if any(substr in filename for substr in match_substrings):
#             print(f"Removing: {filename}")
#             os.remove(file_path)


import os
import shutil

parent_folder = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/rna_rna/a2021_rnaformer25"
combined_db_folder = os.path.join(parent_folder, "combined_database")

# Create the destination folder if it doesn't exist
os.makedirs(combined_db_folder, exist_ok=True)

# Loop through all subdirectories
for folder_name in os.listdir(parent_folder):
    folder_path = os.path.join(parent_folder, folder_name)

    # Proceed only if it's a directory and starts with "config"
    if os.path.isdir(folder_path):# and folder_name.startswith("config"):
        for file in os.listdir(folder_path):
            if file.endswith(".db"):
                source_path = os.path.join(folder_path, file)
                new_file_name = f"{folder_name}.db"
                destination_path = os.path.join(combined_db_folder, new_file_name)

                shutil.copy2(source_path, destination_path)
                print(f"Copied and renamed {file} → {destination_path}")

