import os
import sys
import subprocess
from database.startConfig import StartConfig

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # One level above utils package.

def execute_command(script_path, *args):
    command = f"{script_path} {' '.join(args)}"
    print(f"Executing: {command}")
    subprocess.run(command, shell=True, check=True)

def main():
    if len(sys.argv) != 3:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <parent_folder> <path_to_tsv>")
        sys.exit(1)

    parent_folder = sys.argv[1]
    path_to_tsv = sys.argv[2]

    database_path = "database"
    eval_path = "evaluation"
    init_command = os.path.join(project_root, database_path, "startConfig.py")
    execute_command(f"python {init_command} {parent_folder}")
    config = StartConfig()
    af_path = config.pred_folder
    csv_name = f"{config.pred_table}.csv"
    results_path_csv = os.path.join(project_root, "results", csv_name)

    execute_command(f"python {os.path.join(project_root, database_path, 'insertRNAcompete_selex.py')} {path_to_tsv}")
    execute_command(f"python {os.path.join(project_root, database_path, 'parseAFcif2DB.py')} {af_path}")
    execute_command(f"python {os.path.join(project_root, eval_path, 'extractInteractions.py')} {af_path}")
    execute_command(f"python {os.path.join(project_root, eval_path, 'extractInterfaceArea.py')} {af_path}")
    execute_command(
        f"python {os.path.join(project_root, eval_path, 'estimateBindingAffinity.py')} {config.pred_table}")
    execute_command(f"python {os.path.join(project_root, eval_path, 'af3Metrics_new.py')} {af_path} {'--no-plots'}")
    csv_command = os.path.join(project_root, database_path, "exportDB2CSV.py")
    execute_command(f"python {csv_command} {config.pred_table}")
    # ADD: ALTER TABLE pred_protein_rna ADD COLUMN RNAmotif_score FLOAT;
    # ALTER TABLE pred_protein_rna ADD COLUMN Secondary_Elements TEXT;
    # ALTER TABLE pred_protein_rna ADD COLUMN motif_similarity_score FLOAT;
    # ALTER TABLE pred_protein_rna ADD COLUMN RNAMotif_inExp FLOAT;
    # BEFORE rnaCompete
    execute_command(f"python {os.path.join(project_root, eval_path, 'rnaCompete.py')}") # {'--plot-only'}")

if __name__ == "__main__":
    main()

#DATA-PREP
# 1. EXTRACT SUBFILES TO af folder
# import os
# import shutil
#
# parent_folder = "/path/to/parent_folder"
# output_folder = os.path.join(parent_folder, "af")
#
# os.makedirs(output_folder, exist_ok=True)
#
# # Loop through all subdirectories
# for subdir, _, files in os.walk(parent_folder):
#     if subdir == output_folder:
#         continue  # Skip the output folder itself
#
#     for file in files:
#         source_path = os.path.join(subdir, file)
#         destination_path = os.path.join(output_folder, file)
#
#         # Ensure no overwriting: Add a suffix if filename already exists
#         counter = 1
#         while os.path.exists(destination_path):
#             name, ext = os.path.splitext(file)
#             destination_path = os.path.join(output_folder, f"{name}_{counter}{ext}")
#             counter += 1
#
#         # Move file
#         shutil.move(source_path, destination_path)
#
# print(f"All files have been moved to: {output_folder}")

#######################
# 2. RENAME AF FILES:
# import os
# import re
#
# # Define the folder containing the files
# folder_path = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/af"  # Change this to your actual folder path
#
# # Iterate through all files in the folder
# for file_name in os.listdir(folder_path):
#     file_path = os.path.join(folder_path, file_name)
#
#     if os.path.isfile(file_path):
#         # Extract the base name (x) before "_model", "_confidences", etc.
#         match = re.match(r"(.+?)_(model|confidences|summary_confidences)\.(cif|json)", file_name)
#
#         if match:
#             base_name, file_type, extension = match.groups()
#
#             # Define the new filename based on the pattern
#             if file_type == "model":
#                 new_name = f"{base_name}_s1_model_0.{extension}"
#             elif file_type == "confidences":
#                 new_name = f"{base_name}_s1_full_data_0.{extension}"
#             elif file_type == "summary_confidences":
#                 new_name = f"{base_name}_s1_summary_confidences_0.{extension}"
#             else:
#                 continue  # Skip if the file doesn't match the expected types
#
#             # Rename the file
#             new_path = os.path.join(folder_path, new_name)
#             os.rename(file_path, new_path)
#             print(f"Renamed: {file_name} -> {new_name}")
#
# print("Renaming completed!")

###########
# 3. REMOVE EXCESSIVE UNDERSCORES:
# import os
# import re
#
# # Define the parent directory where the subfolders are located
# parent_dir = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/af"
#
# # Loop through all files in subdirectories
# for root, _, files in os.walk(parent_dir):
#     for filename in files:
#         if "_s1_" in filename:
#             # Split the filename at "_s1_"
#             parts = filename.split("_s1_", 1)
#             if len(parts) == 2:
#                 x = parts[0]  # The part before "_s1_"
#                 rest = parts[1]  # The part after "_s1_"
#
#                 # Remove all underscores from `x`
#                 x_clean = re.sub(r"_", "", x)
#
#                 # Create the new filename
#                 new_filename = f"{x_clean}_s1_{rest}"
#
#                 # Get full paths
#                 old_path = os.path.join(root, filename)
#                 new_path = os.path.join(root, new_filename)
#
#                 # Rename the file
#                 os.rename(old_path, new_path)
#                 print(f"Renamed: {filename} -> {new_filename}")

"3a for SELEX 1st table 30 samples"
# import os
#
# folder_path = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/af"
#
# for filename in os.listdir(folder_path):
#     if filename.startswith("fold") and not filename.startswith("fold_"):
#         new_filename = "fold_" + filename[4:]  # Insert underscore after "fold"
#         old_path = os.path.join(folder_path, filename)
#         new_path = os.path.join(folder_path, new_filename)
#         os.rename(old_path, new_path)
#         print(f"Renamed: {filename} -> {new_filename}")
#
# print("Renaming completed.")


#######
# 4. ADD fold_
# import os
#
# # Define the directory containing the files
# folder_path = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/af"
#
# # Loop through all files in the folder
# for filename in os.listdir(folder_path):
#     old_path = os.path.join(folder_path, filename)
#
#     # Ensure it's a file (not a directory)
#     if os.path.isfile(old_path):
#         new_filename = f"fold_{filename}"  # Add "fold_" prefix
#         new_path = os.path.join(folder_path, new_filename)
#
#         os.rename(old_path, new_path)
#         print(f"Renamed: {filename} -> {new_filename}")

#########
# 5. RNAMotif_inExp
# sqlite3 rbpDatabase.db
# ALTER TABLE pred_protein_rna ADD COLUMN RNAMotif_inExp;