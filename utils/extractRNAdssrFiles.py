import os
import shutil

# Specify the source directory where the files are located
source_directory = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/a2021_rnaOnly/aligned"

# Create the target directory called 'DS' if it doesn't already exist
target_directory = os.path.join(source_directory, "DS")
if not os.path.exists(target_directory):
    os.makedirs(target_directory)

# Loop through all files in the source directory
for filename in os.listdir(source_directory):
    # Check if the file ends with '_rna.pdb'
    if filename.endswith("_rna.pdb"):
        # Define the full file paths
        source_file = os.path.join(source_directory, filename)
        target_file = os.path.join(target_directory, filename)

        # Copy the file to the target directory
        shutil.copy(source_file, target_file)
        print(f"Copied: {filename} to {target_directory}")

print("Extraction completed.")
