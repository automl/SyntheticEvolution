import os
import shutil

# Define the main directory
folder_A = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/af"
processed_folder = os.path.join(folder_A, "processed")

# Create the 'processed' folder if it doesn't exist
os.makedirs(processed_folder, exist_ok=True)

# Loop over files in folder A
for filename in os.listdir(folder_A):
    if filename.endswith("_complex_tm_output.txt"):
        # Extract 'xxxx' from the file name
        xxxx = filename.replace("_complex_tm_output.txt", "")

        # Search for all files that contain 'xxxx'
        for other_file in os.listdir(folder_A):
            if xxxx in other_file and os.path.isfile(os.path.join(folder_A, other_file)):
                # Move the file to 'processed' subfolder
                src = os.path.join(folder_A, other_file)
                dst = os.path.join(processed_folder, other_file)
                shutil.move(src, dst)
                print(f"Moved: {other_file}")

print("Done.")
