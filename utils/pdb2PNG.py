import os
from PIL import Image
from pymol import cmd

input_folder = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/b2021_RF2NA"  # Replace with your PDB folder path
output_images_folder = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/b2021_RF2NA"  # Folder to save PNG images
output_pdf_file = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/b2021_RF2NA/pdb_images.pdf"  # Output PDF file name

# List of PDB files in the input folder
pdb_files = [f for f in os.listdir(input_folder) if f.endswith(".pdb")]

# List to store paths of generated PNG images
image_files = []

# Process each PDB file
for pdb_file in pdb_files:
    pdb_path = os.path.join(input_folder, pdb_file)
    image_name = os.path.splitext(pdb_file)[0] + ".png"  # Image name
    image_path = os.path.join(output_images_folder, image_name)

    # Load and render PDB file in PyMOL
    print(f"Processing {pdb_file}...")
    cmd.reinitialize()  # Clear previous session
    cmd.load(pdb_path)
    # cmd.hide("everything")  # Hide default representation
    # # Show proteins in cartoon representation
    # cmd.select("proteins", "polymer.protein")
    # cmd.show("cartoon", "proteins")
    #
    # # Show nucleotides in both cartoon and sticks
    # cmd.select("nucleotides", "resn A+T+G+C+U")
    # cmd.show("cartoon", "nucleotides")  # Cartoon for nucleotides
    # cmd.show("sticks", "nucleotides")  # Sticks for nucleotides
    # cmd.show("sticks", "organic")
    # cmd.bg_color('white')

    # Apply colors for clarity
    cmd.util.cnc()  # Color chains for clarity
    cmd.zoom("all")  # Zoom to fit the whole structure
    cmd.png(image_path, width=800, height=800, dpi=300)  # Save as PNG
    cmd.delete("all")  # Clear loaded structure

    # Add the image to the list
    image_files.append(image_path)

# Convert all PNG images to a single PDF file
if image_files:
    print("Generating PDF file...")
    images = [Image.open(img).convert("RGB") for img in image_files]
    images[0].save(output_pdf_file, save_all=True, append_images=images[1:])
    print(f"PDF saved to {output_pdf_file}")
else:
    print("No images generated. Please check the input folder.")

print("Processing complete.")
