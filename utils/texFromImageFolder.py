# def write_frame_header(f, title):
#     """Write the frame header with proper spacing"""
#     f.write(f"\\frametitle{{{title}}}\n")
#     f.write("\\hspace{20cm}\n")  # Add space after title
#     f.write("\\begin{minipage}{1.0\\textwidth}\\centering\n")  # Centered block
#
#
# def create_tex_file(image_folder, output_file):
#     """Create LaTeX file from images with proper formatting"""
#     with open(output_file, 'w') as f:
#         # Write document header
#         f.write("\\documentclass{beamer}\n")
#         f.write("\\usepackage{graphicx}\n")
#         f.write("\\begin{document}\n\n")
#
#         # Group images by their prefix
#         image_groups = {}
#         for filename in sorted(os.listdir(image_folder)):
#             if filename.endswith(('.png', '.jpg', '.pdf')):
#                 prefix = filename.split('_')[0]
#                 if prefix not in image_groups:
#                     image_groups[prefix] = []
#                 image_groups[prefix].append(filename)
#
#         # Process each group
#         for prefix, images in image_groups.items():
#             f.write("\\begin{frame}\n")
#
#             # Write frame header
#             write_frame_header(f, prefix)
#
#             # Calculate rows and columns needed
#             n_images = len(images)
#             if n_images <= 2:
#                 cols = n_images
#                 rows = 1
#             else:
#                 cols = 2
#                 rows = (n_images + 1) // 2
#
#             # Calculate image width based on number of columns
#             img_width = 0.48 if cols == 2 else 0.95
#
#             # Write images
#             for i, image in enumerate(images):
#                 if i > 0 and i % cols == 0:
#                     f.write("\\\\[0.5cm]\n")  # Add vertical space between rows
#                 elif i > 0:
#                     f.write("\\hspace{0.2cm}")  # Add horizontal space between images
#
#                 image_path = os.path.join(image_folder, image)
#                 f.write(f"\\includegraphics[width={img_width}\\textwidth]{{{image_path}}}")
#
#             # Close frame
#             f.write("\n\\end{minipage}\n")
#             f.write("\\end{frame}\n\n")
#
#         # Close document
#         f.write("\\end{document}\n")
#
#
# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print(f"Usage: python {os.path.basename(sys.argv[0])} <image_folder> <output_tex_file>")
#         sys.exit(1)
#
#     image_folder = sys.argv[1]
#     output_file = sys.argv[2]
#     create_tex_file(image_folder, output_file)

# import os
# import math

# mainfolder = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/results/plots"
# parent_folders = sorted([d for d in os.listdir(mainfolder) if os.path.isdir(os.path.join(mainfolder, d))])
# print("parent_folders", parent_folders)
# # subfolders = sorted([d for d in os.listdir(mainfolder) if os.path.isdir(os.path.join(mainfolder, d))])

# with open("/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/results/tex/slides.tex", "w") as f:
#     f.write("\\documentclass{beamer}\n")
#     f.write("\\usepackage{graphicx}\n")
#     f.write("\\usepackage{caption}\n")  # Use caption instead of subcaption
#     f.write("\\begin{document}\n")

#     # Iterate over each parent folder
#     for parent_folder in parent_folders:
#         parent_folder_path = os.path.join(mainfolder, parent_folder)
#         subfolders = sorted([d for d in os.listdir(parent_folder_path) if os.path.isdir(os.path.join(parent_folder_path, d))])
#         # Now, process subfolders A1, A2, etc., across all parent folders
#         for sub in subfolders:
#             images = sorted([img for img in os.listdir(os.path.join(parent_folder_path, sub)) if img.endswith(".png")])
#             print(images)

#             num_slides = math.ceil(len(images) / (1 if sub == "domain_metrics" or sub == "rna_metrics" else 4))  # One image per slide for "domain_metrics" and "rna_metrics"

#             # Iterate over the slides
#             for i in range(num_slides):
#                 # Create the new slide title by combining parent folder and subfolder name
#                 slide_title = f"{sub}_{parent_folder}".replace("_", r"\_")  # Escape underscores in the title

#                 f.write(f"\\begin{{frame}}\n\\frametitle{{{slide_title}}}\n")

#                 if sub == "domain_metrics" or sub == "rna_metrics":
#                     # One image per slide for domain_metrics and rna_metrics
#                     f.write("\\vspace{-0.5cm}\n")  # Adjust spacing between title and images
#                     f.write("\\begin{minipage}{1.0\\textwidth}\\centering\n")  # Fixed width, centered block

#                     # Select the next batch of images for this slide
#                     img_subset = images[i * 1:(i + 1) * 1]

#                     for j, img in enumerate(img_subset):
#                         img_escaped = img.replace("_", r"\_")  # Escape underscores in file name

#                         f.write(f"\\includegraphics[width=0.9\\linewidth]{{{mainfolder}/{parent_folder}/{sub}/{img}}}\n")  # Image path
#                         f.write("\\end{minipage}\\\\[1em]\n")  # New row after each image

#                 else:
#                     # Use a fixed-width layout for consistent image placement for other subfolders
#                     f.write("\\hspace{-2cm}\n")  # Adjust spacing between title and images
#                     f.write("\\begin{minipage}{1.0\\textwidth}\\centering\n")  # Fixed width, centered block

#                     # Select the next batch of images for this slide
#                     img_subset = images[i * 4:(i + 1) * 4]

#                     for j, img in enumerate(img_subset):
#                         img_escaped = img.replace("_", r"\_")  # Escape underscores in file name

#                         if j % 2 == 0:
#                             f.write("\\begin{minipage}{0.48\\textwidth}\\centering\n")  # First column

#                         f.write(f"\\includegraphics[width=0.9\\linewidth]{{{mainfolder}/{parent_folder}/{sub}/{img}}}\n")  # Image path

#                         if j % 2 == 1 or j == len(img_subset) - 1:
#                             f.write("\\end{minipage}\n")  # Close minipage for column

#                         if j == 1:  # Start a new row after 2 images
#                             f.write("\\\\[1em]\n")

#                     f.write("\\end{minipage}\n")  # End fixed-width block

#                 f.write("\\end{frame}\n")

#     f.write("\\end{document}\n")

# print("Generated slides.tex! Now compile with pdflatex.")

import os
import math

mainfolder = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/results/plots"
parent_folders = sorted([d for d in os.listdir(mainfolder) if os.path.isdir(os.path.join(mainfolder, d))])
print("parent_folders", parent_folders)

# Get all unique subfolder names across all parent folders
all_subfolders = set()
for parent_folder in parent_folders:
    parent_folder_path = os.path.join(mainfolder, parent_folder)
    subfolders = [d for d in os.listdir(parent_folder_path) if os.path.isdir(os.path.join(parent_folder_path, d))]
    all_subfolders.update(subfolders)
all_subfolders = sorted(list(all_subfolders))

with open("/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/results/tex/slides.tex", "w") as f:
    f.write("\\documentclass{beamer}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{caption}\n")  # Use caption instead of subcaption
    f.write("\\begin{document}\n")

    # Process each subfolder type (A1, A2, etc.) across all parent folders
    for sub in all_subfolders:
        # Process this subfolder type for each parent folder
        for parent_folder in parent_folders:
            parent_folder_path = os.path.join(mainfolder, parent_folder)
            sub_path = os.path.join(parent_folder_path, sub)
            
            # Skip if this subfolder doesn't exist in current parent folder
            if not os.path.exists(sub_path):
                continue
                
            images = sorted([img for img in os.listdir(sub_path) if img.endswith(".png")])
            print(images)

            num_slides = math.ceil(len(images) / (1 if sub == "domain_metrics" or sub == "rna_metrics" or sub == "correlation_matrix" else 4))

            # Rest of the code remains the same...
            for i in range(num_slides):
                slide_title = f"{sub}_{parent_folder}".replace("_", r"\_")

                f.write(f"\\begin{{frame}}\n\\frametitle{{{slide_title}}}\n")

                if sub == "domain_metrics" or sub == "rna_metrics" or sub == "correlation_matrix":
                    f.write("\\vspace{-0.5cm}\n")
                    f.write("\\begin{minipage}{1.0\\textwidth}\\centering\n")

                    img_subset = images[i * 1:(i + 1) * 1]

                    for j, img in enumerate(img_subset):
                        img_escaped = img.replace("_", r"\_")
                        f.write(f"\\includegraphics[width=0.9\\linewidth]{{{mainfolder}/{parent_folder}/{sub}/{img}}}\n")
                        f.write("\\end{minipage}\\\\[1em]\n")

                else:
                    f.write("\\hspace{-2cm}\n")
                    f.write("\\begin{minipage}{1.0\\textwidth}\\centering\n")

                    img_subset = images[i * 4:(i + 1) * 4]

                    for j, img in enumerate(img_subset):
                        img_escaped = img.replace("_", r"\_")

                        if j % 2 == 0:
                            f.write("\\begin{minipage}{0.48\\textwidth}\\centering\n")

                        f.write(f"\\includegraphics[width=0.9\\linewidth]{{{mainfolder}/{parent_folder}/{sub}/{img}}}\n")

                        if j % 2 == 1 or j == len(img_subset) - 1:
                            f.write("\\end{minipage}\n")

                        if j == 1:
                            f.write("\\\\[1em]\n")

                    f.write("\\end{minipage}\n")

                f.write("\\end{frame}\n")

    f.write("\\end{document}\n")

print("Generated slides.tex! Now compile with pdflatex.")