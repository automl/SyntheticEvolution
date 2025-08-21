import os
import sys
from pymol import cmd
from Bio.PDB import MMCIFParser, PDBIO
import tempfile

def get_PNG(file_name_cif_pdb, file_name_cif_af):
    cif_parser = MMCIFParser(QUIET=True)
    structure_reference = cif_parser.get_structure(file_name_cif_pdb, os.path.join(pdb_folder, file_name_cif_pdb))
    # truncate_chain_ids(structure_reference) TO-DO see above
    structure_predicted = cif_parser.get_structure(file_name_cif_af, os.path.join(af_folder, file_name_cif_af))

    # Convert structures to PDB format and save to temporary files
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file_ref:
        pdb_io = PDBIO()
        pdb_io.set_structure(structure_reference)
        pdb_io.save(temp_pdb_file_ref.name)
        temp_pdb_file_ref_path = temp_pdb_file_ref.name

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file_tgt:
        pdb_io.set_structure(structure_predicted)
        pdb_io.save(temp_pdb_file_tgt.name)
        temp_pdb_file_tgt_path = temp_pdb_file_tgt.name

    cmd.remove('solvent')
    cmd.remove('organic')  # Remove all ligands
    cmd.remove('ions')  # Remove all ions - not sufficient
    cmd.remove('resn MG')
    cmd.remove('resn K')
    cmd.remove('resn NA')
    cmd.remove('resn CA')
    cmd.remove('resn ZN')
    cmd.remove('resn CL')
    cmd.remove('resn SO4')

    joined_folder = os.path.join(os.path.dirname(af_folder), "joined")
    if not os.path.exists(joined_folder):
        os.makedirs(joined_folder)
    cmd.load(temp_pdb_file_ref_path, 'reference')
    cmd.load(temp_pdb_file_tgt_path, 'predicted')

    cmd.align(
        "predicted and polymer.nucleic and name P", #protein: polymer and name CA
        "reference and polymer.nucleic and name P",
        reset=1
    )

    cmd.hide('everything')
    cmd.show('cartoon', 'reference')
    cmd.show('cartoon', 'predicted')
    cmd.color('wheat', 'reference')
    cmd.color('green', 'predicted')
    cmd.bg_color('white')
    # cmd.set('ray_color', 'cmyk')  # Set color space for publication (use CMYK)
    cmd.zoom('reference')

    output_png_path = os.path.join(joined_folder, f'{pdb_id}_aligned_structures.png')
    cmd.png(output_png_path, width=800, height=600, dpi=300, ray=True)

    cmd.delete('reference')
    cmd.delete('predicted')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <reference_directory_path> <predicted_directory_path>")
        sys.exit(1)
    pdb_folder = sys.argv[1]
    af_folder = sys.argv[2]

    for file_name_cif_pdb in os.listdir(pdb_folder):
        if file_name_cif_pdb.endswith('.cif'):
            pdb_id = file_name_cif_pdb[:4]
            print(pdb_id)
            file_name_cif_af = None
            for file in os.listdir(af_folder): # Find the corresponding file in the af folder
                if file.startswith(f"fold_{pdb_id}_") and file.endswith('.cif'):
                    file_name_cif_af = file
                    break

            if not file_name_cif_af:
                raise FileNotFoundError(f"No corresponding file found in 'af' folder for PDB-ID {pdb_id}")

            get_PNG(file_name_cif_pdb, file_name_cif_af)

    cmd.quit()