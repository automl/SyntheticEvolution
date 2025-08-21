import os
import sys
from pymol import cmd
from Bio.PDB import MMCIFParser, PDBIO
import tempfile


def get_PNG(file_name_cif_pdb, file_name_cif_af, file_name_cif_informer):
    cif_parser = MMCIFParser(QUIET=True)
    structure_reference = cif_parser.get_structure(file_name_cif_pdb, os.path.join(pdb_folder, file_name_cif_pdb))
    structure_predicted = cif_parser.get_structure(file_name_cif_af, os.path.join(af_folder, file_name_cif_af))
    structure_informer = cif_parser.get_structure(file_name_cif_informer,
                                                  os.path.join(informer_folder, file_name_cif_informer))

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file_ref:
        pdb_io = PDBIO()
        pdb_io.set_structure(structure_reference)
        pdb_io.save(temp_pdb_file_ref.name)
        temp_pdb_file_ref_path = temp_pdb_file_ref.name

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file_tgt:
        pdb_io.set_structure(structure_predicted)
        pdb_io.save(temp_pdb_file_tgt.name)
        temp_pdb_file_tgt_path = temp_pdb_file_tgt.name

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file_inf:
        pdb_io.set_structure(structure_informer)
        pdb_io.save(temp_pdb_file_inf.name)
        temp_pdb_file_inf_path = temp_pdb_file_inf.name

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
    cmd.load(temp_pdb_file_inf_path, 'informer')

    cmd.align(
        "predicted and polymer.nucleic and name P",
        "reference and polymer.nucleic and name P",
        reset=1
    )

    cmd.align(
        "informer and polymer.nucleic and name P",
        "reference and polymer.nucleic and name P",
        reset=0  # Do not reset after this alignment
    )

    cmd.hide('everything')
    cmd.show('cartoon', 'reference')
    cmd.show('cartoon', 'predicted')
    cmd.show('cartoon', 'informer')
    cmd.color('wheat', 'reference')
    cmd.color('green', 'predicted')
    cmd.color('blue', 'informer')
    cmd.bg_color('white')
    cmd.zoom('reference')

    output_png_path = os.path.join(joined_folder, f'{pdb_id}_aligned_structures.png')
    cmd.png(output_png_path, width=800, height=600, dpi=300, ray=True)

    cmd.delete('reference')
    cmd.delete('predicted')
    cmd.delete('informer')


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            f"Usage: python {os.path.basename(sys.argv[0])} <reference_directory_path> <af_directory_path> <informer_directory_path>")
        sys.exit(1)
    pdb_folder = sys.argv[1]
    af_folder = sys.argv[2]
    informer_folder = sys.argv[3]

    for file_name_cif_pdb in os.listdir(pdb_folder):
        if file_name_cif_pdb.endswith('.cif'):
            pdb_id = file_name_cif_pdb[:4]
            print(pdb_id)
            file_name_cif_af = None
            file_name_cif_informer = None

            for file in os.listdir(af_folder):
                if file.startswith(f"{pdb_id}_af") and file.endswith('.cif'):
                    file_name_cif_af = file
                    break

            for file in os.listdir(informer_folder):
                if file.startswith(f"{pdb_id}_rnainformer") and file.endswith('.cif'):
                    file_name_cif_informer = file
                    break

            if not file_name_cif_af:
                raise FileNotFoundError(f"No corresponding file found in 'af' folder for PDB-ID {pdb_id}")
            if not file_name_cif_informer:
                raise FileNotFoundError(f"No corresponding file found in 'informer' folder for PDB-ID {pdb_id}")

            get_PNG(file_name_cif_pdb, file_name_cif_af, file_name_cif_informer)

    cmd.quit()
