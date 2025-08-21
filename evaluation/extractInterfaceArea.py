from Bio.PDB import MMCIFParser, PDBIO
import freesasa
import tempfile
import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from structures.rna import RNA
from structures.dna import DNA
from structures.protein import Protein
# from pymol import cmd, stored
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from pymol import cmd, stored
from database.databaseMethods import DatabaseMethods
# import glob

database_methods = DatabaseMethods()

warnings.simplefilter('ignore', PDBConstructionWarning)

if len(sys.argv) < 2:
    print("Usage: python extractInterfaceArea.py /path/to/cif_files/")
    sys.exit(1)
cif_directory = sys.argv[1]
cif_parser = MMCIFParser()

def convert_to_chain_groups(chain_ids):
    # Remove brackets and quotes from the string
    cleaned_chain_ids = chain_ids.replace('[', '').replace(']', '').replace("'", '').replace('"', '').replace(' ', '')
    # Split into a list by commas and join without commas
    chain_ids_list = cleaned_chain_ids.split(',')
    return ''.join(chain_ids_list)

# https://freesasa.github.io/1.1/Python.html https://freesasa.github.io/doxygen/Geometry.html
def extract_interface_area(file_name_cif, file_path):
    structure = cif_parser.get_structure(file_name_cif, file_path)

    # Write the structure to a PDB format using PDBIO and a temporary file
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_pdb_file:
        # pdb_io = PDBIO()
        # pdb_io.set_structure(structure)
        # pdb_io.save(temp_pdb_file.name)
        # fs_chains = set()

        # Check for invalid chain IDs and handle FreeSASA structure creation
        try:
            pdb_io = PDBIO()
            pdb_io.set_structure(structure)
            pdb_io.save(temp_pdb_file.name)
            fs_chains = set()
            # First check if any chain ID is longer than one character
            for model in structure:
                for chain in model:
                    if len(chain.id) > 1:
                        print(f"Warning: Skipping {file_name_cif} - Invalid chain ID format: {chain.id}")
                        return None

            # If all chain IDs are valid, proceed with FreeSASA
            fs_structure = freesasa.Structure(temp_pdb_file.name)
            for i in range(fs_structure.nAtoms()):
                fs_chains.add(fs_structure.chainLabel(i))
            
            pdb_chains = set(chain.id for model in structure for chain in model)
            missing_chains = pdb_chains - fs_chains

        except Exception as e:
            print(f"Warning: Skipping {file_name_cif} - FreeSASA error: {e}")
            return None

        aap_protein_rna = {}
        aap_protein_dna = {}
        print(file_name_cif)
        file = os.path.splitext(file_name_cif)[0]
        if len(file) in [4, 6]:
            rna = RNA.get_rna_from_db(id=file.upper())
            protein = Protein.get_protein_from_db(id=file.upper())
            dna = DNA.get_dna_from_db(id=file.upper())
        elif len(file) > 6:
            parts = file_name_cif.split('_')
            rna = RNA.get_rna_from_db(id=parts[1].upper(), file_name=file_name_cif)
            protein = Protein.get_protein_from_db(id=parts[1].upper(), file_name=file_name_cif)
            dna = DNA.get_dna_from_db(id=parts[1].upper(), file_name=file_name_cif)

        if not missing_chains:
            # Calculate total SASA
            result = freesasa.calc(fs_structure)
            total_sasa_complex = result.totalArea()
            print(f"Total SASA Complex: {total_sasa_complex:.2f} Å²")

            # chain_ids = []
            if rna.is_rna:
                rna_chain_ids = rna.get_rna_chain_IDs()
                # print("rna", rna_chain_ids, rna.is_single_copy_rna())
                rna_chain_number = rna.get_rna_chain_number()
                selections = freesasa.selectArea({'RNA, resn a+u+g+c'}, fs_structure, result)
                sasa_rna = selections.get('RNA', 0)
                print(f"Complexed RNA SASA: {sasa_rna:.2f} Å²")

            if protein.is_protein:
                protein_chain_ids = protein.get_protein_chain_IDs()
                # print("protein", protein_chain_ids, protein.is_single_copy_protein())
                # print(protein_chain_ids, isinstance(protein_chain_ids, str))
                protein_chain_number = protein.get_protein_chain_number()
                selections = freesasa.selectArea({'Protein, resn ala+arg+asn+asp+cys+gln+glu+gly+his+ile+leu+lys+met+phe+pro+ser+thr+trp+tyr+val'},
                    fs_structure, result)
                sasa_protein = selections.get('Protein', 0)
                print(f"Complexed Protein SASA: {sasa_protein:.2f} Å²")
                # print("protein_chain_ids", protein_chain_ids)
                # if not protein.is_single_copy_protein():
                #     # chain_ids_list = ast.literal_eval(protein_chain_ids)
                #     chain_ids_list = protein_chain_ids.split(',')
                #     protein_chain_ids = chain_ids_list[0] # Get the first item from the list
                # chain_ids.append(protein_chain_ids)
                # print("protein chain_ids", chain_ids)
            if dna.is_dna:
                dna_chain_ids = dna.get_dna_chain_IDs()
                # print("dna", dna_chain_ids, dna.is_single_copy_dna())
                dna_chain_number = dna.get_dna_chain_number()
                selections = freesasa.selectArea({'DNA, resn da+dt+dg+dc'}, fs_structure, result)
                sasa_dna = selections.get('DNA', 0)
                print(f"Complexed DNA SASA: {sasa_dna:.2f} Å²")
                # if not dna.is_single_copy_dna():
                #     # chain_ids_list = ast.literal_eval(dna_chain_ids)
                #     chain_ids_list = dna_chain_ids.split(',')
                #     dna_chain_ids = chain_ids_list[0] # Get the first item from the list
                # chain_ids.append(dna_chain_ids)
            if (dna.is_dna and dna.is_single_copy_dna()) and (rna.is_rna and rna.is_single_copy_rna()):
                dna_chain_groups = convert_to_chain_groups(dna_chain_ids)
                rna_chain_groups = convert_to_chain_groups(rna_chain_ids)
                chain_groups = dna_chain_groups + rna_chain_groups
                structureArray = freesasa.structureArray(temp_pdb_file.name,
                                                         options={"separate-chains": False,
                                                                  "chain-groups": chain_groups})
                for model_index, model in enumerate(structureArray):
                    if model_index == 1:  # Consider only the second model (model_index = 0 is the total area)
                        result = freesasa.calc(model)
                        sasa_value = result.totalArea()
                        complexed_rna_dna_interface_area = sasa_value
                        print(f"Uncomplexed complexed RNA-DNA Chain:{chain_groups}, SASA: {sasa_value:.2f} Å²")

            # https://github.com/freesasa/freesasa-python/issues/7 https://github.com/freesasa/freesasa-python/issues/8
            rna_chain_sasa_dict = {}
            protein_chain_sasa_dict = {}
            dna_chain_sasa_dict = {}

            if (rna.is_rna and rna.is_single_copy_rna()) or (protein.is_protein and protein.is_single_copy_protein())\
                or (dna.is_dna and dna.is_single_copy_dna()):
                structureArray = freesasa.structureArray(temp_pdb_file.name,
                                                         options={'separate-chains': True, 'separate-models': True})
                for model in structureArray:
                    chain_label = model.chainLabel(1)
                    result = freesasa.calc(model)
                    sasa_value = result.totalArea()
                    # print(rna.is_single_copy_rna(), rna.rna_sequence)
                    # print(protein.is_single_copy_protein(), protein.protein_sequence)
                    if (rna.is_rna and rna.is_single_copy_rna()) and (chain_label in rna_chain_ids):
                        print(f"Uncomplexed RNA Chain:{chain_label}, SASA: {sasa_value:.2f} Å²")
                        # rna_sasa_value = sasa_value
                        rna_chain_sasa_dict[chain_label] = sasa_value
                    elif (protein.is_protein and protein.is_single_copy_protein()) and (chain_label in protein_chain_ids):
                        print(f"Uncomplexed protein Chain:{chain_label}, SASA: {sasa_value:.2f} Å²")
                        # protein_sasa_value = sasa_value
                        protein_chain_sasa_dict[chain_label] = sasa_value
                    elif (dna.is_dna and dna.is_single_copy_dna()) and (chain_label in dna_chain_ids):
                        print(f"Uncomplexed DNA Chain:{chain_label}, SASA: {sasa_value:.2f} Å²")
                        # dna_sasa_value = sasa_value
                        dna_chain_sasa_dict[chain_label] = sasa_value

                # http://www.ysbl.york.ac.uk/ccp4bb/2000/msg00654.html
                # calculate the accessible surface area for
                # each of the individual proteins and than calculate the surface of the
                # complex. The buried surface area is S(1) + S(2) - S(1-2) and the
                # interaction surface area is half of this value.
                if (rna.is_rna and rna.is_single_copy_rna()) and (protein.is_protein and protein.is_single_copy_protein()):
                    interface_area = sum(rna_chain_sasa_dict.values()) + sum(protein_chain_sasa_dict.values()) - total_sasa_complex
                    max_chain_number = max(rna_chain_number, protein_chain_number)
                    protein_rna_interface_area = interface_area/max(max_chain_number, 2)
                    print(f"Interface Protein/RNA Area: {protein_rna_interface_area:.2f} Å²")

                if (dna.is_dna and dna.is_single_copy_dna()) and (protein.is_protein and protein.is_single_copy_protein()):
                    interface_area = sum(dna_chain_sasa_dict.values()) + sum(protein_chain_sasa_dict.values()) - total_sasa_complex
                    max_chain_number = max(dna_chain_number, protein_chain_number)
                    protein_dna_interface_area = interface_area / max(max_chain_number, 2)
                    print(f"Interface Protein/DNA Area: {protein_dna_interface_area:.2f} Å²")

                if (rna.is_rna and rna.is_single_copy_rna()) and (dna.is_dna and dna.is_single_copy_dna()):
                    interface_area = sum(dna_chain_sasa_dict.values()) + sum(rna_chain_sasa_dict.values()) - complexed_rna_dna_interface_area
                    max_chain_number = max(rna_chain_number, dna_chain_number)
                    rna_dna_interface_area = interface_area / max(max_chain_number, 2)
                    print(f"Interface RNA/DNA Area: {rna_dna_interface_area:.2f} Å²")

                if (rna.is_rna and rna.is_single_copy_rna()) and not protein.is_protein and not dna.is_dna:
                    interface_area = sum(rna_chain_sasa_dict.values()) - total_sasa_complex
                    rna_rna_interface_area = interface_area/max(rna_chain_number, 2)
                    print(f"Interface RNA/RNA Area: {rna_rna_interface_area:.2f} Å²")

            if (rna.is_rna and not rna.is_single_copy_rna()) or (protein.is_protein and not protein.is_single_copy_protein()) \
                    or (dna.is_dna and not dna.is_single_copy_dna()):
                if (protein.is_protein and not protein.is_single_copy_protein()):
                    # print(protein_chain_ids, protein.is_single_copy_protein())
                    chain_groups = convert_to_chain_groups(protein_chain_ids)
                    structureArray = freesasa.structureArray(temp_pdb_file.name,
                                                             options={"separate-chains": False,
                                                                      "chain-groups": chain_groups})
                    for model_index, model in enumerate(structureArray):
                        if model_index == 1:  # Consider only the second model (model_index = 0 is the total area)
                            result = freesasa.calc(model)
                            sasa_value = result.totalArea()
                            protein_interface_area = sasa_value
                            print(f"Uncomplexed Protein Chain: {chain_groups}, SASA: {sasa_value:.2f} Å²")

                if (rna.is_rna and not rna.is_single_copy_rna()):
                    chain_groups = convert_to_chain_groups(rna_chain_ids)
                    structureArray = freesasa.structureArray(temp_pdb_file.name,
                                                             options={"separate-chains": False,
                                                                      "chain-groups": chain_groups})
                    for model_index, model in enumerate(structureArray):
                        if model_index == 1:  # Consider only the second model (model_index = 0 is the total area)
                            result = freesasa.calc(model)
                            sasa_value = result.totalArea()
                            rna_interface_area = sasa_value
                            print(f"Uncomplexed RNA Chain: {chain_groups}, SASA: {sasa_value:.2f} Å²")
                        # else:
                        # for chain in missing_chains:
                        #     if chain in chain_groups:
                                # Problem only not missing chains are considered: not realistic, too low interface ratio
                                # chain_list = ast.literal_eval(rna_chain_ids)
                                # rna_chain_ids = set(chain_list)
                                # chain_groups = rna_chain_ids - missing_chains
                                # chain_groups = ', '.join(chain_groups)

                if (dna.is_dna and not dna.is_single_copy_dna()):
                    chain_groups = convert_to_chain_groups(dna_chain_ids)
                    structureArray = freesasa.structureArray(temp_pdb_file.name,
                                                             options={"separate-chains": False,
                                                                      "chain-groups": chain_groups})
                    for model_index, model in enumerate(structureArray):
                        if model_index == 1:  # Consider only the second model (model_index = 0 is the total area)
                            result = freesasa.calc(model)
                            sasa_value = result.totalArea()
                            dna_interface_area = sasa_value
                            print(f"Uncomplexed DNA Chain:{chain_groups}, SASA: {sasa_value:.2f} Å²")

                if (dna.is_dna and dna.is_single_copy_dna()) and (rna.is_rna and not rna.is_single_copy_rna()):
                    dna_chain_groups = ''.join(dna_chain_sasa_dict.keys())
                    rna_chain_groups = convert_to_chain_groups(rna_chain_ids)
                    chain_groups = dna_chain_groups + rna_chain_groups
                    structureArray = freesasa.structureArray(temp_pdb_file.name,
                                                             options={"separate-chains": False,
                                                                      "chain-groups": chain_groups})
                    for model_index, model in enumerate(structureArray):
                        if model_index == 1:  # Consider only the second model (model_index = 0 is the total area)
                            result = freesasa.calc(model)
                            sasa_value = result.totalArea()
                            complexed_rna_dna_interface_area = sasa_value
                            print(f"Uncomplexed complexed RNA-DNA Chain:{chain_groups}, SASA: {sasa_value:.2f} Å²")
                if (dna.is_dna and not dna.is_single_copy_dna()) and (rna.is_rna and rna.is_single_copy_rna()):
                    rna_chain_groups = ''.join(rna_chain_sasa_dict.keys())
                    dna_chain_groups = convert_to_chain_groups(dna_chain_ids)
                    chain_groups = dna_chain_groups + rna_chain_groups
                    structureArray = freesasa.structureArray(temp_pdb_file.name,
                                                             options={"separate-chains": False,
                                                                      "chain-groups": chain_groups})
                    for model_index, model in enumerate(structureArray):
                        if model_index == 1:  # Consider only the second model (model_index = 0 is the total area)
                            result = freesasa.calc(model)
                            sasa_value = result.totalArea()
                            complexed_rna_dna_interface_area = sasa_value
                            print(f"Uncomplexed complexed RNA-DNA Chain:{chain_groups}, SASA: {sasa_value:.2f} Å²")
                if (dna.is_dna and not dna.is_single_copy_dna()) and (rna.is_rna and not rna.is_single_copy_rna()):
                    rna_chain_groups = convert_to_chain_groups(rna_chain_ids)
                    dna_chain_groups = convert_to_chain_groups(dna_chain_ids)
                    chain_groups = dna_chain_groups + rna_chain_groups
                    structureArray = freesasa.structureArray(temp_pdb_file.name,
                                                             options={"separate-chains": False,
                                                                      "chain-groups": chain_groups})
                    for model_index, model in enumerate(structureArray):
                        if model_index == 1:  # Consider only the second model (model_index = 0 is the total area)
                            result = freesasa.calc(model)
                            sasa_value = result.totalArea()
                            complexed_rna_dna_interface_area = sasa_value
                            print(f"Uncomplexed complexed RNA-DNA Chain:{chain_groups}, SASA: {sasa_value:.2f} Å²")

                if (rna.is_rna and not rna.is_single_copy_rna()) and (protein.is_protein and not protein.is_single_copy_protein()):
                    interface_area = protein_interface_area + rna_interface_area - total_sasa_complex
                    protein_rna_interface_area = interface_area / 2
                    print(f"Interface Protein/RNA Area: {protein_rna_interface_area:.2f} Å²")

                if (dna.is_dna and not dna.is_single_copy_dna()) and (protein.is_protein and not protein.is_single_copy_protein()):
                    interface_area = dna_interface_area + protein_interface_area - total_sasa_complex
                    protein_dna_interface_area = interface_area / 2
                    print(f"Interface Protein/DNA Area: {protein_dna_interface_area:.2f} Å²")

                if (rna.is_rna and not rna.is_single_copy_rna()) and (dna.is_dna and not dna.is_single_copy_dna()):
                    interface_area = dna_interface_area + rna_interface_area - complexed_rna_dna_interface_area
                    rna_dna_interface_area = interface_area / 2
                    print(f"Interface RNA/DNA Area: {rna_dna_interface_area:.2f} Å²")

                if (rna.is_rna and rna.is_single_copy_rna()) and (dna.is_dna and not dna.is_single_copy_dna()):
                    interface_area = dna_interface_area + sum(rna_chain_sasa_dict.values()) - complexed_rna_dna_interface_area
                    rna_dna_interface_area = interface_area / 2
                    print(f"Interface RNA/DNA Area: {rna_dna_interface_area:.2f} Å²")

                if (rna.is_rna and not rna.is_single_copy_rna()) and (dna.is_dna and dna.is_single_copy_dna()):
                    interface_area = sum(dna_chain_sasa_dict.values()) + rna_interface_area - complexed_rna_dna_interface_area
                    rna_dna_interface_area = interface_area / 2
                    print("rna_dna_interface_area", rna_dna_interface_area)
                    print(f"Interface RNA/DNA Area: {rna_dna_interface_area:.2f} Å²")

                if (rna.is_rna and rna.is_single_copy_rna()) and (protein.is_protein and not protein.is_single_copy_protein()):
                    interface_area = sum(rna_chain_sasa_dict.values()) + protein_interface_area - total_sasa_complex
                    protein_rna_interface_area = interface_area/max(rna_chain_number, 2)
                    print(f"Interface Protein/RNA Area: {protein_rna_interface_area:.2f} Å²")

                if (rna.is_rna and not rna.is_single_copy_rna()) and (protein.is_protein and protein.is_single_copy_protein()):
                    interface_area = rna_interface_area + sum(protein_chain_sasa_dict.values()) - total_sasa_complex
                    protein_rna_interface_area = interface_area/max(protein_chain_number, 2)
                    print(f"Interface Protein/RNA Area: {protein_rna_interface_area:.2f} Å²")

                if (dna.is_dna and not dna.is_single_copy_dna()) and (protein.is_protein and protein.is_single_copy_protein()):
                    interface_area = dna_interface_area + sum(protein_chain_sasa_dict.values()) - total_sasa_complex
                    protein_dna_interface_area = interface_area / max(protein_chain_number, 2)
                    print(f"Interface Protein/DNA Area: {protein_dna_interface_area:.2f} Å²")

                if (dna.is_dna and dna.is_single_copy_dna()) and (protein.is_protein and not protein.is_single_copy_protein()):
                    interface_area = sum(dna_chain_sasa_dict.values()) + protein_interface_area - total_sasa_complex
                    protein_dna_interface_area = interface_area / max(dna_chain_number, 2)
                    print(f"Interface Protein/DNA Area: {protein_dna_interface_area:.2f} Å²")

                if (rna.is_rna and not rna.is_single_copy_rna()) and not protein.is_protein and not dna.is_dna:
                    interface_area = rna_interface_area - total_sasa_complex
                    rna_rna_interface_area = interface_area/2
                    print(f"Interface RNA/RNA Area: {rna_rna_interface_area:.2f} Å²")

            if total_sasa_complex > 0:
                if rna.is_rna and protein.is_protein:
                    protein_rna_ratio = (protein_rna_interface_area / total_sasa_complex) * 100
                    print(f"Protein-RNA Interface Ratio (R_i/s): {protein_rna_ratio:.2f}%")
                    aap_protein_rna = extract_surface_residues(file_name, file_path, protein_rna_interface_area, total_sasa_complex)
                    protein_rna_ratio = round(protein_rna_ratio, 2)
                    protein_rna_interface_area = round(protein_rna_interface_area, 2)
                    # return protein_rna_interface_area, total_sasa_complex, aap
                if dna.is_dna and protein.is_protein:
                    protein_dna_ratio = (protein_dna_interface_area / total_sasa_complex) * 100
                    print(f"Protein-DNA Interface Ratio (R_i/s): {protein_dna_ratio:.2f}%")
                    aap_protein_dna = extract_surface_residues(file_name, file_path, protein_dna_interface_area, total_sasa_complex)
                    protein_dna_ratio = round(protein_dna_ratio, 2)
                    protein_dna_interface_area = round(protein_dna_interface_area, 2)
                    # return protein_dna_interface_area, total_sasa_complex, aap
                if rna.is_rna and dna.is_dna:
                    rna_dna_ratio = (rna_dna_interface_area / total_sasa_complex) * 100
                    rna_dna_ratio = round(rna_dna_ratio, 2)
                    print(f"DNA-RNA Interface Ratio (R_i/s): {rna_dna_ratio:.2f}%")
                    rna_dna_interface_area = round(rna_dna_interface_area, 2)
                    # return rna_dna_interface_area, total_sasa_complex
                if rna.is_rna and not protein.is_protein and not dna.is_dna:
                    rna_rna_ratio = (rna_rna_interface_area / total_sasa_complex) * 100
                    rna_rna_ratio = round(rna_rna_ratio, 2)
                    rna_rna_interface_area = round(rna_rna_interface_area, 2)
                    print(f"RNA-RNA Interface Ratio (R_i/s): {rna_rna_ratio:.2f}%")
                    # return rna_rna_interface_area, total_sasa_complex
            else:
                ratio = 0
            print("-" * 50)
        else:
            if rna.is_rna and protein.is_protein:
                protein_rna_interface_area = None
                protein_rna_ratio = None
                aap_protein_rna = None
            if dna.is_dna and protein.is_protein:
                protein_dna_interface_area = None
                protein_dna_ratio = None
                aap_protein_dna = None
            if rna.is_rna and dna.is_dna:
                protein_rna_interface_area = None
                protein_dna_interface_area = None
                rna_dna_interface_area = None
                protein_rna_ratio = None
                protein_dna_ratio = None
                rna_dna_ratio = None
                aap_protein_rna = None
                aap_protein_dna = None
            if rna.is_rna and not protein.is_protein and not dna.is_dna:
                rna_rna_interface_area = None
                rna_rna_ratio = None

        # Clean up the temporary PDB file
        os.remove(temp_pdb_file.name)

    if aap_protein_rna is None:
        return

    pdb_id = os.path.basename(file_path).split('.')[0][:4].upper()
    file_name_cif = file_name[:-4] # without cif

    if aap_protein_rna:
        aap_protein_rna = ', '.join([f"{key}: {value}" for key, value in aap_protein_rna.items()])
    if aap_protein_dna:
        aap_protein_dna = ', '.join([f"{key}: {value}" for key, value in aap_protein_dna.items()])

    if len(file_name_cif) > 6:
        if protein.is_protein and rna.is_rna and not dna.is_dna:
            database_methods.update_or_insert('pred_protein_rna',
                                              ['ProteinRNAInterfaceArea', 'ProteinRNAInterfaceRatio', 'AAPproteinRNA'],
                                              (protein_rna_interface_area, protein_rna_ratio, aap_protein_rna),
                                              condition=f"FileName = '{file_name}'")

        elif not protein.is_protein and rna.is_rna and not dna.is_dna:
            database_methods.update_or_insert('pred_rna_rna',
                                              ['RNArnaInterfaceArea', 'RNArnaInterfaceRatio'],
                                              (rna_rna_interface_area, rna_rna_ratio),
                                              condition=f"FileName = '{file_name}'")

        elif protein.is_protein and not rna.is_rna and dna.is_dna:
            database_methods.update_or_insert('pred_protein_dna',
                                              ['ProteinDNAInterfaceArea', 'ProteinDNAInterfaceRatio', 'AAPproteinDNA'],
                                              (protein_dna_interface_area, protein_dna_ratio, aap_protein_dna),
                                              condition=f"FileName = '{file_name}'")

        elif protein.is_protein and rna.is_rna and dna.is_dna:
            database_methods.update_or_insert('pred_protein_rna_dna',
                                              ['ProteinRNAInterfaceArea', 'ProteinDNAInterfaceArea', 'RNAdnaInterfaceArea'
                                               'ProteinRNAInterfaceRatio', 'ProteinDNAInterfaceRatio',
                                               'RNAdnaInterfaceRatio', 'AAPproteinRNA', 'AAPproteinDNA'],
                                              (protein_rna_interface_area, protein_dna_interface_area,
                                               rna_dna_interface_area, protein_rna_ratio, protein_dna_ratio,
                                               rna_dna_ratio, aap_protein_rna, aap_protein_dna),
                                              condition=f"FileName = '{file_name}'")

    else:
        if protein.is_protein and rna.is_rna and not dna.is_dna:
            database_methods.update_or_insert('exp_protein_rna',
                                              ['ProteinRNAInterfaceArea', 'ProteinRNAInterfaceRatio', 'AAPproteinRNA'],
                                              (protein_rna_interface_area, protein_rna_ratio, aap_protein_rna),
                                              condition=f"PDBId = '{pdb_id}'")

        elif not protein.is_protein and rna.is_rna and not dna.is_dna:
            database_methods.update_or_insert('exp_rna_rna',
                                              ['RNArnaInterfaceArea', 'RNArnaInterfaceRatio'],
                                              (rna_rna_interface_area, rna_rna_ratio),
                                              condition=f"PDBId = '{pdb_id}'")

        elif protein.is_protein and not rna.is_rna and dna.is_dna:
            database_methods.update_or_insert('exp_protein_dna',
                                              ['ProteinDNAInterfaceArea', 'ProteinDNAInterfaceRatio', 'AAPproteinDNA'],
                                              (protein_dna_interface_area, protein_dna_ratio, aap_protein_dna),
                                              condition=f"PDBId = '{pdb_id}'")

        elif protein.is_protein and rna.is_rna and dna.is_dna:
            database_methods.update_or_insert('exp_protein_rna_dna',
                                              ['ProteinRNAInterfaceArea', 'ProteinDNAInterfaceArea', 'RNAdnaInterfaceArea'
                                               'ProteinRNAInterfaceRatio', 'ProteinDNAInterfaceRatio',
                                               'RNAdnaInterfaceRatio', 'AAPproteinRNA', 'AAPproteinDNA'],
                                              (protein_rna_interface_area, protein_dna_interface_area,
                                               rna_dna_interface_area, protein_rna_ratio, protein_dna_ratio,
                                               rna_dna_ratio, aap_protein_rna, aap_protein_dna),
                                              condition=f"PDBId = '{pdb_id}'")

def extract_surface_residues(file_name, file_path, interface_area, total_sasa_complex, sasa_cutoff=50):
    cmd.reinitialize()
    try:
        cmd.load(file_path, object='protein_structure')
    except Exception as e:
        print(f"Error loading mmCIF file: {e}")
        sys.exit(1)

    # Remove all solvents and other heteroatoms to focus on protein residues
    cmd.remove('solvent')
    cmd.remove('not polymer.protein')
    cmd.set('dot_solvent', 1)
    cmd.set('dot_density', 4)
    cmd.get_area(selection='all', load_b=1)
    model = cmd.get_model('all')
    surface_residue_sasa = {}

    for atom in model.atom:
        key = (atom.chain, atom.resn, atom.resi)
        sasa = atom.b
        if key in surface_residue_sasa:
            surface_residue_sasa[key] += sasa
        else:
            surface_residue_sasa[key] = sasa

    sasa_sums = {}
    for (chain, resn, resi), sasa in surface_residue_sasa.items():
        if sasa > sasa_cutoff:
            if resn in sasa_sums:
                sasa_sums[resn] += sasa # Sum up the SASA values for each amino acid (resn)
            else:
                sasa_sums[resn] = sasa
    sasa_sums = {resn: round(sasa, 2) for resn, sasa in sasa_sums.items()}
    # print(sasa_sums)

    file_name_wo_cif = os.path.splitext(file_name)[0]
    if len(file_name_wo_cif) in [4, 6]:
        protein = Protein.get_protein_from_db(id=file_name_wo_cif.upper())
    elif len(file_name_wo_cif) > 6:
        parts = file_name_wo_cif.split('_')
        protein = Protein.get_protein_from_db(id=parts[1].upper(), file_name=file_name)

    aa_motif = protein.get_aa_motif()
    if not aa_motif:
        print(f"Error: Add {file_name_wo_cif} to the database.")
        pass
    else:
        residues = aa_motif.split(", ")
        resn_list = []
        resi_list = []

        for residue in residues:
            try:
                resn, resi = residue.split("(")
                resn_list.append(resn)
                resi_list.append(resi.rstrip(")"))
            except ValueError:
                # If split fails, log the PDB and skip this residue
                pdb_id = os.path.basename(file_path).split('.')[0][:4].upper()
                if pdb_id not in skipped_pdbs:
                    skipped_pdbs.append(pdb_id)
                print(f"Warning: Malformed residue format in {pdb_id}: {residue}")
                continue

        # aa_chain = protein.get_protein_chain_IDs()
        interface_sasa = {}
        for atom in model.atom:
            chain = atom.chain
            atom_resn = atom.resn
            atom_resi = str(atom.resi)  # Convert to string to match with resi_list

            # Check if the atom's resn and resi match to any in the aa_motif
            if atom_resn in resn_list and atom_resi in resi_list:
                key = (chain, atom_resn, atom_resi)
                sasa = atom.b  # atom.b stores the SASA value

                if key in interface_sasa:
                    interface_sasa[key] += sasa
                else:
                    interface_sasa[key] = sasa

        interface_sums = {}
        for (chain, resn, resi), sasa in interface_sasa.items():
            if resn in interface_sums:
                interface_sums[resn] += sasa
            else:
                interface_sums[resn] = sasa
        interface_sums = {resn: round(sasa, 2) for resn, sasa in interface_sums.items()}
        # print(interface_sums)

        # aap = (sum(interface_sums) / interface_area) / (sum(sasa_sums) / total_sasa_complex)
        aap_values = {}
        for resn in interface_sums.keys():
            # Ratio of amino acid SASAs in the interface to total interface area
            numerator = interface_sums[resn] / interface_area
            # Ratio of amino acid SASAs on the surface to total surface SASA
            if resn in sasa_sums:  # Ensure the residue is in both interface and surface lists
                denominator = sasa_sums[resn] / total_sasa_complex
            else:
                denominator = 0
            if denominator > 0:
                aap = numerator / denominator
                aap_values[resn] = aap
            else:
                aap_values[resn] = float('inf')
        if aap_values == {}: # No interactions
            sorted_aap_values = None
        else:
            aap_values = {resn: round(sasa, 2) for resn, sasa in aap_values.items()}
            sorted_aap_values = dict(sorted(aap_values.items()))

        # print("sasa_sums", sasa_sums)
        print("aap", sorted_aap_values)

        return sorted_aap_values

if __name__ == "__main__":
    # Initialize global list for skipped PDBs
    skipped_pdbs = []

    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <directory_path>")
        sys.exit(1)
    cif_directory = sys.argv[1]
    if not os.path.isdir(cif_directory):
        print(f"Error: {cif_directory} is not a valid directory.")
        sys.exit(1)
    cif_parser = MMCIFParser()

    for file_name in os.listdir(cif_directory):
        if file_name.endswith(".cif"):
            file_path = os.path.join(cif_directory, file_name)
            # extract_interface_area(file_name, file_path)
            try:
                extract_interface_area(file_name, file_path)
            except Exception as e:
                print(f"Error processing {file_name}: {str(e)}")
                continue

    # Print skipped PDBs at the end
    if skipped_pdbs:
        print("\nThe following PDB IDs were skipped due to errors:")
        print(", ".join(skipped_pdbs))
        print(f"Total skipped: {len(skipped_pdbs)}")
    else:
        print("\nNo PDB files were skipped.")

    database_methods.close_connection()