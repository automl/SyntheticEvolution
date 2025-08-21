
import os
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import sqlite3
import json
import re
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from utils.annotate_rna_family import annotate_families_from_sequence
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqUtils import seq1
from collections import defaultdict

unmodified_residues = {'C', 'A', 'G', 'U', 'T', 'DG', 'DT', 'DC', 'DA', 'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'CYS', 'MET',
                       'PHE', 'TRP', 'PRO', 'SER', 'THR', 'TYR', 'ASN', 'GLN', 'ASP', 'GLU', 'HIS', 'LYS', 'ARG'}
modified_residues = {'6MZ', 'PSU', '5MC', 'OMC', '4OC', '5MU', 'OMU', 'UR3', 'A2M', 'MA6', '6MZ', '2MG', 'OMG', '7MG',
                     'RSQ', '5CM', 'C34', '5HC', '6OG', '6MA', '1CC', '8OG', '5FC', '3DR', 'SEP', 'TPO', 'PTR',
                     'NEP', 'HIP', 'ALY', 'MLY', 'M3L', 'MLZ', '2MR', 'AGM', 'MCS', 'HYP', 'HY3', 'LYZ', 'AHB',
                     'P1L', 'SNN', 'SNC', 'TRF', 'KCR', 'CIR', 'YHA'}

def calculate_length(sequence):
    mod_pattern = re.compile(r'\(.*?\)') # Regex pattern to match modifications in parentheses
    cleaned_sequence = mod_pattern.sub('X', sequence) # Remove all parentheses and count them as 1 character 'X'
    return len(cleaned_sequence.replace(' ', '').replace(';', '').replace('\n', ''))

def parse_cif_sequence(cif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_file)
    model = structure[0]
    sequences = defaultdict(list)
    pdb_info = MMCIF2Dict(cif_file)
    # entity_poly_seq = pdb_info.get('_entity_poly.pdbx_seq_one_letter_code', [])
    poly_types = pdb_info.get('_entity_poly.type', [])
    strand_ids = pdb_info.get('_entity_poly.pdbx_strand_id', [])
    protein_sequences = []
    rna_sequences = []
    protein_chain_ids = []
    rna_chain_ids = []
    protein_seq_lengths = []
    rna_seq_lengths = []
    number_proteins = 0
    number_RNAs = 0
    protein_modified_aas = []
    rna_modified_nucleotides = []

    for poly_type, strand_id in zip(poly_types, strand_ids):
        chain_ids = strand_id.split(',') # Split strand IDs in case there are multiple strands for an entity
        for chain in model:
            chain_id = chain.get_id()
            if chain_id not in chain_ids: # Only process chains listed in the strand_ids
                continue
            for residue in chain: # Process each residue in the chain
                resname = residue.get_resname()

                if poly_type == 'polypeptide(L)':
                    if chain_id not in protein_chain_ids:
                        protein_chain_ids.append(chain_id)
                    if resname in modified_residues:
                        sequences[chain_id].append(f"({resname})")
                        if resname not in protein_modified_aas:
                            protein_modified_aas.append(resname)
                    elif resname in unmodified_residues:
                        sequences[chain_id].append(seq1(resname))

                elif poly_type == 'polyribonucleotide':
                    if chain_id not in rna_chain_ids:
                        rna_chain_ids.append(chain_id)
                    if resname in modified_residues:
                        sequences[chain_id].append(f"({resname})")
                        if resname not in rna_modified_nucleotides:
                            rna_modified_nucleotides.append(resname)
                    elif resname in unmodified_residues:
                        sequences[chain_id].append(resname)

    for chain_id, seq_list in sequences.items():
        sequence = ''.join(seq_list)
        if any(poly_type == 'polypeptide(L)' for poly_type, chain_ids in zip(poly_types, strand_ids) if
               chain_id in chain_ids.split(',')):
            protein_sequences.append(sequence)
            length = calculate_length(sequence)
            protein_seq_lengths.append(length)
            number_proteins += 1
        elif any(poly_type == 'polyribonucleotide' for poly_type, chain_ids in zip(poly_types, strand_ids) if
                 chain_id in chain_ids.split(',')):
            rna_sequences.append(sequence)
            length = calculate_length(sequence)
            rna_seq_lengths.append(length)
            number_RNAs += 1

    return protein_sequences, rna_sequences, protein_chain_ids, rna_chain_ids, \
        protein_seq_lengths, rna_seq_lengths, number_proteins, number_RNAs, \
        protein_modified_aas, rna_modified_nucleotides

# Extract data from the PDB CIF file only
def parse_cif_file(cif_file):
    pdb_info = MMCIF2Dict(cif_file)

    pdb_id = pdb_info.get('_entry.id', [''])[0]
    print(f"Processing {pdb_id}...")
    deposited_date = pdb_info.get('_pdbx_audit_revision_history.revision_date', [''])[0]
    experiment = pdb_info.get('_exptl.method', [''])[0]
    xray_resolution = pdb_info.get('_refine.ls_d_res_high', [''])[0]

    is_protein = False
    is_rna = False
    protein_names = []
    uniprot_ids = []
    protein_descriptions = []
    protein_sources = []
    protein_organisms = []
    rna_descriptions = []
    rna_sources = []
    rna_families = []

    protein_sequences, rna_sequences, protein_chain_ids, rna_chain_ids, \
    protein_lengths, rna_lengths,number_proteins, number_RNAs, \
    protein_modified_aas, rna_modified_nucleotides= parse_cif_sequence(cif_file)

    entity_info_map = {}
    for i in range(len(pdb_info.get('_entity.id', []))):
        entity_id = pdb_info['_entity.id'][i]
        entity_info_map[entity_id] = {
            'type': pdb_info['_entity.type'][i],
            'src_method': pdb_info['_entity.src_method'][i],
            'description': pdb_info['_entity.pdbx_description'][i],
        }
    struct_ref_entity_ids = pdb_info.get('_struct_ref.entity_id', [])
    for i in range(len(struct_ref_entity_ids)):
        entity_id = struct_ref_entity_ids[i]
        if entity_id in entity_info_map:
            entity_info_map[entity_id]['id_name'] = pdb_info['_struct_ref.db_code'][i]
            entity_info_map[entity_id]['uniprot_id'] = pdb_info['_struct_ref.pdbx_db_accession'][i]

    # Map entity ID to info from _entity_src_gen
    entity_src_gen_map = {}
    for i in range(len(pdb_info.get('_entity_src_gen.entity_id', []))):
        entity_id = pdb_info['_entity_src_gen.entity_id'][i]
        entity_src_gen_map[entity_id] = {
            'organism_scientific': pdb_info['_entity_src_gen.pdbx_gene_src_scientific_name'][i],
        }

    entity_poly_type = pdb_info.get('_entity_poly.type', [])
    entity_poly_seq = pdb_info.get('_entity_poly.pdbx_seq_one_letter_code', [])
    entity_ids = pdb_info.get('_entity_poly.entity_id', [])
    strand_ids = pdb_info.get('_entity_poly.pdbx_strand_id', [])
    for entity_id, entity_type, seq, strand_id in zip(entity_ids, entity_poly_type, entity_poly_seq,
                                                                           strand_ids):
        seq = seq.replace(';', '').replace('\n', '')
        rna_sequence = ""
        if entity_type == 'polypeptide(L)':
            is_protein = True
            if entity_id in entity_info_map:
                protein_description = entity_info_map[entity_id]['description']
                protein_source = entity_info_map[entity_id]['src_method']
                uniprot_id = entity_info_map[entity_id]['uniprot_id']
                protein_name = entity_info_map[entity_id]['id_name']
                organism = entity_src_gen_map.get(entity_id, {}).get('organism_scientific')
                protein_names.append(protein_name)
                uniprot_ids.append(uniprot_id)
                protein_descriptions.append(protein_description)
                protein_sources.append(protein_source)
                protein_organisms.append(organism)

        elif entity_type == 'polyribonucleotide':
            rna_sequence += seq
            if '(D' in rna_sequence and ')' in rna_sequence:
                print("ERROR: RNA Sequence such as 'AUCCAGGUGCAC(DA)(DA)(DA)(DG)(DA)(DA)' not supported. File: ",
                      cif_file)
                break
            is_rna = True
            if entity_id in entity_info_map:
                # rna_name = entity_info_map[entity_id]['id_name']
                rna_description = entity_info_map[entity_id]['description']
                rna_source = entity_info_map[entity_id]['src_method']
                rna_descriptions.append(rna_description)
                rna_sources.append(rna_source)
                project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                rfam_path = os.path.join(project_root, "evaluation/data/rfam")
                rna_family = annotate_families_from_sequence(rna_sequence, rfam_path)
                if rna_family is None:
                    rna_family = "None"
                else:
                    rna_family = rna_family['families'][0]
                rna_families.append(rna_family)
                print("rna_families", rna_families)

    if is_rna and is_protein:
        classification = "protein_rna"
        return classification, (pdb_id, uniprot_ids, protein_names, protein_descriptions, protein_sources,
                                protein_modified_aas, protein_organisms,protein_sequences, protein_chain_ids,
                                protein_lengths, number_proteins, rna_descriptions, rna_sources, rna_families,
                                rna_modified_nucleotides, rna_sequences, rna_chain_ids, rna_lengths, number_RNAs,
                                deposited_date, experiment, xray_resolution)

    elif is_rna and not is_protein:
        classification = "rna_rna"
        return classification, (pdb_id, rna_descriptions, rna_sources, rna_families, rna_modified_nucleotides,
            rna_sequences, rna_chain_ids, rna_lengths, number_RNAs, deposited_date, experiment, xray_resolution)

def insert_data_to_db(db_path, classification, data):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    def handle_list_or_value(item):
        if isinstance(item, list):
            if len(item) == 1:
                return item[0]
            return json.dumps(item)
        return item

    record_values = [handle_list_or_value(data[i]) for i in range(len(data))]

    if classification == "protein_rna":
        cursor.execute('''CREATE TABLE IF NOT EXISTS exp_protein_rna (
                            Id INTEGER PRIMARY KEY AUTOINCREMENT,
                            PDBId TEXT,
                            UniprotId TEXT,
                            ProteinName TEXT,
                            ProteinDescription TEXT,
                            ProteinSource TEXT,
                            ProteinModified TEXT,
                            ProteinOrganism TEXT,
                            ProteinSequence TEXT,
                            ProteinChainIDs TEXT,
                            ProteinLength INTEGER,
                            NumberProteins INTEGER,
                            RNADescription TEXT,
                            RNASource TEXT,
                            RNAFamily TEXT,
                            RNAModified TEXT,
                            RNASequence TEXT,
                            RNAChainIDs TEXT,
                            RNALength INTEGER,
                            NumberRNAs INTEGER,
                            Deposited_date TEXT,
                            Experiment TEXT,
                            XRayResolution FLOAT
                          )''')

        # Check if PDBId already exists in the database
        cursor.execute("SELECT 1 FROM exp_protein_rna WHERE PDBId = ?", (data[0],))
        if cursor.fetchone() is None:
            # Insert new record
            cursor.execute('''INSERT INTO exp_protein_rna (PDBId, UniprotId, ProteinName, ProteinDescription, ProteinSource, 
                              ProteinModified, ProteinOrganism, ProteinSequence, ProteinChainIDs, ProteinLength, NumberProteins, 
                              RNADescription, RNASource, RNAFamily, RNAModified, RNASequence, RNAChainIDs, RNALength, NumberRNAs, 
                              Deposited_date, Experiment, XRayResolution)
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', record_values)
        else:
            # Update existing record
            cursor.execute('''UPDATE exp_protein_rna SET 
                              UniprotId=?, ProteinName=?, ProteinDescription=?, ProteinSource=?, 
                              ProteinModified=?, ProteinOrganism=?, ProteinSequence=?, ProteinChainIDs=?, ProteinLength=?, 
                              NumberProteins=?,
                              RNADescription=?, RNASource=?, RNAFamily=?, RNAModified=?, RNASequence=?, RNAChainIDs=?, RNALength=?, NumberRNAs=?, 
                              Deposited_date=?, Experiment=?, XRayResolution=?
                              WHERE PDBId=?''',
                              record_values[1:] + [record_values[0]])  # Move PDBId to the end for WHERE clause

    if classification == "rna_rna":
        cursor.execute('''CREATE TABLE IF NOT EXISTS exp_rna_rna (
                            Id INTEGER PRIMARY KEY AUTOINCREMENT,
                            PDBId TEXT,
                            RNADescription TEXT,
                            RNASource TEXT,
                            RNAFamily TEXT,
                            RNAModified TEXT,
                            RNASequence TEXT,
                            RNAChainIDs TEXT,
                            RNALength INTEGER,
                            NumberRNAs INTEGER,
                            Deposited_date TEXT,
                            Experiment TEXT,
                            XRayResolution FLOAT
                          )''')

        # Check if PDBId already exists in the database
        cursor.execute("SELECT 1 FROM exp_rna_rna WHERE PDBId = ?", (data[0],))
        if cursor.fetchone() is None:
            # Insert new record
            cursor.execute('''INSERT INTO exp_rna_rna (PDBId, 
                              RNADescription, RNASource, RNAFamily, RNAModified, RNASequence, RNAChainIDs, RNALength, NumberRNAs, 
                              Deposited_date, Experiment, XRayResolution)
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', record_values)
        else:
            # Update existing record
            cursor.execute('''UPDATE exp_rna_rna SET 
                              RNADescription=?, RNASource=?, RNAFamily=?, RNAModified=?, RNASequence=?, RNAChainIDs=?, RNALength=?, NumberRNAs=?, 
                              Deposited_date=?, Experiment=?, XRayResolution=?
                              WHERE PDBId=?''',
                              record_values[1:] + [record_values[0]])  # Move PDBId to the end for WHERE clause

    conn.commit()
    conn.close()

def process_cif_files(directory_path, db_path):
    check_folder = os.path.join(os.path.dirname(directory_path), 'check_pdb')
    if not os.path.exists(check_folder):
        os.makedirs(check_folder)

    cif_files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.cif')]
    for cif_file in cif_files:
        try:
            classification, data = parse_cif_file(cif_file)
            if classification and data:  # Only insert if we got valid data
                insert_data_to_db(db_path, classification, data)
            else:
                # Move problematic file to check folder
                file_name = os.path.basename(cif_file)
                dst_path = os.path.join(check_folder, file_name)
                print(f"Moving {file_name} to {check_folder}")
                import shutil
                shutil.move(cif_file, dst_path)
        except Exception as e:
            # Move file to check folder on any error
            file_name = os.path.basename(cif_file)
            dst_path = os.path.join(check_folder, file_name)
            print(f"Error processing {file_name}: {e}")
            print(f"Moving {file_name} to {check_folder}")
            import shutil
            shutil.move(cif_file, dst_path)

# Example usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <directory_path>")
    else:
        directory_path = sys.argv[1]
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # One level above utils package.
        db_path = os.path.join(project_root, "database/rbpDatabase.db")
        process_cif_files(directory_path, db_path)