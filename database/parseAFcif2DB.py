import os
import re
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import sqlite3
import json
import shutil

unmodified_residues = {'C', 'A', 'G', 'U', 'T', 'DG', 'DT', 'DC', 'DA', 'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'CYS', 'MET',
                       'PHE', 'TRP', 'PRO', 'SER', 'THR', 'TYR', 'ASN', 'GLN', 'ASP', 'GLU', 'HIS', 'LYS', 'ARG'}
modified_residues = {'6MZ', 'PSU', '5MC', 'OMC', '4OC', '5MU', 'OMU', 'UR3', 'A2M', 'MA6', '6MZ', '2MG', 'OMG', '7MG',
                     'RSQ', '5CM', 'C34', '5HC', '6OG', '6MA', '1CC', '8OG', '5FC', '3DR', 'SEP', 'TPO', 'PTR',
                     'NEP', 'HIP', 'ALY', 'MLY', 'M3L', 'MLZ', '2MR', 'AGM', 'MCS', 'HYP', 'HY3', 'LYZ', 'AHB',
                     'P1L', 'SNN', 'SNC', 'TRF', 'KCR', 'CIR', 'YHA'}

def calculate_length(sequence):
    mod_pattern = re.compile(r'\(.*?\)') # Regex pattern to match modifications in parentheses
    cleaned_sequence = mod_pattern.sub('X', sequence) # Remove all parentheses and count them as 1 character 'X'
    # print("cleaned_sequence", cleaned_sequence)
    return len(cleaned_sequence.replace(' ', '').replace(';', '').replace('\n', ''))

def get_strand_type_from_id(strand_id, mapping):
    # Find strand_type for given strand_id
    for (strand_type, id_key), entity_type in mapping.items():
        if id_key == strand_id:
            return id_key
    return None

def parse_cif_file(file_path): # AF3's cif file does not contain '_entity_poly.pdbx_seq_one_letter_code'
    # Parse the cif file using MMCIF2Dict to get metadata
    af_info = MMCIF2Dict(file_path)
    strand_ids = af_info.get('_entity_poly.pdbx_strand_id', []) # Chains A, B...
    entity_types = af_info.get('_entity_poly.type', [])
    # strand_types = af_info.get('_entity_poly.entity_id', []) # Chains 1, 2...

    # Mapping of strand IDs to entity types
    strand_id_to_type = {strand_id: entity_type for strand_id, entity_type in zip(strand_ids, entity_types)}
    # print("strand_id_to_type", strand_id_to_type)

    # {('1', 'A'): 'polyribonucleotide', ('2', 'B'): 'polyribonucleotide'} = (strand_types, strand_ids: entity_types)
    # strand_id_to_entity_type = {
    #     (strand_type, strand_id): entity_type
    #     for strand_type, strand_id, entity_type in zip(strand_types, strand_ids, entity_types)
    # }
    # print(strand_id_to_entity_type)

    protein_lengths = []
    protein_sequences = []
    rna_lengths = []
    rna_sequences = []
    dna_lengths = []
    dna_sequences = []
    protein_chain_ids = []
    rna_chain_ids = []
    dna_chain_ids = []
    number_proteins = 0
    number_RNAs = 0
    number_DNAs = 0

    # Parse the structure using Bio.PDB
    parser = MMCIFParser()
    structure = parser.get_structure('structure', file_path)

    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            if chain_id in strand_id_to_type:
                entity_type = strand_id_to_type[chain_id]
                # print("chain_id, entity_type", chain_id, entity_type)
                protein_sequence = ""
                rna_sequence = ""
                dna_sequence = ""
                for residue in chain:
                    hetero_flag, res_seq_nr, insertion_code = residue.get_id()
                    resname = residue.resname
                    if entity_type == 'polypeptide(L)':
                        if hetero_flag == ' ':
                            protein_sequence += resname
                        else:
                            protein_sequence += f"({resname})"
                    elif entity_type == 'polyribonucleotide':
                        if hetero_flag == ' ':
                            rna_sequence += resname  # Regular residue
                        else:
                            rna_sequence += f"({resname})"  # Modified residue
                    elif entity_type == 'polydeoxyribonucleotide':
                        if hetero_flag == ' ':
                            dna_sequence += resname
                        else:
                            dna_sequence += f"({resname})"
                if entity_type == 'polypeptide(L)':
                    number_proteins += 1
                    protein_sequence = aa_transformation_three2one_letter(protein_sequence)
                    protein_sequences.append(protein_sequence)
                    protein_chain_ids.append(chain_id)
                    protein_length = calculate_length(protein_sequence)
                    protein_lengths.append(protein_length)
                elif entity_type == 'polyribonucleotide':
                    number_RNAs += 1
                    rna_sequence = rna_sequence.replace(' ', '')
                    rna_sequences.append(rna_sequence)
                    rna_chain_ids.append(chain_id)
                    print("rna_sequence", rna_sequence)
                    rna_length = calculate_length(rna_sequence)
                    print("rna_length", rna_length)
                    rna_lengths.append(rna_length)
                elif entity_type == 'polydeoxyribonucleotide':
                    number_DNAs += 1
                    dna_sequence = dna_sequence.replace(' ', '')
                    dna_sequences.append(dna_sequence)
                    dna_chain_ids.append(chain_id)
                    dna_length = calculate_length(dna_sequence)
                    dna_lengths.append(dna_length)

    # print("protein_sequences", protein_sequences)
    # print("rna_sequences", rna_sequences)
    # print("dna_sequences", dna_sequences)
    file_name = os.path.basename(file_path)
    parts = file_name.split('_')
    if len(parts) >= 5:
        exp_db_id = parts[1].upper()
        match = re.search(r'_s(\d+)_model_(\d+)', file_name)
        if match:
            seed = int(match.group(1))
            model = int(match.group(2))

        # print(f"exp_db_id Name: {exp_db_id}")
        # print(f"seed: {seed}")
        # print(f"model: {model}")
        # print(f"Protein Sequence: {protein_sequence}")
        # print(f"RNA Sequence: {rna_sequence}")
        # print(f"Protein Length: {protein_length}")
        # print(f"RNA Length: {rna_length}")

    if rna_sequences and protein_sequences and not dna_sequences:
        classification = "protein_rna"
        return classification, (exp_db_id, file_name, protein_sequences, protein_chain_ids, protein_lengths, number_proteins,
                                "", "", "", rna_sequences, rna_chain_ids, rna_lengths, number_RNAs, "", "", "", "", "",
                                "", "", "", "", "", "", seed, model, "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
                                "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
                                "", "", "", "", "")

    elif rna_sequences and protein_sequences and dna_sequences:
        classification = "protein_rna_dna"
        return classification, (exp_db_id, file_name, protein_sequences, protein_chain_ids, protein_lengths, number_proteins,
                                "", "", "", "", rna_sequences, rna_chain_ids, rna_lengths, number_RNAs, "", "", "",
                                dna_sequences, dna_chain_ids, dna_lengths, number_DNAs, "", "", "", "", "", "", "", "", "",
                                "", "", "", "", "", "",  seed, model, "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
                                "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
                                "", "", "", "", "", "", "")

    elif not rna_sequences and protein_sequences and dna_sequences:
        classification = "protein_dna"
        return classification, (exp_db_id, file_name, protein_sequences, protein_chain_ids, protein_lengths, number_proteins, "",
                                "", "", dna_sequences, dna_chain_ids, dna_lengths, number_DNAs, "", "", "", "", "",
                                "", "", "", "", "", seed, model, "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
                                "", "", "", "", "", "")

    elif rna_sequences and not protein_sequences and not dna_sequences:
        classification = "rna_rna"
        return classification, (exp_db_id, file_name, rna_sequences, rna_chain_ids, rna_lengths, number_RNAs, "", "", "", "",
                                "", "", "", "", "", "", "", seed, model, "", "", "", "", "", "", "", "", "", "", "",
                                "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")

    else:
        print(f"Error: The file name {file_name} format does not match the expected pattern for seed and model. Correct file name: fold_id_s1_model_0.cif")
        return None

def aa_transformation_three2one_letter(three_letter_sequence):
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    one_letter_sequence = ''
    i = 0
    while i < len(three_letter_sequence):
        # Check if we're at a modified residue (in parentheses)
        if three_letter_sequence[i] == '(':
            # Find the closing parenthesis
            end = three_letter_sequence.find(')', i)
            if end != -1:
                # Keep the entire modified residue as is
                one_letter_sequence += three_letter_sequence[i:end+1]
                i = end + 1
                continue
        
        # Handle regular three-letter codes
        if i + 2 < len(three_letter_sequence):
            three_letter_code = three_letter_sequence[i:i+3]
            if three_letter_code in three_to_one:
                one_letter_sequence += three_to_one[three_letter_code]
            else:
                one_letter_sequence += 'X'
            i += 3
        else:
            # Handle any remaining characters
            one_letter_sequence += three_letter_sequence[i:]
            break

    return one_letter_sequence

def insert_pred_data(classification, data):
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # One level above utils package.
    db_path = os.path.join(project_root, "database/rbpDatabase.db")
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()

    def handle_list_or_value(item):
        if isinstance(item, list):
            if len(item) == 1:
                return item[0]
            return json.dumps(item)
        return item

    record_values = [handle_list_or_value(data[i]) for i in range(len(data))]

    if classification == "protein_rna":
        cursor.execute('''CREATE TABLE IF NOT EXISTS pred_protein_rna (
                            Id INTEGER PRIMARY KEY AUTOINCREMENT,
                            exp_db_id TEXT,
                            FileName TEXT,
                            ProteinSequence TEXT,
                            ProteinChainIDs TEXT,
                            ProteinLength INTEGER,
                            NumberProteins INTEGER,
                            AAMotif TEXT,
                            AAC TEXT,
                            AAPproteinRNA TEXT,
                            RNASequence TEXT,
                            RNAChainIDs TEXT,
                            RNALength INTEGER,
                            NumberRNAs INTEGER,
                            RNAMotif TEXT,
                            RNAMotifLength TEXT,
                            RNAModified TEXT,
                            ContactList TEXT,
                            ChainIDpairList_proteinRNA TEXT,
                            Hbond_proteinRNA TEXT,
                            vdWbond_proteinRNA TEXT,
                            ProteinRNAInterfaceArea FLOAT,
                            ProteinRNAInterfaceRatio FLOAT,
                            Free_energy REAL,
                            Binding_affinity_kd TEXT,
                            Seed INTEGER,
                            Model INTEGER,
                            af3_protein_pTM TEXT,
                            af3_rna_pTM TEXT,
                            af3_protein_ipTM TEXT,
                            af3_rna_ipTM TEXT,
                            af3_protein_pLDDT_avg FLOAT,
                            af3_rna_pLDDT_avg FLOAT,
                            af3_global_pae_avg FLOAT,
                            af3_chain_pair_pae_min TEXT,
                            af3_fraction_disordered FLOAT,
                            af3_has_clash FLOAT,
                            af3_ranking_score FLOAT,
                            af3_protein_MSA INTEGER,
                            af3_rna_MSA INTEGER,
                            rna_prot_interface_atom_ids TEXT, 
                            interface_rna_atom_ids TEXT,
                            Complex_RMSD TEXT,
                            Protein_RMSD TEXT, 
                            RNA_RMSD TEXT,
                            Protein_LDDT FLOAT,
                            RNA_LDDT FLOAT, 
                            Protein_TM FLOAT,
                            RNA_TM FLOAT, 
                            Complex_TM FLOAT, 
                            Complex_LDDT FLOAT,  
                            RNA_INF TEXT, 
                            RNA_DI FLOAT,
                            RNA_f1_score FLOAT,
                            RNA_precision FLOAT,
                            RNA_recall FLOAT,
                            RNA_mcc FLOAT,
                            RNA_wl FLOAT,
                            RNA_GlycosidicBond TEXT, 
                            RNA_SugarPucker TEXT, 
                            RNA_GammaAngle TEXT,
                            RNAmotif_score FLOAT,
                            RNA_Stems TEXT,
                            RNA_HairpinLoops TEXT,
                            RNA_InternalLoops TEXT,
                            RNA_MultibranchLoops TEXT,
                            RNA_DanglingEnds TEXT,
                            RNA_Pseudoknots TEXT,
                            RNA_isDuplex INTEGER,
                            RNAMotif_isDuplex INTEGER,
                            FOREIGN KEY (exp_db_id) REFERENCES exp_protein_rna(UniprotId),
                            FOREIGN KEY (exp_db_id) REFERENCES exp_protein_rna(PDBId)
                          )''')

        # Check if exp_db_id already exists in the database
        cursor.execute("SELECT 1 FROM pred_protein_rna WHERE exp_db_id = ?", (data[0],))
        if cursor.fetchone() is None:
            # Insert new record
            cursor.execute('''INSERT INTO pred_protein_rna (exp_db_id, FileName, ProteinSequence, ProteinChainIDs, 
                              ProteinLength, NumberProteins, AAMotif, AAC, AAPproteinRNA, RNASequence, RNAChainIDs, RNALength, NumberRNAs, 
                              RNAMotif, RNAMotifLength, RNAModified, ContactList, ChainIDpairList_proteinRNA, 
                              Hbond_proteinRNA, vdWbond_proteinRNA, ProteinRNAInterfaceArea, 
                              ProteinRNAInterfaceRatio, Free_energy, Binding_affinity_kd, Seed, Model, 
                              af3_protein_pTM, af3_rna_pTM, af3_protein_ipTM, af3_rna_ipTM, af3_protein_pLDDT_avg,
                              af3_rna_pLDDT_avg, af3_global_pae_avg, af3_chain_pair_pae_min, af3_fraction_disordered,
                              af3_has_clash, af3_ranking_score, af3_protein_MSA, af3_rna_MSA,
                              rna_prot_interface_atom_ids, interface_rna_atom_ids,
                              Complex_RMSD, Protein_RMSD, RNA_RMSD, Protein_LDDT, RNA_LDDT, Protein_TM, RNA_TM, 
                              Complex_TM, Complex_LDDT, RNA_INF, RNA_DI,
                              RNA_f1_score, RNA_precision, RNA_recall, RNA_mcc, RNA_wl,
                              RNA_GlycosidicBond, RNA_SugarPucker, RNA_GammaAngle, RNAmotif_score, RNA_Stems, 
                              RNA_HairpinLoops, RNA_InternalLoops, 
                              RNA_MultibranchLoops, RNA_DanglingEnds, RNA_Pseudoknots, RNA_isDuplex, RNAMotif_isDuplex)
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                              ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 
                              ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                              )''',
                              record_values)
        else:
            # Update existing record
            cursor.execute('''UPDATE pred_protein_rna SET 
                              FileName=?, ProteinSequence=?, ProteinChainIDs=?, 
                              ProteinLength=?, NumberProteins=?, AAMotif=?, AAC=?, AAPproteinRNA=?, RNASequence=?, RNAChainIDs=?, RNALength=?, NumberRNAs=?, 
                              RNAMotif=?, RNAMotifLength=?, RNAModified=?, ContactList=?, ChainIDpairList_proteinRNA=?, 
                              Hbond_proteinRNA=?, vdWbond_proteinRNA=?, ProteinRNAInterfaceArea=?, 
                              ProteinRNAInterfaceRatio=?, Free_energy=?, Binding_affinity_kd=?, Seed=?, Model=?, 
                              af3_protein_pTM=?, af3_rna_pTM=?, af3_protein_ipTM=?, af3_rna_ipTM=?, af3_protein_pLDDT_avg=?,
                              af3_rna_pLDDT_avg=?, af3_global_pae_avg=?, af3_chain_pair_pae_min=?, af3_fraction_disordered=?,
                              af3_has_clash=?, af3_ranking_score=?, af3_protein_MSA=?, af3_rna_MSA=?,
                              rna_prot_interface_atom_ids=?, interface_rna_atom_ids=?,
                              Complex_RMSD=?, Protein_RMSD=?, RNA_RMSD=?, Protein_LDDT=?, RNA_LDDT=?, Protein_TM=?, RNA_TM=?, 
                              Complex_TM=?, Complex_LDDT=?, RNA_INF=?, RNA_DI=?,
                              RNA_f1_score=?, RNA_precision=?, RNA_recall=?, RNA_mcc=?, RNA_wl=?,
                              RNA_GlycosidicBond=?, RNA_SugarPucker=?, RNA_GammaAngle=?, RNAmotif_score=?, RNA_Stems=?, 
                              RNA_HairpinLoops=?, RNA_InternalLoops=?, 
                              RNA_MultibranchLoops=?, RNA_DanglingEnds=?, RNA_Pseudoknots=?, RNA_isDuplex=?, RNAMotif_isDuplex=?
                              WHERE exp_db_id=?''',
                              record_values[1:] + [record_values[0]])  # Move exp_db_id to the end for WHERE clause

    if classification == "protein_rna_dna":
        cursor.execute('''CREATE TABLE IF NOT EXISTS pred_protein_rna_dna (
                            Id INTEGER PRIMARY KEY AUTOINCREMENT,
                            exp_db_id TEXT,
                            FileName TEXT,
                            ProteinSequence TEXT,
                            ProteinChainIDs TEXT,
                            ProteinLength INTEGER,
                            NumberProteins INTEGER,
                            AAMotif TEXT,
                            AAC TEXT,
                            AAPproteinRNA TEXT,
                            AAPproteinDNA TEXT,
                            RNASequence TEXT,
                            RNAChainIDs TEXT,
                            RNALength INTEGER,
                            NumberRNAs INTEGER,
                            RNAMotif TEXT,
                            RNAMotifLength TEXT,
                            RNAModified TEXT,
                            DNASequence TEXT,
                            DNAChainIDs TEXT,
                            DNALength INTEGER,
                            NumberDNAs INTEGER,
                            DNAMotif TEXT,
                            DNAModified TEXT,
                            ContactList TEXT,
                            ChainIDpairList_proteinRNA TEXT,
                            Hbond_proteinRNA TEXT,
                            vdWbond_proteinRNA TEXT,
                            ChainIDpairList_proteinDNA TEXT,
                            Hbond_proteinDNA TEXT,
                            vdWbond_proteinDNA TEXT,
                            ProteinRNAInterfaceArea FLOAT,
                            ProteinDNAInterfaceArea FLOAT,
                            RNAdnaInterfaceArea FLOAT,
                            ProteinRNAInterfaceRatio FLOAT,
                            ProteinDNAInterfaceRatio FLOAT,
                            RNAdnaInterfaceRatio FLOAT,
                            Seed INTEGER,
                            Model INTEGER,
                            af3_protein_pTM TEXT,
                            af3_rna_pTM TEXT,
                            af3_dna_pTM TEXT,
                            af3_protein_ipTM TEXT,
                            af3_rna_ipTM TEXT,
                            af3_dna_ipTM TEXT,
                            af3_protein_pLDDT_avg FLOAT,
                            af3_rna_pLDDT_avg FLOAT,
                            af3_dna_pLDDT_avg FLOAT,
                            af3_global_pae_avg FLOAT,
                            af3_chain_pair_pae_min TEXT,
                            af3_fraction_disordered FLOAT,
                            af3_has_clash FLOAT,
                            af3_ranking_score FLOAT,
                            af3_protein_MSA INTEGER,
                            af3_rna_MSA INTEGER,
                            af3_dna_MSA INTEGER,
                            rna_prot_interface_atom_ids TEXT, 
                            dna_prot_interface_atom_ids TEXT,
                            interface_rna_atom_ids TEXT,
                            interface_dna_atom_ids TEXT,
                            All_RMSD TEXT,
                            Protein_RMSD TEXT, 
                            RNA_RMSD TEXT,
                            DNA_RMSD TEXT,
                            Protein_LDDT FLOAT,
                            RNA_LDDT FLOAT, 
                            DNA_LDDT FLOAT,
                            Protein_TM FLOAT,
                            RNA_TM FLOAT, 
                            DNA_TM FLOAT,
                            Protein_GDT_TS FLOAT, 
                            RNA_GDT_TS FLOAT,  
                            DNA_GDT_TS FLOAT,  
                            RNA_INF TEXT, 
                            RNA_DI FLOAT,
                            RNA_f1_score FLOAT,
                            RNA_precision FLOAT,
                            RNA_recall FLOAT,
                            RNA_mcc FLOAT,
                            RNA_wl FLOAT,
                            RNA_GlycosidicBond TEXT, 
                            RNA_SugarPucker TEXT, 
                            RNA_GammaAngle TEXT,
                            FOREIGN KEY (exp_db_id) REFERENCES exp_protein_rna(UniprotId),
                            FOREIGN KEY (exp_db_id) REFERENCES exp_protein_rna(PDBId)
                          )''')

        # Check if exp_db_id already exists in the database
        cursor.execute("SELECT 1 FROM pred_protein_rna_dna WHERE exp_db_id = ?", (data[0],))
        if cursor.fetchone() is None:
            # Insert new record
            cursor.execute('''INSERT INTO pred_protein_rna_dna (exp_db_id, FileName, ProteinSequence, ProteinChainIDs, 
                              ProteinLength, NumberProteins, AAMotif, AAC, AAPproteinRNA, AAPproteinDNA, RNASequence, 
                              RNAChainIDs, RNALength, NumberRNAs, 
                              RNAMotif, RNAMotifLength, RNAModified, DNASequence, DNAChainIDs, DNALength, NumberDNAs, 
                              DNAMotif, DNAModified, ContactList, ChainIDpairList_proteinRNA, Hbond_proteinRNA, 
                              vdWbond_proteinRNA, ChainIDpairList_proteinDNA, Hbond_proteinDNA, vdWbond_proteinDNA, 
                              ProteinRNAInterfaceArea, ProteinDNAInterfaceArea, RNAdnaInterfaceArea,
                              ProteinRNAInterfaceRatio, ProteinDNAInterfaceRatio, RNAdnaInterfaceRatio,
                              Seed, Model, af3_protein_pTM, af3_rna_pTM, af3_dna_pTM, af3_protein_ipTM, af3_rna_ipTM,
                              af3_dna_ipTM, af3_protein_pLDDT_avg, af3_rna_pLDDT_avg, af3_dna_pLDDT_avg,
                              af3_global_pae_avg, af3_chain_pair_pae_min, af3_fraction_disordered, af3_has_clash,
                              af3_ranking_score, af3_protein_MSA, af3_rna_MSA, af3_dna_MSA, 
                              rna_prot_interface_atom_ids, dna_prot_interface_atom_ids, 
                              interface_rna_atom_ids, interface_dna_atom_ids, 
                              All_RMSD, Protein_RMSD, RNA_RMSD, DNA_RMSD, Protein_LDDT, RNA_LDDT, DNA_LDDT, Protein_TM, 
                              RNA_TM, DNA_TM, Protein_GDT_TS, RNA_GDT_TS, DNA_GDT_TS, RNA_INF, RNA_DI,
                              RNA_f1_score, RNA_precision, RNA_recall, RNA_mcc, RNA_wl,
                              RNA_GlycosidicBond, RNA_SugarPucker, RNA_GammaAngle)
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 
                              ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                              ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                              record_values)
        else:
            # Update existing record
            cursor.execute('''UPDATE pred_protein_rna_dna SET 
                              FileName=?, ProteinSequence=?, ProteinChainIDs=?, 
                              ProteinLength=?, NumberProteins=?, AAMotif=?, AAC=?, AAPproteinRNA=?, AAPproteinDNA=?, RNASequence=?, 
                              RNAChainIDs=?, RNALength=?, NumberRNAs=?, 
                              RNAMotif=?, RNAMotifLength=?, RNAModified=?, DNASequence=?, DNAChainIDs=?, DNALength=?, NumberDNAs=?, 
                              DNAMotif=?, DNAModified=?, ContactList=?, ChainIDpairList_proteinRNA=?, Hbond_proteinRNA=?, 
                              vdWbond_proteinRNA=?, ChainIDpairList_proteinDNA=?, Hbond_proteinDNA=?, vdWbond_proteinDNA=?, 
                              ProteinRNAInterfaceArea=?, ProteinDNAInterfaceArea=?, RNAdnaInterfaceArea=?,
                              ProteinRNAInterfaceRatio=?, ProteinDNAInterfaceRatio=?, RNAdnaInterfaceRatio=?,
                              Seed=?, Model=?, af3_protein_pTM=?, af3_rna_pTM=?, af3_dna_pTM=?, af3_protein_ipTM=?, af3_rna_ipTM=?,
                              af3_dna_ipTM=?, af3_protein_pLDDT_avg=?, af3_rna_pLDDT_avg=?, af3_dna_pLDDT_avg=?,
                              af3_global_pae_avg=?, af3_chain_pair_pae_min=?, af3_fraction_disordered=?, af3_has_clash=?,
                              af3_ranking_score=?, af3_protein_MSA=?, af3_rna_MSA=?, af3_dna_MSA=?, 
                              rna_prot_interface_atom_ids=?, dna_prot_interface_atom_ids=?, 
                              interface_rna_atom_ids=?, interface_dna_atom_ids=?, 
                              All_RMSD=?, Protein_RMSD=?, RNA_RMSD=?, DNA_RMSD=?, Protein_LDDT=?, RNA_LDDT=?, DNA_LDDT=?, Protein_TM=?, 
                              RNA_TM=?, DNA_TM=?, Protein_GDT_TS=?, RNA_GDT_TS=?, DNA_GDT_TS=?, RNA_INF=?, RNA_DI=?,
                              RNA_f1_score=?, RNA_precision=?, RNA_recall=?, RNA_mcc=?, RNA_wl=?,
                              RNA_GlycosidicBond=?, RNA_SugarPucker=?, RNA_GammaAngle=?
                              WHERE exp_db_id=?''',
                              record_values[1:] + [record_values[0]])  # Move exp_db_id to the end for WHERE clause

    if classification == "protein_dna":
        cursor.execute('''CREATE TABLE IF NOT EXISTS pred_protein_dna (
                            Id INTEGER PRIMARY KEY AUTOINCREMENT,
                            exp_db_id TEXT,
                            FileName TEXT,
                            ProteinSequence TEXT,
                            ProteinChainIDs TEXT,
                            ProteinLength INTEGER,
                            NumberProteins INTEGER,
                            AAMotif TEXT,
                            AAC TEXT,
                            AAPproteinDNA TEXT,
                            DNASequence TEXT,
                            DNAChainIDs TEXT,
                            DNALength INTEGER,
                            NumberDNAs INTEGER,
                            DNAMotif TEXT,
                            DNAModified TEXT,
                            ContactList TEXT,
                            ChainIDpairList_proteinDNA TEXT,
                            Hbond_proteinDNA TEXT,
                            vdWbond_proteinDNA TEXT,
                            ProteinDNAInterfaceArea FLOAT,
                            ProteinDNAInterfaceRatio FLOAT,
                            Free_energy REAL,
                            Binding_affinity_kd TEXT,
                            Seed INTEGER,
                            Model INTEGER,
                            af3_protein_pTM TEXT,
                            af3_dna_pTM TEXT,
                            af3_protein_ipTM TEXT,
                            af3_dna_ipTM TEXT,
                            af3_protein_pLDDT_avg FLOAT,
                            af3_dna_pLDDT_avg FLOAT,
                            af3_global_pae_avg FLOAT,
                            af3_chain_pair_pae_min TEXT,
                            af3_fraction_disordered FLOAT,
                            af3_has_clash FLOAT,
                            af3_ranking_score FLOAT,
                            af3_protein_MSA INTEGER,
                            af3_dna_MSA INTEGER,
                            dna_prot_interface_atom_ids TEXT,  
                            interface_dna_atom_ids TEXT, 
                            All_RMSD TEXT,
                            Protein_RMSD TEXT, 
                            DNA_RMSD TEXT,
                            Protein_LDDT FLOAT,
                            DNA_LDDT FLOAT,
                            Protein_TM FLOAT,
                            DNA_TM FLOAT,
                            Protein_GDT_TS FLOAT, 
                            DNA_GDT_TS FLOAT, 
                            FOREIGN KEY (exp_db_id) REFERENCES exp_protein_rna(UniprotId),
                            FOREIGN KEY (exp_db_id) REFERENCES exp_protein_rna(PDBId)
                          )''')

        # Check if exp_db_id already exists in the database
        cursor.execute("SELECT 1 FROM pred_protein_dna WHERE exp_db_id = ?", (data[0],))
        if cursor.fetchone() is None:
            # Insert new record
            cursor.execute('''INSERT INTO pred_protein_dna (exp_db_id, FileName, ProteinSequence, ProteinChainIDs,
                              ProteinLength, NumberProteins, AAMotif, AAC, AAPproteinDNA, DNASequence, DNAChainIDs, 
                              DNALength, NumberDNAs, 
                              DNAMotif, DNAModified, ContactList, ChainIDpairList_proteinDNA,Hbond_proteinDNA, vdWbond_proteinDNA,
                              ProteinDNAInterfaceArea, ProteinDNAInterfaceRatio, Free_energy, Binding_affinity_kd,
                              Seed, Model, af3_protein_pTM, af3_rna_pTM, af3_protein_ipTM, af3_rna_ipTM, af3_protein_pLDDT_avg,
                              af3_rna_pLDDT_avg, af3_global_pae_avg, af3_chain_pair_pae_min, af3_fraction_disordered,
                              af3_has_clash, af3_ranking_score, af3_protein_MSA,
                              af3_dna_MSA, dna_prot_interface_atom_ids, 
                              interface_dna_atom_ids, All_RMSD, Protein_RMSD, DNA_RMSD, Protein_LDDT, 
                              DNA_LDDT, Protein_TM, DNA_TM, Protein_GDT_TS, DNA_GDT_TS)
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                              ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', record_values)
        else:
            # Update existing record
            cursor.execute('''UPDATE pred_protein_dna SET 
                              FileName=?, ProteinSequence=?, ProteinChainIDs=?,
                              ProteinLength=?, NumberProteins=?, AAMotif=?, AAC=?, AAPproteinDNA=?, DNASequence=?, DNAChainIDs=?, 
                              DNALength=?, NumberDNAs=?, 
                              DNAMotif=?, DNAModified=?, ContactList=?, ChainIDpairList_proteinDNA=?, Hbond_proteinDNA=?, vdWbond_proteinDNA=?,
                              ProteinDNAInterfaceArea=?, ProteinDNAInterfaceRatio=?, Free_energy=?, Binding_affinity_kd=?,
                              Seed=?, Model=?, af3_protein_pTM=?, af3_rna_pTM=?, af3_protein_ipTM=?, af3_rna_ipTM=?, af3_protein_pLDDT_avg=?,
                              af3_rna_pLDDT_avg=?, af3_global_pae_avg=?, af3_chain_pair_pae_min=?, af3_fraction_disordered=?,
                              af3_has_clash=?, af3_ranking_score=?, af3_protein_MSA=?,
                              af3_dna_MSA=?, dna_prot_interface_atom_ids=?, 
                              interface_dna_atom_ids=?, All_RMSD=?, Protein_RMSD=?, DNA_RMSD=?, Protein_LDDT=?, 
                              DNA_LDDT=?, Protein_TM=?, DNA_TM=?, Protein_GDT_TS=?, DNA_GDT_TS=?
                              WHERE exp_db_id=?''',
                              record_values[1:] + [record_values[0]])  # Move exp_db_id to the end for WHERE clause

    if classification == "rna_rna":
        cursor.execute('''CREATE TABLE IF NOT EXISTS pred_rna_rna (
                            Id INTEGER PRIMARY KEY AUTOINCREMENT,
                            exp_db_id TEXT,
                            FileName TEXT,
                            RNASequence TEXT,
                            RNAChainIDs TEXT,
                            RNALength INTEGER,
                            NumberRNAs INTEGER,
                            RNAMotif TEXT,
                            RNAMotifLength TEXT,
                            RNAModified TEXT,
                            ContactList TEXT,
                            ChainIDpairList_rnaRNA TEXT,
                            Hbond_rnaRNA TEXT,
                            vdWbond_rnaRNA TEXT,
                            RNArnaInterfaceArea FLOAT,
                            RNArnaInterfaceRatio FLOAT,
                            Free_energy REAL,
                            Binding_affinity_kd TEXT,
                            Seed INTEGER,
                            Model INTEGER, 
                            af3_rna_pTM TEXT, 
                            af3_rna_ipTM TEXT, 
                            af3_rna_pLDDT_avg FLOAT,
                            af3_global_pae_avg FLOAT,
                            af3_chain_pair_pae_min TEXT,
                            af3_fraction_disordered FLOAT,
                            af3_has_clash FLOAT,
                            af3_ranking_score FLOAT,
                            af3_rna_MSA INTEGER,
                            rna_rna_interface_atom_ids TEXT,  
                            interface_rna_atom_ids TEXT, 
                            Complex_RMSD TEXT,
                            RNA_RMSD TEXT,
                            RNA_LDDT FLOAT, 
                            Complex_LDDT FLOAT,
                            RNA_TM FLOAT, 
                            Complex_TM FLOAT, 
                            RNA_INF TEXT, 
                            RNA_DI FLOAT,
                            RNA_f1_score FLOAT,
                            RNA_precision FLOAT,
                            RNA_recall FLOAT,
                            RNA_mcc FLOAT,
                            RNA_wl FLOAT,
                            RNA_GlycosidicBond TEXT, 
                            RNA_SugarPucker TEXT, 
                            RNA_GammaAngle TEXT,
                            FOREIGN KEY (exp_db_id) REFERENCES exp_protein_rna(UniprotId),
                            FOREIGN KEY (exp_db_id) REFERENCES exp_protein_rna(PDBId)
                          )''')

        # Check if exp_db_id already exists in the database
        cursor.execute("SELECT 1 FROM pred_rna_rna WHERE exp_db_id = ?", (data[0],))
        if cursor.fetchone() is None:
            # Insert new record
            cursor.execute('''INSERT INTO pred_rna_rna (exp_db_id, FileName, RNASequence, RNAChainIDs, RNALength, 
                              NumberRNAs, RNAMotif, RNAMotifLength, RNAModified, ContactList,
                              ChainIDpairList_rnaRNA, Hbond_rnaRNA, vdWbond_rnaRNA, 
                              RNArnaInterfaceArea, RNArnaInterfaceRatio, Free_energy, Binding_affinity_kd, Seed, Model, 
                              af3_rna_pTM, af3_rna_ipTM, 
                              af3_rna_pLDDT_avg, af3_global_pae_avg, af3_chain_pair_pae_min, af3_fraction_disordered,
                              af3_has_clash, af3_ranking_score, af3_rna_MSA,
                              rna_rna_interface_atom_ids, interface_rna_atom_ids, 
                              Complex_RMSD, RNA_RMSD, RNA_LDDT, Complex_LDDT, RNA_TM, Complex_TM, RNA_INF, RNA_DI, 
                              RNA_f1_score, RNA_precision, RNA_recall, RNA_mcc, RNA_wl,
                              RNA_GlycosidicBond, RNA_SugarPucker, RNA_GammaAngle)
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                              ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', record_values)
        else:
            # Update existing record
            cursor.execute('''UPDATE pred_rna_rna SET 
                              FileName=?, RNASequence=?, RNAChainIDs=?, RNALength=?, 
                              NumberRNAs=?, RNAMotif=?, RNAMotifLength=?, RNAModified=?, ContactList=?,
                              ChainIDpairList_rnaRNA=?, Hbond_rnaRNA=?, vdWbond_rnaRNA=?, 
                              RNArnaInterfaceArea=?, RNArnaInterfaceRatio=?, Free_energy=?, Binding_affinity_kd=?, Seed=?, Model=?, 
                              af3_rna_pTM=?, af3_rna_ipTM=?, 
                              af3_rna_pLDDT_avg=?, af3_global_pae_avg=?, af3_chain_pair_pae_min=?, af3_fraction_disordered=?,
                              af3_has_clash=?, af3_ranking_score=?, af3_rna_MSA=?,
                              rna_rna_interface_atom_ids=?, interface_rna_atom_ids=?, 
                              Complex_RMSD=?, RNA_RMSD=?, RNA_LDDT=?, Complex_LDDT=?, RNA_TM=?, Complex_TM=?, RNA_INF=?, RNA_DI=?, 
                              RNA_f1_score=?, RNA_precision=?, RNA_recall=?, RNA_mcc=?, RNA_wl=?,
                              RNA_GlycosidicBond=?, RNA_SugarPucker=?, RNA_GammaAngle=?
                              WHERE exp_db_id=?''',
                              record_values[1:] + [record_values[0]])  # Move exp_db_id to the end for WHERE clause

    connection.commit()
    connection.close()

# If all or most AF3 cif files are processable:
# def process_cif_files(folder_path):
#     for filename in os.listdir(folder_path):
#         if filename.endswith(".cif"):
#             print(filename)
#             file_path = os.path.join(folder_path, filename)
#             classification, data = parse_cif_file(file_path)
#             # print(data)
#             # print(classification)
#             insert_pred_data(classification, data)

# if __name__ == "__main__":
#     import sys
#
#     if len(sys.argv) != 2:
#         print("Usage: python parseAFcif2DB.py <folder_path>")
#     else:
#         folder_path = sys.argv[1]
#         process_cif_files(folder_path)


def process_cif_files(af_folder):
    """Process AF files and move failed ones to check folder"""

    # Create check folder if it doesn't exist
    check_folder = os.path.join(os.path.dirname(af_folder), 'af_check_cif')
    if not os.path.exists(check_folder):
        os.makedirs(check_folder)

    cif_files = [f for f in os.listdir(af_folder) if f.endswith('.cif')]

    for cif_file in cif_files:
        try:
            # Get PDB ID from filename
            pdb_id = cif_file.split('_s1_model_')[0].split('fold_')[1]
            print(f"Processing {pdb_id}")

            # Get all related files
            related_files = [
                cif_file,  # The CIF file
                f"fold_{pdb_id}_s1_full_data_0.json",  # Full data JSON
                f"fold_{pdb_id}_s1_summary_confidences_0.json"  # Summary JSON
            ]

            # Try to parse CIF file
            cif_path = os.path.join(af_folder, cif_file)
            try:
                classification, data = parse_cif_file(cif_path)
                insert_pred_data(classification, data)
            except Exception as e:
                print(f"Error parsing {cif_file}: {e}")

                # Move all related files to check folder
                for related_file in related_files:
                    src_path = os.path.join(af_folder, related_file)
                    if os.path.exists(src_path):
                        dst_path = os.path.join(check_folder, related_file)
                        print(f"Moving {related_file} to {check_folder}")
                        shutil.move(src_path, dst_path)
                continue

        except Exception as e:
            print(f"Error processing {cif_file}: {e}")
            continue



if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python parseAFcif2DB.py <af_folder>")
        sys.exit(1)

    af_folder = sys.argv[1]
    if not os.path.isdir(af_folder):
        print(f"Error: {af_folder} is not a valid directory")
        sys.exit(1)

    process_cif_files(af_folder)