from pymol import cmd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os
import sys
ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import glob
from structures.rna import RNA
from structures.dna import DNA
import json
import ast
from collections import Counter
import re
from itertools import product
from database.databaseMethods import DatabaseMethods
import glob
import shutil

database_methods = DatabaseMethods()
print("Database connection initialized.")
# https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/#MLformula
# https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
# https://www.ebi.ac.uk/pdbe/entry/pdb/5szc/modified/MLZ
amino_acid_atoms = {
    "ALA": {"CB": "side_chain"},
    "ARG": {"CB": "side_chain", "CG": "side_chain", "CD": "side_chain", "NE": "side_chain", "NH1": "side_chain", "NH2": "side_chain", "CZ": "side_chain"},
    "ASN": {"CB": "side_chain", "CG": "side_chain", "OD1": "side_chain", "ND2": "side_chain"},
    "ASP": {"CB": "side_chain", "CG": "side_chain", "OD1": "side_chain", "OD2": "side_chain"},
    "CYS": {"CB": "side_chain", "SG": "side_chain"},
    "GLN": {"CB": "side_chain", "CG": "side_chain", "CD": "side_chain", "OE1": "side_chain", "NE2": "side_chain"},
    "GLU": {"CB": "side_chain", "CG": "side_chain", "CD": "side_chain", "OE1": "side_chain", "OE2": "side_chain"},
    "GLY": {},
    "HIS": {"CB": "side_chain", "CG": "side_chain", "ND1": "side_chain", "NE2": "side_chain", "CD2": "side_chain", "CE1": "side_chain"},
    "ILE": {"CB": "side_chain", "CG1": "side_chain", "CG2": "side_chain", "CD1": "side_chain"},
    "LEU": {"CB": "side_chain", "CG": "side_chain", "CD1": "side_chain", "CD2": "side_chain"},
    "LYS": {"CB": "side_chain", "CG": "side_chain", "CD": "side_chain", "CE": "side_chain", "NZ": "side_chain"},
    "MET": {"CB": "side_chain", "CG": "side_chain", "SD": "side_chain", "CE": "side_chain"},
    "PHE": {"CB": "side_chain", "CG": "side_chain", "CD1": "side_chain", "CD2": "side_chain", "CE1": "side_chain", "CE2": "side_chain", "CZ": "side_chain"},
    "PRO": {"CB": "side_chain", "CG": "side_chain", "CD": "side_chain"},
    "SER": {"CB": "side_chain", "OG": "side_chain"},
    "THR": {"CB": "side_chain", "OG1": "side_chain", "CG2": "side_chain"},
    "TRP": {"CB": "side_chain", "CG": "side_chain", "CD1": "side_chain", "CD2": "side_chain", "NE1": "side_chain", "CE2": "side_chain", "CE3": "side_chain", "CZ2": "side_chain", "CZ3": "side_chain", "CH2": "side_chain"},
    "TYR": {"CB": "side_chain", "CG": "side_chain", "CD1": "side_chain", "CD2": "side_chain", "CE1": "side_chain", "CE2": "side_chain", "CZ": "side_chain", "OH": "side_chain"},
    "VAL": {"CB": "side_chain", "CG1": "side_chain", "CG2": "side_chain"},
    "MLZ": {"CB": "side_chain", "CG": "side_chain", "CD": "side_chain", "CE": "side_chain", "NZ": "side_chain", "CM": "side_chain"},
    "SEP": {"CB": "side_chain", "OG": "side_chain", "OXT": "side_chain", "P": "side_chain", "O1P": "side_chain",
            "O2P": "side_chain", "O3P": "side_chain"},
    "N": "backbone",
    "CA": "backbone",
    "C": "backbone",
    "O": "backbone"
}

rna_sugar = ["C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'", "O4'", "O5'", "CM'", "CM2"] # "CM'" of A2M, "CM2" of OMU/C/G
common_phosphate = ["P", "OP1", "OP2"]

# Combined dictionary for both single-letter and two-letter nucleotide representations https://www.nakb.org/modifiednt.html?
nucleotide_atoms = {
    "A": {"base": ["N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"]},
    "A2M": {"base": ["N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"]},
    "G": {"base": ["N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9"]},
    "C": {"base": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]},
    "DA": {"base": ["N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"]},
    "DG": {"base": ["N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9"]},
    "DC": {"base": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]},
    "5MC": { # 5-METHYLCYTIDINE-5'-MONOPHOSPHATE: 5MC: P, OP1, OP2, O5', C5', C4', O4', C3', O3', C2', O2', C1', N1, C2, O2, N3, C4, N4, C5, C6, CM5
        "base": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6", "CM5"]},
    "6MZ": { # N6-METHYLADENOSINE-5'-MONOPHOSPHATE
        "base": ["N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9", "C9"]},
    "DT": {"base": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6"]},
    "U": {"base": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]},
    "OMU": {"base": ["N1", "C2", "N3", "C4", "C5", "C6", "O2", "O4"]},
    "OMC": {"base": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]},
    "OMG": {"base": ["N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9"]},
    "5MU": {"base": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "C5M"]},
    "PSU": {"base": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]},
    "UR3": {"base": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "C3U"]}
}

# AF-supported modified nucleotides:  Biologically common chemical modifications of the nucleic acids:
#     DNA: Methylation of cytosine, guanine, and adenine
#         Carboxylation of cytosine
#         Oxidation of guanine
#         Formylation of cytosine
#     RNA:
#         Methylation of cytosine, guanine, adenine, and uracil
#         Isomerisation of uridine into pseudouridine
#         Formylation of cytosine

def classify_atom(resn, name):
    if resn in amino_acid_atoms and len(resn) == 3:
        if name in amino_acid_atoms[resn]:
            return amino_acid_atoms[resn][name]
        elif name in amino_acid_atoms:
            return amino_acid_atoms[name]
        else:
            return "unknown"
    elif resn in nucleotide_atoms:
        for key, atom_list in nucleotide_atoms[resn].items():
            if name in atom_list:
                return key
            elif name in rna_sugar:
                return "sugar"
            elif name in common_phosphate:
                return "phosphate"
        return "unknown"
    else:
        return "unknown"

def calculate_highest_contacts(chain_pairs):
    max_contacts = 0
    best_peptide_chain = None
    for polynucleotide_chain, peptide_chain in chain_pairs:
        hb = cmd.find_pairs(f"chain {peptide_chain} & e. n+o", f"chain {polynucleotide_chain} & e. n+o", cutoff=3.2, mode=1, angle=50)
        num_contacts = len(hb)
        print(f"Number of contacts between {peptide_chain} and {polynucleotide_chain}: {num_contacts}")
        if num_contacts > max_contacts:
            max_contacts = num_contacts
            best_peptide_chain = peptide_chain
    if best_peptide_chain:
        print(f"The peptide chain with the highest contacts is {best_peptide_chain}, with {max_contacts} contacts.")
    else:
        print("No contacts found.")
    return best_peptide_chain


def get_resi_from_polynucleotide_chain(polynucleotide_chain, file_path):  # Get full resi list to avoid motif: A(-3), U(-2), U(-1), U(0)...
    # Create check folder if it doesn't exist
    check_folder = os.path.join(os.path.dirname(os.path.dirname(file_path)), 'check_pdb')
    if not os.path.exists(check_folder):
        os.makedirs(check_folder)

    resi_set = set()

    try:
        # Select all residues from the specific polynucleotide chain
        selection = f"chain {polynucleotide_chain} and polymer"
        cmd.iterate(selection, 'resi_set.add(resi)', space={'resi_set': resi_set})

        # Try to sort numerically
        resi_list = sorted(resi_set, key=int)
        return resi_list
    except ValueError:
        # If sorting fails due to non-numeric residue IDs
        print(f"Warning: Found non-numeric residue IDs {resi_set}. Skipping this chain.")
        # Move file to check folder
        file_name = os.path.basename(file_path)
        dst_path = os.path.join(check_folder, file_name)
        print(f"Moving {file_name} to {check_folder}")
        import shutil
        shutil.move(file_path, dst_path)
        return None
    except Exception as e:
        print(f"Error processing chain {polynucleotide_chain}: {e}")
        # Move file to check folder
        file_name = os.path.basename(file_path)
        dst_path = os.path.join(check_folder, file_name)
        print(f"Moving {file_name} to {check_folder}")
        import shutil
        shutil.move(file_path, dst_path)
        return None

# Based on Robert L. Campbell, 2010, https://static.igem.org/mediawiki/2015/a/a1/Space-p_list_hbonds.txt
# List interactions between selections specified in a .cif file and save the results in the same directory as the input file.
def list_interactions_from_cif(file_path, list_name):
    # Create check folder if it doesn't exist
    check_folder = os.path.join(os.path.dirname(os.path.dirname(file_path)), 'check_pdb')
    if not os.path.exists(check_folder):
        os.makedirs(check_folder)

    try:
        global skipped_pdbs  # Make it accessible globally
        print(file_path)
        pdb_info = MMCIF2Dict(file_path)
        strand_ids = pdb_info.get('_entity_poly.pdbx_strand_id', [])
        entity_ids = pdb_info.get('_entity_poly.entity_id', [])
        entity_types = pdb_info.get('_entity_poly.type', [])
        is_protein = False
        is_rna = False
        is_dna = False
        rna_motif_stringS = []
        dna_motif_stringS = []
        aa_motif_stringS = []
        contact_listS = []
        rna_prot_interface_atom_idsS = []
        interface_rna_atom_idsS = []
        dna_prot_interface_atom_idsS = []
        interface_dna_atom_idsS = []
        length_rna_motif = []
        aac_list = []
        chain_id_pairs_rna = []
        h_bonds_rna = []
        vdW_bonds_rna = []
        chain_id_pairs_dna = []
        h_bonds_dna = []
        vdW_bonds_dna = []

        if not entity_ids or not strand_ids or not entity_types:
            raise ValueError("Entity IDs, strand IDs, or entity types not found in the CIF file.")

        polynucleotide_chains = [strand for strand, etype in zip(strand_ids, entity_types) if
                                 'polyribonucleotide' in etype or 'polydeoxyribonucleotide' in etype]

        peptide_chains = [strand for strand, etype in zip(strand_ids, entity_types) if 'polypeptide(L)' in etype]

        if len(strand_ids) < 2: # Only for some RNA-RNA cif files from pdb
            adjusted_strand_ids = []
            for strand in strand_ids:
                if ',' in strand: # _entity_poly.type: polyribonucleotide, _entity_poly.pdbx_strand_id: ['A,B'] instead of ['A','B']
                    adjusted_strand_ids.extend(strand.split(','))
                else:
                    raise ValueError("Not enough entity types found in the CIF file.")
            strand_ids = adjusted_strand_ids
            polynucleotide_chains = strand_ids

        if not polynucleotide_chains:
            raise ValueError("No polynucleotide chains found in the CIF file.")

            # If RNA-RNA also accounted, then:
            # Determine the loops to use based on the presence of peptide chains
        def split_and_pair(chain_pairs):
            final_pairs = []
            for pair in chain_pairs:
                polynucleotide_chain, peptide_chain = pair
                poly_chains = polynucleotide_chain.split(',')
                pep_chains = peptide_chain.split(',')
                # Cartesian product between poly_chains and pep_chains
                final_pairs.extend(product(poly_chains, pep_chains))

            return final_pairs

        if peptide_chains:
            chain_pairs = [(polynucleotide_chain, pep_chain) for polynucleotide_chain in polynucleotide_chains for pep_chain
                           in
                           peptide_chains]
        else:
            chain_pairs = [(polynucleotide_chain, next_poly_chain) for polynucleotide_chain in polynucleotide_chains for
                           next_poly_chain in
                           polynucleotide_chains if polynucleotide_chain != next_poly_chain]

        chain_pairs = split_and_pair(chain_pairs)

        # if peptide_chain and polynucleotide_chain:
        for polynucleotide_chain, peptide_chain in chain_pairs:
            # print("polynucleotide_chain, peptide_chain polynucleotide_chain, peptide_chain ", polynucleotide_chain, peptide_chain )
        # if polynucleotide_chain and peptide_chain:
            # print("for polynucleotide_chain, peptide_chain in chain_pairs", "polynucleotide_chain", polynucleotide_chain,
            #       "peptide_chain", peptide_chain)
            if 'polypeptide(L)' in entity_types:
                is_protein = True
            # # Ensure each pair is processed only once
            # if (polynucleotide_chain, peptide_chain) in processed_pairs or (peptide_chain, polynucleotide_chain) in processed_pairs:
            #     break
            # processed_pairs.add((polynucleotide_chain, peptide_chain))

        # If only peptide and polynucleotide chains involved, the following suffices (instead of the above starting from if peptide_chains):
        # Don't forget to indent the block
        # for polynucleotide_chain in polynucleotide_chains:
        #     for peptide_chain in peptide_chains:
            # print(f"Polynucleotide chain: {polynucleotide_chain}, Type: {entity_types[strand_ids.index(polynucleotide_chain)]}")

            # Find hydrogen bonds involving nitrogen and oxygen, not fluorine.
            hb = cmd.find_pairs(f"chain {peptide_chain} & e. n+o", f"chain {polynucleotide_chain} & e. n+o", cutoff=3.2, mode=1, angle=50)
            # DOCUMENTAION: angle = float: hydrogen bond angle cutoff, only if mode=1 {default: 45.0}
            # cutoff = DA distance (not HA distance)
            # angle = angle range of the donor–acceptor–acceptor antecedent angle, which is smaller than D-H-A???
            # e.g. DA distance, DHA angle, H-A distance, H-A-AA angle, D-A-AA angle in HBplus
            # 3.04 157.5 2.09 115.7 113.3 -> VDW @ 40, 3.2 A; HB @ 50, 3.2 A (cutoff: 50 -> 130º?)
            # 3.04 159.9 2.08 144.4 144.2 -> WDW @ 30, 3.2 A; HB @ 40, 3.2 A & 50, 3.2 A (cutoff: 40 -> 140º?)
            # To qualify as a hydrogen bond the hydrogen bond donor and acceptor atoms should be separated by 3.5 A
            # or less and the donor-acceptor-acceptor antecedent angle (Ĥ) should be larger than 100°.
            # The GROMACS convention is the angle corresponds to H-D-A rather than D-H-A like some other
            # > >> programs. H-D-A cutoff of 30 degrees is equivalent to the "conventional" 150-degree D-H-A cutoff.

            # Get residue pairs involved in hydrogen bonds
            hb_residue_pairs = set()
            for pair in hb:
                chain1 = cmd.get_model(f"{pair[0][0]} and index {pair[0][1]}").atom[0].resi
                chain2 = cmd.get_model(f"{pair[1][0]} and index {pair[1][1]}").atom[0].resi
                hb_residue_pairs.add((chain1, chain2))

            # Find van der Waals interactions excluding hydrogen atoms
            vdw = cmd.find_pairs(f"chain {peptide_chain} and not elem H", f"chain {polynucleotide_chain} and not elem H", cutoff= 3.9, mode = 0)

            shortest_interactions = {}
            combined_interactions = {}
            number_hbonds = 0
            number_vdWBonds = 0

            # Add hydrogen bonds to the combined list
            for pairs in hb:
                chain1_info = cmd.get_model(f"{pairs[0][0]} and index {pairs[0][1]}").atom[0]
                chain2_info = cmd.get_model(f"{pairs[1][0]} and index {pairs[1][1]}").atom[0]
                distance = cmd.distance(list_name, f"{pairs[0][0]} and index {pairs[0][1]}",
                                        f"{pairs[1][0]} and index {pairs[1][1]}")

                if chain1_info.resn != "HOH" and chain2_info.resn != "HOH":  # Exclude water molecules
                    interaction_type = "HB"
                    chain1_class = classify_atom(chain1_info.resn, chain1_info.name)
                    chain2_class = classify_atom(chain2_info.resn, chain2_info.name)
                    residue_pair = (chain1_info.resi, chain2_info.resi)
                    combined_interactions[residue_pair] = (chain1_info, chain2_info, distance, interaction_type, chain1_class, chain2_class, pairs[0][1], pairs[1][1])
                    # print(", pairs[0][1], pairs[1][1]", pairs[0][1], pairs[1][1])
                    number_hbonds += 1

            # Add van der Waals interactions to the combined list
            for pairs in vdw:
                chain1_info = cmd.get_model(f"{pairs[0][0]} and index {pairs[0][1]}").atom[0]
                chain2_info = cmd.get_model(f"{pairs[1][0]} and index {pairs[1][1]}").atom[0]
                distance = cmd.distance(list_name, f"{pairs[0][0]} and index {pairs[0][1]}",
                                        f"{pairs[1][0]} and index {pairs[1][1]}")

                if chain1_info.resn != "HOH" and chain2_info.resn != "HOH":  # Exclude water molecules
                    residue_pair = (chain1_info.resi, chain2_info.resi)
                    hb_residue_pair = (chain1_info.resi, chain2_info.resi)
                    if hb_residue_pair not in hb_residue_pairs:
                        if residue_pair not in shortest_interactions or distance < shortest_interactions[residue_pair][2]:
                            shortest_interactions[residue_pair] = (chain1_info, chain2_info, distance)
                            interaction_type = "VDW"
                            chain1_class = classify_atom(chain1_info.resn, chain1_info.name)
                            chain2_class = classify_atom(chain2_info.resn, chain2_info.name)
                            combined_interactions[residue_pair] = (chain1_info, chain2_info, distance, interaction_type, chain1_class, chain2_class, pairs[0][1], pairs[1][1])
                            number_vdWBonds += 1

            # print(file_path)
            try:
                sorted_combined_interactions = sorted(combined_interactions.values(),
                                              key=lambda x: (int(x[0].resi), int(x[1].resi)))
            except ValueError as e:
                # Move file to check folder on residue numbering error
                file_name = os.path.basename(file_path)
                dst_path = os.path.join(check_folder, file_name)
                print(f"Error processing residue numbers in {file_name}: {e}")
                print(f"Moving {file_name} to {check_folder}")
                import shutil
                shutil.move(file_path, dst_path)
                return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None

            def find_entity_type(chain, strand_ids, entity_types):
                for i, strand_id in enumerate(strand_ids):
                    split_strands = strand_id.split(',')
                    if chain in split_strands:
                        return entity_types[i]

            if sorted_combined_interactions:
                output_file = os.path.splitext(file_path)[0] + "_" + peptide_chain + polynucleotide_chain + "_interactions.txt"
                with open(output_file, 'w') as f:
                    if (find_entity_type(polynucleotide_chain, strand_ids, entity_types) == 'polyribonucleotide'):
                    # if (entity_types[strand_ids.index(polynucleotide_chain)] == 'polyribonucleotide'):
                        is_rna = True
                        pair = peptide_chain, polynucleotide_chain
                        chain_id_pairs_rna.append(pair)
                        h_bonds_rna.append(number_hbonds)
                        vdW_bonds_rna.append(number_vdWBonds)
                        # if (entity_types[strand_ids.index(peptide_chain)] != 'polypeptide(L)'):
                        if (find_entity_type(peptide_chain, strand_ids, entity_types) != 'polypeptide(L)'):
                            f.write("RNA Nucleotide Motif:\n")
                            rna_motif_string, length_rna_motif_int = extract_base_to_base_nucleotides(sorted_combined_interactions)
                            if rna_motif_string:
                                rna_motif_stringS.append(rna_motif_string)
                                length_rna_motif.append(length_rna_motif_int)
                            f.write(rna_motif_string + "\n")
                        else:
                            polynucleotide_chainA = get_resi_from_polynucleotide_chain(polynucleotide_chain, file_path)
                            if polynucleotide_chainA == None:
                                return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None
                            rna_motif_string, rna_motif_visual, rna_motif_indices, chain_id, length_rna_motif_int = extract_rna_motif(
                                sorted_combined_interactions, polynucleotide_chainA, file_path)
                            rna_motif_stringS.append(rna_motif_string)
                            length_rna_motif.append(length_rna_motif_int)
                            f.write("RNA Nucleotide Motif:\n")
                            f.write(rna_motif_string + "\n")
                            f.write("RNA Extracted Nucleotide Motif:\n")
                            f.write(rna_motif_visual + "\n")
                            full_rna_sequence, chain = get_full_rna_string_from_db(file_path, True)
                            # If a list of sequences:
                            if full_rna_sequence.startswith('['):
                                chain = ast.literal_eval(chain)
                                index = next((i for i, pair in enumerate(chain) if chain_id in pair.split(',')), None)
                                full_rna_sequence = ast.literal_eval(full_rna_sequence)
                                full_rna_sequence = full_rna_sequence[index]
                            highlighted_sequence = highlight_sequence(full_rna_sequence, rna_motif_indices)
                            f.write("Full RNA Sequence with Highlights:\n")
                            f.write(highlighted_sequence + "\n")
                    # elif (entity_types[strand_ids.index(polynucleotide_chain)] == 'polydeoxyribonucleotide'):
                    elif (find_entity_type(polynucleotide_chain, strand_ids, entity_types) == 'polydeoxyribonucleotide'):
                        polynucleotide_chainA = get_resi_from_polynucleotide_chain(polynucleotide_chain, file_path)
                        if polynucleotide_chainA is None:
                            return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None
                        dna_motif_string, dna_motif_visual, dna_motif_indices, chain_id = extract_dna_motif(sorted_combined_interactions, polynucleotide_chainA, file_path)
                        dna_motif_stringS.append(dna_motif_string)
                        is_dna = True
                        pair = peptide_chain, polynucleotide_chain
                        chain_id_pairs_dna.append(pair)
                        h_bonds_dna.append(number_hbonds)
                        vdW_bonds_dna.append(number_vdWBonds)
                        f.write("DNA Nucleotide Motif:\n")
                        f.write(dna_motif_string + "\n")
                        f.write("DNA Extracted Nucleotide Motif:\n")
                        f.write(dna_motif_visual + "\n")
                        full_na_sequence, chain = get_full_rna_string_from_db(file_path, False)
                        if full_na_sequence.startswith('['):
                            chain = ast.literal_eval(chain)
                            index = next((i for i, pair in enumerate(chain) if chain_id in pair.split(',')), None)
                            full_na_sequence = ast.literal_eval(full_na_sequence)
                            full_na_sequence = full_na_sequence[index]
                        full_dna_sequence = transform_dna_sequence(full_na_sequence)
                        highlighted_sequence = highlight_sequence(full_dna_sequence, dna_motif_indices)
                        f.write("Full DNA Sequence with Highlights:\n")
                        f.write(highlighted_sequence + "\n")

                    # if not (polynucleotide_chain in polynucleotide_chains and peptide_chain in polynucleotide_chains):
                    if not all(chain in polynucleotide_chains for chain in (polynucleotide_chain, peptide_chain)):
                        f.write("\n")
                        f.write("Peptide Binding Pocket:\n")
                        aa_motif_string = extract_aa_motif(sorted_combined_interactions)
                        aa_motif_stringS.append(aa_motif_string)
                        aac = calculate_aac(aa_motif_string)
                        aac_list.append(aac)
                        f.write(aa_motif_string + "\n")
                    else:
                        aa_motif_string = ""

                    f.write("\n")
                    f.write("Hydrogen Bond and van der Waals Interactions:\n")
                    for chain1_info, chain2_info, distance, interaction_type, chain1_class, chain2_class, chain1_id, chain2_id in sorted_combined_interactions:
                        f.write(f"{interaction_type}Chain1: {chain1_info.chain}/{chain1_info.resn}{chain1_info.resi}/{chain1_info.name} ({chain1_class}) ")
                        f.write(f"{interaction_type}Chain2: {chain2_info.chain}/{chain2_info.resn}{chain2_info.resi}/{chain2_info.name} ({chain2_class}) ")
                        f.write(f"{interaction_type}Distance: {distance:.2f}\n")

                    f.write("\n")
                    f.write("Number of Hydrogen Bonds:\n")
                    f.write(str(number_hbonds) + "\n")
                    f.write("Number of van der Waals Bonds:\n")
                    f.write(str(number_vdWBonds) + "\n")

                    f.write("\n")
                    f.write("Contact List:\n")
                    contact_list = []
                    for chain1_info, chain2_info, distance, interaction_type, chain1_class, chain2_class, chain1_id, chain2_id in sorted_combined_interactions:
                        contact_list.append(
                            f"{chain1_info.resn}({chain1_info.resi})-{chain2_info.resn}({chain2_info.resi}): {distance:.2f}")
                    contact_list_string = ', '.join(contact_list)
                    f.write(contact_list_string)

                    # if (entity_types[strand_ids.index(polynucleotide_chain)] == 'polyribonucleotide'):
                    if (find_entity_type(polynucleotide_chain, strand_ids, entity_types) == 'polyribonucleotide'):
                        rna_prot_interface_atom_ids = extract_atom_id(sorted_combined_interactions)
                        rna_prot_interface_atom_idsS.append(rna_prot_interface_atom_ids)
                        interface_rna_atom_ids = extract_rna_atom_id(sorted_combined_interactions, polynucleotide_chains)
                        interface_rna_atom_idsS.append(interface_rna_atom_ids)
                    # elif (entity_types[strand_ids.index(polynucleotide_chain)] == 'polydeoxyribonucleotide'):
                    elif (find_entity_type(polynucleotide_chain, strand_ids, entity_types) == 'polydeoxyribonucleotide'):
                        dna_prot_interface_atom_ids = extract_atom_id(sorted_combined_interactions)
                        dna_prot_interface_atom_idsS.append(dna_prot_interface_atom_ids)
                        interface_dna_atom_ids = extract_rna_atom_id(sorted_combined_interactions, polynucleotide_chains)
                        interface_dna_atom_idsS.append(interface_dna_atom_ids)
                    contact_listS.append(contact_list)

        return rna_motif_stringS, length_rna_motif, dna_motif_stringS, aa_motif_stringS, contact_listS, rna_prot_interface_atom_idsS, \
            interface_rna_atom_idsS, dna_prot_interface_atom_idsS, interface_dna_atom_idsS, aac_list, sorted_combined_interactions, \
            is_protein, is_dna, is_rna, chain_id_pairs_rna, h_bonds_rna, vdW_bonds_rna, chain_id_pairs_dna, h_bonds_dna, vdW_bonds_dna

    except Exception as e:
        # Move file to check folder on any other error
        file_name = os.path.basename(file_path)
        dst_path = os.path.join(check_folder, file_name)
        print(f"Error processing {file_name}: {e}")
        print(f"Moving {file_name} to {check_folder}")
        import shutil
        shutil.move(file_path, dst_path)
        return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None

def adjust_motif_sequence(full_sequence_list, motif_sequence):
    # Mapping from full_sequence_list index to its actual value
    sequence_map = {int(value): idx + 1 for idx, value in enumerate(full_sequence_list)}
    adjusted_motif_sequence = []
    for motif in motif_sequence:
        if isinstance(motif, tuple):  # Handling case with nucleotides
            nucleotide, motif_index = motif
            motif_index = int(motif_index)
        else: # In case there are no nucleotides, just an index
            nucleotide = None
            motif_index = int(motif)
        if motif_index in sequence_map:
            # Adjust the motif index based on its position in the full_sequence_list
            adjusted_index = sequence_map[motif_index]
            if nucleotide:
                adjusted_motif_sequence.append((nucleotide, adjusted_index))
            else:
                adjusted_motif_sequence.append(adjusted_index)
    return adjusted_motif_sequence

def extract_rna_motif(sorted_combined_interactions, polynucleotide_chain, file_path):
    file_name = os.path.basename(file_path)
    file = os.path.splitext(file_name)[0]

    nucleotides = {'A', 'U', 'G', 'C'}
    modified_nucleotides = {'6MZ', '5MC', 'OMU', 'A2M', 'OMC', 'OMG', 'UR3'}# Add more as needed
    rna_motif_set = set()
    chain_id = None  # Probe what chain id it is: e.g. A, B, C
    for interaction in sorted_combined_interactions:
        chain1_info, chain2_info, _, interaction_type, chain1_class, chain2_class, _, _ = interaction
        if chain1_info.resn in modified_nucleotides: # Is nucleotide
            rna_motif_set.add((f"({chain1_info.resn})", chain1_info.resi))
        if chain1_info.resn in nucleotides: # Is nucleotide
            rna_motif_set.add((chain1_info.resn, chain1_info.resi))
            if chain_id == None:
                chain_id = chain1_info.chain
        elif chain2_info.resn in modified_nucleotides: # Is nucleotide
            rna_motif_set.add((f"({chain2_info.resn})", chain2_info.resi))
        if chain2_info.resn in nucleotides:
            rna_motif_set.add((chain2_info.resn, chain2_info.resi))
            if chain_id == None:
                chain_id = chain2_info.chain

    length_rna_motif = len(rna_motif_set)
    sorted_rna_motif = sorted(rna_motif_set, key=lambda x: int(x[1]))
    rna_motif_indices = [int(resi) for resn, resi in sorted_rna_motif]

    if len(file) in [4, 6]:
        rna_motif_indices = adjust_motif_sequence(polynucleotide_chain, rna_motif_indices)
        updated_rna_motif_set = set()
        for i, (nucleotide, _) in enumerate(sorted_rna_motif):
            updated_rna_motif_set.add((nucleotide, rna_motif_indices[i]))
        rna_motif_set = updated_rna_motif_set
        sorted_rna_motif = sorted(rna_motif_set, key=lambda x: int(x[1]))
        # print("sorted_rna_motif after", sorted_rna_motif)
        rna_motif_indices = [int(resi) for resn, resi in sorted_rna_motif]

    rna_motif_string = ', '.join([f"{resn}({resi})" for resn, resi in sorted_rna_motif])

    # To separate non-consecutive motifs by "|":
    if not rna_motif_indices:
        rna_motif_visual = ''
    else:
        rna_motif_visual_parts = []
        prev_index = rna_motif_indices[0]
        current_part = [sorted_rna_motif[0][0]]

        for i in range(1, len(rna_motif_indices)):
            current_index = rna_motif_indices[i]
            if current_index == prev_index + 1:
                current_part.append(sorted_rna_motif[i][0])
            else:
                rna_motif_visual_parts.append(''.join(current_part))
                rna_motif_visual_parts.append('|')
                current_part = [sorted_rna_motif[i][0]]
            prev_index = current_index

        # Append the last part
        rna_motif_visual_parts.append(''.join(current_part))
        rna_motif_visual = ''.join(rna_motif_visual_parts)

        if '|' in rna_motif_visual:
            motifs = rna_motif_visual.split('|')
            longest_motif = max(motifs, key=len)
            length_rna_motif = len(longest_motif)

    return rna_motif_string, rna_motif_visual, rna_motif_indices, chain_id, length_rna_motif

def extract_base_to_base_nucleotides(sorted_combined_interactions):
    base_to_base_interactions = []

    for interaction in sorted_combined_interactions:
        chain1_info, chain2_info, distance, interaction_type, chain1_class, chain2_class, chain1_id, chain2_id = interaction
        # Check if both interacting elements are bases
        if chain1_class == "base" and chain2_class == "base":
            base_to_base_interactions.append(
                (chain1_info.resn, chain1_info.resi, chain1_info.chain, chain2_info.resn, chain2_info.resi))
            # base_to_base_interactions.append(
            #     (chain1_info.resn, chain1_info.resi, chain2_info.resn, chain2_info.resi, interaction_type, round(distance, 2)))

    def filter_interactions(interactions):
        # Count occurrences of (chain1_resn, chain1_resi) pairs
        chain1_pair_count = Counter((interaction[0], interaction[1]) for interaction in interactions)
        filtered = []
        added_pairs = set()

        for interaction in interactions:
            chain1_pair = (interaction[0], interaction[1])
            if chain1_pair_count[chain1_pair] >= 2 and chain1_pair not in added_pairs: # Bp with at least 2 bonds
                filtered.append(interaction[:3])  # Only include chain1_resn and chain1_resi
                added_pairs.add(chain1_pair)
        return filtered

    tuples_list = filter_interactions(base_to_base_interactions)
    length_rnaRNA_motif = len(tuples_list)
    rna_rna_string = ', '.join([f"{resn}({chain}{resi})" for resn, resi, chain in tuples_list])
    return rna_rna_string, length_rnaRNA_motif

def extract_dna_motif(sorted_combined_interactions, polynucleotide_chain, file_path):
    file_name = os.path.basename(file_path)
    file = os.path.splitext(file_name)[0]
    def transform_nucleotide(two_letter_code):
        two_to_one_letter_mapping = {"DA": "A","DT": "T","DG": "G","DC": "C"}
        modified_nucleotide = {"3DR"}
        if two_letter_code in two_to_one_letter_mapping:
            return two_to_one_letter_mapping.get(two_letter_code, '')
        # If it's a special case like (3DR), just return it as is
        if two_letter_code in modified_nucleotide:
            return two_letter_code

    chain_id = None  # Probe what chain id it is: e.g. A, B, C
    dna_motif_set = set()
    for interaction in sorted_combined_interactions:
        chain1_info, chain2_info, _, interaction_type, chain1_class, chain2_class, _, _ = interaction
        # Two-letter nucleotides
        if len(chain1_info.resn) == 2:
            dna_motif_set.add((transform_nucleotide(chain1_info.resn), chain1_info.resi))
            if chain_id == None:
                chain_id = chain1_info.chain
        elif len(chain2_info.resn) == 2:
            dna_motif_set.add((transform_nucleotide(chain2_info.resn), chain2_info.resi))
            if chain_id == None:
                chain_id = chain2_info.chain

    sorted_dna_motif = sorted(dna_motif_set, key=lambda x: int(x[1]))
    dna_motif_indices = [int(resi) for resn, resi in sorted_dna_motif]

    if len(file) in [4, 6] and not all(index in polynucleotide_chain for index in map(str, dna_motif_indices)):
        # Filter out invalid indices from dna_motif_indices and update sorted_dna_motif accordingly
        # e.g. 701 in 5xpa.cif, which represents the Mg2+ instead of polynucleotide_chain
        filtered_motif_data = [(nucleotide, idx) for (nucleotide, idx) in sorted_dna_motif if
                               idx in polynucleotide_chain]
        # Extract the filtered indices and nucleotides separately
        dna_motif_indices = [int(idx) for _, idx in filtered_motif_data]
        sorted_dna_motif = [(nucleotide, idx) for nucleotide, idx in filtered_motif_data]
        updated_dna_motif_set = set()
        for i, (nucleotide, _) in enumerate(sorted_dna_motif):
            updated_dna_motif_set.add((nucleotide, dna_motif_indices[i]))
        dna_motif_set = updated_dna_motif_set
        sorted_dna_motif = sorted(dna_motif_set, key=lambda x: int(x[1]))
        dna_motif_indices = [int(resi) for resn, resi in sorted_dna_motif]

    first_chain_index = int(polynucleotide_chain[0])
    if first_chain_index != 1:
        offset = first_chain_index - 1
        # Adjust dna_motif_indices by subtracting the same offset
        dna_motif_indices = [int(resi) - offset for resi in dna_motif_indices]
        updated_dna_motif_set = set()
        for i, (nucleotide, _) in enumerate(sorted_dna_motif):
            updated_dna_motif_set.add((nucleotide, dna_motif_indices[i]))
        dna_motif_set = updated_dna_motif_set
        sorted_dna_motif = sorted(dna_motif_set, key=lambda x: int(x[1]))

    dna_motif_string = ', '.join([f"{resn}({resi})" for resn, resi in sorted_dna_motif])

    # To separate non-consecutive motifs by "|" use:
    if not dna_motif_indices:
        dna_motif_visual = ''
    else:
        dna_motif_visual_parts = []
        prev_index = dna_motif_indices[0]
        current_part = [sorted_dna_motif[0][0]]

        for i in range(1, len(dna_motif_indices)):
            current_index = dna_motif_indices[i]
            if current_index == prev_index + 1:
                current_part.append(sorted_dna_motif[i][0])
            else:
                dna_motif_visual_parts.append(''.join(current_part))
                dna_motif_visual_parts.append('|')
                current_part = [sorted_dna_motif[i][0]]
            prev_index = current_index

        # Append the last part
        dna_motif_visual_parts.append(''.join(current_part))
        dna_motif_visual = ''.join(dna_motif_visual_parts)

    return dna_motif_string, dna_motif_visual, dna_motif_indices, chain_id

def extract_aa_motif(sorted_combined_interactions):
    allowed_residues = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
        'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V',
        'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'NEP': 'H', 'HIP': 'H', 'ALY': 'K', 'MLY': 'K', 'M3L': 'K', 'MLZ': 'K',
        '2MR': 'R', 'AGM': 'R', 'MCS': 'C', 'HYP': 'P', 'HY3': 'H', 'LYZ': 'K', 'AHB': 'A', 'P1L': 'P', 'SNN': 'S',
        'SNC': 'C', 'TRF': 'W', 'KCR': 'K', 'CIR': 'R', 'YHA': 'Y'
    }
    aa_motif_set = set()
    for interaction in sorted_combined_interactions:
        chain1_info, chain2_info, _, interaction_type, chain1_class, chain2_class, _, _ = interaction

        if len(chain1_info.resn) == 3 and chain1_info.resn in allowed_residues:
            # print("chain1_info.resn", chain1_info.resn, chain1_info.resi)
            aa_motif_set.add((chain1_info.resn, chain1_info.resi))
        elif len(chain2_info.resn) == 3 and chain2_info.resn in allowed_residues:
            # print("chain2_info.resn", chain2_info.resn, chain2_info.resi)
            aa_motif_set.add((chain2_info.resn, chain2_info.resi))

    sorted_aa_motif = sorted(aa_motif_set, key=lambda x: int(x[1]))
    aa_motif_string = ', '.join([f"{resn}({resi})" for resn, resi in sorted_aa_motif])

    return aa_motif_string


def calculate_aac(amino_acid_motif):
    amino_acids = [
        "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"
    ]
    motif_counts = {aa: 0 for aa in amino_acids}
    amino_acid_motif = re.findall(r'[A-Z]{3}', amino_acid_motif)
    total_motif_amino_acids = len(amino_acid_motif)
    for aa in amino_acid_motif: # Count the occurrences of each amino acid in the motif
        if aa in motif_counts:
            motif_counts[aa] += 1

    aac_dict = {aa: round((count / total_motif_amino_acids), 2) for aa, count in motif_counts.items() if count > 0}
    aac_list = [f"{key}: {value}" for key, value in aac_dict.items()]
    aac = f"[{', '.join(aac_list)}]"

    return aac

def extract_atom_id(sorted_combined_interactions):
    atom_ids = set()
    for interaction in sorted_combined_interactions:
        chain1_info, chain2_info, _, _, _, _, chain1_id, chain2_id = interaction
        atom_ids.add(chain1_id)
        atom_ids.add(chain2_id)
    sorted_atom_ids = sorted(atom_ids)
    return sorted_atom_ids

def extract_rna_atom_id(sorted_combined_interactions, polynucleotide_chains):
    atom_ids = set()
    for interaction in sorted_combined_interactions:
        chain1_info, chain2_info, _, _, _, _, chain1_id, chain2_id = interaction
        # Check if the chain IDs belong to polynucleotide chains
        if chain1_info.chain in polynucleotide_chains:
            atom_ids.add(chain1_id)
        if chain2_info.chain in polynucleotide_chains:
            atom_ids.add(chain2_id)
    sorted_rna_atom_ids = sorted(atom_ids)
    return sorted_rna_atom_ids

def get_full_rna_string_from_db(file_path, is_rna):
    def query_full_sequence(file, file_name, is_rna):
        if len(file) in [4, 6] and is_rna:
            # RNA instance with PDBId
            rna = RNA.get_rna_from_db(id=file.upper())
            sequence = rna.get_rna_sequence()
            chain = rna.get_rna_chain_IDs()
        elif len(file) > 6 and is_rna:
            # RNA instance with FileName
            parts = file_name.split('_')
            rna = RNA.get_rna_from_db(id=parts[1].upper(), file_name=file_name)
            sequence = rna.get_rna_sequence()
            chain = rna.get_rna_chain_IDs()
        if len(file) in [4, 6] and not is_rna:
            dna = DNA.get_dna_from_db(id=file.upper())
            sequence = dna.get_dna_sequence()
            chain = dna.get_dna_chain_IDs()
        elif len(file) > 6 and not is_rna:
            parts = file_name.split('_')
            dna = DNA.get_dna_from_db(id=parts[1].upper(), file_name=file_name)
            sequence = dna.get_dna_sequence()
            chain = dna.get_dna_chain_IDs()

        return sequence, chain

    file_name = os.path.basename(file_path)
    file = os.path.splitext(file_name)[0]
    # Query the RNA sequence
    full_rna_sequence = query_full_sequence(file, file_name, is_rna)

    if not full_rna_sequence: # TO-DO: use cifParser.py
        print("Parse .cif file first before extracting interactions.")

    return full_rna_sequence

def get_longest_chain_id_from_db(file_path, is_rna):
    file_name = os.path.basename(file_path)
    file = os.path.splitext(file_name)[0]
    if len(file) in [4, 6] and is_rna:
        rna = RNA.get_rna_from_db(id=file.upper())
        chain_id = rna.get_longest_chain_id()
    elif len(file) > 6 and is_rna:
        parts = file_name.split('_')
        rna = RNA.get_rna_from_db(id=parts[1].upper(), file_name=file_name)
        chain_id = rna.get_longest_chain_id()
    if len(file) in [4, 6] and not is_rna:
        dna = DNA.get_dna_from_db(id=file.upper())
        chain_id = dna.get_longest_chain_id()
    elif len(file) > 6 and not is_rna:
        parts = file_name.split('_')
        dna = DNA.get_dna_from_db(id=parts[1].upper(), file_name=file_name)
        chain_id = dna.get_longest_chain_id()
        # length = dna.get_rna_length()
        # chain = dna.get_dna_chain_IDs()

    return chain_id

def highlight_sequence(full_rna_string, indices):
    def handle_modified_nucleotides(sequence):
        # Identify modified nucleotides like (5MC)
        return re.findall(r'\([A-Za-z0-9]+\)|[A-Za-z]', sequence)

    rna_nucleotides = handle_modified_nucleotides(full_rna_string)

    highlighted_sequence = ""
    index_set = set(indices)
    nucleotide_index = 1  # Start index for nucleotides

    for nucleotide in rna_nucleotides:
        # Check if the current nucleotide's index is in the indices to be highlighted
        if nucleotide_index in index_set:
            highlighted_sequence += f"[{nucleotide}]"
        else:
            highlighted_sequence += nucleotide
        nucleotide_index += 1

    return highlighted_sequence

def transform_dna_sequence(two_letter_sequence):
    two_to_one_letter_mapping = {"DA": "A", "DT": "T", "DG": "G", "DC": "C"}
    modified_nucleotides = {"(3DR)"}
    # print(two_letter_sequence)
    pattern = r"\(\w+\)|\w{2}"
    matches = re.findall(pattern, two_letter_sequence)
    one_letter_sequence = ""
    for match in matches:
        if match in two_to_one_letter_mapping:
            one_letter_sequence += two_to_one_letter_mapping[match]
        elif match in modified_nucleotides:
            one_letter_sequence += match
        else:
            raise ValueError(f"Invalid code: {match}")

    return one_letter_sequence

def update_exp_pdb_with_motifs(file_path, rna_motif, length_rna_motif, dna_motif, aa_motif, contact_list, rna_prot_interface_atom_ids,
                               interface_rna_atom_idsS, dna_prot_interface_atom_ids, interface_dna_atom_ids, aac_list, is_protein,
                               is_dna, is_rna, chain_id_pairs_rna, h_bonds_rna, vdW_bonds_rna, chain_id_pairs_dna, h_bonds_dna, vdW_bonds_dna):
    pdb_id = os.path.basename(file_path).split('.')[0][:4].upper()
    base_name = os.path.basename(file_path)
    file_name = base_name[:-4]
    data_motif = rna_motif, length_rna_motif, dna_motif, aa_motif, aac_list, contact_list
    data_interface = rna_prot_interface_atom_ids, interface_rna_atom_idsS, dna_prot_interface_atom_ids, interface_dna_atom_ids

    def handle_list_or_value(item):
        if isinstance(item, list):
            if len(item) == 1:
                return item[0]
            return json.dumps(item)
        return item

    motif_values = [handle_list_or_value(data_motif[i]) for i in range(len(data_motif))]
    interface_values = [handle_list_or_value(data_interface[i]) for i in range(len(data_interface))]
    # list -> string:
    motif_values = tuple(
        json.dumps(motif_values) if isinstance(motif_values, list) else motif_values for motif_values in
        motif_values)
    interface_values = tuple(
        json.dumps(interface_values) if isinstance(interface_values, list) else interface_values for interface_values in
        interface_values)
    rna_motif, length_rna_motif, dna_motif, aa_motif, aac_list, contact_list = motif_values
    rna_prot_interface_atom_ids, interface_rna_atom_ids, dna_prot_interface_atom_ids, interface_dna_atom_ids = interface_values

    def handle(value):
        return value if not isinstance(value, list) else json.dumps(value)
    bond_values = chain_id_pairs_rna, h_bonds_rna, vdW_bonds_rna, chain_id_pairs_dna, h_bonds_dna, vdW_bonds_dna
    bond_values = [handle(bond_values[i]) for i in range(len(bond_values))]
    bond_values = tuple(bond_values)
    # for i, bond_value in enumerate(bond_values, start=1):
    #     print(f"Parameter {i}:(Type: {type(bond_value)})")
    chain_id_pairs_rna, h_bonds_rna, vdW_bonds_rna, chain_id_pairs_dna, h_bonds_dna, vdW_bonds_dna = bond_values

    if len(file_name) > 6:
        if is_protein and is_rna and not is_dna:
            database_methods.update_or_insert('pred_protein_rna',
                                              ['RNAMotif', 'RNAMotifLength', 'AAMotif', 'AAC', 'ContactList',
                                               'ChainIDpairList_proteinRNA', 'Hbond_proteinRNA', 'vdWbond_proteinRNA',
                                               'rna_prot_interface_atom_ids', 'interface_rna_atom_ids'],
                                              (rna_motif, length_rna_motif, aa_motif, aac_list, contact_list,
                                               chain_id_pairs_rna, h_bonds_rna, vdW_bonds_rna,
                                               rna_prot_interface_atom_ids, interface_rna_atom_ids),
                                              condition=f"FileName = '{base_name}'")

        if not is_protein and is_rna and not is_dna:
            database_methods.update_or_insert('pred_rna_rna',
                                              ['RNAMotif', 'RNAMotifLength', 'ContactList',
                                               'ChainIDpairList_proteinRNA', 'Hbond_proteinRNA', 'vdWbond_proteinRNA',
                                               'rna_rna_interface_atom_ids', 'interface_rna_atom_ids'],
                                              (rna_motif, length_rna_motif, contact_list, chain_id_pairs_rna,
                                               h_bonds_rna, vdW_bonds_rna, rna_prot_interface_atom_ids,
                                               interface_rna_atom_ids),
                                              condition=f"FileName = '{base_name}'")

        if is_protein and not is_rna and is_dna:
            database_methods.update_or_insert('pred_protein_dna',
                                              ['DNAMotif', 'AAMotif', 'AAC', 'ContactList',
                                               'ChainIDpairList_proteinDNA', 'Hbond_proteinDNA', 'vdWbond_proteinDNA',
                                               'dna_prot_interface_atom_ids', 'interface_dna_atom_ids'],
                                              (dna_motif, aa_motif, aac_list, contact_list,
                                               chain_id_pairs_dna, h_bonds_dna, vdW_bonds_dna,
                                               dna_prot_interface_atom_ids, interface_dna_atom_ids),
                                              condition=f"FileName = '{base_name}'")

        if is_protein and is_rna and is_dna:
            database_methods.update_or_insert('pred_protein_rna_dna',
                                              ['RNAMotif', 'RNAMotifLength', 'DNAMotif', 'AAMotif', 'AAC', 'ContactList',
                                               'ChainIDpairList_proteinRNA', 'Hbond_proteinRNA', 'vdWbond_proteinRNA',
                                               'ChainIDpairList_proteinDNA', 'Hbond_proteinDNA', 'vdWbond_proteinDNA',
                                               'rna_prot_interface_atom_ids', 'interface_rna_atom_ids',
                                               'dna_prot_interface_atom_ids', 'interface_dna_atom_ids'],
                                              (rna_motif, length_rna_motif, dna_motif, aa_motif, aac_list, contact_list,
                                               chain_id_pairs_rna, h_bonds_rna, vdW_bonds_rna, chain_id_pairs_dna,
                                               h_bonds_dna, vdW_bonds_dna, rna_prot_interface_atom_ids,
                                               interface_rna_atom_ids, dna_prot_interface_atom_ids, interface_dna_atom_ids),
                                              condition=f"FileName = '{base_name}'")

    else:
        if is_protein and is_rna and not is_dna:
            database_methods.update_or_insert('exp_protein_rna',
                                              ['RNAMotif', 'RNAMotifLength', 'AAMotif', 'AAC', 'ContactList',
                                               'ChainIDpairList_proteinRNA', 'Hbond_proteinRNA', 'vdWbond_proteinRNA'],
                                              (rna_motif, length_rna_motif, aa_motif, aac_list, contact_list,
                                               chain_id_pairs_rna, h_bonds_rna, vdW_bonds_rna),
                                              condition=f"PDBId = '{pdb_id}'")

        if not is_protein and is_rna and not is_dna:
            database_methods.update_or_insert('exp_rna_rna',
                                              ['RNAMotif', 'RNAMotifLength', 'ContactList',
                                               'ChainIDpairList_proteinRNA', 'Hbond_proteinRNA', 'vdWbond_proteinRNA'],
                                              (rna_motif, length_rna_motif, contact_list, chain_id_pairs_rna,
                                               h_bonds_rna, vdW_bonds_rna),
                                              condition=f"PDBId = '{pdb_id}'")

        if is_protein and not is_rna and is_dna:
            database_methods.update_or_insert('exp_protein_dna',
                                              ['DNAMotif', 'AAMotif', 'AAC', 'ContactList',
                                               'ChainIDpairList_proteinDNA', 'Hbond_proteinDNA', 'vdWbond_proteinDNA'],
                                              (dna_motif, aa_motif, aac_list, contact_list,
                                               chain_id_pairs_dna, h_bonds_dna, vdW_bonds_dna),
                                              condition=f"PDBId = '{pdb_id}'")

        if is_protein and is_rna and is_dna:
            database_methods.update_or_insert('pred_protein_rna_dna',
                                              ['RNAMotif', 'RNAMotifLength', 'DNAMotif', 'AAMotif', 'AAC',
                                               'ContactList',
                                               'ChainIDpairList_proteinRNA', 'Hbond_proteinRNA', 'vdWbond_proteinRNA',
                                               'ChainIDpairList_proteinDNA', 'Hbond_proteinDNA', 'vdWbond_proteinDNA'],
                                              (rna_motif, length_rna_motif, dna_motif, aa_motif, aac_list, contact_list,
                                               chain_id_pairs_rna, h_bonds_rna, vdW_bonds_rna, chain_id_pairs_dna,
                                               h_bonds_dna, vdW_bonds_dna),
                                              condition=f"PDBId = '{pdb_id}'")

if __name__ == "__main__":
    # Initialize global list for skipped PDBs
    skipped_pdbs = []

    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <directory_path>")
        sys.exit(1)
    directory_path = sys.argv[1]
    if not os.path.isdir(directory_path):
        print(f"Error: {directory_path} is not a valid directory.")
        sys.exit(1)

    cif_files = glob.glob(os.path.join(directory_path, '*.cif'))
    for file_path in cif_files:
        file_name = os.path.basename(file_path)
        file_base = os.path.splitext(file_name)[0]
        output_file = file_base + "_interactions.txt"
        cmd.load(file_path, 'structure') # Load the structure into PyMOL
        rna_motif_string, length_rna_motif, dna_motif_string, aa_motif_string, contact_list, rna_prot_interface_atom_ids, \
        interface_rna_atom_ids, interface_dna_atom_ids, dna_prot_interface_atom_ids, aac_list, sorted_combined_interactions, \
        is_protein, is_dna, is_rna, chain_id_pairs_rna, h_bonds_rna, vdW_bonds_rna, chain_id_pairs_dna, h_bonds_dna, vdW_bonds_dna = list_interactions_from_cif(file_path, output_file)
        # print(rna_motif_string, aa_motif_string, contact_list, interface_atom_ids, interface_rna_atom_ids)
        update_exp_pdb_with_motifs(file_path, rna_motif_string, length_rna_motif, dna_motif_string, aa_motif_string, contact_list,
                                   rna_prot_interface_atom_ids, interface_rna_atom_ids, interface_dna_atom_ids,
                                   dna_prot_interface_atom_ids, aac_list, is_protein, is_dna, is_rna, chain_id_pairs_rna,
                                   h_bonds_rna, vdW_bonds_rna, chain_id_pairs_dna, h_bonds_dna, vdW_bonds_dna)

        cmd.delete('structure') # Delete loaded structure to free memory (optional)

    # Print skipped PDBs at the end
    if skipped_pdbs:
        print("\nThe following PDB IDs were skipped due to non-numeric residue IDs:")
        print(", ".join(skipped_pdbs))
        print(f"Total skipped: {len(skipped_pdbs)}")
    else:
        print("\nNo PDB files were skipped.")

    database_methods.close_connection()
    cmd.quit()