import os
import sqlite3
import json
import ast
from database.protein_family_request import get_protein_family_from_sequence
from database.protein_family_request import extract_interpro_entries
from structures.protein import Protein
import startConfig

db_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'database', 'rbpDatabase.db')

def get_pdb_ids(table_name):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    query = f"SELECT PDBId FROM {table_name}"
    cursor.execute(query)
    pdb_ids = cursor.fetchall()
    conn.close()
    return [pdb_id[0] for pdb_id in pdb_ids]

def get_uniprot_ids(table_name):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    query = f"SELECT UniprotId FROM {table_name}"
    cursor.execute(query)
    pdb_ids = cursor.fetchall()
    conn.close()
    return [pdb_id[0] for pdb_id in pdb_ids]

# def insert_family_into_db(pdb_ids):
#     conn = sqlite3.connect(db_path)
#     cursor = conn.cursor()
#
#     def handle_list_or_value(item):
#         if isinstance(item, list):
#             if len(set(item)) == 1:
#                 return str(item[0])
#             # Return JSON-formatted string if elements are not identical
#             return json.dumps(item)
#         return str(item)
#
#     for pdb_id in pdb_ids:
#         protein = Protein.get_protein_from_db(pdb_id)
#         sequence = protein.get_protein_sequence()
#         if sequence.startswith("["):
#             sequence_list = ast.literal_eval(sequence)
#             sequence = max(sequence_list, key=len)
#         print(pdb_id, sequence)
#         protein_families = get_protein_family_from_sequence(sequence, pdb_id)
#         print("protein_families", protein_families)
#         if protein_families:
#             protein_familyH = [handle_list_or_value(protein_families[i]) for i in range(len(protein_families))]
#             protein_family = json.dumps(protein_familyH)
#         else:
#             protein_family = 'None'
#         if table_name == 'exp_protein_rna':
#             cursor.execute('''UPDATE exp_protein_rna
#                                   SET ProteinFamily = ?
#                                   WHERE PDBId = ?''', (protein_family, pdb_id))
#             print("success", protein_family, pdb_id)
#             conn.commit()
#         if table_name == 'exp_rna_rna':
#             cursor.execute('''UPDATE exp_rna_rna
#                                   SET ProteinFamily = ?
#                                   WHERE PDBId = ?''', (protein_family, pdb_id))
#             conn.commit()
#         if table_name == 'exp_protein_dna':
#             cursor.execute('''UPDATE exp_protein_dna
#                                   SET ProteinFamily = ?
#                                   WHERE PDBId = ?''', (protein_family, pdb_id))
#             conn.commit()
#         if table_name == 'exp_protein_rna_dna':
#             cursor.execute('''UPDATE exp_protein_rna_dna
#                                   SET ProteinFamily = ?
#                                   WHERE PDBId = ?''', (protein_family, pdb_id))
#             conn.commit()
#
#     conn.close()

def insert_family_into_db(pdb_ids):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    def handle_list_or_value(item):
        if isinstance(item, list):
            if len(set(item)) == 1:
                return str(item[0])
            # Return JSON-formatted string if elements are not identical
            return json.dumps(item)
        return str(item)

    for pdb_id in pdb_ids:
        protein = Protein.get_protein_from_db(pdb_id)
        sequence = protein.get_protein_sequence()
        if sequence.startswith("["):
            sequence_list = ast.literal_eval(sequence)
            sequence = max(sequence_list, key=len)
        print(pdb_id, sequence)
        protein_families = extract_interpro_entries(pdb_id)
        print("protein_families", protein_families)
        if protein_families:
            protein_familyH = [handle_list_or_value(protein_families[i]) for i in range(len(protein_families))]
            protein_family = json.dumps(protein_familyH)
        else:
            protein_family = 'None'
        if table_name == 'exp_protein_rna':
            cursor.execute('''UPDATE exp_protein_rna
                                  SET ProteinFamily = ?
                                  WHERE PDBId = ?''', (protein_family, pdb_id))
            print("success", protein_family, pdb_id)
            conn.commit()
        if table_name == 'exp_rna_rna':
            cursor.execute('''UPDATE exp_rna_rna
                                  SET ProteinFamily = ?
                                  WHERE PDBId = ?''', (protein_family, pdb_id))
            conn.commit()
        if table_name == 'exp_protein_dna':
            cursor.execute('''UPDATE exp_protein_dna
                                  SET ProteinFamily = ?
                                  WHERE PDBId = ?''', (protein_family, pdb_id))
            conn.commit()
        if table_name == 'exp_protein_rna_dna':
            cursor.execute('''UPDATE exp_protein_rna_dna
                                  SET ProteinFamily = ?
                                  WHERE PDBId = ?''', (protein_family, pdb_id))
            conn.commit()

    conn.close()

if __name__ == "__main__":
    table_name = 'exp_protein_rna'
    pdb_ids = get_pdb_ids(table_name)
    uniprot_ids = get_uniprot_ids(table_name)
    print("pdb_ids", pdb_ids)
    print("uniprot_ids", uniprot_ids)
    # input_string = "'7OFW', '7VSJ', '7PMQ', '7VKL', '7S02', '7SOV', '7D8O', '7UJ1', '7TUV', '8H0S', '7F36', '7M50', '7PMM', '7ENI', '7OZQ', '7E8O', '7M2V', '7PDV', '7KI3', '7W9S', '7M3T', '7DWH', '7WNU', '7RZZ', '7ENR', '8IDF', '8T2B'"
    # pdb_ids = [x.strip().strip("'") for x in input_string.split(",")]
    print("uniprot_ids", uniprot_ids)
    # insert_family_into_db(uniprot_ids)
    def handle_list_or_value(value):
        """Handles JSON lists or single values."""
        if isinstance(value, list):
            return value
        try:
            return json.loads(value)
        except (ValueError, TypeError):
            return value

    db_file = db_path
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Retrieve PDBId and UniprotId from the database
    cursor.execute(f"SELECT PDBId, UniprotId, ProteinFamily FROM {table_name}")
    rows = cursor.fetchall()

    for row in rows:
        pdb_id = row[0]
        uniprot_id_raw = row[1]
        protein_family_value = row[2]

        # Skip if ProteinFamily already has a value (i.e., not an empty string or whitespace)
        if protein_family_value and protein_family_value.strip():  # Check if it's not empty or just spaces
            print(f"Skipping PDBId {pdb_id} because ProteinFamily already has a value {protein_family_value}.")
            continue

        # Parse UniprotId (could be JSON array or a single string)
        uniprot_ids = handle_list_or_value(uniprot_id_raw)
        print("uniprot_ids", uniprot_ids)

        if isinstance(uniprot_ids, list): # CASE 1: more than 1 uniprot_ids
            # Check if any of the Uniprot IDs match the PDBId
            if pdb_id in uniprot_ids:
                # Case: At least one Uniprot ID equals PDBId
                if all(uid == pdb_id for uid in uniprot_ids): # CASE 1a: more than 1 uniprot_ids & all are pdb_id
                    protein = Protein.get_protein_from_db(pdb_id)
                    sequence = protein.get_protein_sequence()
                    if sequence.startswith("["):
                        sequence_list = ast.literal_eval(sequence)
                        sequence = max(sequence_list, key=len)
                    print(pdb_id, sequence)
                    protein_families = get_protein_family_from_sequence(sequence, pdb_id)
                    protein_familyH = [handle_list_or_value(protein_families[i]) for i in range(len(protein_families))]
                    protein_family = json.dumps(protein_familyH)
                    # Update the database
                    cursor.execute(f'''UPDATE {table_name}
                                    SET ProteinFamily = ?
                                    WHERE PDBId = ?''', (protein_family, pdb_id))
                    conn.commit()
                else: # CASE 1a: more than 1 uniprot_ids & not all are pdb_id
                    # Case: Not all Uniprot IDs match PDBId
                    # interpro_entries = extract_interpro_entries(uniprot_ids)
                    all_protein_families = []
                    # Store all protein_families_per_entry for comparison
                    protein_families_collection = []
                    for uid in uniprot_ids:
                        if uid != pdb_id:
                            # Extract InterPro entries for each Uniprot ID
                            protein_families_per_entry = extract_interpro_entries(uid)
                            protein_families_collection.append(protein_families_per_entry)
                    seen_families = []
                    unique_protein_families = []

                    # Iterate over all entries and keep only unique ones
                    for entry in all_protein_families:
                        if entry["FamilyName"] not in seen_families:
                            unique_protein_families.append(entry)
                            seen_families.append(entry["FamilyName"])  # Add to seen_families list

                    if not unique_protein_families:
                        protein = Protein.get_protein_from_db(pdb_id)
                        sequence = protein.get_protein_sequence()
                        if sequence.startswith("["):
                            sequence_list = ast.literal_eval(sequence)
                            sequence = max(sequence_list, key=len)
                        protein_family = get_protein_family_from_sequence(sequence, pdb_id)
                        if protein_family:
                            unique_protein_families = [handle_list_or_value(protein_family[i]) for i in
                                                range(len(protein_family))]
                        else:
                            unique_protein_families = None
                    # Convert unique_protein_families to JSON format
                    protein_family_json = json.dumps(unique_protein_families)
                    cursor.execute(f'''UPDATE {table_name}
                                       SET ProteinFamily = ?
                                       WHERE PDBId = ?''', (protein_family_json, pdb_id))
                    conn.commit()
            else:
                # Case: None of the Uniprot IDs match PDBId
                # Initialize an empty list to collect unique protein families
                all_protein_families = []
                # Store all protein_families_per_entry for comparison
                protein_families_collection = []
                for uid in uniprot_ids:
                    # Extract InterPro entries for each Uniprot ID
                    protein_families_per_entry = extract_interpro_entries(uid)
                    protein_families_collection.append(protein_families_per_entry)
                seen_families = []
                unique_protein_families = []

                # Iterate over all entries and keep only unique ones
                for entry in all_protein_families:
                    if entry["FamilyName"] not in seen_families:
                        unique_protein_families.append(entry)
                        seen_families.append(entry["FamilyName"])  # Add to seen_families list

                if not unique_protein_families:
                    protein = Protein.get_protein_from_db(pdb_id)
                    sequence = protein.get_protein_sequence()
                    if sequence.startswith("["):
                        sequence_list = ast.literal_eval(sequence)
                        sequence = max(sequence_list, key=len)
                    protein_family = get_protein_family_from_sequence(sequence, pdb_id)
                    if protein_family:
                        unique_protein_families = [handle_list_or_value(protein_family[i]) for i in
                                            range(len(protein_family))]
                    else:
                        unique_protein_families = None

                # Convert unique_protein_families to JSON format
                protein_family_json = json.dumps(unique_protein_families)
                # Update the database
                cursor.execute(f'''UPDATE {table_name}
                                   SET ProteinFamily = ?
                                   WHERE PDBId = ?''', (protein_family_json, pdb_id))
                conn.commit()

        else:
            # Case: Single Uniprot ID
            if pdb_id == uniprot_ids:
                protein = Protein.get_protein_from_db(pdb_id)
                sequence = protein.get_protein_sequence()
                if sequence.startswith("["):
                    sequence_list = ast.literal_eval(sequence)
                    sequence = max(sequence_list, key=len)
                print(pdb_id, sequence)
                protein_families = get_protein_family_from_sequence(sequence, pdb_id)
                if protein_families:
                    protein_familyH = [handle_list_or_value(protein_families[i]) for i in range(len(protein_families))]
                    protein_family = json.dumps(protein_familyH)
                else:
                    protein_family = "None"
                # Update the database
                cursor.execute(f'''UPDATE {table_name}
                                                      SET ProteinFamily = ?
                                                      WHERE PDBId = ?''', (protein_family, pdb_id))
                conn.commit()

            else:
                interpro_entries = extract_interpro_entries(uniprot_ids)
                if not interpro_entries:
                    protein = Protein.get_protein_from_db(pdb_id)
                    sequence = protein.get_protein_sequence()
                    if sequence.startswith("["):
                        sequence_list = ast.literal_eval(sequence)
                        sequence = max(sequence_list, key=len)
                    protein_family = get_protein_family_from_sequence(sequence, pdb_id)
                    if protein_family:
                        protein_families = [handle_list_or_value(protein_family[i]) for i in range(len(protein_family))]
                    else:
                        protein_families = None
                else:
                    protein_families = [entry['FamilyName'] for entry in interpro_entries]
                protein_family_json = json.dumps(protein_families)
                cursor.execute(f'''UPDATE {table_name}
                                   SET ProteinFamily = ?
                                   WHERE PDBId = ?''', (protein_family_json, pdb_id))
                conn.commit()

    rna_binding_domains = ["CSD", "Cold_shock", "RRM", "La", "PAZ", "Piwi", "KH", "RBD", "dsRBD", "DEAD", "DEAH", "PUM",
                           "SAM", "S1", "R3H", "PUA", "PUM", "THUMP", "YTH", "Zn_finger", "Znf", "GRP", "RNA_binding",
                           "disordered", "pentatricopeptide"]


    conn.close()
