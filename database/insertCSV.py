import sqlite3
import csv
import json

conn = sqlite3.connect('rbpDatabase.db')
cursor = conn.cursor()

# Open and read the CSV file
with open('ValidatedRNA1.csv', 'r', newline='', encoding='utf-8') as file:
    reader = csv.DictReader(file, delimiter=';')

    # # Print fieldnames to verify headers
    # print(reader.fieldnames)
    # # Print a sample row to verify data structure
    # sample_row = next(reader)
    # print(sample_row)
    # # Reset file pointer to beginning for actual data insertion
    # file.seek(0)
    # # Skip the header row
    # next(file)

    csv_data = []

    for row in reader:
        # Combine 'Comment' and 'Paper' fields
        comment = row.get('Comment', '')
        paper = row.get('Paper', '')
        combined_list = [comment, paper]
        combined_string = json.dumps(combined_list)
        # print(isinstance(combined_string, str), combined_string)

        csv_data.append(
            ('0', row['UniprotId'], row['\ufeffProteinName'], '', '', '', row['Organism'], row['ProteinSequence'], '',
             row['ProteinLength'], '1', '', row['RNAName'], '', '', '0', row['RNASequence'], '1', row['RNALength'],
             '1', row['RNAMotif'], '', '', row['Experiment'], '', combined_string))

cursor.executemany('''
    INSERT INTO exp_protein_rna (PDBId, UniprotId, ProteinName, ProteinDescription, ProteinSource, 
                              ProteinFamily, ProteinOrganism, ProteinSequence, ProteinChainIDs, ProteinLength, NumberProteins, AAMotif, 
                              RNADescription, RNASource, RNAFamily, RNAModified, RNASequence, RNAChainIDs, RNALength, NumberRNAs, 
                              RNAMotif, ContactList, Deposited_date, Experiment, XRayResolution, Comment)
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', csv_data)

conn.commit()
conn.close()