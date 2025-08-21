import sqlite3
import os

class rbpDatabase:
    def __init__(self, 'rbpDatabase.db': str):
        self.connection = sqlite3.connect('rbpDatabase.db')
        self.cursor = self.connection.cursor()
        self.initialize_db()

    def initialize_db(self):
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS exp (
            Id INTEGER PRIMARY KEY AUTOINCREMENT,
            ProteinName TEXT,
            RNAName TEXT,
            UniprotId TEXT,
            PDBId TEXT,
            Organism TEXT,
            ProteinSequence TEXT,
            RNASequence TEXT,
            ProteinLength INTEGER,
            RNALength INTEGER,
            RNAMotif TEXT,
            AAMotif TEXT,
            Experiment TEXT,
            Comment TEXT,
            Distance FLOAT,
            ContactList TEXT
        );
        ''')

        cursor.execute('''
        CREATE TABLE IF NOT EXISTS pred (
            Id INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_db_id INTEGER,
            FileName TEXT,
            ProteinSequence TEXT,
            RNASequence TEXT,
            ProteinLength INTEGER,
            RNALength INTEGER,
            Seed INTEGER,
            Model INTEGER,
            RNAMotif TEXT,
            AAMotif TEXT,
            Distance FLOAT,
            Comment TEXT,
            ContactList TEXT,
            FOREIGN KEY (exp_db_id) REFERENCES exp(UniprotId),
            FOREIGN KEY (exp_db_id) REFERENCES exp(PDBId)
        );
        ''')
        self.connection.commit()

    def insert_data(self, table_name: str, data: dict): # keys are column names and the values are the data to be inserted
        placeholders = ', '.join(['?' for _ in data])
        columns = ', '.join(data.keys())
        sql = f'INSERT INTO {table_name} ({columns}) VALUES ({placeholders})'
        self.cursor.execute(sql, list(data.values()))
        self.connection.commit()

    def get_data(self, table_name: str, conditions: dict = None):
        if conditions:
            condition_string = ' AND '.join([f"{k} = ?" for k in conditions])
            sql = f'SELECT * FROM {table_name} WHERE {condition_string}'
            self.cursor.execute(sql, list(conditions.values()))
        else:
            sql = f'SELECT * FROM {table_name}'
            self.cursor.execute(sql)
        return self.cursor.fetchall()

    def update_motifs(self, identifier, rna_motif, aa_motif, distance, contacts, is_pred=False):
        rna_motif_string = ', '.join(rna_motif)
        aa_motif_string = ', '.join(aa_motif)
        formatted_contacts = [
            f"{contact[2]}({contact[1]})-{contact[5]}({contact[4]}): {round(contact[6], 2)}"
            for contact in contacts
        ]
        contact_list = ', '.join(formatted_contacts)

        if is_pred:
            self.cursor.execute('''
                UPDATE pred
                SET RNAMotif = ?, AAMotif = ?, Distance = ?, ContactList = ?
                WHERE FileName = ?
            ''', (rna_motif_string, aa_motif_string, distance, contact_list, identifier))
        else:
            self.cursor.execute('''
                UPDATE exp
                SET RNAMotif = ?, AAMotif = ?, Distance = ?, ContactList = ?
                WHERE PDBId = ?
            ''', (rna_motif_string, aa_motif_string, distance, contact_list, identifier))

        self.connection.commit()

    def close(self):
        self.connection.close()

# db = rbpDatabase('path/to/database.db')
# db.insert_data('exp', {'PDBId': '1XYZ', 'RNAMotif': 'motif1, motif2', 'AAMotif': 'motifA, motifB', 'Distance': 2.5, 'ContactList': 'contact1, contact2'})
# data = db.fetch_data('exp', {'PDBId': '1XYZ'})
# print(data)
# db.close()
