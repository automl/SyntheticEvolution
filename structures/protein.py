import sqlite3
import os
import ast
from database.startConfig import StartConfig
from database.databaseMethods import DatabaseMethods

class Protein:
    config = StartConfig()
    pred_table_name = config.pred_table
    exp_table_name = config.ref_table
    db_path = config.db_path
    def __init__(self, id: str, is_pred: bool = None, protein_length: int = None, file_name: str = None, protein_sequence: str = "",
                 protein_chain_IDs: str = "", protein_chain_number: int = None,
                 is_protein: bool = None):
        self.id = id
        self.is_pred = is_pred
        self.file_name = file_name if is_pred else None
        self.protein_sequence = protein_sequence
        self.protein_chain_IDs = protein_chain_IDs
        self.protein_length = protein_length
        self.protein_chain_number = protein_chain_number
        self.is_protein = is_protein

    def get_protein_sequence(self):
        return self.protein_sequence

    def get_protein_chain_IDs(self):
        return self.protein_chain_IDs

    def get_protein_length(self):
        return self.protein_length

    def get_protein_chain_number(self):
        if self.protein_sequence.startswith('[') and self.protein_sequence.endswith(']'):
            count = self.protein_sequence.count(',') if len(self.protein_sequence) > 2 else 0
            self.protein_chain_number = count + 1
        elif ',' in self.get_protein_chain_IDs():
            self.protein_chain_number = self.get_protein_chain_IDs().count(',') + 1
        else:
            self.protein_chain_number = 1
        return self.protein_chain_number

    def is_single_copy_protein(self):
        if self.protein_sequence.startswith('[') and self.protein_sequence.endswith(']'):
            count = self.get_protein_chain_number()
            protein_sequence_list = ast.literal_eval(self.protein_sequence)
            if not all(protein_sequence_list[0] == protein_sequence_list[i] for i in range(1, count)):
                return True
            else:
                return False
        elif ',' in self.get_protein_chain_IDs():
            return False
        else:
            return True

    @classmethod
    def get_protein_from_db(cls, id: str, file_name=None):
        database_methods = DatabaseMethods()
        id = id.upper()
        if file_name != None:
            row = database_methods.get_table_row(table_name=cls.pred_table_name,
                                                 condition=f"exp_db_id = '{id}' AND FileName = '{file_name}'")
        elif len(id) == 6 and id[0].isalpha():
            row = database_methods.get_table_row(table_name=cls.exp_table_name,
                                                 condition=f"UniprotId = '{id}'")
        elif len(id) == 4 and id[0].isdigit():
            row = database_methods.get_table_row(table_name=cls.exp_table_name,
                                                 condition=f"PDBId = '{id}'")
        # database_methods.close_connection()

        if row is not None and 'ProteinSequence' in row.keys() and row['ProteinSequence']:
            if 'UniprotId' in row:
                id = row['UniprotId']
            elif 'PDBId' in row:
                id = row['PDBId']
            elif 'exp_db_id' in row:
                id = row['exp_db_id']
            protein = Protein(
                id=id,
                protein_sequence=row['ProteinSequence'],
                protein_chain_IDs=row['ProteinChainIDs'],
                protein_length=row['ProteinLength'],
                is_protein=True
            )
            return protein
        else:
            protein = Protein(id=id, is_protein=False)
            return protein

    def get_longest_chain_id(self):
        chain_lengths = self.get_protein_length()
        if isinstance(chain_lengths, str):
            chain_lengths = ast.literal_eval(chain_lengths)
        elif isinstance(chain_lengths, int):
            chain_id = self.get_protein_chain_IDs()
            return chain_id
        print("chain_lengths", chain_lengths)
        max_length = max(chain_lengths)
        print("max_length", max_length)
        max_index = chain_lengths.index(max_length)
        chain_ids = self.get_protein_chain_IDs()
        if isinstance(chain_ids, str):
            chain_ids = ast.literal_eval(chain_ids)
        # print("chain_ids", chain_ids)
        # print("max_index", max_index)
        longest_chain = chain_ids[max_index]
        # print("longest_chain", longest_chain)
        return longest_chain

# protein = Protein.get_protein_from_db(id="8t2r")
# print(protein.get_protein_length())
# print(protein.get_protein_sequence())
# print(protein.is_single_copy_protein())
# print(protein.get_protein_chain_number(), protein.get_longest_chain_id())