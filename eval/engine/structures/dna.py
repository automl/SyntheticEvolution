import sqlite3
import os
import ast
from database.startConfig import StartConfig
from database.databaseMethods import DatabaseMethods

class DNA:
    config = StartConfig()
    pred_table_name = config.pred_table
    exp_table_name = config.ref_table
    db_path = config.db_path
    def __init__(self,  id: str, is_pred: bool = None, dna_length: int = None, file_name: str = None, dna_sequence: str = "",
                 dna_chain_IDs: str = "", dna_metrics=None, dna_motif: str = "", dna_chain_number: int = None, is_dna: bool = None):
        self.id = id
        self.is_pred = is_pred
        self.file_name = file_name if is_pred else None
        self.dna_sequence = dna_sequence
        self.dna_chain_IDs = dna_chain_IDs
        self.dna_metrics = dna_metrics if dna_metrics is not None else []
        self.dna_motif = dna_motif
        # self.class_type = class_type
        self.dna_length = dna_length
        # print(f"Initialized RNA with id: {id}, is_pred: {is_pred}, file_name: {file_name}")
        self.dna_chain_number = dna_chain_number
        self.is_dna = is_dna

    @classmethod
    def get_dna_from_db(cls, id: str, file_name=None):
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

        if row is not None and 'DNASequence' in row.keys() and row['DNASequence']:
            if 'UniprotId' in row:
                id = row['UniprotId']
            elif 'PDBId' in row:
                id = row['PDBId']
            elif 'exp_db_id' in row:
                id = row['exp_db_id']
            dna = DNA(
                id=id,
                dna_sequence=row['DNASequence'],
                dna_chain_IDs=row['DNAChainIDs'],
                dna_length=row['DNALength'],
                dna_motif=row['DNAMotif'],
                is_dna=True
            )
            return dna
        else:
            dna = DNA(id=id, is_dna=False)
            return dna

    def get_dna_sequence(self):
        return self.dna_sequence

    def get_dna_chain_IDs(self):
        return self.dna_chain_IDs

    def get_dna_length(self):
        return self.dna_length

    def get_dna_motif(self):
        return self.dna_motif

    def get_dna_chain_number(self):
        if self.dna_sequence.startswith('[') and self.dna_sequence.endswith(']'):
            count = self.dna_sequence.count(',') if len(self.dna_sequence) > 2 else 0
            self.dna_chain_number = count + 1
        elif ',' in self.get_dna_chain_IDs():
            self.dna_chain_number = self.get_dna_chain_IDs().count(',') + 1
        else:
            self.dna_chain_number = 1
        return self.dna_chain_number

    def is_single_copy_dna(self):
        if self.dna_sequence.startswith('[') and self.dna_sequence.endswith(']'):
            count = self.get_dna_chain_number()
            dna_sequence_list = ast.literal_eval(self.dna_sequence)
            if not all(dna_sequence_list[0] == dna_sequence_list[i] for i in range(1, count)):
                return True
            else:
                return False
        elif ',' in self.get_dna_chain_IDs():
            return False
        else:
            return True

    def get_longest_chain_id(self):
        chain_lengths = self.get_dna_length()
        if isinstance(chain_lengths, str):
            chain_lengths = ast.literal_eval(chain_lengths)
        elif isinstance(chain_lengths, (int, float)):
            chain_id = self.get_dna_chain_IDs()
            # print("chain_lengths",chain_lengths)
            return chain_id
        # print("chain_lengths", chain_lengths)
        max_length = max(chain_lengths)
        # print("max_length", max_length)
        max_index = chain_lengths.index(max_length)
        chain_ids = self.get_dna_chain_IDs()
        if isinstance(chain_ids, str):
            chain_ids = ast.literal_eval(chain_ids)
        longest_chain = chain_ids[max_index]
        return longest_chain

#
# dna = DNA.get_dna_from_unifiedTable(id="7YEY")
# print(dna.get_dna_sequence())
# print(dna.is_single_copy_dna())
# print(dna.get_dna_chain_number())