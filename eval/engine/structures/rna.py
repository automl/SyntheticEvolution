import sqlite3
import os
import ast
from database.startConfig import StartConfig
from database.databaseMethods import DatabaseMethods
import numpy as np

class RNA:
    config = StartConfig()
    pred_table_name = config.pred_table
    exp_table_name = config.ref_table
    db_path = config.db_path
    def __init__(self,  id: str, is_pred: bool = None, rna_length: int = None, file_name: str = None, rna_sequence: str = "",
                 rna_chain_IDs: str = "", rna_metrics=None, rna_motif: str = "", rna_chain_number: int = None,
                 is_rna: bool = None, chain_pairs: str = "", motif_lengths: str = ""):
        self.id = id
        self.is_pred = is_pred
        self.file_name = file_name if is_pred else None
        self.rna_sequence = rna_sequence
        self.rna_chain_IDs = rna_chain_IDs
        self.rna_metrics = rna_metrics if rna_metrics is not None else []
        self.rna_motif = rna_motif
        # self.class_type = class_type
        self.rna_length = rna_length
        # print(f"Initialized RNA with id: {id}, is_pred: {is_pred}, file_name: {file_name}")
        self.rna_chain_number = rna_chain_number
        self.is_rna = is_rna
        self.chain_pairs = chain_pairs
        self.motif_lengths = motif_lengths

    @classmethod
    def get_rna_from_db(cls, id: str, file_name=None):
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

        # database_methods.close_connection() # remove?
        if row is not None and 'RNASequence' in row.keys() and row['RNASequence']:
            if 'UniprotId' in row:
                id = row['UniprotId']
            elif 'PDBId' in row:
                id = row['PDBId']
            elif 'exp_db_id' in row:
                id = row['exp_db_id']

            if cls.exp_table_name == "exp_protein_rna":
                rna = RNA(
                    id=id,
                    rna_sequence=row['RNASequence'],
                    rna_chain_IDs=row['RNAChainIDs'],
                    rna_length=row['RNALength'],
                    rna_motif=row['RNAMotif'],
                    is_rna=True,
                    chain_pairs=row['ChainIDpairList_proteinRNA'],
                    motif_lengths=row['RNAMotifLength']
                )
            elif cls.exp_table_name == "exp_rna_rna":
                rna = RNA(
                    id=id,
                    rna_sequence=row['RNASequence'],
                    rna_chain_IDs=row['RNAChainIDs'],
                    rna_length=row['RNALength'],
                    rna_motif=row['RNAMotif'],
                    is_rna=True,
                    chain_pairs=row['ChainIDpairList_rnaRNA'],
                    motif_lengths=row['RNAMotifLength']
                )

            return rna
        else:
            rna = RNA(id=id, is_rna=False)
            return rna

    def get_rna_sequence(self):
        return self.rna_sequence

    def get_rna_chain_IDs(self):
        return self.rna_chain_IDs

    def get_rna_length(self):
        return self.rna_length

    def get_rna_motif(self):
        return self.rna_motif

    # def get_rna_motif_length(self):
    #     return self.motif_lengths

    def get_rna_chain_number(self):
        if self.rna_sequence.startswith('[') and self.rna_sequence.endswith(']'):
            count = self.rna_sequence.count(',') if len(self.rna_sequence) > 2 else 0
            self.rna_chain_number = count + 1
        elif ',' in self.get_rna_chain_IDs():
            self.rna_chain_number = self.get_rna_chain_IDs().count(',') + 1
        else:
            self.rna_chain_number = 1
        return self.rna_chain_number

    def is_single_copy_rna(self):
        if self.rna_sequence.startswith('[') and self.rna_sequence.endswith(']'):
            count = self.get_rna_chain_number()
            rna_sequence_list = ast.literal_eval(self.rna_sequence)
            if not all(rna_sequence_list[0] == rna_sequence_list[i] for i in range(1, count)):
                return True
            else:
                return False
        elif ',' in self.get_rna_chain_IDs():
            return False
        else:
            return True

    @classmethod
    def get_filename_from_id(cls, select_table, pdb_id: str):
        connection = sqlite3.connect(cls.db_path)
        connection.row_factory = sqlite3.Row
        cursor = connection.cursor()
        query = f"SELECT filename FROM {select_table} WHERE exp_db_id = ?"
        cursor.execute(query, (pdb_id,))
        file_name = cursor.fetchone()
        if file_name:
            return file_name[0]
        else:
            return None
            connection.close()

    def get_longest_chain_id(self):
        chain_lengths = self.get_rna_length()
        if isinstance(chain_lengths, str):
            chain_lengths = ast.literal_eval(chain_lengths)
        elif isinstance(chain_lengths, int):
            chain_id = self.get_rna_chain_IDs()
            # print("chain_lengths", chain_lengths)
            return chain_id
        # print("chain_lengths", chain_lengths)
        max_length = max(chain_lengths)
        # print("max_length", max_length)
        max_index = chain_lengths.index(max_length)
        chain_ids = self.get_rna_chain_IDs()
        if isinstance(chain_ids, str):
            chain_ids = ast.literal_eval(chain_ids)
        # print("chain_ids", chain_ids)
        # print("max_index", max_index)
        longest_chain = chain_ids[max_index]
        # print("longest_chain", longest_chain)
        return longest_chain

    # def is_rna_modified(self):
    #     if not self.rna_modified or self.rna_modified == ["[]"]:
    #         return False
    #     return True

    def get_chain_pair_with_longest_motif(self):
        """Get the chain pair that corresponds to the longest RNA motif"""
        try:
            # print("self.chain_pairs", self.chain_pairs)
            # print("self.motif_lengths", self.motif_lengths)
            # Convert string representations to lists
            chain_pairs = eval(self.chain_pairs) if isinstance(self.chain_pairs, str) else self.chain_pairs
            motif_lengths = eval(self.motif_lengths) if isinstance(self.motif_lengths, str) else self.motif_lengths
            # print("chain_pairs",chain_pairs)
            # print("motif_lengths", motif_lengths)

            if isinstance(motif_lengths, int):
                max_length_idx = 0
            else:
                # Find index of longest motif
                max_length_idx = motif_lengths.index(max(motif_lengths))
                # print("chain_pairs[max_length_idx]", chain_pairs[max_length_idx])
            
            # Return corresponding chain pair
            return chain_pairs[max_length_idx]
            
        except (ValueError, SyntaxError, IndexError) as e:
            print(f"Error processing chain pairs and motif lengths: {e}")
            return None

# rna_instance = RNA.get_rna_from_db(id="8hb1", file_name='fold_8hb1_s1_model_0.cif')
# print("rna_instance", rna_instance.is_pred)
# print(rna_instance.get_rna_chain_number())
# print(rna_instance.get_rna_chain_IDs(), rna_instance.get_rna_length())
# print(rna_instance.get_longest_chain_id())
# rna_instance2 = RNA.get_rna_from_unifiedTable(id="7XWZ", file_name='fold_7xwz_s1_model_0.cif')
# print("rna_instance2", rna_instance2.is_pred)
# print(rna_instance2.get_rna_sequence())
# print(rna_instance2.get_rna_chain_number(), rna_instance2.is_single_copy_rna(), rna_instance2.get_rna_chain_IDs())

# rna_instance2 = RNA.get_rna_from_db(id="8t2r")
# print("rna_instance2", rna_instance2.get_chain_pair_with_longest_motif())
# print(rna_instance2.is_single_copy_rna())
# print(rna_instance2.get_rna_chain_IDs(), rna_instance2.is_single_copy_rna())
# print(rna_instance2.chain_pairs)

# rna_instance2 = RNA.get_rna_from_db(id="2ded")
# print(rna_instance2.is_rna)