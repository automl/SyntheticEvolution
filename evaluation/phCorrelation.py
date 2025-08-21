import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import re
import ast
import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from plots.plotCreator import PlotCreator
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig
from Bio.Align import PairwiseAligner
from collections import defaultdict
from collections import Counter


class PHCorrelationAnalyzer(DatabaseMethods):
    def __init__(self, db_path='database/interface.db', single_chain_only=False, msa_option=None):
        DatabaseMethods.__init__(self)
        self.plot_creator = PlotCreator('pH_correlation', msa_option, single_chain_only)
        self.config = StartConfig()

    def calculate_residue_match(self, pred_contacts, exp_contacts):
        """Calculate how well predicted residues match experimental residues"""
        if not pred_contacts or not exp_contacts:
            return 0.0

        pred_set = set((c['res_name'], c['res_num']) for c in pred_contacts)
        exp_set = set((c['res_name'], c['res_num']) for c in exp_contacts)

        # Flexible matching: name matches and index is within ±1
        matched = set()
        for pname, pnum in pred_set:
            for ename, enum in exp_set:
                if pname == ename and abs(pnum - enum) <= 1:
                    matched.add((pname, pnum))
                    break  # Only count each pred residue once

        # Union for Jaccard: all unique (name, num) in both sets
        union = pred_set | exp_set
        return len(matched) / len(union) if union else 0.0

    def align_motif_to_sequence(self, contacts, protein_sequence_data):
        """Aligns a motif to a protein sequence using PairwiseAligner."""
        try:
            if not contacts or not protein_sequence_data:
                return contacts

            sequences = []
            try:
                if isinstance(protein_sequence_data, str) and protein_sequence_data.strip().startswith('['):
                    sequences = ast.literal_eval(protein_sequence_data)
                elif isinstance(protein_sequence_data, list):
                    sequences = protein_sequence_data
                else:
                    sequences = [str(protein_sequence_data)]
            except (ValueError, SyntaxError):
                sequences = [str(protein_sequence_data)]

            if not sequences:
                return contacts

            target_protein_sequence = max(sequences, key=len)

            contacts_sorted = sorted(contacts, key=lambda c: c['res_num'])
            aa3_to_1 = {
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
                'HIS': 'H',
                'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
                'TRP': 'W',
                'TYR': 'Y', 'VAL': 'V', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'NEP': 'H', 'HIP': 'H', 'ALY': 'K',
                'MLY': 'K',
                'M3L': 'K', 'MLZ': 'K', '2MR': 'R', 'AGM': 'R', 'MCS': 'C', 'HYP': 'P', 'HY3': 'H', 'LYZ': 'K',
                'AHB': 'A',
                'P1L': 'P', 'SNN': 'S', 'SNC': 'C', 'TRF': 'W', 'KCR': 'K', 'CIR': 'R', 'YHA': 'Y'
            }

            unique_contacts_for_alignment = []
            seen_contacts = set()
            for c in contacts_sorted:
                contact_tuple = (c['res_name'], c['res_num'])
                if contact_tuple not in seen_contacts:
                    if c['res_name'] in aa3_to_1:
                        unique_contacts_for_alignment.append(c)
                        seen_contacts.add(contact_tuple)

            if not unique_contacts_for_alignment: return contacts

            motif_letters = [aa3_to_1[c['res_name']] for c in unique_contacts_for_alignment]
            motif_indices = [c['res_num'] for c in unique_contacts_for_alignment]

            gapped_motif = motif_letters[0]
            for i in range(1, len(motif_letters)):
                gap = motif_indices[i] - motif_indices[i - 1] - 1
                gapped_motif += '-' * max(0, gap) + motif_letters[i]

            aligner = PairwiseAligner()
            aligner.mode = 'local'
            aligner.open_gap_score = -10
            aligner.extend_gap_score = -0.5
            alignments = aligner.align(target_protein_sequence, gapped_motif)

            if not alignments: return contacts

            alignment = alignments[0]
            print(alignment)
            prot_start = alignment.aligned[0][0][0]

            correction_map = {}
            motif_char_idx = 0
            for i, align_char in enumerate(alignment.query):
                if align_char != '-':
                    original_contact = unique_contacts_for_alignment[motif_char_idx]
                    old_res_name = original_contact['res_name']
                    old_res_num = original_contact['res_num']
                    new_res_num = prot_start + i + 1
                    correction_map[(old_res_name, old_res_num)] = new_res_num
                    motif_char_idx += 1

            final_corrected_contacts = []
            for c in contacts:
                contact_key = (c['res_name'], c['res_num'])
                if contact_key in correction_map:
                    new_c = c.copy()
                    new_c['res_num'] = correction_map[contact_key]
                    final_corrected_contacts.append(new_c)

            return final_corrected_contacts
        except Exception:
            return contacts

    def align_rna_motif_to_sequence(self, contacts, rna_sequence_data):
        """Aligns an RNA motif to an RNA sequence using PairwiseAligner."""
        try:
            if not contacts or not rna_sequence_data:
                return contacts

            sequences = []
            try:
                if isinstance(rna_sequence_data, str) and rna_sequence_data.strip().startswith('['):
                    sequences = ast.literal_eval(rna_sequence_data)
                elif isinstance(rna_sequence_data, list):
                    sequences = rna_sequence_data
                else:
                    sequences = [str(rna_sequence_data)]
            except (ValueError, SyntaxError):
                sequences = [str(rna_sequence_data)]

            if not sequences:
                return contacts

            target_rna_sequence = max(sequences, key=len)

            contacts_sorted = sorted(contacts, key=lambda c: c['res_num'])
            
            # RNA nucleotides are already single letters
            unique_contacts_for_alignment = []
            seen_contacts = set()
            for c in contacts_sorted:
                contact_tuple = (c['res_name'], c['res_num'])
                if contact_tuple not in seen_contacts:
                    unique_contacts_for_alignment.append(c)
                    seen_contacts.add(contact_tuple)

            if not unique_contacts_for_alignment: return contacts

            motif_letters = [c['res_name'] for c in unique_contacts_for_alignment]
            motif_indices = [c['res_num'] for c in unique_contacts_for_alignment]

            gapped_motif = motif_letters[0]
            for i in range(1, len(motif_letters)):
                gap = motif_indices[i] - motif_indices[i - 1] - 1
                gapped_motif += '-' * max(0, gap) + motif_letters[i]

            aligner = PairwiseAligner()
            aligner.mode = 'local'
            aligner.open_gap_score = -10
            aligner.extend_gap_score = -0.5
            alignments = aligner.align(target_rna_sequence, gapped_motif)

            if not alignments: return contacts

            alignment = alignments[0]
            rna_start = alignment.aligned[0][0][0]

            correction_map = {}
            motif_char_idx = 0
            for i, align_char in enumerate(alignment.query):
                if align_char != '-':
                    original_contact = unique_contacts_for_alignment[motif_char_idx]
                    old_res_name = original_contact['res_name']
                    old_res_num = original_contact['res_num']
                    new_res_num = rna_start + i + 1
                    correction_map[(old_res_name, old_res_num)] = new_res_num
                    motif_char_idx += 1

            final_corrected_contacts = []
            for c in contacts:
                contact_key = (c['res_name'], c['res_num'])
                if contact_key in correction_map:
                    new_c = c.copy()
                    new_c['res_num'] = correction_map[contact_key]
                    final_corrected_contacts.append(new_c)

            return final_corrected_contacts
        except Exception:
            return contacts

    def safe_literal_eval(self, val):
        if isinstance(val, str) and val.strip() not in ["", "None", "nan"]:
            try:
                return ast.literal_eval(val)
            except Exception:
                # If literal_eval fails, return the original string as a single item
                return [val.strip()]
        return []

    def get_sequence_to_motifs(self, row):
        motifs = self.parse_motifs(row['AAMotif'])
        chain_pairs = self.safe_literal_eval(row.get('ChainIDpairList_proteinRNA', ''))
        protein_chain_ids = self.safe_literal_eval(row.get('ProteinChainIDs', ''))
        protein_sequences = self.safe_literal_eval(row.get('ProteinSequence', ''))
        seq_to_motifs = {}
        
        # Handle case where there's only one protein sequence
        if len(protein_sequences) == 1:
            # All motifs should be mapped to this single sequence
            seq_to_motifs[protein_sequences[0]] = motifs
            return seq_to_motifs
        
        # Handle multiple protein sequences
        for i, motif in enumerate(motifs):
            # Find the protein chain for this motif
            protein_seq = None
            if i < len(chain_pairs):
                protein_chain = chain_pairs[i][0]
                if protein_chain in protein_chain_ids:
                    idx = protein_chain_ids.index(protein_chain)
                    protein_seq = protein_sequences[idx]
            elif len(protein_chain_ids) == 1:
                protein_seq = protein_sequences[0]
            if protein_seq:
                seq_to_motifs.setdefault(protein_seq, []).append(motif)
        return seq_to_motifs

    def get_rna_sequence_to_motifs(self, row, modified_residues_nucleotide):
        """Get RNA sequence to motifs mapping, handling modified nucleotides"""
        motifs = self.parse_rna_motifs(row['RNAMotif'], modified_residues_nucleotide)
        chain_pairs = self.safe_literal_eval(row.get('ChainIDpairList_proteinRNA', ''))
        rna_chain_ids = self.safe_literal_eval(row.get('RNAChainIDs', ''))
        rna_sequences = self.safe_literal_eval(row.get('RNASequence', ''))
        seq_to_motifs = {}
        
        # Handle case where there's only one RNA sequence
        if len(rna_sequences) == 1:
            # All motifs should be mapped to this single sequence
            seq_to_motifs[rna_sequences[0]] = motifs
            return seq_to_motifs
        
        # Handle multiple RNA sequences
        for i, motif in enumerate(motifs):
            # Find the RNA chain for this motif
            rna_seq = None
            if i < len(chain_pairs):
                rna_chain = chain_pairs[i][1]  # RNA chain is second in the pair
                if rna_chain in rna_chain_ids:
                    idx = rna_chain_ids.index(rna_chain)
                    rna_seq = rna_sequences[idx]
            elif len(rna_chain_ids) == 1:
                rna_seq = rna_sequences[0]
            if rna_seq:
                seq_to_motifs.setdefault(rna_seq, []).append(motif)
        return seq_to_motifs

    def update_aa_match_scores(self):
        """Calculate and update AAmatch_score in the database"""
        print("Calculating and updating AAmatch_score in database...")
        pred_table = self.config.get_predicted_table_name()
        try:
            check_query = f"SELECT COUNT(*) FROM pragma_table_info('{pred_table}') WHERE name='AAmatch_score'"
            cursor = self.connection.cursor()
            cursor.execute(check_query)
            if not bool(cursor.fetchone()[0]):
                print("Creating AAmatch_score column...")
                cursor.execute(f"ALTER TABLE {pred_table} ADD COLUMN AAmatch_score REAL")
                self.connection.commit()
            else:
                print("AAmatch_score column already exists.")
        except Exception as e:
            print(f"Error checking/creating AAmatch_score column: {e}")
            return
        query = f"""
        SELECT 
            p.exp_db_id, p.AAMotif as pred_AAMotif, p.ChainIDpairList_proteinRNA as pred_ChainIDpairList_proteinRNA, p.ProteinChainIDs as pred_ProteinChainIDs, p.ProteinSequence as pred_ProteinSequence,
            e.AAMotif as exp_AAMotif, e.ChainIDpairList_proteinRNA as exp_ChainIDpairList_proteinRNA, e.ProteinChainIDs as exp_ProteinChainIDs, e.ProteinSequence as exp_ProteinSequence
        FROM {pred_table} p
        LEFT JOIN exp_protein_rna e ON p.exp_db_id = e.PDBId
        WHERE p.AAMotif IS NOT NULL AND e.AAMotif IS NOT NULL
        """
        try:
            df = pd.read_sql_query(query, self.connection)
            print(f"Retrieved {len(df)} samples for AAmatch_score calculation")
            final_pred_list, final_exp_list, final_score_list = [], [], []
            for idx, row in df.iterrows():
                # Build sequence-to-motifs mapping for predicted and experimental
                pred_row = {
                    'AAMotif': row['pred_AAMotif'],
                    'ChainIDpairList_proteinRNA': row['pred_ChainIDpairList_proteinRNA'],
                    'ProteinChainIDs': row['pred_ProteinChainIDs'],
                    'ProteinSequence': row['pred_ProteinSequence']
                }
                exp_row = {
                    'AAMotif': row['exp_AAMotif'],
                    'ChainIDpairList_proteinRNA': row['exp_ChainIDpairList_proteinRNA'],
                    'ProteinChainIDs': row['exp_ProteinChainIDs'],
                    'ProteinSequence': row['exp_ProteinSequence']
                }
                pred_seq_to_motifs = self.get_sequence_to_motifs(pred_row)
                exp_seq_to_motifs = self.get_sequence_to_motifs(exp_row)
                
                # Debug information
                if idx < 5:  # Only print for first 5 samples to avoid spam
                    print(f"\nSample {idx}: {row['exp_db_id']}")
                    print(f"Pred sequences: {len(pred_seq_to_motifs)}")
                    print(f"Exp sequences: {len(exp_seq_to_motifs)}")
                    if pred_seq_to_motifs:
                        print(f"Pred motifs per seq: {[len(motifs) for motifs in pred_seq_to_motifs.values()]}")
                    if exp_seq_to_motifs:
                        print(f"Exp motifs per seq: {[len(motifs) for motifs in exp_seq_to_motifs.values()]}")
                
                best_score = 0.0
                best_pred_motif = []
                best_exp_motif = []
                all_motif_pairs = []  # Store all motif pairs and their scores
                
                for seq in pred_seq_to_motifs:
                    if seq in exp_seq_to_motifs:
                        for pred_motif in pred_seq_to_motifs[seq]:
                            for exp_motif in exp_seq_to_motifs[seq]:
                                # Debug information for first few comparisons
                                if idx < 3:  # Only for first 3 samples
                                    print(f"  Comparing motifs:")
                                    print(f"    Pred: {pred_motif[:3]}...")  # Show first 3 residues
                                    print(f"    Exp:  {exp_motif[:3]}...")   # Show first 3 residues
                                # Align both motifs to the sequence before scoring
                                pred_aligned = self.align_motif_to_sequence(
                                    [{'res_name': n, 'res_num': num} for n, num in pred_motif], seq)
                                exp_aligned = self.align_motif_to_sequence(
                                    [{'res_name': n, 'res_num': num} for n, num in exp_motif], seq)
                                # Convert back to (res_name, res_num) tuples for scoring
                                pred_aligned_tuples = [(c['res_name'], c['res_num']) for c in pred_aligned]
                                exp_aligned_tuples = [(c['res_name'], c['res_num']) for c in exp_aligned]
                                score = self.calculate_residue_match_simple(pred_aligned_tuples, exp_aligned_tuples)
                                if score == 0.0:
                                    triplet_score, triplet_exp = self.motif_triplet_offset_match(pred_motif, exp_motif)
                                    if triplet_score > score:
                                        score = triplet_score
                                        best_pred_motif = pred_aligned_tuples
                                        best_exp_motif = triplet_exp
                                # Store this motif pair and its score
                                all_motif_pairs.append((pred_aligned_tuples, exp_aligned_tuples, score))
                
                # Calculate weighted average score based on motif lengths
                if all_motif_pairs:
                    total_weight = 0
                    weighted_sum = 0
                    
                    # Debug information for first few samples
                    if idx < 3:
                        print(f"  Motif pairs for {row['exp_db_id']}:")
                        for i, (pred_motif, exp_motif, score) in enumerate(all_motif_pairs):
                            print(f"    Pair {i+1}: Pred({len(pred_motif)} residues) vs Exp({len(exp_motif)} residues) = {score:.3f}")
                    
                    for pred_motif, exp_motif, score in all_motif_pairs:
                        # Weight by the length of the longer motif
                        weight = max(len(pred_motif), len(exp_motif))
                        weighted_sum += score * weight
                        total_weight += weight
                    
                    best_score = weighted_sum / total_weight if total_weight > 0 else 0.0
                    
                    # Debug information
                    if idx < 3:
                        print(f"  Weighted average score: {best_score:.3f}")
                    
                    # Use the motif pair with the highest individual score for storage
                    best_pair = max(all_motif_pairs, key=lambda x: x[2])
                    best_pred_motif = best_pair[0]
                    best_exp_motif = best_pair[1]
                final_pred_list.append(best_pred_motif)
                final_exp_list.append(best_exp_motif)
                final_score_list.append(best_score)
            df['pred_contact_data'] = final_pred_list
            df['exp_contact_data'] = final_exp_list
            df['AAmatch_score'] = final_score_list
            print("Updating database with AAmatch_score values...")
            cursor = self.connection.cursor()
            for idx, row in df.iterrows():
                try:
                    cursor.execute(f"UPDATE {pred_table} SET AAmatch_score = ? WHERE exp_db_id = ?", (row['AAmatch_score'] if pd.notna(row['AAmatch_score']) else None, row['exp_db_id']))
                except Exception as e:
                    print(f"Error updating for {row['exp_db_id']}: {e}")
            self.connection.commit()
            valid_scores = pd.to_numeric(df['AAmatch_score'], errors='coerce').dropna()
            print(f"\nAAmatch_score Statistics:\nTotal samples processed: {len(df)}\nValid scores: {len(valid_scores)}")
            if not valid_scores.empty:
                print(f"Mean: {valid_scores.mean():.3f}, Std: {valid_scores.std():.3f}, Min: {valid_scores.min():.3f}, Max: {valid_scores.max():.3f}")
        except Exception as e:
            print(f"Error in update_aa_match_scores: {e}")
            self.connection.rollback()

    def update_rna_match_scores(self):
        """Calculate and update RNAmatch_score in the database"""
        print("Calculating and updating RNAmatch_score in database...")
        pred_table = self.config.get_predicted_table_name()
        try:
            check_query = f"SELECT COUNT(*) FROM pragma_table_info('{pred_table}') WHERE name='RNAmatch_score'"
            cursor = self.connection.cursor()
            cursor.execute(check_query)
            if not bool(cursor.fetchone()[0]):
                print("Creating RNAmatch_score column...")
                cursor.execute(f"ALTER TABLE {pred_table} ADD COLUMN RNAmatch_score REAL")
                self.connection.commit()
            else:
                print("RNAmatch_score column already exists.")
        except Exception as e:
            print(f"Error checking/creating RNAmatch_score column: {e}")
            return
        
        # Modified nucleotide mapping
        modified_residues_nucleotide = {
            '6MZ': 'A', 'PSU': 'U', '5MC': 'C', 'OMC': 'C', '4OC': 'C', '5MU': 'U',
            'OMU': 'U', 'UR3': 'U', 'A2M': 'A', 'MA6': 'A', '2MG': 'G', 'OMG': 'G',
            '7MG': 'G', 'RSQ': 'G', '5CM': 'C', 'C34': 'C', '5HC': 'C', '6OG': 'G',
            '6MA': 'A', '1CC': 'C', '8OG': 'G', '5FC': 'C', '3DR': 'T'
        }
        
        query = f"""
        SELECT 
            p.exp_db_id, p.RNAMotif as pred_RNAMotif, p.ChainIDpairList_proteinRNA as pred_ChainIDpairList_proteinRNA, p.RNAChainIDs as pred_RNAChainIDs, p.RNASequence as pred_RNASequence,
            e.RNAMotif as exp_RNAMotif, e.ChainIDpairList_proteinRNA as exp_ChainIDpairList_proteinRNA, e.RNAChainIDs as exp_RNAChainIDs, e.RNASequence as exp_RNASequence
        FROM {pred_table} p
        LEFT JOIN exp_protein_rna e ON p.exp_db_id = e.PDBId
        WHERE p.RNAMotif IS NOT NULL AND e.RNAMotif IS NOT NULL
        """
        try:
            df = pd.read_sql_query(query, self.connection)
            print(f"Retrieved {len(df)} samples for RNAmatch_score calculation")
            final_pred_list, final_exp_list, final_score_list = [], [], []
            for idx, row in df.iterrows():
                # Build sequence-to-motifs mapping for predicted and experimental
                pred_row = {
                    'RNAMotif': row['pred_RNAMotif'],
                    'ChainIDpairList_proteinRNA': row['pred_ChainIDpairList_proteinRNA'],
                    'RNAChainIDs': row['pred_RNAChainIDs'],
                    'RNASequence': row['pred_RNASequence']
                }
                exp_row = {
                    'RNAMotif': row['exp_RNAMotif'],
                    'ChainIDpairList_proteinRNA': row['exp_ChainIDpairList_proteinRNA'],
                    'RNAChainIDs': row['exp_RNAChainIDs'],
                    'RNASequence': row['exp_RNASequence']
                }
                pred_seq_to_motifs = self.get_rna_sequence_to_motifs(pred_row, modified_residues_nucleotide)
                exp_seq_to_motifs = self.get_rna_sequence_to_motifs(exp_row, modified_residues_nucleotide)
                best_score = 0.0
                best_pred_motif = []
                best_exp_motif = []
                all_motif_pairs = []  # Store all motif pairs and their scores
                
                for seq in pred_seq_to_motifs:
                    if seq in exp_seq_to_motifs:
                        for pred_motif in pred_seq_to_motifs[seq]:
                            for exp_motif in exp_seq_to_motifs[seq]:
                                # Align both motifs to the sequence before scoring
                                pred_aligned = self.align_rna_motif_to_sequence(
                                    [{'res_name': n, 'res_num': num} for n, num in pred_motif], seq)
                                exp_aligned = self.align_rna_motif_to_sequence(
                                    [{'res_name': n, 'res_num': num} for n, num in exp_motif], seq)
                                # Convert back to (res_name, res_num) tuples for scoring
                                pred_aligned_tuples = [(c['res_name'], c['res_num']) for c in pred_aligned]
                                exp_aligned_tuples = [(c['res_name'], c['res_num']) for c in exp_aligned]
                                score = self.calculate_residue_match_simple(pred_aligned_tuples, exp_aligned_tuples)
                                if score == 0.0:
                                    triplet_score, triplet_exp = self.rna_motif_triplet_offset_match(pred_motif, exp_motif)
                                    if triplet_score > score:
                                        score = triplet_score
                                        pred_aligned_tuples = pred_motif
                                        exp_aligned_tuples = triplet_exp
                                # Store this motif pair and its score
                                all_motif_pairs.append((pred_aligned_tuples, exp_aligned_tuples, score))
                
                # Calculate weighted average score based on motif lengths
                if all_motif_pairs:
                    total_weight = 0
                    weighted_sum = 0
                    
                    # Debug information for first few samples
                    if idx < 3:
                        print(f"  RNA Motif pairs for {row['exp_db_id']}:")
                        for i, (pred_motif, exp_motif, score) in enumerate(all_motif_pairs):
                            print(f"    Pair {i+1}: Pred({len(pred_motif)} residues) vs Exp({len(exp_motif)} residues) = {score:.3f}")
                    
                    for pred_motif, exp_motif, score in all_motif_pairs:
                        # Weight by the length of the longer motif
                        weight = max(len(pred_motif), len(exp_motif))
                        weighted_sum += score * weight
                        total_weight += weight
                    
                    best_score = weighted_sum / total_weight if total_weight > 0 else 0.0
                    
                    # Debug information
                    if idx < 3:
                        print(f"  RNA Weighted average score: {best_score:.3f}")
                    
                    # Use the motif pair with the highest individual score for storage
                    best_pair = max(all_motif_pairs, key=lambda x: x[2])
                    best_pred_motif = best_pair[0]
                    best_exp_motif = best_pair[1]
                
                final_pred_list.append(best_pred_motif)
                final_exp_list.append(best_exp_motif)
                final_score_list.append(best_score)
            df['pred_rna_contact_data'] = final_pred_list
            df['exp_rna_contact_data'] = final_exp_list
            df['RNAmatch_score'] = final_score_list
            print("Updating database with RNAmatch_score values...")
            cursor = self.connection.cursor()
            for idx, row in df.iterrows():
                try:
                    cursor.execute(f"UPDATE {pred_table} SET RNAmatch_score = ? WHERE exp_db_id = ?", (row['RNAmatch_score'] if pd.notna(row['RNAmatch_score']) else None, row['exp_db_id']))
                except Exception as e:
                    print(f"Error updating for {row['exp_db_id']}: {e}")
            self.connection.commit()
            valid_scores = pd.to_numeric(df['RNAmatch_score'], errors='coerce').dropna()
            print(f"\nRNAmatch_score Statistics:\nTotal samples processed: {len(df)}\nValid scores: {len(valid_scores)}")
            if not valid_scores.empty:
                print(f"Mean: {valid_scores.mean():.3f}, Std: {valid_scores.std():.3f}, Min: {valid_scores.min():.3f}, Max: {valid_scores.max():.3f}")
        except Exception as e:
            print(f"Error in update_rna_match_scores: {e}")
            self.connection.rollback()

    def extract_first_value(self, text_list):
        """Extract first value from a text list like '[0.46, 16.56]'"""
        if not text_list or text_list == 'None':
            return None
        try:
            numbers = re.findall(r'[-+]?\d*\.\d+|\d+', text_list)
            if numbers:
                return float(numbers[0])
            return None
        except:
            return None

    def parse_contact_list(self, contact_list_str):
        """Parse contact list string into structured data"""
        if not contact_list_str or contact_list_str == 'None':
            return []
        try:
            if contact_list_str.strip().startswith('['):
                contacts = ast.literal_eval(contact_list_str)
            else:
                contacts = [contact.strip() for contact in contact_list_str.split(',')]

            parsed_contacts = []
            for contact in contacts:
                match = re.match(r'([A-Z]+)\((\d+)\)-([A-Z])\((\d+)\): (\d+\.\d+)', contact)
                if match:
                    res_name, res_num, rna_chain, rna_pos, distance = match.groups()
                    parsed_contacts.append({
                        'res_name': res_name,
                        'res_num': int(res_num),
                        'rna_chain': rna_chain,
                        'rna_pos': int(rna_pos),
                        'distance': float(distance)
                    })
            return parsed_contacts
        except:
            return []

    def correct_exp_contacts(self, row):
        """Correct experimental contact residue indices by aligning to protein sequence"""
        exp_contacts = row['exp_contact_data']
        pred_contacts = row.get('pred_contact_data', [])

        if not exp_contacts or len(exp_contacts) <= 2:
            return exp_contacts

        exp_indices = set(c['res_num'] for c in exp_contacts)
        pred_indices = set(c['res_num'] for c in pred_contacts)
        if exp_indices & pred_indices:
            return exp_contacts

        exp_names = [c['res_name'] for c in exp_contacts]
        exp_nums = [c['res_num'] for c in exp_contacts]
        pred_names = [c['res_name'] for c in pred_contacts]
        pred_nums = [c['res_num'] for c in pred_contacts]

        def find_matching_triplets(names1, nums1, names2, nums2):
            if len(names1) < 3 or len(names2) < 3:
                return []
            for i in range(len(names1) - 2):
                triplet1 = names1[i:i+3]
                for j in range(len(names2) - 2):
                    triplet2 = names2[j:j+3]
                    if triplet1 == triplet2:
                        d1_1, d1_2 = nums1[i+1] - nums1[i], nums1[i+2] - nums1[i+1]
                        d2_1, d2_2 = nums2[j+1] - nums2[j], nums2[j+2] - nums2[j+1]
                        if d1_1 == d2_1 and d1_2 == d2_2:
                            return [(i, j)]
            return []

        matches = find_matching_triplets(exp_names, exp_nums, pred_names, pred_nums)
        if matches:
            offset = pred_nums[matches[0][1]] - exp_nums[matches[0][0]]
            
            temp_corrected = []
            has_negative_index = False
            for c in exp_contacts:
                new_num = c['res_num'] + offset
                if new_num < 1:
                    has_negative_index = True
                    break
                new_c = c.copy()
                new_c['res_num'] = new_num
                temp_corrected.append(new_c)

            if not has_negative_index:
                print(f"[Remap] Used triplet match for {row.get('exp_db_id', 'unknown')}, offset {offset}")
                return temp_corrected
            else:
                print(f"[Remap] Triplet match for {row.get('exp_db_id', 'unknown')} failed (negative index). Falling back.")

        return exp_contacts

    def calculate_contact_metrics(self, contacts):
        """Calculate various contact-based metrics"""
        if not contacts:
            return {
                'total_contacts': 0,
                'his_contacts': 0,
                'acidic_contacts': 0,
                'basic_contacts': 0,
                'contact_score': 0,
                'avg_distance': 0
            }
        
        his_contacts = [c for c in contacts if c['res_name'] == 'HIS']
        acidic_contacts = [c for c in contacts if c['res_name'] in ['ASP', 'GLU']]
        basic_contacts = [c for c in contacts if c['res_name'] in ['LYS', 'ARG']]
        
        # Calculate contact score (sum of 1/distance for contacts < 3.5Å)
        close_contacts = [c for c in contacts if c['distance'] < 3.5]
        contact_score = sum(1/c['distance'] for c in close_contacts)
        
        return {
            'total_contacts': len(contacts),
            'his_contacts': len(his_contacts),
            'acidic_contacts': len(acidic_contacts),
            'basic_contacts': len(basic_contacts),
            'contact_score': contact_score,
            'avg_distance': sum(c['distance'] for c in contacts) / len(contacts) if contacts else 0
        }

    def get_data(self, single_chain_only=False, msa_option=None):
        """Fetch and process data from database"""
        # Get predicted data
        pred_query = """
        SELECT 
            p.exp_db_id,
            p.Complex_RMSD,
            p.RNA_RMSD,
            p.RNA_LDDT,
            p.RNA_TM,
            p.Complex_TM,
            p.Complex_LDDT,
            p.af3_ranking_score,
            p.af3_rna_pLDDT_avg,
            p.af3_chain_pair_pae_min,
            p.RNA_DI,
            p.RNAmotif_score,
            p.ContactList as pred_contacts,
            p.NumberRNAs,
            p.RNALength,
            p.RNAChainIDs,
            p.ChainIDpairList_proteinRNA,
            p.RNAMotifLength
        FROM pred_protein_rna p
        """
        
        # Get experimental data
        exp_query = """
        SELECT 
            e.PDBId,
            e.exp_pH,
            e.ContactList as exp_contacts,
            e.RNALength,
            e.RNAChainIDs,
            e.ChainIDpairList_proteinRNA,
            e.RNAMotifLength,
            e.ProteinSequence,
            e.ProteinChainIDs
        FROM exp_protein_rna e
        WHERE e.exp_pH IS NOT NULL
        """
        
        pred_df = pd.read_sql_query(pred_query, self.connection)
        exp_df = pd.read_sql_query(exp_query, self.connection)

        # Print initial data info
        print("\nInitial data info:")
        print(f"Predicted samples: {len(pred_df)}")
        print(f"Experimental samples: {len(exp_df)}")
        print(f"Unique PDBIds in predicted: {pred_df['exp_db_id'].nunique()}")
        print(f"Unique PDBIds in experimental: {exp_df['PDBId'].nunique()}")

        # Check for any PDBIds that might be causing issues
        print("\nSample of experimental PDBIds:")
        print(exp_df['PDBId'].head(10))
        
        # Process contact lists for multiple RNA chains
        def process_contacts(row):
            if not row['pred_contacts'] or row['pred_contacts'] == 'None':
                return None
            # print("1row['RNAChainIDs']", row['RNAChainIDs'])

            try:
                # Check if the string starts with double brackets (multiple lists)
                if row['pred_contacts'].strip().startswith('[['):
                    # print("2row['RNAChainIDs']", row['RNAChainIDs'])
                    # Parse the outer list
                    contacts = ast.literal_eval(row['pred_contacts'])
                    if not isinstance(contacts, list) or len(contacts) == 0:
                        return row['pred_contacts']

                    # Handle RNAChainIDs which could be a single string or a list
                    try:
                        rna_chains = ast.literal_eval(row['RNAChainIDs'])
                        # Parse RNA and chain information
                        rna_lengths = ast.literal_eval(row['RNALength'])
                    except (SyntaxError, ValueError):
                        rna_chains = row['RNAChainIDs']  # Use as is if it's a single character
                        # print("@£$!£$£$", rna_chains)
                        
                    chain_pairs = ast.literal_eval(row['ChainIDpairList_proteinRNA'])
                    
                    # Check if we have a single RNA chain
                    if isinstance(rna_chains, str):
                        # print("3row['RNAChainIDs']", row['RNAChainIDs'])
                        # Single RNA chain, use RNAMotifLength to find the best protein chain
                        rna_motif_lengths = ast.literal_eval(row['RNAMotifLength'])
                        max_motif_idx = rna_motif_lengths.index(max(rna_motif_lengths))
                        return str(contacts[max_motif_idx])
                    else:
                        # Multiple RNA chains, use the first longest one
                        max_length = max(rna_lengths)
                        longest_idx = rna_lengths.index(max_length)
                        target_rna_chain = rna_chains[longest_idx]
                        
                        # Find the chain pair containing the target RNA chain
                        for i, pair in enumerate(chain_pairs):
                            if target_rna_chain in pair:
                                return str(contacts[i])

                        return str(contacts[0])  # Fallback to first list if no match found
                else:
                    # Single contact list, return as is
                    return row['pred_contacts']
                
            except Exception as e:
                print(f"Error processing contacts for {row['exp_db_id']}: {str(e)}")
                return row['pred_contacts']  # Return original string on error

        # Apply the same processing to experimental contacts
        def process_exp_contacts(row):
            if not row['exp_contacts'] or row['exp_contacts'] == 'None':
                return None
            # print("1row['RNAChainIDs']", row['RNAChainIDs'])

            try:
                # Check if the string starts with double brackets (multiple lists)
                if row['exp_contacts'].strip().startswith('[['):
                    # Parse the outer list
                    contacts = ast.literal_eval(row['exp_contacts'])
                    if not isinstance(contacts, list) or len(contacts) == 0:
                        return row['exp_contacts']

                    # Handle RNAChainIDs which could be a single string or a list
                    try:
                        rna_chains = ast.literal_eval(row['RNAChainIDs'])
                        # Parse RNA and chain information
                        rna_lengths = ast.literal_eval(row['RNALength'])
                    except (SyntaxError, ValueError):
                        rna_chains = row['RNAChainIDs']  # Use as is if it's a single character
                        # print("2@£$!£$£$", rna_chains)

                    chain_pairs = ast.literal_eval(row['ChainIDpairList_proteinRNA'])

                    if isinstance(rna_chains, str):
                        # print("3row['RNAChainIDs']", row['RNAChainIDs'])
                        # Single RNA chain, use RNAMotifLength to find the best protein chain
                        rna_motif_lengths = ast.literal_eval(row['RNAMotifLength'])
                        max_motif_idx = rna_motif_lengths.index(max(rna_motif_lengths))
                        # print("best protein chain", row['RNAChainIDs'], row['RNAMotifLength'], max_motif_idx, str(contacts[max_motif_idx]))
                        return str(contacts[max_motif_idx])
                    else:
                        # Multiple RNA chains, use the first longest one
                        max_length = max(rna_lengths)
                        longest_idx = rna_lengths.index(max_length)
                        target_rna_chain = rna_chains[longest_idx]

                        # Find the chain pair containing the target RNA chain
                        for i, pair in enumerate(chain_pairs):
                            if target_rna_chain in pair:
                                return str(contacts[i])

                        return str(contacts[0])  # Fallback to first list if no match found
                else:
                    # Single contact list, return as is
                    return row['exp_contacts']

            except Exception as e:
                print(f"Error processing contacts for {row['PDBId']}: {str(e)}")
                return row['exp_contacts']  # Return original string on error
        
        # Process contacts for both predicted and experimental data
        pred_df['pred_contacts'] = pred_df.apply(process_contacts, axis=1)
        # print("@£$%$£%$£%£$", len(pred_df))

        exp_df['exp_contacts'] = exp_df.apply(process_exp_contacts, axis=1)
        # print("2435@£$%$£%$£%£$", len(exp_df))

        # Apply single chain filter if specified
        if single_chain_only:
            pred_df = pred_df[pred_df['NumberRNAs'] == 1]
            exp_df = exp_df[exp_df['PDBId'].isin(pred_df['exp_db_id'])]
        
        # Apply MSA filter if specified
        if msa_option:
            if msa_option == '+MSA':
                pred_df = pred_df[pred_df['exp_db_id'].str.contains('_msa')]
            elif msa_option == '-MSA':
                pred_df = pred_df[~pred_df['exp_db_id'].str.contains('_msa')]
            exp_df = exp_df[exp_df['PDBId'].isin(pred_df['exp_db_id'])]
        
        # Merge predicted and experimental data
        df = pd.merge(pred_df, exp_df, left_on='exp_db_id', right_on='PDBId', how='inner')
        
        # Process Complex_RMSD and RNA_RMSD to get first value
        df['Complex_RMSD'] = df['Complex_RMSD'].apply(self.extract_first_value)
        df['RNA_RMSD'] = df['RNA_RMSD'].apply(self.extract_first_value)
        
        # Convert pH to numeric
        df['exp_pH'] = pd.to_numeric(df['exp_pH'], errors='coerce')
        
        # Convert all metric columns to numeric
        metrics = [
            'Complex_RMSD', 'RNA_RMSD', 'RNA_LDDT', 'RNA_TM', 'Complex_TM',
            'Complex_LDDT', 'af3_ranking_score', 'af3_rna_pLDDT_avg',
            'af3_chain_pair_pae_min', 'RNA_DI', 'RNAmotif_score'
        ]
        
        for metric in metrics:
            df[metric] = pd.to_numeric(df[metric], errors='coerce')
        
        # Add match scores to dataframe
        print("Adding match scores to dataframe...")
        match_scores_query = f"""
        SELECT exp_db_id, AAmatch_score, RNAmatch_score
        FROM {self.config.get_predicted_table_name()}
        WHERE AAmatch_score IS NOT NULL OR RNAmatch_score IS NOT NULL
        """
        try:
            match_scores_df = pd.read_sql_query(match_scores_query, self.connection)
            df = pd.merge(df, match_scores_df, on='exp_db_id', how='left')
            print(f"Added match scores for {len(match_scores_df)} samples")
            print(f"Match scores columns: {list(match_scores_df.columns)}")
        except Exception as e:
            print(f"Warning: Could not add match scores: {e}")
            # Add empty columns if match scores don't exist
            df['AAmatch_score'] = np.nan
            df['RNAmatch_score'] = np.nan
        
        # Print available columns for debugging
        print(f"Available columns in dataframe: {list(df.columns)}")
        
        # Process contact lists and calculate metrics
        df['pred_contact_data'] = df['pred_contacts'].apply(self.parse_contact_list)
        df['exp_contact_data'] = df['exp_contacts'].apply(self.parse_contact_list)
        
        # Correct experimental contact residue indices
        df['exp_contact_data'] = df.apply(self.correct_exp_contacts, axis=1)
        
        # Calculate metrics for both predicted and experimental contacts
        pred_metrics = df['pred_contact_data'].apply(self.calculate_contact_metrics)
        exp_metrics = df['exp_contact_data'].apply(self.calculate_contact_metrics)
        
        # Add contact metrics to dataframe with prefixes
        for metric in pred_metrics.iloc[0].keys():
            df[f'pred_{metric}'] = pred_metrics.apply(lambda x: x[metric])
            df[f'exp_{metric}'] = exp_metrics.apply(lambda x: x[metric])
        
        # Add pH range category
        df['pH_category'] = pd.cut(df['exp_pH'], 
                                 bins=[0, 6, 8, 14],
                                 labels=['Acidic', 'Neutral', 'Basic'])
        
        # Replace inf values with NaN
        df = df.replace([np.inf, -np.inf], np.nan)
        
        # Print data quality information
        print("\nData Quality Information:")
        print("=" * 80)
        print("\nNumber of samples:", len(df))
        print("\nMissing values per column:")
        print(df.isna().sum())
        print("\nValue ranges:")
        print(df.describe())
        print("\nSamples per pH category:")
        print(df['pH_category'].value_counts())
        
        return df

    def calculate_correlations(self, df):
        """Calculate correlations between pH and metrics"""
        metrics = [
            'Complex_RMSD', 'RNA_RMSD', 'RNA_LDDT', 'RNA_TM', 'Complex_TM',
            'Complex_LDDT', 'af3_ranking_score', 'af3_rna_pLDDT_avg',
            'af3_chain_pair_pae_min', 'RNA_DI', 'RNAmotif_score',
            'pred_total_contacts', 'pred_his_contacts', 'pred_acidic_contacts', 
            'pred_basic_contacts', 'pred_contact_score', 'pred_avg_distance',
            'exp_total_contacts', 'exp_his_contacts', 'exp_acidic_contacts',
            'exp_basic_contacts', 'exp_contact_score', 'exp_avg_distance'
        ]
        
        correlations = {}
        for metric in metrics:
            # Create a clean subset of data for this metric
            clean_df = df[['exp_pH', metric]].dropna()
            
            if len(clean_df) < 2:
                print(f"\nWarning: Not enough valid data points for {metric}")
                correlations[metric] = {
                    'pearson_corr': np.nan,
                    'pearson_p': np.nan,
                    'spearman_corr': np.nan,
                    'spearman_p': np.nan,
                    'n_samples': len(clean_df)
                }
                continue
            
            try:
                # Calculate Pearson correlation
                pearson_corr, pearson_p = stats.pearsonr(clean_df['exp_pH'], clean_df[metric])
                
                # Calculate Spearman correlation (non-linear)
                spearman_corr, spearman_p = stats.spearmanr(clean_df['exp_pH'], clean_df[metric])
                
                correlations[metric] = {
                    'pearson_corr': pearson_corr,
                    'pearson_p': pearson_p,
                    'spearman_corr': spearman_corr,
                    'spearman_p': spearman_p,
                    'n_samples': len(clean_df)
                }
            except Exception as e:
                print(f"\nError calculating correlations for {metric}: {str(e)}")
                correlations[metric] = {
                    'pearson_corr': np.nan,
                    'pearson_p': np.nan,
                    'spearman_corr': np.nan,
                    'spearman_p': np.nan,
                    'n_samples': len(clean_df)
                }
            
        return correlations

    def plot_correlations(self, df, correlations, single_chain_only=False, msa_option=None):
        """Create plots for each metric vs pH"""
        # Update plot_creator with current parameters
        self.plot_creator = PlotCreator('pH_correlation', msa_option, single_chain_only)
        metrics = [
            'Complex_RMSD', 'RNA_RMSD', 'RNA_LDDT', 'RNA_TM', 'Complex_TM',
            'Complex_LDDT', 'af3_ranking_score', 'af3_rna_pLDDT_avg',
            'af3_chain_pair_pae_min', 'RNA_DI', 'RNAmotif_score',
            'pred_total_contacts', 'pred_his_contacts', 'pred_acidic_contacts', 
            'pred_basic_contacts', 'pred_contact_score', 'pred_avg_distance',
            'exp_total_contacts', 'exp_his_contacts', 'exp_acidic_contacts',
            'exp_basic_contacts', 'exp_contact_score', 'exp_avg_distance'
        ]
        
        # Set up the plotting style
        plt.style.use('seaborn-v0_8-whitegrid')
        
        for metric in metrics:
            # Create clean subset for plotting
            plot_df = df[['exp_pH', metric]].dropna()
            
            if len(plot_df) < 2:
                print(f"\nSkipping plot for {metric} - insufficient data")
                continue
                
            plt.figure(figsize=(10, 6))
            
            # Create scatter plot
            sns.scatterplot(data=plot_df, x='exp_pH', y=metric, alpha=0.6)
            
            # Add trend line
            sns.regplot(data=plot_df, x='exp_pH', y=metric, scatter=False, 
                       line_kws={'color': 'red', 'alpha': 0.5})
            
            # Add correlation information
            corr_info = correlations[metric]
            plt.title(f'{metric} vs pH\n'
                     f'Pearson: {corr_info["pearson_corr"] if not np.isnan(corr_info["pearson_corr"]) else "N/A"} (p={corr_info["pearson_p"] if not np.isnan(corr_info["pearson_p"]) else "N/A"})')
            
            plt.xlabel('Experimental pH')
            plt.ylabel(metric)
            
            # Save plot using PlotCreator
            self.plot_creator.get_scatterplot(
                table_source='ph_correlation',
                xAxis_score=plot_df['exp_pH'].tolist(),
                xAxis_label="Experimental pH",
                name=f"pH_{metric}",
                yAxis_label=metric,
                yAxis_score=plot_df[metric].tolist()
            )
            
            plt.close()

        # Create box plots for contact metrics by pH category
        contact_metrics = [
            'pred_total_contacts', 'pred_his_contacts', 'pred_acidic_contacts', 
            'pred_basic_contacts', 'pred_contact_score', 'pred_avg_distance',
            'exp_total_contacts', 'exp_his_contacts', 'exp_acidic_contacts',
            'exp_basic_contacts', 'exp_contact_score', 'exp_avg_distance'
        ]
        
        for metric in contact_metrics:
            plt.figure(figsize=(10, 6))
            sns.boxplot(data=df, x='pH_category', y=metric)
            plt.title(f'{metric} by pH Category')
            plt.xlabel('pH Category')
            plt.ylabel(metric)
            
            # Save plot using PlotCreator
            self.plot_creator.get_scatterplot(
                table_source='ph_correlation',
                xAxis_score=df['pH_category'].tolist(),
                xAxis_label="pH Category",
                name=f"pH_category_{metric}",
                yAxis_label=metric,
                yAxis_score=df[metric].tolist()
            )
            
            plt.close()

    def plot_pred_vs_exp(self, df, single_chain_only=False, msa_option=None):
        """Create scatter plots comparing predicted vs experimental metrics"""
        # Update plot_creator with current parameters
        self.plot_creator = PlotCreator('pH_correlation', msa_option, single_chain_only)
        # Define pH ranges for coloring
        df['pH_range'] = pd.cut(df['exp_pH'], 
                               bins=[0, 5, 6, 7, 8, 14],
                               labels=['<5', '5-6', '6-7', '7-8', '>8'])
        
        # Define metrics to compare
        metrics = [
            ('total_contacts', 'Total Contacts'),
            ('his_contacts', 'Histidine Contacts'),
            ('acidic_contacts', 'Acidic Contacts'),
            ('basic_contacts', 'Basic Contacts'),
            ('contact_score', 'Contact Score'),
            ('avg_distance', 'Average Distance')
        ]
        
        # Set up the plotting style
        plt.style.use('seaborn-v0_8-whitegrid')
        
        # Create plots directory if it doesn't exist
        os.makedirs('plots', exist_ok=True)
        
        for metric, title in metrics:
            pred_col = f'pred_{metric}'
            exp_col = f'exp_{metric}'
            
            # Create clean subset for plotting
            plot_df = df[[pred_col, exp_col, 'exp_pH']].dropna()
            
            if len(plot_df) < 2:
                print(f"\nSkipping plot for {metric} - insufficient data")
                continue
            
            # Print some debug information
            print(f"\nCreating plot for {metric}")
            print(f"Number of data points: {len(plot_df)}")
            print(f"Experimental range: {plot_df[exp_col].min():.2f} to {plot_df[exp_col].max():.2f}")
            print(f"Predicted range: {plot_df[pred_col].min():.2f} to {plot_df[pred_col].max():.2f}")
            
            # Create figure for continuous pH coloring
            plt.figure(figsize=(10, 8))

            # Create scatter plot with pH-based coloring
            scatter = plt.scatter(plot_df[exp_col], plot_df[pred_col],
                                c=plot_df['exp_pH'], cmap='viridis',
                                alpha=0.6, s=50)

            # Add colorbar
            cbar = plt.colorbar(scatter)
            cbar.set_label('pH')

            # Add perfect correlation line
            min_val = min(plot_df[exp_col].min(), plot_df[pred_col].min())
            max_val = max(plot_df[exp_col].max(), plot_df[pred_col].max())
            plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5,
                    label='Perfect Correlation')

            # Calculate correlation
            corr = stats.pearsonr(plot_df[exp_col], plot_df[pred_col])

            # Create figure for discrete pH ranges
            plt.figure(figsize=(10, 8))

            # Create scatter plot with discrete pH ranges
            for pH_range in sorted(df['pH_range'].unique()):
                mask = df['pH_range'] == pH_range
                if mask.any():  # Only plot if there are points in this pH range
                    plt.scatter(df.loc[mask, exp_col], df.loc[mask, pred_col],
                              label=f'pH {pH_range}', alpha=0.6, s=50)

            # Add perfect correlation line
            plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5,
                    label='Perfect Correlation')

            plt.title(f'Predicted vs Experimental {title}\n'
                     f'Pearson correlation: {corr[0]:.3f} (p={corr[1]:.3f})\n'
                     f'N={len(plot_df)}')

            plt.xlabel(f'Experimental {title}')
            plt.ylabel(f'Predicted {title}')

            # Add legend
            plt.legend()

            # Make plot square
            plt.axis('square')

            # Save plot using PlotCreator
            self.plot_creator.get_scatterplot(
                table_source='ph_correlation',
                xAxis_score=df.loc[mask, exp_col].tolist(),
                xAxis_label=f"Experimental {title}",
                name=f"pred_vs_exp_{metric}_discrete",
                yAxis_label=f"Predicted {title}",
                yAxis_score=df.loc[mask, pred_col].tolist()
            )

            plt.close()
            
            print(f"Created plots for {metric}")

    def analyze_region_matching(self, df, single_chain_only=False, msa_option=None):
        """Analyze how well predicted residues match experimental regions"""
        # Update plot_creator with current parameters
        self.plot_creator = PlotCreator('pH_correlation', msa_option, single_chain_only)
        def calculate_residue_match(pred_contacts, exp_contacts):
            """Calculate how well predicted residues match experimental residues"""
            if not pred_contacts or not exp_contacts:
                return 0.0
            
            pred_set = set((c['res_name'], c['res_num']) for c in pred_contacts)
            exp_set = set((c['res_name'], c['res_num']) for c in exp_contacts)

            # Flexible matching: name matches and index is within ±1
            matched = set()
            for pname, pnum in pred_set:
                for ename, enum in exp_set:
                    if pname == ename and abs(pnum - enum) <= 1:
                        matched.add((pname, pnum))
                        break  # Only count each pred residue once

            # Union for Jaccard: all unique (name, num) in both sets
            union = pred_set | exp_set
            return len(matched) / len(union) if union else 0.0

        def calculate_spatial_proximity(pred_contacts, exp_contacts):
            """Calculate spatial proximity between predicted and experimental residues"""
            if not pred_contacts or not exp_contacts:
                return 0.0
            
            # Extract residue numbers and positions
            pred_res = [(c['res_num'], c['rna_pos']) for c in pred_contacts]
            exp_res = [(c['res_num'], c['rna_pos']) for c in exp_contacts]
            
            if not pred_res or not exp_res:
                return 0.0
            
            # Calculate distances between each predicted residue and its nearest experimental residue
            total_score = 0
            for pred_num, pred_rna_pos in pred_res:
                min_distance = float('inf')
                for exp_num, exp_rna_pos in exp_res:
                    # Calculate Euclidean distance in 2D space (residue number, RNA position)
                    distance = ((pred_num - exp_num) ** 2 + (pred_rna_pos - exp_rna_pos) ** 2) ** 0.5
                    min_distance = min(min_distance, distance)
                
                # Convert distance to a score between 0 and 1
                # Using exponential decay: score = e^(-distance/scale)
                # scale parameter controls how quickly the score decays with distance
                scale = 10.0  # Adjust this value to control the rate of decay
                score = np.exp(-min_distance / scale)
                total_score += score
            
            # Normalize by number of predicted residues
            return total_score / len(pred_res)

        def calculate_rna_pos_match(pred_contacts, exp_contacts):
            """Calculate how well predicted RNA positions match experimental RNA positions"""
            if not pred_contacts or not exp_contacts:
                return 0.0
            
            # Create sets of (rna_chain, rna_pos) tuples for both predicted and experimental
            pred_pos = set((c['rna_chain'], c['rna_pos']) for c in pred_contacts)
            exp_pos = set((c['rna_chain'], c['rna_pos']) for c in exp_contacts)
            
            # Calculate match percentage
            if not pred_pos or not exp_pos:
                return 0.0
            
            # Calculate intersection and union
            intersection = len(pred_pos & exp_pos)
            union = len(pred_pos | exp_pos)
            
            # Calculate Jaccard similarity
            return intersection / union if union > 0 else 0.0

        def calculate_res_name_match(pred_contacts, exp_contacts):
            """Calculate how well predicted residue names match experimental residue names"""
            if not pred_contacts or not exp_contacts:
                return 0.0
            
            # Create sets of just res_name for both predicted and experimental
            pred_names = set(c['res_name'] for c in pred_contacts)
            exp_names = set(c['res_name'] for c in exp_contacts)
            
            # Calculate match percentage
            if not pred_names or not exp_names:
                return 0.0
            
            # Calculate intersection and union
            intersection = len(pred_names & exp_names)
            union = len(pred_names | exp_names)
            
            # Calculate Jaccard similarity
            return intersection / union if union > 0 else 0.0

        # Calculate all metrics for each sample
        df['residue_match'] = df.apply(
            lambda row: calculate_residue_match(row['pred_contact_data'], row['exp_contact_data']),
            axis=1
        )
        
        df['spatial_proximity'] = df.apply(
            lambda row: calculate_spatial_proximity(row['pred_contact_data'], row['exp_contact_data']),
            axis=1
        )
        
        df['rna_pos_match'] = df.apply(
            lambda row: calculate_rna_pos_match(row['pred_contact_data'], row['exp_contact_data']),
            axis=1
        )

        df['res_name_match'] = df.apply(
            lambda row: calculate_res_name_match(row['pred_contact_data'], row['exp_contact_data']),
            axis=1
        )
        
        # Create pH ranges
        df['pH_range'] = pd.cut(df['exp_pH'], 
                               bins=[0, 5, 6, 7, 8, 14],
                               labels=['<5', '5-6', '6-7', '7-8', '>8'])
        
        # Calculate percentage of samples with at least one residue match
        df['has_any_match'] = df['residue_match'] > 0
        any_match_stats = df.groupby('pH_range', observed=True)['has_any_match'].agg(['mean', 'std', 'count']).reset_index()
        
        # Create bar plot for any match using PlotCreator
        data_dict = dict(zip(any_match_stats['pH_range'], any_match_stats['mean'] * 100))  # Convert to percentage
        self.plot_creator.get_barplot(
            table_source='ph_correlation',
            data=data_dict,
            title='Percentage of Samples with Any Residue Match by pH Range',
            xlabel='pH Range',
            ylabel='Percentage of Samples with Any Match'
        )
        
        # Create scatter plots for region match vs RMSD
        rmsd_metrics = [
            ('Complex_RMSD', 'Complex RMSD'),
            ('RNA_RMSD', 'RNA RMSD')
        ]
        
        match_metrics = [
            ('residue_match', 'Residue Match'),
            ('rna_pos_match', 'RNA Position Match'),
            ('res_name_match', 'Residue Name Match')
        ]
        
        for rmsd_metric, rmsd_title in rmsd_metrics:
            for match_metric, match_title in match_metrics:
                # Create clean subset for plotting
                plot_df = df[[rmsd_metric, match_metric, 'exp_pH']].dropna()
                
                if len(plot_df) < 2:
                    print(f"\nSkipping plot for {match_metric} vs {rmsd_metric} - insufficient data")
                    continue
                
                # Define pH ranges
                plot_df['pH_range'] = pd.cut(plot_df['exp_pH'], 
                                            bins=[0, 5, 6, 7, 8, 14],
                                            labels=['<5', '5-6', '6-7', '7-8', '>8'])
                
                # Calculate correlation
                corr = stats.pearsonr(plot_df[rmsd_metric], plot_df[match_metric])
                
                # Create plot using new method
                self.plot_creator.get_colour_scatterplot(
                    table_source='ph_correlation',
                    xAxis_score=plot_df[rmsd_metric].tolist(),
                    xAxis_label=rmsd_title,
                    name=f"{match_metric}_vs_{rmsd_metric.lower()}",
                    yAxis_label=f"{match_title} Score",
                    yAxis_score=plot_df[match_metric].tolist(),
                    pH_ranges=plot_df['pH_range'].tolist()
                )
                
                # Print some debug information
                print(f"\nCorrelation analysis for {match_title} vs {rmsd_title}:")
                print(f"Number of samples: {len(plot_df)}")
                print(f"Range of {rmsd_title}: {plot_df[rmsd_metric].min():.2f} to {plot_df[rmsd_metric].max():.2f}")
                print(f"Range of {match_title}: {plot_df[match_metric].min():.2f} to {plot_df[match_metric].max():.2f}")
                print(f"Pearson correlation: {corr[0]:.3f} (p={corr[1]:.3f})")
        
        # Create plots for all metrics vs pH
        metrics = [
            ('residue_match', 'Residue Match'),
            ('spatial_proximity', 'Spatial Proximity'),
            ('rna_pos_match', 'RNA Position Match'),
            ('res_name_match', 'Residue Name Match')
        ]
        
        for metric, title in metrics:
            # Define pH ranges
            df['pH_range'] = pd.cut(df['exp_pH'], 
                                   bins=[0, 5, 6, 7, 8, 14],
                                   labels=['<5', '5-6', '6-7', '7-8', '>8'])
            
            # Calculate correlation
            corr = stats.pearsonr(df['exp_pH'], df[metric])
            
            # Create plot using new method
            self.plot_creator.get_colour_scatterplot(
                table_source='ph_correlation',
                xAxis_score=df['exp_pH'].tolist(),
                xAxis_label="Experimental pH",
                name=f"{metric}_vs_ph",
                yAxis_label=f"{title} Score",
                yAxis_score=df[metric].tolist(),
                pH_ranges=df['pH_range'].tolist()
            )
            
            # Print some debug information
            print(f"\nCorrelation analysis for {title}:")
            print(f"Number of samples: {len(df)}")
            print(f"Range of pH: {df['exp_pH'].min():.2f} to {df['exp_pH'].max():.2f}")
            print(f"Range of {title}: {df[metric].min():.2f} to {df[metric].max():.2f}")
            print(f"Pearson correlation: {corr[0]:.3f} (p={corr[1]:.3f})")
        
        # Create pH-color coded scatter plots for match scores vs various metrics
        print("\nCreating pH-color coded scatter plots for match scores...")
        
        # Define pH ranges for coloring
        df['pH_range'] = pd.cut(df['exp_pH'], 
                               bins=[0, 5, 6, 7, 8, 14],
                               labels=['<5', '5-6', '6-7', '7-8', '>8'])
        
        # Match score vs metric combinations
        match_score_combinations = [
            ('RNAmatch_score', 'Complex_RMSD', 'RNA Match Score vs Complex RMSD'),
            ('RNAmatch_score', 'Complex_TM', 'RNA Match Score vs Complex TM'),
            ('RNAmatch_score', 'Complex_LDDT', 'RNA Match Score vs Complex LDDT'),
            ('AAmatch_score', 'Complex_TM', 'Residue Match vs Complex TM'),
            ('AAmatch_score', 'Complex_LDDT', 'Residue Match vs Complex LDDT'),
            ('AAmatch_score', 'RNAmatch_score', 'AA Residue Match Score vs RNA Match Score')
        ]
        
        for match_score, metric, title in match_score_combinations:
            # Create clean subset for plotting
            plot_df = df[[match_score, metric, 'exp_pH']].dropna()
            
            if len(plot_df) < 2:
                print(f"\nSkipping plot for {title} - insufficient data")
                continue
            
            # Create pH_range for this subset
            plot_df['pH_range'] = pd.cut(plot_df['exp_pH'], 
                                        bins=[0, 5, 6, 7, 8, 14],
                                        labels=['<5', '5-6', '6-7', '7-8', '>8'])
            
            # Calculate correlation
            corr = stats.pearsonr(plot_df[match_score], plot_df[metric])
            
            # Create pH-color coded scatter plot
            self.plot_creator.get_colour_scatterplot(
                table_source='ph_correlation',
                xAxis_score=plot_df[metric].tolist(),
                xAxis_label=metric,
                name=f"{match_score}_vs_{metric.lower()}",
                yAxis_label=f"{match_score}",
                yAxis_score=plot_df[match_score].tolist(),
                pH_ranges=plot_df['pH_range'].tolist(),
                lim=(-0.02, 1.02)
            )
            
            # Print debug information
            print(f"\nCorrelation analysis for {title}:")
            print(f"Number of samples: {len(plot_df)}")
            print(f"Range of {metric}: {plot_df[metric].min():.2f} to {plot_df[metric].max():.2f}")
            print(f"Range of {match_score}: {plot_df[match_score].min():.2f} to {plot_df[match_score].max():.2f}")
            print(f"Pearson correlation: {corr[0]:.3f} (p={corr[1]:.3f})")
        
        return any_match_stats

    def analyze_residue_occurrences(self, df, single_chain_only=False, msa_option=None):
        """Analyze and plot residue name occurrences in predicted vs experimental contacts"""
        # Update plot_creator with current parameters
        self.plot_creator = PlotCreator('pH_correlation', msa_option, single_chain_only)
        def get_residue_counts(contacts):
            """Count occurrences of each residue type in contacts"""
            if not contacts:
                return {}
            counts = {}
            for contact in contacts:
                res_name = contact['res_name']
                counts[res_name] = counts.get(res_name, 0) + 1
            return counts

        # Calculate residue counts for each sample
        df['pred_res_counts'] = df['pred_contact_data'].apply(get_residue_counts)
        df['exp_res_counts'] = df['exp_contact_data'].apply(get_residue_counts)

        # Get all unique residue types
        all_res_types = set()
        for counts in df['pred_res_counts']:
            all_res_types.update(counts.keys())
        for counts in df['exp_res_counts']:
            all_res_types.update(counts.keys())

        # Calculate average counts per pH range
        pH_ranges = ['<5', '5-6', '6-7', '7-8', '>8']
        stats = {}

        for res_type in sorted(all_res_types):
            stats[res_type] = {
                'pred': {},
                'exp': {}
            }
            
            # Calculate average counts for each pH range
            for pH_range in pH_ranges:
                mask = df['pH_range'] == pH_range
                if mask.any():
                    # Predicted counts
                    pred_counts = [counts.get(res_type, 0) for counts in df.loc[mask, 'pred_res_counts']]
                    stats[res_type]['pred'][pH_range] = {
                        'mean': np.mean(pred_counts),
                        'std': np.std(pred_counts),
                        'count': len(pred_counts)
                    }
                    
                    # Experimental counts
                    exp_counts = [counts.get(res_type, 0) for counts in df.loc[mask, 'exp_res_counts']]
                    stats[res_type]['exp'][pH_range] = {
                        'mean': np.mean(exp_counts),
                        'std': np.std(exp_counts),
                        'count': len(exp_counts)
                    }

        # Print overall statistics
        print("\nResidue Occurrence Statistics:")
        print("=" * 80)
        for res_type in sorted(all_res_types):
            print(f"\n{res_type}:")
            print("Predicted vs Experimental (overall):")
            pred_total = sum(counts.get(res_type, 0) for counts in df['pred_res_counts'])
            exp_total = sum(counts.get(res_type, 0) for counts in df['exp_res_counts'])
            print(f"Predicted: {pred_total} occurrences")
            print(f"Experimental: {exp_total} occurrences")
            print(f"Difference: {pred_total - exp_total} (predicted - experimental)")

        # Create focused comparison for specific residues
        residue_groups = {
            'Charged': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU'],
            'Non-polar': ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'MET', 'PHE', 'TRP', 'CYS'],
            'Polar-uncharged': ['SER', 'THR', 'ASN', 'GLN', 'TYR']
        }
        
        # Prepare data for grouped bar plot
        pH_ranges = ['<5', '5-6', '6-7', '7-8', '>8']
        
        # Create plots for each residue group
        for group_name, focus_residues in residue_groups.items():
            data_dict = {}
            
            # Collect data for each pH range
            for pH_range in pH_ranges:
                pred_means = []
                exp_means = []
                for res_type in focus_residues:
                    if res_type in stats and pH_range in stats[res_type]['pred']:
                        pred_means.append(stats[res_type]['pred'][pH_range]['mean'])
                        exp_means.append(stats[res_type]['exp'][pH_range]['mean'])
                    else:
                        pred_means.append(0)
                        exp_means.append(0)
                
                data_dict[pH_range] = {
                    'pred': pred_means,
                    'exp': exp_means
                }

            # Create and save the grouped bar plot
            self.plot_creator.get_range_barplot(
                table_source='ph_correlation',
                data_dict=data_dict,
                labels=focus_residues,
                title=f'Occurrences of {group_name} Residues: Predicted vs Experimental by pH Range',
                xlabel='Residue Type',
                ylabel='Average Number of Occurrences',
                name=f'{group_name.lower().replace("-", "_")}_residue_occurrences_by_ph'
            )
            
            # Print statistics for this group
            print(f"\n{group_name} Residue Statistics:")
            print("=" * 80)
            for res_type in focus_residues:
                if res_type in stats:
                    print(f"\n{res_type}:")
                    pred_total = sum(counts.get(res_type, 0) for counts in df['pred_res_counts'])
                    exp_total = sum(counts.get(res_type, 0) for counts in df['exp_res_counts'])
                    print(f"Predicted: {pred_total} occurrences")
                    print(f"Experimental: {exp_total} occurrences")
                    print(f"Difference: {pred_total - exp_total} (predicted - experimental)")

    def analyze_mispredictions(self, df, single_chain_only=False, msa_option=None):
        """Analyze the top mispredicted residues and their relationships"""
        # Update plot_creator with current parameters
        self.plot_creator = PlotCreator('pH_correlation', msa_option, single_chain_only)
        def get_residue_counts(contacts):
            if not contacts:
                return {}
            counts = {}
            for contact in contacts:
                res_name = contact['res_name']
                counts[res_name] = counts.get(res_name, 0) + 1
            return counts

        # Filter data based on options
        if single_chain_only:
            df = df[df['exp_db_id'].str.contains('_') == False]  # Only single chain entries
        if msa_option:
            if msa_option == '+MSA':
                df = df[df['exp_db_id'].str.contains('_msa')]  # Only MSA entries
            elif msa_option == '-MSA':
                df = df[~df['exp_db_id'].str.contains('_msa')]  # Exclude MSA entries

        # Calculate residue counts for each sample
        df['pred_res_counts'] = df['pred_contact_data'].apply(get_residue_counts)
        df['exp_res_counts'] = df['exp_contact_data'].apply(get_residue_counts)

        # Calculate differences for each residue type by pH range
        pH_ranges = ['<5', '5-6', '6-7', '7-8', '>8']
        all_res_types = set()
        for counts in df['pred_res_counts']:
            all_res_types.update(counts.keys())
        for counts in df['exp_res_counts']:
            all_res_types.update(counts.keys())

        # Calculate differences for each residue type by pH range
        diff_stats = {}
        for pH_range in pH_ranges:
            mask = df['pH_range'] == pH_range
            if mask.any():
                diff_stats[pH_range] = {}
                for res_type in all_res_types:
                    pred_counts = [counts.get(res_type, 0) for counts in df.loc[mask, 'pred_res_counts']]
                    exp_counts = [counts.get(res_type, 0) for counts in df.loc[mask, 'exp_res_counts']]
                    mean_diff = np.mean([p - e for p, e in zip(pred_counts, exp_counts)])
                    diff_stats[pH_range][res_type] = mean_diff

        # Find top 3 mispredicted residues for each pH range
        top_mispredicted = {}
        for pH_range in pH_ranges:
            if pH_range in diff_stats:
                # Sort by absolute difference
                sorted_diffs = sorted(diff_stats[pH_range].items(), 
                                   key=lambda x: abs(x[1]), 
                                   reverse=True)
                top_mispredicted[pH_range] = sorted_diffs[:3]

        # Create plot for top mispredicted residues
        # Find pH ranges that have data
        available_pH_ranges = [ph for ph in pH_ranges if ph in top_mispredicted and top_mispredicted[ph]]
        
        if available_pH_ranges:
            plt.figure(figsize=(12, 6))
            x = np.arange(len(available_pH_ranges))
            width = 0.2

            # Use the first available pH range as reference
            reference_ph = available_pH_ranges[0]
            for i, (res_type, diff) in enumerate(top_mispredicted[reference_ph]):
                diffs = []
                for ph in available_pH_ranges:
                    if ph in top_mispredicted and i < len(top_mispredicted[ph]):
                        diffs.append(top_mispredicted[ph][i][1])
                    else:
                        diffs.append(0)  # No data for this pH range
                plt.bar(x + i*width, diffs, width, label=f'{res_type}')

            plt.title('Top 3 Mispredicted Residues by pH Range')
            plt.xlabel('pH Range')
            plt.ylabel('Mean Difference (Predicted - Experimental)')
            plt.xticks(x + width, available_pH_ranges)
            plt.legend()
            plt.grid(True, alpha=0.3)

            # Save the plot
            self.plot_creator.save_plot('top_mispredicted_residues_by_ph')
            plt.close()
        else:
            print("No data available for any pH range in misprediction analysis")

        # Analyze relationship between residue match and RNA position match
        plt.figure(figsize=(10, 6))
        plt.scatter(df['residue_match'], df['rna_pos_match'], alpha=0.6)
        plt.title('Relationship between Residue Match and RNA Position Match')
        plt.xlabel('Residue Match Score')
        plt.ylabel('RNA Position Match Score')
        
        # Add trend line
        z = np.polyfit(df['residue_match'], df['rna_pos_match'], 1)
        p = np.poly1d(z)
        plt.plot(df['residue_match'], p(df['residue_match']), "r--", alpha=0.8)
        
        # Calculate correlation
        corr = stats.pearsonr(df['residue_match'], df['rna_pos_match'])
        plt.text(0.05, 0.95, f'Pearson correlation: {corr[0]:.3f} (p={corr[1]:.3f})',
                transform=plt.gca().transAxes)

        # Save the plot
        self.plot_creator.save_plot('residue_vs_rna_position_match')
        plt.close()

        # Analyze relationship between bad residue match and mispredicted residues
        plt.figure(figsize=(12, 6))

        # Calculate average difference for each residue type
        overall_diffs = {}
        for res_type in all_res_types:
            pred_counts = [counts.get(res_type, 0) for counts in df['pred_res_counts']]
            exp_counts = [counts.get(res_type, 0) for counts in df['exp_res_counts']]
            overall_diffs[res_type] = np.mean([abs(p - e) for p, e in zip(pred_counts, exp_counts)])

        # Sort residues by their average difference
        sorted_diffs = sorted(overall_diffs.items(), key=lambda x: x[1], reverse=True)
        top_residues = [res for res, _ in sorted_diffs[:5]]  # Top 5 mispredicted residues

        # Create box plot of residue match scores for samples with high misprediction
        data_to_plot = []
        labels = []

        for res_type in top_residues:
            # Find samples with high misprediction for this residue
            high_mispred = df.apply(lambda row:
                                    abs(row['pred_res_counts'].get(res_type, 0) -
                                        row['exp_res_counts'].get(res_type, 0)) > 1, axis=1)

            if high_mispred.any():
                data_to_plot.append(df.loc[high_mispred, 'residue_match'])
                labels.append(f'{res_type}\n(n={high_mispred.sum()})')

        plt.boxplot(data_to_plot, labels=labels)
        plt.title('Residue Match Scores for Samples with High Misprediction')
        plt.xlabel('Residue Type (number of samples)')
        plt.ylabel('Residue Match Score')
        plt.xticks(rotation=45)
        plt.grid(True, alpha=0.3)

        # Save the plot
        self.plot_creator.save_plot('residue_match_vs_misprediction')
        plt.close()

        # Print statistics
        print("\nMisprediction Analysis:")
        print("=" * 80)
        print(f"\nTotal samples analyzed: {len(df)}")
        for pH_range in pH_ranges:
            print(f"\npH Range {pH_range}:")
            if pH_range in top_mispredicted and top_mispredicted[pH_range]:
                print("Top 3 mispredicted residues:")
                for res_type, diff in top_mispredicted[pH_range]:
                    print(f"{res_type}: {diff:.2f} (predicted - experimental)")
            else:
                print("No samples in this pH range")
        
        print("\nOverall misprediction statistics:")
        for res_type, diff in sorted_diffs[:5]:
            print(f"\n{res_type}:")
            print(f"Average absolute difference: {diff:.2f}")
            print(f"Number of samples with misprediction: {sum(1 for counts in df['pred_res_counts'] if counts.get(res_type, 0) > 0)}")

    def compare_lengths(self, pred_length, exp_length):
        """Compare lengths between predicted and experimental data"""
        try:
            # Try to parse as lists
            pred_list = ast.literal_eval(pred_length)
            exp_list = ast.literal_eval(exp_length)
            
            # If both are lists, compare sets (order doesn't matter)
            if isinstance(pred_list, list) and isinstance(exp_list, list):
                return set(pred_list) == set(exp_list)
            
            # If one is list and other isn't, they don't match
            return False
            
        except (SyntaxError, ValueError):
            # If parsing fails, compare as single values
            try:
                return int(pred_length) == int(exp_length)
            except (ValueError, TypeError):
                return False

    def get_protein_sequence_for_contacts(self, row, contact_set_idx, is_predicted=True):
        """
        Given a row, contact set index, and whether it's predicted or experimental,
        return the correct protein sequence for alignment.
        """
        # Get the relevant columns
        protein_chain_ids = row['ProteinChainIDs']
        protein_sequences = row['ProteinSequence']
        chain_id_pair_list = row['ChainIDpairList_proteinRNA']
        
        # Parse if necessary
        if isinstance(protein_chain_ids, str) and protein_chain_ids.strip().startswith('['):
            protein_chain_ids = ast.literal_eval(protein_chain_ids)
        if isinstance(protein_sequences, str) and protein_sequences.strip().startswith('['):
            protein_sequences = ast.literal_eval(protein_sequences)
        if isinstance(chain_id_pair_list, str) and chain_id_pair_list.strip().startswith('['):
            chain_id_pair_list = ast.literal_eval(chain_id_pair_list)
        
        # If only one protein chain, return the only sequence
        if isinstance(protein_chain_ids, str):
            return protein_sequences[0] if isinstance(protein_sequences, list) else protein_sequences
        if isinstance(protein_chain_ids, list) and len(protein_chain_ids) == 1:
            return protein_sequences[0] if isinstance(protein_sequences, list) else protein_sequences
        
        # Otherwise, find the protein chain for this contact set
        if isinstance(chain_id_pair_list, list) and len(chain_id_pair_list) > contact_set_idx:
            protein_chain = chain_id_pair_list[contact_set_idx][0]  # protein chain is always first in the pair
            if protein_chain in protein_chain_ids:
                idx = protein_chain_ids.index(protein_chain)
                return protein_sequences[idx]
        # Fallback: return the longest sequence
        if isinstance(protein_sequences, list):
            return max(protein_sequences, key=len)
        return protein_sequences

    def parse_motifs(self, motif_str):
        """Parse AAMotif into a list of motifs, each as a list of (res_name, res_num)."""
        if not motif_str or motif_str == 'None':
            return []
        try:
            motif_str = motif_str.strip()
            if motif_str.startswith('['):
                motifs = ast.literal_eval(motif_str)
                parsed = []
                for motif in motifs:
                    motif_list = []
                    for part in motif.split(','):
                        match = re.match(r'([A-Z]+)\((\d+)\)', part.strip())
                        if match:
                            res_name, res_num = match.groups()
                            motif_list.append((res_name, int(res_num)))
                    parsed.append(motif_list)
                return parsed
            else:
                # Single motif string, not a list
                motif_list = []
                for part in motif_str.split(','):
                    match = re.match(r'([A-Z]+)\((\d+)\)', part.strip())
                    if match:
                        res_name, res_num = match.groups()
                        motif_list.append((res_name, int(res_num)))
                return [motif_list]
        except Exception as e:
            print(f"Error parsing motif: {motif_str}: {e}")
            return []

    def parse_rna_motifs(self, motif_str, modified_residues_nucleotide):
        """Parse RNAMotif into a list of motifs, each as a list of (res_name, res_num), handling modified nucleotides."""
        if not motif_str or motif_str == 'None':
            return []
        try:
            motif_str = motif_str.strip()
            if motif_str.startswith('['):
                motifs = ast.literal_eval(motif_str)
                parsed = []
                for motif in motifs:
                    motif_list = []
                    for part in motif.split(','):
                        # Handle both standard nucleotides and modified nucleotides
                        match = re.match(r'([A-Z0-9]+)\((\d+)\)', part.strip())
                        if match:
                            res_name, res_num = match.groups()
                            # Replace modified nucleotides with standard ones
                            if res_name in modified_residues_nucleotide:
                                res_name = modified_residues_nucleotide[res_name]
                            motif_list.append((res_name, int(res_num)))
                    parsed.append(motif_list)
                return parsed
            else:
                # Single motif string, not a list
                motif_list = []
                for part in motif_str.split(','):
                    match = re.match(r'([A-Z0-9]+)\((\d+)\)', part.strip())
                    if match:
                        res_name, res_num = match.groups()
                        # Replace modified nucleotides with standard ones
                        if res_name in modified_residues_nucleotide:
                            res_name = modified_residues_nucleotide[res_name]
                        motif_list.append((res_name, int(res_num)))
                return [motif_list]
        except Exception as e:
            print(f"Error parsing RNA motif: {motif_str}: {e}")
            return []

    def get_aamotifs(self, row, is_predicted=True):
        """Return a list of motifs (list of (res_name, res_num)) for each motif in the row, using AAMotif only."""
        motif_col = 'AAMotif'
        motif_str = row[motif_col] if motif_col in row and row[motif_col] not in [None, 'None', ''] else None
        return self.parse_motifs(motif_str)

    def calculate_residue_match_simple(self, pred_motif, exp_motif):
        """Calculate match score for two motifs (lists of (res_name, res_num))."""
        if not pred_motif or not exp_motif:
            return 0.0
        pred_set = set(pred_motif)
        exp_set = set(exp_motif)
        matched = set()
        for pname, pnum in pred_set:
            for ename, enum in exp_set:
                if pname == ename and abs(pnum - enum) <= 1:
                    matched.add((pname, pnum))
                    break
        union = pred_set | exp_set
        return len(matched) / len(union) if union else 0.0

    def motif_triplet_offset_match(self, pred_motif, exp_motif):
        # pred_motif, exp_motif: list of (res_name, res_num)
        if len(pred_motif) < 3 or len(exp_motif) < 3:
            return 0.0, exp_motif  # can't do triplet match
        pred_names = [n for n, _ in pred_motif]
        pred_nums = [num for _, num in pred_motif]
        exp_names = [n for n, _ in exp_motif]
        exp_nums = [num for _, num in exp_motif]
        for i in range(len(exp_names) - 2):
            triplet1 = exp_names[i:i+3]
            for j in range(len(pred_names) - 2):
                triplet2 = pred_names[j:j+3]
                if triplet1 == triplet2:
                    d1_1, d1_2 = exp_nums[i+1] - exp_nums[i], exp_nums[i+2] - exp_nums[i+1]
                    d2_1, d2_2 = pred_nums[j+1] - pred_nums[j], pred_nums[j+2] - pred_nums[j+1]
                    if d1_1 == d2_1 and d1_2 == d2_2:
                        offset = pred_nums[j] - exp_nums[i]
                        # Remap exp_motif
                        new_exp = [(n, num + offset) for n, num in exp_motif]
                        # Only keep if all new indices are positive
                        if all(num > 0 for _, num in new_exp):
                            score = self.calculate_residue_match_simple(pred_motif, new_exp)
                            return score, new_exp
        return 0.0, exp_motif

    def rna_motif_triplet_offset_match(self, pred_motif, exp_motif):
        # pred_motif, exp_motif: list of (res_name, res_num) for RNA
        if len(pred_motif) < 3 or len(exp_motif) < 3:
            return 0.0, exp_motif  # can't do triplet match
        pred_names = [n for n, _ in pred_motif]
        pred_nums = [num for _, num in pred_motif]
        exp_names = [n for n, _ in exp_motif]
        exp_nums = [num for _, num in exp_motif]
        for i in range(len(exp_names) - 2):
            triplet1 = exp_names[i:i+3]
            for j in range(len(pred_names) - 2):
                triplet2 = pred_names[j:j+3]
                if triplet1 == triplet2:
                    d1_1, d1_2 = exp_nums[i+1] - exp_nums[i], exp_nums[i+2] - exp_nums[i+1]
                    d2_1, d2_2 = pred_nums[j+1] - pred_nums[j], pred_nums[j+2] - pred_nums[j+1]
                    if d1_1 == d2_1 and d1_2 == d2_2:
                        offset = pred_nums[j] - exp_nums[i]
                        # Remap exp_motif
                        new_exp = [(n, num + offset) for n, num in exp_motif]
                        # Only keep if all new indices are positive
                        if all(num > 0 for _, num in new_exp):
                            score = self.calculate_residue_match_simple(pred_motif, new_exp)
                            return score, new_exp
        return 0.0, exp_motif

    def analyze(self, single_chain_only=False, msa_option=None):
        """Main analysis function"""
        # Update AAmatch_score in database first
        print("Updating AAmatch_score in database...")
        self.update_aa_match_scores()

        # Update RNAmatch_score in database
        print("Updating RNAmatch_score in database...")
        self.update_rna_match_scores()
        
        # Get data with filters applied
        df = self.get_data(single_chain_only, msa_option)
        
        # Print some debug information about the data
        print("\nData Overview:")
        print("=" * 80)
        print(f"Total number of samples: {len(df)}")
        print("\nColumns in dataframe:")
        print(df.columns.tolist())
        print("\nSample of data:")
        print(df[['exp_pH', 'pred_his_contacts', 'exp_his_contacts']].head())
        
        # Calculate correlations
        correlations = self.calculate_correlations(df)
        
        # Print correlation results
        print("\nCorrelation Analysis Results:")
        print("=" * 80)
        for metric, corr_data in correlations.items():
            print(f"\n{metric}:")
            print(f"Number of samples: {corr_data['n_samples']}")
            print(f"Pearson correlation: {corr_data['pearson_corr']:.3f} (p={corr_data['pearson_p']:.3f})")
            print(f"Spearman correlation: {corr_data['spearman_corr']:.3f} (p={corr_data['spearman_p']:.3f})")
        
        # Create plots
        self.plot_correlations(df, correlations, single_chain_only, msa_option)
        
        # Create predicted vs experimental comparison plots
        print("\nCreating predicted vs experimental comparison plots...")
        self.plot_pred_vs_exp(df, single_chain_only, msa_option)
        
        # Filter samples based on matching lengths
        print("\nFiltering samples based on matching lengths...")
        print(f"Initial number of samples: {len(df)}")
        
        try:
            # Get protein and RNA lengths from database
            pred_query = """
            SELECT exp_db_id, ProteinLength, RNALength
            FROM pred_protein_rna
            """
            exp_query = """
            SELECT PDBId, ProteinLength, RNALength
            FROM exp_protein_rna
            """
            
            pred_lengths = pd.read_sql_query(pred_query, self.connection)
            exp_lengths = pd.read_sql_query(exp_query, self.connection)
            
            # Merge length information
            df = pd.merge(df, pred_lengths, on='exp_db_id', how='left')
            df = pd.merge(df, exp_lengths, on='PDBId', how='left', suffixes=('_pred', '_exp'))
            
            # Track failed comparisons
            failed_protein_lengths = []
            failed_rna_lengths = []
            
            # Filter based on matching lengths
            length_matches = []
            for idx, row in df.iterrows():
                protein_match = self.compare_lengths(row['ProteinLength_pred'], row['ProteinLength_exp'])
                rna_match = self.compare_lengths(row['RNALength_pred'], row['RNALength_exp'])
                
                if not protein_match:
                    failed_protein_lengths.append((row['exp_db_id'], row['ProteinLength_pred'], row['ProteinLength_exp']))
                if not rna_match:
                    failed_rna_lengths.append((row['exp_db_id'], row['RNALength_pred'], row['RNALength_exp']))
                
                length_matches.append(protein_match and rna_match)
            
            df = df[length_matches]
            print(f"Number of samples after length filtering: {len(df)}")
            
            # Print failed comparisons
            if failed_protein_lengths:
                print("\nFailed Protein Length Comparisons:")
                print("PDB ID | Predicted Length | Experimental Length")
                print("-" * 50)
                for pdb_id, pred_len, exp_len in failed_protein_lengths:
                    print(f"{pdb_id} | {pred_len} | {exp_len}")
            
            if failed_rna_lengths:
                print("\nFailed RNA Length Comparisons:")
                print("PDB ID | Predicted Length | Experimental Length")
                print("-" * 50)
                for pdb_id, pred_len, exp_len in failed_rna_lengths:
                    print(f"{pdb_id} | {pred_len} | {exp_len}")
            
        except Exception as e:
            print(f"Warning: Length filtering failed: {str(e)}")
            print("Continuing with unfiltered dataset...")
        
        # Analyze region matching
        print("\nAnalyzing region matching...")
        self.analyze_region_matching(df, single_chain_only, msa_option)
        
        # Analyze residue occurrences
        print("\nAnalyzing residue occurrences...")
        self.analyze_residue_occurrences(df, single_chain_only, msa_option)
        
        # Analyze mispredictions
        print("\nAnalyzing mispredictions...")
        self.analyze_mispredictions(df, single_chain_only, msa_option)
        
        # Close database connection
        self.close_connection()

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) not in [1, 2, 3]:
        print("Usage: python phCorrelation.py [single_chain|all] [+MSA|-MSA]")
        sys.exit(1)

    single_chain_only = sys.argv[1] == 'single_chain' if len(sys.argv) > 1 else False
    msa_option = sys.argv[2] if len(sys.argv) > 2 else None

    if msa_option and msa_option.upper() not in ['+MSA', '-MSA']:
        print("MSA option must be either +MSA or -MSA")
        sys.exit(1)

    analyzer = PHCorrelationAnalyzer(single_chain_only=single_chain_only, msa_option=msa_option)
    analyzer.analyze(single_chain_only, msa_option)