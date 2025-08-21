import os
import sys
import ast
import re
from Bio.Align import PairwiseAligner

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from plots.plotCreator import PlotCreator
from database.databaseMethods import DatabaseMethods

class ProteinRNAInterface(DatabaseMethods):
    def __init__(self, plot_creator=None, id: str = "", protein_rna_interface_area: str = "", protein_rna_r_is: str = "", protein_rna_aap: str = "",
                 aac: str = ""):
        super().__init__()
        self.plot_creator = plot_creator
        self.id = id
        self.protein_rna_interface_area = protein_rna_interface_area
        self.protein_rna_r_is = protein_rna_r_is
        self.protein_rna_aap = protein_rna_aap
        self.aac = aac

    def get_proteinRNA_metrics(self, pred_table_name, exp_table_name):
        common_pdb = self.get_intersect_values(table1=f"{pred_table_name}",
                                                   column1="exp_db_id",
                                                   table2=f"{exp_table_name}",
                                                   column2="PDBId"
                                                   )
        
        filtered_pdb_ids = self.plot_creator.get_filtered_pdb_ids(pred_table_name)

        common_pdb_ids = list(set(filtered_pdb_ids) & set(common_pdb))

        results = {}
        columns = [
            "ProteinSequence", "ProteinChainIDs", "RNASequence", "RNAChainIDs",
            "AAC", "AAPproteinRNA", "ChainIDpairList_proteinRNA", "Hbond_proteinRNA",
            "vdWbond_proteinRNA", "ProteinRNAInterfaceArea", "ProteinRNAInterfaceRatio",
            "RNALength", "RNAMotif", "RNAMotifLength"
        ]
        # print("£$%£@$%%£$%£ common_pdb_ids", common_pdb_ids)
        for pdb_id in common_pdb_ids:
            pred_data = self.get_table_columns(pred_table_name, columns, condition=f"exp_db_id = '{pdb_id}'")
            exp_data = self.get_table_columns(exp_table_name, columns, condition=f"PDBId = '{pdb_id}'")

            results[pdb_id] = {
                'predicted': pred_data,
                'experimental': exp_data
            }

        return results


class RnaRNAInterface(DatabaseMethods):
    def __init__(self, id: str = "", rna_rna_interface_area: str = "", rna_rna_r_is: str = ""):
        super().__init__()
        self.id = id
        self.rna_rna_interface_area = rna_rna_interface_area
        self.rna_rna_r_is = rna_rna_r_is

    def get_rnaRNA_metrics(self, pred_table_name, exp_table_name):
        common_pdb_ids = self.get_intersect_values(table1=f"{pred_table_name}",
                                                   column1="exp_db_id",
                                                   table2=f"{exp_table_name}",
                                                   column2="PDBId"
                                                   )
        results = {}
        columns = [
            "RNArnaInterfaceArea", "rnaRNAInterfaceRatio", "RNALength"
        ]

        for pdb_id in common_pdb_ids:
            pred_data = self.get_table_columns(pred_table_name, columns, condition=f"exp_db_id = '{pdb_id}'")
            exp_data = self.get_table_columns(exp_table_name, columns, condition=f"PDBId = '{pdb_id}'")

            results[pdb_id] = {
                'predicted': pred_data,
                'experimental': exp_data
            }

        return results

def get_rna_motif_scores(result, plot_creator=None):
    database = DatabaseMethods()  # Create database connection
    score_length_mapping = {}
    protein_chain_interaction_differences = {
        'interaction_diff': [],
        'max_rna_length': []
    }
    for pdb_id, data in result.items():
        print(f"PDB ID: {pdb_id}")
        (
            pred_protein_seq, pred_protein_chains, pred_rna_seq, pred_rna_chains, 
            pred_aac, pred_aap, pred_chain_id_pairs, pred_h_bonds, pred_vdW_bonds, 
            pred_interface_area, pred_interface_ratio, rna_length, pred_rna_motif, 
            pred_rna_motif_length_list
        ) = data['predicted']
        
        (
            ref_protein_seq, ref_protein_chains, ref_rna_seq, ref_rna_chains, 
            ref_aac, ref_aap, ref_chain_id_pairs, ref_h_bonds, ref_vdW_bonds, 
            ref_interface_area, ref_interface_ratio, ref_rna_length, ref_rna_motif, 
            ref_rna_motif_length_list
        ) = data['experimental']
        
        # print(pred_chain_id_pairs, pred_protein_chains, pred_rna_chains, pred_rna_motif,
        #     ref_chain_id_pairs, ref_protein_chains, ref_rna_chains, ref_rna_motif, rna_length)

        def convert_parameters(pred_chain_id_pairs, pred_protein_chains, pred_rna_chains, pred_rna_motif,
                               ref_chain_id_pairs, ref_protein_chains, ref_rna_chains, ref_rna_motif, rna_length):
            try:
                pred_chain_id_pairs = ast.literal_eval(pred_chain_id_pairs) if pred_chain_id_pairs else []
            except (SyntaxError, ValueError) as e:
                print(f"Error parsing pred_chain_id_pairs: {pred_chain_id_pairs}")
                print(f"Error: {e}")
                pred_chain_id_pairs = []

            if len(pred_protein_chains) == 1:
                pred_protein_chains = [pred_protein_chains]
            else:
                try:
                    pred_protein_chains = ast.literal_eval(pred_protein_chains) if pred_protein_chains else []
                except (SyntaxError, ValueError) as e:
                    print(f"Error parsing pred_protein_chains: {pred_protein_chains}")
                    print(f"Error: {e}")
                    pred_protein_chains = []

            if len(pred_rna_chains) == 1:
                pred_rna_chains = [pred_rna_chains]
            else:
                try:
                    pred_rna_chains = ast.literal_eval(pred_rna_chains) if pred_rna_chains else []
                except (SyntaxError, ValueError) as e:
                    print(f"Error parsing pred_rna_chains: {pred_rna_chains}")
                    print(f"Error: {e}")
                    pred_rna_chains = []

            if pred_rna_motif.startswith('['):
                try:
                    pred_rna_motif = ast.literal_eval(pred_rna_motif) if pred_rna_motif else []
                except (SyntaxError, ValueError) as e:
                    print(f"Error parsing pred_rna_motif: {pred_rna_motif}")
                    print(f"Error: {e}")
                    pred_rna_motif = []
            else:
                pred_rna_motif = [pred_rna_motif]

            # Add debugging for ref_chain_id_pairs
            print(f"Debug - ref_chain_id_pairs before parsing: '{ref_chain_id_pairs}'")
            print(f"Debug - type of ref_chain_id_pairs: {type(ref_chain_id_pairs)}")
            
            try:
                if not ref_chain_id_pairs or ref_chain_id_pairs.strip() == '':
                    print("Warning: ref_chain_id_pairs is empty or whitespace")
                    ref_chain_id_pairs = []
                else:
                    ref_chain_id_pairs = ast.literal_eval(ref_chain_id_pairs)
            except (SyntaxError, ValueError) as e:
                print(f"Error parsing ref_chain_id_pairs: '{ref_chain_id_pairs}'")
                print(f"Error: {e}")
                ref_chain_id_pairs = []

            if len(ref_protein_chains) == 1:
                ref_protein_chains = [ref_protein_chains]
            else:
                try:
                    ref_protein_chains = ast.literal_eval(ref_protein_chains) if ref_protein_chains else []
                except (SyntaxError, ValueError) as e:
                    print(f"Error parsing ref_protein_chains: {ref_protein_chains}")
                    print(f"Error: {e}")
                    ref_protein_chains = []

            if len(ref_rna_chains) == 1:
                ref_rna_chains = [ref_rna_chains]
            else:
                try:
                    ref_rna_chains = ast.literal_eval(ref_rna_chains) if ref_rna_chains else []
                except (SyntaxError, ValueError) as e:
                    print(f"Error parsing ref_rna_chains: {ref_rna_chains}")
                    print(f"Error: {e}")
                    ref_rna_chains = []

            if ref_rna_motif.startswith('['):
                try:
                    ref_rna_motif = ast.literal_eval(ref_rna_motif) if ref_rna_motif else []
                except (SyntaxError, ValueError) as e:
                    print(f"Error parsing ref_rna_motif: {ref_rna_motif}")
                    print(f"Error: {e}")
                    ref_rna_motif = []
            else:
                ref_rna_motif = [ref_rna_motif]

            if isinstance(rna_length, int):
                rna_length = [rna_length]
            else:
                try:
                    rna_length = ast.literal_eval(rna_length) if rna_length else []
                except (SyntaxError, ValueError) as e:
                    print(f"Error parsing rna_length: {rna_length}")
                    print(f"Error: {e}")
                    rna_length = []

            return pred_chain_id_pairs, pred_protein_chains, pred_rna_chains, pred_rna_motif, \
                ref_chain_id_pairs, ref_protein_chains, ref_rna_chains, ref_rna_motif, rna_length

        pred_chain_id_pairs, pred_protein_chains, pred_rna_chains, pred_rna_motif, \
            ref_chain_id_pairs, ref_protein_chains, ref_rna_chains, ref_rna_motif, rna_length \
            = convert_parameters(pred_chain_id_pairs, pred_protein_chains, pred_rna_chains, pred_rna_motif,
                                 ref_chain_id_pairs, ref_protein_chains, ref_rna_chains, ref_rna_motif, rna_length)

        def get_combined_rna_motif(rna_chain, rna_motif, chain_id_pairs):
            if len(rna_motif) < len(chain_id_pairs):
                return rna_motif[0]
            # Find the chain pairs related to the specified RNA chain
            associated_motifs = []
            for i, pair in enumerate(chain_id_pairs):
                if pair[1] == rna_chain:  # The second element in pair is the RNA chain
                    associated_motifs.append(rna_motif[i])
            
            nucleotide_positions = []
            # Updated regex to handle both regular and modified nucleotides
            motif_pattern = re.compile(r'\(?([AUGC]|OMU|OMC|OMG|PSU|2MG|5MC|5MU|UR3|A2M|MA6|6MZ|7MG|RSQ|5CM|C34|5HC|6OG|6MA|1CC|8OG|5FC|3DR)\)?\((\d+)\)')
            
            for motif in associated_motifs:
                for match in motif_pattern.findall(motif):
                    nucleotide, position = match[0], int(match[1])
                    nucleotide_positions.append((nucleotide, position))

            unique_nucleotides = list(set(nucleotide_positions))
            sorted_nucleotides = sorted(unique_nucleotides, key=lambda x: x[1])
            
            # Format the output to match the input format
            combined_motif = ', '.join([f'({nt})({pos})' if nt in ['OMU', 'OMC', 'OMG', 'PSU', '2MG', '5MC', '5MU', 'UR3', 'A2M', 'MA6', '6MZ', '7MG', 'RSQ', '5CM', 'C34', '5HC', '6OG', '6MA', '1CC', '8OG', '5FC', '3DR'] else f'{nt}({pos})' for nt, pos in sorted_nucleotides])
            
            return combined_motif

        # Group Protein-RNA Chain Pairs by RNA Chain
        ref_rna_chain_pairs = {rna_chain: [] for rna_chain in ref_rna_chains}
        for pair in ref_chain_id_pairs:
            protein_chain, rna_chain = pair
            if rna_chain in ref_rna_chain_pairs:
                ref_rna_chain_pairs[rna_chain].append(pair)
        # print("!!!!ref rna_chain_pairs", ref_rna_chain_pairs)

        pred_rna_chain_pairs = {rna_chain: [] for rna_chain in pred_rna_chains}
        for pair in pred_chain_id_pairs:
            protein_chain, rna_chain = pair
            if rna_chain in pred_rna_chain_pairs:
                pred_rna_chain_pairs[rna_chain].append(pair)
        # print("!!!!pred rna_chain_pairs", pred_rna_chain_pairs)

        # Count RNA-Protein Interactions
        ref_interactions = len(ref_chain_id_pairs)
        pred_interactions = len(pred_chain_id_pairs)
        interaction_difference = pred_interactions - ref_interactions
        # print(f"Reference interactions: {ref_interactions}")
        # print(f"Predicted interactions: {pred_interactions}")
        # print(f"Interaction difference: {interaction_difference}")
        max_rna_length = max(rna_length)
        protein_chain_interaction_differences['interaction_diff'].append(interaction_difference)
        protein_chain_interaction_differences['max_rna_length'].append(max_rna_length)
        # print("protein_chain_interaction_differences", protein_chain_interaction_differences)

        # if (pred_protein_chains != ref_protein_chains or pred_rna_chains != ref_rna_chains) and len(pred_rna_chains) > 1:
        #     pred_protein_chains, pred_rna_chains = reorder_pred_chains(ref_protein_chains, ref_rna_chains, pred_protein_chains, pred_rna_chains)
        #     pred_chain_id_pairs = reorder_pred_pairs(pred_chain_id_pairs, ref_protein_chains, ref_rna_chains)
        #     print("REORDER pred_chain_id_pairs", pred_chain_id_pairs)
        #     print("SAME REF_chain_id_pairs", ref_chain_id_pairs)
        #
        # # Group Protein-RNA Chain Pairs by RNA Chain
        # ref_rna_chain_pairs = {rna_chain: [] for rna_chain in ref_rna_chains}
        # for pair in ref_chain_id_pairs:
        #     protein_chain, rna_chain = pair
        #     if rna_chain in ref_rna_chain_pairs:
        #         ref_rna_chain_pairs[rna_chain].append(pair)
        # print("ref rna_chain_pairs", ref_rna_chain_pairs)
        #
        # pred_rna_chain_pairs = {rna_chain: [] for rna_chain in pred_rna_chains}
        # for pair in pred_chain_id_pairs:
        #     protein_chain, rna_chain = pair
        #     if rna_chain in pred_rna_chain_pairs:
        #         pred_rna_chain_pairs[rna_chain].append(pair)
        # print("pred rna_chain_pairs", pred_rna_chain_pairs)

        # for rna_chain, associated_pairs in rna_chain_pairs.items():
        #     print(f"RNA Chain '{rna_chain}': {associated_pairs}")

        def fill_gaps_motif(motif):
            if not motif:
                return ""
        
            # Extract nucleotides and their positions using regular expressions
            # Updated regex to handle both regular and modified nucleotides
            matches = re.findall(r'\(?([AUGC]|OMU|OMC|OMG|PSU|2MG|5MC|5MU|UR3|A2M|MA6|6MZ|7MG|RSQ|5CM|C34|5HC|6OG|6MA|1CC|8OG|5FC|3DR)\)?\((\d+)\)', motif)

            if not matches:
                return ""

            # Initialize the resulting sequence and track the previous position
            result_sequence = ''
            previous_position = None

            # Loop over the nucleotides and their positions
            for nt, pos in matches:
                pos = int(pos)

                # If this isn't the first nucleotide and there is a gap, insert 'N' for each missing position
                if previous_position is not None and pos > previous_position + 1:
                    gap_size = pos - previous_position - 1
                    result_sequence += 'N' * gap_size

                # Add the current nucleotide to the result
                result_sequence += nt

                # Update the previous position
                previous_position = pos

            return result_sequence

        def align_motifs(pred_motif, ref_motif):
            # Fill gaps in both motifs
            pred_str = fill_gaps_motif(pred_motif)
            ref_str = fill_gaps_motif(ref_motif)
            
            # Check if either string is empty
            if not pred_str or not ref_str:
                print(f"Warning: Empty sequence detected. pred_str: '{pred_str}', ref_str: '{ref_str}'")
                return None, 0
        
            max_motif_length = max(len(pred_str), len(ref_str))
            
            # Configure the aligner
            aligner = PairwiseAligner()
            aligner.mode = 'global'  # global alignment (Needleman-Wunsch), (local: Smith-Waterman)
            aligner.match_score = 1  # Score for matching nucleotides
            aligner.mismatch_score = -2  # Penalty for mismatch # -2.0
            aligner.gap_score = -2 # General penalty for gaps inside the sequence.
            aligner.open_gap_score = -2  # Penalty for opening a gap
            aligner.extend_gap_score = -2
            aligner.end_gap_score = -2
            max_alignments = 1000000

            try:
                # Perform the alignment
                alignments = aligner.align(ref_str, pred_str)

                # Check if the number of alignments is unreasonably high
                if len(alignments) > max_alignments:
                    print(f"Warning: Too many alignments generated ({len(alignments)}).")
                    return None, 0

                # Choose the best alignment (first one is usually the best in the sorted list)
                best_alignment = alignments[0] if alignments else None

                return best_alignment, max_motif_length

            except OverflowError:
                print("Error: Number of alignments is too large to handle. Consider changing alignment parameters or input sequences.")
                return None, 0
            except ValueError as e:
                print(f"Error during alignment: {e}")
                print(f"pred_str: '{pred_str}'")
                print(f"ref_str: '{ref_str}'")
                return None, 0

        def score_alignment(alignment):
            aligned_ref = alignment.aligned[0]
            aligned_pred = alignment.aligned[1]

            match_count = 0
            total_count = 0

            for ref_seg, pred_seg in zip(aligned_ref, aligned_pred):
                # Iterate through aligned segments
                ref_start, ref_end = ref_seg
                pred_start, pred_end = pred_seg

                ref_aligned_part = alignment.target[ref_start:ref_end]
                pred_aligned_part = alignment.query[pred_start:pred_end]

                # Count the matches and total aligned nucleotides
                for r, p in zip(ref_aligned_part, pred_aligned_part):
                    if r == p:
                        match_count += 1
                    total_count += 1

            if total_count == 0:
                return 0

            score = match_count / total_count
            return score

        def normalize_score(alignment, max_motif_length):
            # Calculate the length of the longer predicted (query) or reference (target) motifs
            # aligned_pred_str = alignment.aligned[0]
            # aligned_ref_str = alignment.aligned[1]
            # print("aligned_ref_str", aligned_ref_str)
            # print("aligned_pred_str", aligned_pred_str)
            # aligned_pred_length = sum(end - start for start, end in aligned_pred_str)
            # aligned_ref_length = sum(end - start for start, end in aligned_ref_str)
            # alignment_length = max(aligned_pred_length, aligned_ref_length)
            # print("&&&&&alignment_length", aligned_pred_length, aligned_ref_length,alignment_length)
            # Maximum possible score (all matches)
            max_score = max_motif_length * 1
            # Minimum possible score (all gaps)
            min_score = max_motif_length * (-2)
            actual_score = alignment.score
            normalized_score = (actual_score - min_score) / (max_score - min_score)

            return normalized_score

        # ref_longest_rna_motifs = []
        # pred_longest_rna_motifs = []
        # # Extract the Longest Motif for Each RNA Chain
        # for rna_chain in ref_rna_chains:
        #     ref_longest_motif = get_longest_rna_motif(rna_chain, ref_rna_motif, ref_chain_id_pairs)
        #     ref_longest_rna_motifs.append(ref_longest_motif)
        #     print(f"Longest motif for RNA chain {rna_chain}: {ref_longest_motif}")
        #     pred_longest_motif = get_longest_rna_motif(rna_chain, pred_rna_motif, pred_chain_id_pairs)
        #     pred_longest_rna_motifs.append(pred_longest_motif)
        #     print(f"PRED Longest motif for RNA chain {rna_chain}: {pred_longest_motif}")
        #     alignment = align_motifs(pred_longest_motif, ref_longest_motif)
        #     # print("alignment", alignment)
        #     # Check if alignment was found
        #     if alignment:
        #         formatted_alignment, score = score_alignment(pred_longest_motif, ref_longest_motif)  # , alignment)
        #         print(f"PRED Aligned Motif for Chain {rna_chain}:\n{formatted_alignment}\nMotif score: {score}")
        #     else:
        #         # No alignment found, assigning -100% score
        #         print(f"No alignment found for Chain {rna_chain}, assigning -1")
        #         score = -1

        # (rna_chain, ref_rna_motif, ref_chain_id_pairs)
        # print("ref_rna_motif", ref_rna_motif)
        # print("ref_chain_id_pairs", ref_chain_id_pairs)


        def process_rna_motifs(ref_rna_chain_pairs, pred_rna_chain_pairs, ref_rna_motif, pred_rna_motif):
            motif_scores = []
            ref_rna_keys = list(ref_rna_chain_pairs.keys())
            pred_rna_keys = list(pred_rna_chain_pairs.keys())
            num_keys = min(len(ref_rna_keys), len(pred_rna_keys))
            
            print(f"Number of RNA chains to process: {num_keys}")
            print(f"Reference RNA chains: {ref_rna_keys}")
            print(f"Predicted RNA chains: {pred_rna_keys}")
            
            # Handle case where we have only one chain pair
            if len(ref_rna_keys) == 1 and len(pred_rna_keys) == 1:
                ref_rna_chain = ref_rna_keys[0]
                pred_rna_chain = pred_rna_keys[0]
                ref_interactions = ref_rna_chain_pairs[ref_rna_chain]
                pred_interactions = pred_rna_chain_pairs[pred_rna_chain]
                
                print(f"Processing single RNA Chain pair: {ref_rna_chain} (Ref) and {pred_rna_chain} (Pred)")
                print(f"Reference interactions: {ref_interactions}")
                print(f"Predicted interactions: {pred_interactions}")
                
                if not ref_interactions or not pred_interactions:
                    print(f"No interactions found for single chain pair, returning 0")
                    return 0
                    
                ref_combined_motif = get_combined_rna_motif(ref_rna_chain, ref_rna_motif, ref_interactions)
                pred_combined_motif = get_combined_rna_motif(pred_rna_chain, pred_rna_motif, pred_interactions)
                
                print(f"Reference combined motif: {ref_combined_motif}")
                print(f"Predicted combined motif: {pred_combined_motif}")
                
                alignment, max_motif_length = align_motifs(pred_combined_motif, ref_combined_motif)
                print(f"Alignment result: {alignment}")
                
                if alignment:
                    score = normalize_score(alignment, max_motif_length)
                    score = round(score, 2)
                    print(f"Alignment score: {score}")
                    return score
                else:
                    print(f"No alignment found for single chain pair, returning 0")
                    return 0
            
            # Process multiple chain pairs
            for idx in range(num_keys):
                ref_rna_chain = ref_rna_keys[idx]
                pred_rna_chain = pred_rna_keys[idx]
                ref_interactions = ref_rna_chain_pairs[ref_rna_chain]
                pred_interactions = pred_rna_chain_pairs[pred_rna_chain]

                print(f"Processing RNA Chain {ref_rna_chain} (Ref) and {pred_rna_chain} (Pred)")
                print(f"Reference interactions: {ref_interactions}")
                print(f"Predicted interactions: {pred_interactions}")
                
                if not ref_interactions or not pred_interactions:
                    print(f"No comparison can be made between RNA chain {ref_rna_chain} (Ref) and {pred_rna_chain} (Pred)")
                    continue
                    
                ref_combined_motif = get_combined_rna_motif(ref_rna_chain, ref_rna_motif, ref_interactions)
                pred_combined_motif = get_combined_rna_motif(pred_rna_chain, pred_rna_motif, pred_interactions)
                
                print(f"Reference combined motif: {ref_combined_motif}")
                print(f"Predicted combined motif: {pred_combined_motif}")
                
                alignment, max_motif_length = align_motifs(pred_combined_motif, ref_combined_motif)
                print(f"Alignment result: {alignment}")
                
                if alignment:
                    score = normalize_score(alignment, max_motif_length)
                    score = round(score, 2)
                    print(f"Alignment score: {score}")
                    motif_scores.append(score)
                else:
                    print(f"No alignment found for Chain {ref_rna_chain}, assigning 0")
                    score = 0
                    motif_scores.append(score)
            
            print(f"Final motif scores: {motif_scores}")
            
            if not motif_scores:
                print("No valid motif comparisons were made, returning 0")
                return 0
                
            return sum(motif_scores) / len(motif_scores)


        # def process_rna_motifs(pred_rna_motif, ref_rna_motif, ref_rna_chains, pred_rna_chains):
        #     ref_combined_rna_motifs = []
        #     pred_combined_rna_motifs = []
        #     motif_scores = []
        #     # Extract the Longest Motif for Each RNA Chain
        #     for rna_chain in ref_rna_chains:
        #         ref_combined_motif =  (rna_chain, ref_rna_motif, ref_chain_id_pairs)
        #         ref_combined_rna_motifs.append(ref_combined_motif)
        #         print(f"REF combined motif for RNA chain {rna_chain}: {ref_combined_motif}")
        #         for rna_chain in pred_rna_chains:
        #             pred_combined_motif = get_combined_rna_motif(rna_chain, pred_rna_motif, pred_chain_id_pairs)
        #             pred_combined_rna_motifs.append(pred_combined_motif)
        #             print(f"PRED combined motif for RNA chain {rna_chain}: {pred_combined_motif}")
        #             if pred_combined_motif is None:
        #                 break
        #     # print("pred_rna_motif", isinstance(pred_rna_motif, str))
        #     # print("pred_rna_motif", pred_rna_motif)
        #     # print("pred_rna_motif", isinstance(longest_rna_motifs, str))
        #     # print("longest_rna_motifs", longest_rna_motifs)
        #     # for i, (pred_motif, ref_motif) in enumerate(zip(pred_combined_rna_motifs, ref_combined_rna_motifs)):
        #         # print("pred_motif", pred_motif)
        #         # print("longest_rna_motif", longest_rna_motif)
        #         # Alignment using the longest ref RNA motif
        #         alignment = align_motifs(pred_combined_motif, ref_combined_motif)
        #         print("alignment", alignment)
        #         # Check if alignment was found
        #         if alignment:
        #             # Format the alignment and score it
        #             score = score_alignment(alignment)
        #             score = round(score, 2)
        #             # print(f"Aligned Motif for Chain {rna_chain}:\n{formatted_alignment}\nMotif score: {score}")
        #             print(f"Aligned Motif for Chain {rna_chain}:\nMotif score: {score}")
        #             motif_scores.append(score)
        #         else:
        #             print(f"No alignment found for Chain {rna_chain}, assigning -1")
        #             score = -1
        #             motif_scores.append(score)
        #     return sum(motif_scores) / len(motif_scores)

        total_score = process_rna_motifs(ref_rna_chain_pairs, pred_rna_chain_pairs, ref_rna_motif, pred_rna_motif)#(pred_rna_motif, ref_rna_motif, ref_rna_chains, pred_rna_chains)
        print("total_score", total_score)
        max_rna_length = max(rna_length)
        print("max_rna_length", max_rna_length)
        if not isinstance(total_score, list):
            total_score = [total_score]  # Wrap in a list if it's a single float

        if not isinstance(max_rna_length, list):
            max_rna_length = [max_rna_length]
        for score, length in zip(total_score, max_rna_length):
            if length in score_length_mapping:
                # If the length already exists, append the new score to the list
                score_length_mapping[length].append(score)
            else:
                score_length_mapping[length] = [score]
        print("score_length_mapping", score_length_mapping)

        try:
            update_query = """
            UPDATE pred_protein_rna 
            SET RNAmotif_score = ? 
            WHERE exp_db_id = ?
            """
            database.connection.execute(update_query, (score, pdb_id))
            database.connection.commit()
        except Exception as e:
            print(f"Error saving RNAmotif score for {pdb_id}: {e}")

    # Ensure we have a plot_creator
    if plot_creator is None:
        raise ValueError("plot_creator must be provided to get_rna_motif_scores")
        
    sorted_mapping = sorted(score_length_mapping.items()) # (rna_length, total_score)
    rna_lengths = []
    total_scores = []
    for rna_length, scores in sorted_mapping:
        rna_lengths.extend([rna_length] * len(scores))
        total_scores.extend(scores)
        
    # Use the passed plot_creator which has the correct filtering options
    plot_creator.get_scatterplot(
        table_source='motif_metrics',  # Use consistent table_source
        xAxis_score=rna_lengths,
        xAxis_label="RNA Sequence Length",
        name="RNAMotif_Length",
        yAxis_label="RNA Motif Similarity",
        yAxis_score=total_scores
    )
    
    plot_creator.get_scatterplot(
        table_source='motif_metrics',  # Use consistent table_source
        xAxis_score=protein_chain_interaction_differences['max_rna_length'],
        xAxis_label="RNA Sequence Length",
        name="Interacting_ProteinChains",
        yAxis_label="Interaction Difference",
        yAxis_score=protein_chain_interaction_differences['interaction_diff']
    )

if __name__ == "__main__":
    if len(sys.argv) not in [2, 3, 4]:
        print("Usage: python motifMetrics.py <table_name> [single_chain|all] [+MSA|-MSA]")
        sys.exit(1)

    table = sys.argv[1]
    single_chain_only = sys.argv[2] == 'single_chain' if len(sys.argv) > 2 else False
    msa_option = sys.argv[3] if len(sys.argv) > 3 else None

    # Create plot creator first
    plot_creator = PlotCreator('motif_metrics', msa_option, single_chain_only)

    if table == "pred_protein_rna":
        # Pass plot_creator to metrics class
        protein_rna_metrics = ProteinRNAInterface(plot_creator=plot_creator)
        result = protein_rna_metrics.get_proteinRNA_metrics("pred_protein_rna", "exp_protein_rna")
        get_rna_motif_scores(result, plot_creator)
        protein_rna_metrics.close_connection()
    elif table == "pred_protein_rna_dna":
        # Pass plot_creator to metrics class
        protein_rna_metrics = ProteinRNAInterface(plot_creator=plot_creator)
        result = protein_rna_metrics.get_proteinRNA_metrics("pred_protein_rna_dna", "exp_protein_rna_dna")
        get_rna_motif_scores(result, plot_creator)
        protein_rna_metrics.close_connection()
    # elif table == "pred_rna_rna": #to-do
    #     rna_rna_metrics = RnaRNAInterface()
    #     result = rna_rna_metrics.get_rnaRNA_metrics("pred_rna_rna", "exp_rna_rna")
    #     get_rna_motif_scores(result)
    #     rna_rna_metrics.close_conection()
