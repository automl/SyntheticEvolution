import ast
import sys
import numpy as np
import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from database.databaseMethods import DatabaseMethods
from plots.plotCreator import PlotCreator
from database.startConfig import StartConfig

def get_longest_rna_info(rna_lengths, rna_chain_ids):
    # print("rna_lengths, rna_chain_ids", rna_lengths, rna_chain_ids)
    """Get info for the longest RNA chain"""
    try:
        # Handle single RNA case (string of integer)
        if isinstance(rna_lengths, str):
            if '[' not in rna_lengths:  # Single RNA case
                return int(rna_lengths), rna_chain_ids[0] if isinstance(rna_chain_ids, list) else rna_chain_ids
            # Multiple RNA case
            rna_lengths = ast.literal_eval(rna_lengths)
            rna_chain_ids = ast.literal_eval(rna_chain_ids) if isinstance(rna_chain_ids, str) else rna_chain_ids

            # Get index of longest RNA
            max_length_idx = rna_lengths.index(max(rna_lengths))
            return rna_lengths[max_length_idx], rna_chain_ids[max_length_idx]
    except Exception as e:
        print(f"Error processing RNA lengths/chains for lengths '{rna_lengths}' and chains '{rna_chain_ids}': {e}")
        return None, None


def extract_motif_indices(rna_motifs, chain_pairs, target_chain):
    """Extract indices from RNA motifs for the target chain"""
    # try:
    # print("rna_motifs", rna_motifs, isinstance(rna_motifs, str))

    # Handle single RNA case (either single motif or list of motifs for one chain)
    if isinstance(rna_motifs, str):
        if '[' not in rna_motifs:
            # Single motif
            return [int(idx) for idx in rna_motifs.replace('(', ' ').replace(')', ' ').split()
                    if idx.isdigit()]
        else:
            # List of motifs for single chain
            motifs_list = ast.literal_eval(rna_motifs)
            all_indices = set()
            for motif in motifs_list:
                indices = [int(idx) for idx in motif.replace('(', ' ').replace(')', ' ').split()
                           if idx.isdigit()]
                all_indices.update(indices)
            return sorted(list(all_indices))

    # print("chain_pairs", chain_pairs)
    # Multiple RNA case with chain pairs
    if isinstance(chain_pairs, str):
        chain_pairs = ast.literal_eval(chain_pairs)
    # print("target_chain", target_chain)

    # Get indices of pairs containing target chain
    target_indices = [i for i, pair in enumerate(chain_pairs)
                      if target_chain in pair]
    # print("target_indices", target_indices)

    # Extract and combine indices from relevant motifs
    all_indices = set()
    for idx in target_indices:
        if idx < len(rna_motifs):
            motif = rna_motifs[idx]
            indices = [int(idx) for idx in motif.replace('(', ' ').replace(')', ' ').split()
                       if idx.isdigit()]
            all_indices.update(indices)

    return sorted(list(all_indices))
    # except Exception as e:
    #     print(f"Error extracting motif indices: {e}")
    #     return []


def check_indices_in_features(indices, features_tuple, is_duplex):
    """Check which indices appear in which secondary structure features"""
    results = {
        'stems': 0,
        'hairpin_loops': 0,
        'internal_loops': 0,
        'multibranch_loops': 0,
        'dangling_ends': 0,
        'single_stranded': 0,
        'pseudoknots': 0
    }

    try:
        indices_set = set(indices)

        # Unpack features from tuple
        stems, hairpins, internal, multibranch, dangling, pseudoknots = features_tuple

        # Check stems
        if stems:
            stems_list = ast.literal_eval(stems)
            for stem_residues in stems_list:
                stem_indices = set(stem_residues)
                results['stems'] += len(indices_set.intersection(stem_indices))

        # Check hairpin loops
        if hairpins:
            hairpins_list = ast.literal_eval(hairpins)
            for loop_indices in hairpins_list:
                loop_indices_set = set(loop_indices)
                results['hairpin_loops'] += len(indices_set.intersection(loop_indices_set))

        # Check internal loops
        if internal:
            internal_list = ast.literal_eval(internal)
            for loop_indices in internal_list:
                loop_indices_set = set(loop_indices)
                results['internal_loops'] += len(indices_set.intersection(loop_indices_set))

        # Check multibranch loops
        if multibranch:
            multibranch_list = ast.literal_eval(multibranch)
            for loop_indices in multibranch_list:
                loop_indices_set = set(loop_indices)
                results['multibranch_loops'] += len(indices_set.intersection(loop_indices_set))

        # Check dangling ends
        if dangling:
            dangling_list = ast.literal_eval(dangling)
            if dangling_list == [[], []] or all(not end_ranges for end_ranges in dangling_list):
                results['single_stranded'] = 1
            else:
                for end_ranges in dangling_list:
                    for start_end in end_ranges:
                        if start_end:
                            start, end = start_end
                            range_indices = set(range(start, end + 1))
                            intersection = len(indices_set.intersection(range_indices))
                            if is_duplex:
                                results['dangling_ends'] += intersection
                            else:
                                results['single_stranded'] += intersection

        # Check pseudoknots
        if pseudoknots:
            pseudoknots_list = ast.literal_eval(pseudoknots)
            if pseudoknots_list:
                for start, end in pseudoknots_list:
                    range_indices = set(range(start, end + 1))
                    results['pseudoknots'] += len(indices_set.intersection(range_indices))

    except Exception as e:
        print(f"Error checking indices in features: {e}")

    return results


def is_structured_motif(results):
    """Check if motif appears in structured regions"""
    structured_features = ['stems', 'internal_loops', 'multibranch_loops', 'pseudoknots']
    return any(results[feature] > 0 for feature in structured_features)


def update_motif_duplex_status(db, results_dict, table, filtered_ids):
    """Update RNAMotif_isDuplex column based on motif locations"""
    # Update each entry
    for pdb_id, indices in results_dict.items():
        if pdb_id not in filtered_ids:  # Skip if not in filtered set
            continue
            
        # Get current features
        if table == 'exp_protein_rna':
            query = """
            SELECT RNA_Stems, RNA_HairpinLoops, RNA_InternalLoops, RNA_MultibranchLoops,
                   RNA_DanglingEnds, RNA_Pseudoknots
            FROM exp_protein_rna WHERE PDBId = ?
            """
        else:
            query = """
            SELECT RNA_Stems, RNA_HairpinLoops, RNA_InternalLoops, RNA_MultibranchLoops,
                   RNA_DanglingEnds, RNA_Pseudoknots
            FROM pred_protein_rna WHERE exp_db_id = ?
            """

        features = db.execute_query(query, (pdb_id,))[0]
        results = check_indices_in_features(indices, features, None)  # is_duplex not needed here

        # Update RNAMotif_isDuplex
        is_duplex = 1 if is_structured_motif(results) else 0
        if table == 'exp_protein_rna':
            update_query = """
            UPDATE exp_protein_rna
            SET RNAMotif_isDuplex = ?
            WHERE PDBId = ?
            """
        else:
            update_query = """
            UPDATE pred_protein_rna
            SET RNAMotif_isDuplex = ?
            WHERE exp_db_id = ?
            """
        db.execute_query(update_query, (is_duplex, pdb_id))


def create_duplex_comparison_plots(db, plot_creator, exp_df, pred_df):
    """Create plots comparing duplex status between experimental and predicted structures"""
    # Use filtered DataFrames instead of new query
    merged_data = pd.merge(
        exp_df[['PDBId', 'RNA_isDuplex', 'RNAMotif_isDuplex']],
        pred_df[['exp_db_id', 'RNA_isDuplex', 'RNAMotif_isDuplex']],
        left_on='PDBId',
        right_on='exp_db_id'
    )

    # Prepare data for RNA_isDuplex comparison
    rna_duplex_diff = {
        'pdb_ids': [],
        'exp_values': [],
        'pred_values': []
    }

    # Prepare data for RNAMotif_isDuplex comparison
    motif_duplex_diff = {
        'pdb_ids': [],
        'exp_values': [],
        'pred_values': []
    }

    # Process each row
    for _, row in merged_data.iterrows():
        # Convert values to numeric, defaulting to 0 if conversion fails
        try:
            exp_rna = float(row['RNA_isDuplex_x'])
            pred_rna = float(row['RNA_isDuplex_y'])
            exp_motif = float(row['RNAMotif_isDuplex_x'])
            pred_motif = float(row['RNAMotif_isDuplex_y'])
        except (ValueError, TypeError):
            continue

        # Only include if values are different
        if exp_rna != pred_rna:
            rna_duplex_diff['pdb_ids'].append(row['PDBId'])
            rna_duplex_diff['exp_values'].append(exp_rna)
            rna_duplex_diff['pred_values'].append(pred_rna)

        if exp_motif != pred_motif:
            motif_duplex_diff['pdb_ids'].append(row['PDBId'])
            motif_duplex_diff['exp_values'].append(exp_motif)
            motif_duplex_diff['pred_values'].append(pred_motif)

    # Create plots only if there are differences to show
    if rna_duplex_diff['pdb_ids']:
        plot_creator.get_grouped_barplot(
            table_source='rna_duplex_comparison',
            groups=[rna_duplex_diff['exp_values'], rna_duplex_diff['pred_values']],
            group_labels=rna_duplex_diff['pdb_ids'],
            legend_labels=['Experimental', 'Predicted'],
            ylabel='RNA_isDuplex',
            title='Differences in RNA_isDuplex between Experimental and Predicted Structures'
        )

    if motif_duplex_diff['pdb_ids']:
        plot_creator.get_grouped_barplot(
            table_source='rna_motif_duplex_comparison',
            groups=[motif_duplex_diff['exp_values'], motif_duplex_diff['pred_values']],
            group_labels=motif_duplex_diff['pdb_ids'],
            legend_labels=['Experimental', 'Predicted'],
            ylabel='RNAMotif_isDuplex',
            title='Differences in RNAMotif_isDuplex between Experimental and Predicted Structures'
        )


def create_duplex_status_distribution(db, plot_creator, exp_df, pred_df):
    """Create plot showing distribution of duplex statuses"""
    # Use filtered DataFrames instead of new query
    merged_data = pd.merge(
        exp_df[['PDBId', 'RNA_isDuplex', 'RNAMotif_isDuplex']],
        pred_df[['exp_db_id', 'RNA_isDuplex', 'RNAMotif_isDuplex']],
        left_on='PDBId',
        right_on='exp_db_id'
    )
    data = merged_data[[
        'RNA_isDuplex_x', 'RNAMotif_isDuplex_x',
        'RNA_isDuplex_y', 'RNAMotif_isDuplex_y'
    ]].values

    # Initialize counters for different categories
    categories = {
        'Exp_1_Pred_1': {'RNA': 0, 'Motif': 0},
        'Exp_0_Pred_0': {'RNA': 0, 'Motif': 0},
        'Exp_1_Pred_0': {'RNA': 0, 'Motif': 0},
        'Exp_0_Pred_1': {'RNA': 0, 'Motif': 0}
    }

    total = len(data)

    # Count occurrences
    for exp_rna, exp_motif, pred_rna, pred_motif in data:
        # RNA_isDuplex categories
        if exp_rna == 1 and pred_rna == 1:
            categories['Exp_1_Pred_1']['RNA'] += 1
        elif exp_rna == 0 and pred_rna == 0:
            categories['Exp_0_Pred_0']['RNA'] += 1
        elif exp_rna == 1 and pred_rna == 0:
            categories['Exp_1_Pred_0']['RNA'] += 1
        elif exp_rna == 0 and pred_rna == 1:
            categories['Exp_0_Pred_1']['RNA'] += 1

        # RNAMotif_isDuplex categories
        if exp_motif == 1 and pred_motif == 1:
            categories['Exp_1_Pred_1']['Motif'] += 1
        elif exp_motif == 0 and pred_motif == 0:
            categories['Exp_0_Pred_0']['Motif'] += 1
        elif exp_motif == 1 and pred_motif == 0:
            categories['Exp_1_Pred_0']['Motif'] += 1
        elif exp_motif == 0 and pred_motif == 1:
            categories['Exp_0_Pred_1']['Motif'] += 1

    # Calculate percentages
    category_labels = ['Exp 1/Pred 1', 'Exp 0/Pred 0', 'Exp 1/Pred 0', 'Exp 0/Pred 1']
    rna_percentages = [
        (categories['Exp_1_Pred_1']['RNA'] / total * 100),
        (categories['Exp_0_Pred_0']['RNA'] / total * 100),
        (categories['Exp_1_Pred_0']['RNA'] / total * 100),
        (categories['Exp_0_Pred_1']['RNA'] / total * 100)
    ]
    motif_percentages = [
        (categories['Exp_1_Pred_1']['Motif'] / total * 100),
        (categories['Exp_0_Pred_0']['Motif'] / total * 100),
        (categories['Exp_1_Pred_0']['Motif'] / total * 100),
        (categories['Exp_0_Pred_1']['Motif'] / total * 100)
    ]

    # Create grouped bar plot
    plot_creator.get_grouped_barplot(
        table_source='duplex_distribution',
        groups=[rna_percentages, motif_percentages],
        group_labels=category_labels,
        legend_labels=['RNA_isDuplex', 'RNAMotif_isDuplex'],
        ylabel='Percentage of samples [%]',
        title='Distribution of Duplex Agreement between Experimental and Predicted Structures'
    )


def create_motif_score_correlation_plot(db, plot_creator, exp_df, pred_df):
    """Create violin plot comparing RNA motif scores across duplex status categories"""
    # Use filtered DataFrames instead of new query
    merged_data = pd.merge(
        exp_df[['PDBId', 'RNAMotif_isDuplex']],
        pred_df[['exp_db_id', 'RNAMotif_isDuplex', 'RNAmotif_score', 'af3_rna_MSA']],
        left_on='PDBId',
        right_on='exp_db_id'
    )
    data = merged_data[[
        'PDBId', 'RNAMotif_isDuplex_x', 'RNAMotif_isDuplex_y',
        'RNAmotif_score', 'af3_rna_MSA'
    ]].values

    # Original duplex categories
    scores_by_category = {
        'Exp ds/Pred ds': [],
        'Exp ss/Pred ss': [],
        'Exp ds/Pred ss': [],
        'Exp ss/Pred ds': []
    }

    # Extended categories with MSA info
    scores_by_category_msa = {
        '+MSA_Exp ds/Pred ds': [],
        '+MSA_Exp ss/Pred ss': [],
        '+MSA_Exp ds/Pred ss': [],
        '+MSA_Exp ss/Pred ds': [],
        '-MSA_Exp ds/Pred ds': [],
        '-MSA_Exp ss/Pred ss': [],
        '-MSA_Exp ds/Pred ss': [],
        '-MSA_Exp ss/Pred ds': []
    }

    # Collect scores for each category
    print(data)
    for pdb_id, exp_motif, pred_motif, score, msa_list in data:
        # Determine MSA status
        try:
            print(msa_list)
            msa_values = ast.literal_eval(msa_list)
            is_msa = msa_values[0] == 1 if len(msa_values) == 1 else False
        except Exception as e:
            print(f"Error processing MSA data for PDB {pdb_id}: {e}")
            continue

        # Add to original categories
        if exp_motif == 1 and pred_motif == 1:
            scores_by_category['Exp ds/Pred ds'].append(score)
            scores_by_category_msa[f"{'+MSA' if is_msa else '-MSA'}_Exp ds/Pred ds"].append(score)
        elif exp_motif == 0 and pred_motif == 0:
            scores_by_category['Exp ss/Pred ss'].append(score)
            scores_by_category_msa[f"{'+MSA' if is_msa else '-MSA'}_Exp ss/Pred ss"].append(score)
        elif exp_motif == 1 and pred_motif == 0:
            scores_by_category['Exp ds/Pred ss'].append(score)
            scores_by_category_msa[f"{'+MSA' if is_msa else '-MSA'}_Exp ds/Pred ss"].append(score)
        elif exp_motif == 0 and pred_motif == 1:
            scores_by_category['Exp ss/Pred ds'].append(score)
            scores_by_category_msa[f"{'+MSA' if is_msa else '-MSA'}_Exp ss/Pred ds"].append(score)

    # Create original violin plot
    labels = ['Exp ds/Pred ds', 'Exp ss/Pred ss', 'Exp ds/Pred ss', 'Exp ss/Pred ds']
    data_values = [scores_by_category[label] for label in labels]

    plot_creator.get_violin_plot(
        table_source='motif_score_distribution',
        data_values=data_values,
        labels=labels,
        name='RNAMotif_isDuplex_Distribution',
        yAxis_label='RNA motif Score',
        rotation=90
    )

    # Create MSA-aware boxplot
    msa_labels = [
        '+MSA\nExp ds/Pred ds', '+MSA\nExp ss/Pred ss', '+MSA\nExp ds/Pred ss', '+MSA\nExp ss/Pred ds',
        '-MSA\nExp ds/Pred ds', '-MSA\nExp ss/Pred ss', '-MSA\nExp ds/Pred ss', '-MSA\nExp ss/Pred ds'
    ]
    msa_data_values = [scores_by_category_msa[label.replace('\n', '_')] for label in msa_labels]

    plot_creator.get_boxplot(
        table_source='motif_score_distribution',
        data_values=msa_data_values,
        labels=msa_labels,
        name='RNAMotif_Score_MSA_Duplex_Distribution1',
        yAxis_label='RNA motif Score',
        rotation = 90
    )

    plot_creator.get_violin_plot(
        table_source='motif_score_distribution',
        data_values=msa_data_values,
        labels=msa_labels,
        name='RNAMotif_Score_MSA_Duplex_Distribution2',
        yAxis_label='RNA motif Score',
        rotation=90# Angled labels for better readability
    )

    # Save statistics for both plots
    summary_text = f"\nRNAmotif Score Statistics by Category:\n" + "=" * 30 + "\n"

    # Original categories statistics
    summary_text += "\nOriginal Categories:\n" + "-" * 20 + "\n"
    for category, scores in scores_by_category.items():
        if scores:
            avg_score = sum(scores) / len(scores)
            summary_text += f"{category}:\n"
            summary_text += f"  Count: {len(scores)}\n"
            summary_text += f"  Average score: {avg_score:.2f}\n\n"

    # MSA categories statistics
    summary_text += "\nMSA Categories:\n" + "-" * 20 + "\n"
    for category, scores in scores_by_category_msa.items():
        if scores:
            avg_score = sum(scores) / len(scores)
            summary_text += f"{category}:\n"
            summary_text += f"  Count: {len(scores)}\n"
            summary_text += f"  Average score: {avg_score:.2f}\n\n"

    save_analysis_data(summary_text)


def is_single_stranded_dangling(dangling_list, stems=None, hairpins=None, internal=None, multibranch=None,
                                pseudoknots=None):
    """
    Check if dangling ends indicate single-stranded region or dangling ends
    [[(1, 6)], []] -> single-stranded (only 5' has a long range AND no other features)
    [[], [(25, 35)]] -> dangling ends (only 3' has a range)
    [[(1, 2)], [(22, 24)]] -> dangling ends (both ends have ranges)
    [[], []] -> neither
    """
    if not dangling_list or dangling_list == [[], []]:
        return None  # No dangling ends and no single stranded

    # Check if structure has any other features
    has_other_features = False
    if stems and ast.literal_eval(stems):
        has_other_features = True
    if hairpins and ast.literal_eval(hairpins):
        has_other_features = True
    if internal and ast.literal_eval(internal):
        has_other_features = True
    if multibranch and ast.literal_eval(multibranch):
        has_other_features = True
    if pseudoknots and ast.literal_eval(pseudoknots):
        has_other_features = True

    # Single stranded: only 5' end has a significant range AND no other features
    if dangling_list[0] and not dangling_list[1] and not has_other_features:
        for start, end in dangling_list[0]:
            if end - start > 0:  # Range spans multiple residues
                return 'single_stranded'

    # All other cases with ranges are dangling ends
    if any(ranges for ranges in dangling_list):
        return 'dangling_ends'

    return None

def analyze_and_plot_motifs():
    db = DatabaseMethods()
    
    # Get initial data as DataFrames for filtering
    pred_df = pd.read_sql_query("SELECT * FROM pred_protein_rna", db.connection)
    exp_df = pd.read_sql_query("SELECT * FROM exp_protein_rna", db.connection)
    
    # Apply filters if specified in command line args
    single_chain_only = sys.argv[2] == 'single_chain' if len(sys.argv) > 2 else False
    msa_option = sys.argv[3].lower() if len(sys.argv) > 3 else None  # Convert to lowercase
    
    if single_chain_only:
        exp_df, pred_df = db.filter_by_singleChain(exp_df, pred_df)
        
    if msa_option is not None:
        # Handle both +MSA/-MSA and +msa/-msa formats
        msa_option = '+MSA' if msa_option in ['+msa', '+msa'] else '-MSA' if msa_option in ['-msa', '-msa'] else None
        if msa_option:
            pred_df = db.filter_by_msa(pred_df, msa_option)
            # Filter experimental data to match
            exp_df = exp_df[exp_df['PDBId'].isin(pred_df['exp_db_id'])]
    
    plot_creator = PlotCreator('RNA_secondary_feature', msa_option, single_chain_only)
    
    # Convert filtered DataFrames to results format, only including filtered entries
    results = []
    filtered_pdb_ids = set(exp_df['PDBId'].values)  # Create a set of filtered PDB IDs
    filtered_pred_ids = set(pred_df['exp_db_id'].values)  # Create a set of filtered predicted IDs
    
    # Only include entries that exist in both filtered DataFrames
    for _, exp_row in exp_df.iterrows():
        if exp_row['PDBId'] in filtered_pred_ids:  # Check if this PDB ID is in the filtered predicted data
            pred_row = pred_df[pred_df['exp_db_id'] == exp_row['PDBId']].iloc[0]
            results.append((
                exp_row['PDBId'], exp_row['RNALength'], exp_row['RNAChainIDs'],
                exp_row['ChainIDpairList_proteinRNA'], exp_row['RNAMotif'],
                exp_row['RNA_Stems'], exp_row['RNA_HairpinLoops'], exp_row['RNA_InternalLoops'],
                exp_row['RNA_MultibranchLoops'], exp_row['RNA_Pseudoknots'], exp_row['RNA_DanglingEnds'],
                pred_row['RNAMotif'], pred_row['RNA_Stems'], pred_row['RNA_HairpinLoops'],
                pred_row['RNA_InternalLoops'], pred_row['RNA_MultibranchLoops'],
                pred_row['RNA_Pseudoknots'], pred_row['RNA_DanglingEnds'],
                pred_row['RNAmotif_score'], pred_row['af3_rna_MSA']
            ))
    
    # Print debug information about filtering
    print(f"Total entries in exp_df: {len(exp_df)}")
    print(f"Total entries in pred_df: {len(pred_df)}")
    print(f"Total entries in filtered results: {len(results)}")
    
    # Initialize containers for results
    exp_results = {}
    pred_results = {}

    # Initialize counters for motif locations
    feature_counts = {
        'Stems': {'exp': 0, 'pred': 0},
        'Hairpin Loops': {'exp': 0, 'pred': 0},
        'Internal Loops': {'exp': 0, 'pred': 0},
        'Multibranch Loops': {'exp': 0, 'pred': 0},
        'Pseudoknots': {'exp': 0, 'pred': 0},
        'Dangling Ends': {'exp': 0, 'pred': 0},
        'Single Stranded': {'exp': 0, 'pred': 0}
    }
    total_motifs = {'exp': 0, 'pred': 0}

    # Initialize MSA score lists
    msa_scores = []
    non_msa_scores = []

    for row in results:
        pdb_id = row[0]
        rna_lengths = row[1]
        rna_chain_ids = row[2]
        exp_motif = row[4]
        pred_motif = row[11]

        # Process motif residues
        exp_residues = set()
        pred_residues = set()
        # print("pdb_id", pdb_id)
        # Single RNA case is determined by RNAMotif format
        if exp_motif:
            indices = extract_motif_indices(exp_motif, row[3], row[2])
            if indices:
                exp_residues.update(indices)
                exp_results[pdb_id] = indices
                # print("indices",indices)

        if pred_motif:
            indices = extract_motif_indices(pred_motif, row[3], row[2])
            if indices:
                pred_residues.update(indices)
                pred_results[pdb_id] = indices
                # print("indices", indices)

        # Check experimental features
        if exp_residues:
            total_motifs['exp'] += 1

            # Check stems
            if row[5]:  # RNA_Stems
                stem_residues = set([r for stem in ast.literal_eval(row[5]) for r in stem])
                if exp_residues.intersection(stem_residues):
                    feature_counts['Stems']['exp'] += 1

            # Check hairpin loops
            if row[6]:  # RNA_HairpinLoops
                hairpin_residues = set([r for loop in ast.literal_eval(row[6]) for r in loop])
                if exp_residues.intersection(hairpin_residues):
                    feature_counts['Hairpin Loops']['exp'] += 1

            # Check internal loops
            if row[7]:  # RNA_InternalLoops
                internal_residues = set([r for loop in ast.literal_eval(row[7]) for r in loop])
                if exp_residues.intersection(internal_residues):
                    feature_counts['Internal Loops']['exp'] += 1

            # Check multibranch loops
            if row[8]:  # RNA_MultibranchLoops
                multibranch_residues = set([r for loop in ast.literal_eval(row[8]) for r in loop])
                if exp_residues.intersection(multibranch_residues):
                    feature_counts['Multibranch Loops']['exp'] += 1

            # Check pseudoknots
            if row[9]:  # RNA_Pseudoknots
                pseudoknot_residues = set([r for pair in ast.literal_eval(row[9]) for r in pair])
                if exp_residues.intersection(pseudoknot_residues):
                    feature_counts['Pseudoknots']['exp'] += 1

            # Check dangling ends and single stranded regions
            if row[10]:  # RNA_DanglingEnds
                dangling_list = ast.literal_eval(row[10])
                region_type = is_single_stranded_dangling(
                    dangling_list,
                    stems=row[5],  # RNA_Stems
                    hairpins=row[6],  # RNA_HairpinLoops
                    internal=row[7],  # RNA_InternalLoops
                    multibranch=row[8],  # RNA_MultibranchLoops
                    pseudoknots=row[9]  # RNA_Pseudoknots
                )

                if region_type == 'single_stranded':
                    feature_counts['Single Stranded']['exp'] += 1
                elif region_type == 'dangling_ends':
                    for end_ranges in dangling_list:
                        for start_end in end_ranges:
                            if start_end:
                                start, end = start_end
                                range_indices = set(range(start, end + 1))
                                if exp_residues.intersection(range_indices):
                                    feature_counts['Dangling Ends']['exp'] += 1
                                    break

        # Check predicted features
        if pred_residues:
            total_motifs['pred'] += 1

            # Check predicted stems
            if row[12]:  # pred_stems
                stem_residues = set([r for stem in ast.literal_eval(row[12]) for r in stem])
                if pred_residues.intersection(stem_residues):
                    feature_counts['Stems']['pred'] += 1

            # Check predicted hairpin loops
            if row[13]:  # pred_hairpins
                hairpin_residues = set([r for loop in ast.literal_eval(row[13]) for r in loop])
                if pred_residues.intersection(hairpin_residues):
                    feature_counts['Hairpin Loops']['pred'] += 1

            # Check predicted internal loops
            if row[14]:  # pred_internal
                internal_residues = set([r for loop in ast.literal_eval(row[14]) for r in loop])
                if pred_residues.intersection(internal_residues):
                    feature_counts['Internal Loops']['pred'] += 1

            # Check predicted multibranch loops
            if row[15]:  # pred_multibranch
                multibranch_residues = set([r for loop in ast.literal_eval(row[15]) for r in loop])
                if pred_residues.intersection(multibranch_residues):
                    feature_counts['Multibranch Loops']['pred'] += 1

            # Check predicted pseudoknots
            if row[16]:  # pred_pseudoknots
                pseudoknot_residues = set([r for pair in ast.literal_eval(row[16]) for r in pair])
                if pred_residues.intersection(pseudoknot_residues):
                    feature_counts['Pseudoknots']['pred'] += 1

            # Check predicted dangling ends and single stranded regions
            if row[17]:  # pred_dangling
                dangling_list = ast.literal_eval(row[17])
                region_type = is_single_stranded_dangling(
                    dangling_list,
                    stems=row[12],  # pred_stems
                    hairpins=row[13],  # pred_hairpins
                    internal=row[14],  # pred_internal
                    multibranch=row[15],  # pred_multibranch
                    pseudoknots=row[16]  # pred_pseudoknots
                )

                if region_type == 'single_stranded':
                    feature_counts['Single Stranded']['pred'] += 1
                elif region_type == 'dangling_ends':
                    for end_ranges in dangling_list:
                        for start_end in end_ranges:
                            if start_end:
                                start, end = start_end
                                range_indices = set(range(start, end + 1))
                                if pred_residues.intersection(range_indices):
                                    feature_counts['Dangling Ends']['pred'] += 1
                                    break

        # Process MSA scores
        score = row[18]  # RNAmotif_score
        msa_list = row[19]  # af3_rna_MSA

        if score is not None and msa_list is not None:
            try:
                msa_values = ast.literal_eval(msa_list)
                # Skip MSA processing if there's only one element
                if len(msa_values) == 1:
                    if msa_values[0] == 1:
                        msa_scores.append(score)
                    else:
                        non_msa_scores.append(score)
                else:
                    length_idx, _ = get_longest_rna_info(rna_lengths, rna_chain_ids)
                    if length_idx is not None:
                        if isinstance(length_idx, int) and length_idx < len(msa_values):
                            if msa_values[length_idx] == 1:
                                msa_scores.append(score)
                            else:
                                non_msa_scores.append(score)
            except Exception as e:
                # print(f"Error processing MSA data for PDB {pdb_id}: {e}")
                # print(f"RNA lengths: {rna_lengths}, Chain IDs: {rna_chain_ids}")
                # print(f"MSA list: {msa_list}")
                continue

        # Store results for duplex analysis
        if exp_residues:
            exp_results[pdb_id] = list(exp_residues)
        if pred_residues:
            pred_results[pdb_id] = list(pred_residues)

    if total_motifs['exp'] > 0 and total_motifs['pred'] > 0:
        exp_percentages = [count['exp'] / total_motifs['exp'] * 100 for count in feature_counts.values()]
        pred_percentages = [count['pred'] / total_motifs['pred'] * 100 for count in feature_counts.values()]

        plot_creator.get_grouped_barplot(
            table_source='rna_motif_locations',
            groups=[exp_percentages, pred_percentages],
            group_labels=list(feature_counts.keys()),
            legend_labels=['Experimental', 'Predicted'],
            ylabel='Percentage of Motifs (%)',
            title='Distribution of RNA Motif Residues in Secondary Structure Features'
        )

    # Save feature location statistics
    summary_text = "\nRNA Motif Location Statistics:\n" + "=" * 30 + "\n"
    summary_text += f"Total motifs analyzed - Experimental: {total_motifs['exp']}, Predicted: {total_motifs['pred']}\n\n"

    for feature, counts in feature_counts.items():
        exp_pct = counts['exp'] / total_motifs['exp'] * 100 if total_motifs['exp'] > 0 else 0
        pred_pct = counts['pred'] / total_motifs['pred'] * 100 if total_motifs['pred'] > 0 else 0
        summary_text += f"{feature}:\n"
        summary_text += f"  Experimental: {counts['exp']} ({exp_pct:.1f}%)\n"
        summary_text += f"  Predicted: {counts['pred']} ({pred_pct:.1f}%)\n\n"
    save_analysis_data(summary_text)

    # Create MSA comparison plots
    if msa_scores and non_msa_scores:
        # Separate feature counts by MSA status for both exp and pred
        msa_feature_counts = {
            'exp': {
                'Stems': 0,
                'Hairpin Loops': 0,
                'Internal Loops': 0,
                'Multibranch Loops': 0,
                'Pseudoknots': 0,
                'Dangling Ends': 0,
                'Single Stranded': 0
            },
            'pred': {
                'Stems': 0,
                'Hairpin Loops': 0,
                'Internal Loops': 0,
                'Multibranch Loops': 0,
                'Pseudoknots': 0,
                'Dangling Ends': 0,
                'Single Stranded': 0
            }
        }
        non_msa_feature_counts = {
            'exp': msa_feature_counts['exp'].copy(),
            'pred': msa_feature_counts['pred'].copy()
        }

        total_motifs = {
            'msa': {'exp': 0, 'pred': 0},
            'non_msa': {'exp': 0, 'pred': 0}
        }

        # Reprocess the data to count features by MSA status
        for row in results:
            score = row[18]  # RNAmotif_score
            msa_list = row[19]  # af3_rna_MSA
            exp_motif = row[4]  # exp RNAMotif
            pred_motif = row[11]  # pred RNAMotif

            if score is not None and msa_list is not None and exp_motif and pred_motif:
                try:
                    # Determine MSA status
                    msa_values = ast.literal_eval(msa_list)
                    if len(msa_values) == 1:
                        is_msa = msa_values[0] == 1
                    else:
                        length_idx, _ = get_longest_rna_info(row[1], row[2])
                        if length_idx is not None and isinstance(length_idx, int) and length_idx < len(msa_values):
                            is_msa = msa_values[length_idx] == 1
                        else:
                            continue

                    # Get motif residues for both exp and pred
                    exp_residues = set()
                    pred_residues = set()

                    exp_indices = extract_motif_indices(exp_motif, row[3], row[2])
                    pred_indices = extract_motif_indices(pred_motif, row[3], row[2])

                    if exp_indices:
                        exp_residues.update(exp_indices)
                    if pred_indices:
                        pred_residues.update(pred_indices)

                    if not exp_residues and not pred_residues:
                        continue

                    # Count features for both structures based on MSA status
                    feature_dict = msa_feature_counts if is_msa else non_msa_feature_counts
                    total_dict = 'msa' if is_msa else 'non_msa'

                    # Process experimental features
                    if exp_residues:
                        total_motifs[total_dict]['exp'] += 1

                        # Check stems
                        if row[5] and exp_residues.intersection(
                                set([r for stem in ast.literal_eval(row[5]) for r in stem])):
                            feature_dict['exp']['Stems'] += 1

                        # Check hairpin loops
                        if row[6] and exp_residues.intersection(
                                set([r for loop in ast.literal_eval(row[6]) for r in loop])):
                            feature_dict['exp']['Hairpin Loops'] += 1

                        # Check internal loops
                        if row[7] and exp_residues.intersection(
                                set([r for loop in ast.literal_eval(row[7]) for r in loop])):
                            feature_dict['exp']['Internal Loops'] += 1

                        # Check multibranch loops
                        if row[8] and exp_residues.intersection(
                                set([r for loop in ast.literal_eval(row[8]) for r in loop])):
                            feature_dict['exp']['Multibranch Loops'] += 1

                        # Check pseudoknots
                        if row[9] and exp_residues.intersection(
                                set([r for pair in ast.literal_eval(row[9]) for r in pair])):
                            feature_dict['exp']['Pseudoknots'] += 1

                        # Check dangling ends and single stranded
                        if row[10]:
                            dangling_list = ast.literal_eval(row[10])
                            region_type = is_single_stranded_dangling(
                                dangling_list,
                                stems=row[5],
                                hairpins=row[6],
                                internal=row[7],
                                multibranch=row[8],
                                pseudoknots=row[9]
                            )
                            if region_type == 'single_stranded':
                                feature_dict['exp']['Single Stranded'] += 1
                            elif region_type == 'dangling_ends':
                                feature_dict['exp']['Dangling Ends'] += 1

                    # Process predicted features
                    if pred_residues:
                        total_motifs[total_dict]['pred'] += 1

                        # Check stems
                        if row[12] and pred_residues.intersection(
                                set([r for stem in ast.literal_eval(row[12]) for r in stem])):
                            feature_dict['pred']['Stems'] += 1

                        # Check hairpin loops
                        if row[13] and pred_residues.intersection(
                                set([r for loop in ast.literal_eval(row[13]) for r in loop])):
                            feature_dict['pred']['Hairpin Loops'] += 1

                        # Check internal loops
                        if row[14] and pred_residues.intersection(
                                set([r for loop in ast.literal_eval(row[14]) for r in loop])):
                            feature_dict['pred']['Internal Loops'] += 1

                        # Check multibranch loops
                        if row[15] and pred_residues.intersection(
                                set([r for loop in ast.literal_eval(row[15]) for r in loop])):
                            feature_dict['pred']['Multibranch Loops'] += 1

                        # Check pseudoknots
                        if row[16] and pred_residues.intersection(
                                set([r for pair in ast.literal_eval(row[16]) for r in pair])):
                            feature_dict['pred']['Pseudoknots'] += 1

                        # Check dangling ends and single stranded
                        if row[17]:
                            dangling_list = ast.literal_eval(row[17])
                            region_type = is_single_stranded_dangling(
                                dangling_list,
                                stems=row[12],
                                hairpins=row[13],
                                internal=row[14],
                                multibranch=row[15],
                                pseudoknots=row[16]
                            )
                            if region_type == 'single_stranded':
                                feature_dict['pred']['Single Stranded'] += 1
                            elif region_type == 'dangling_ends':
                                feature_dict['pred']['Dangling Ends'] += 1

                except Exception as e:
                    print(f"Error processing MSA data for PDB {row[0]}: {e}")
                    continue

        # Create plots for experimental and predicted structures
        for struct_type in ['exp', 'pred']:
            # Calculate percentages
            msa_percentages = [count / total_motifs['msa'][struct_type] * 100
                               if total_motifs['msa'][struct_type] > 0 else 0
                               for count in msa_feature_counts[struct_type].values()]
            non_msa_percentages = [count / total_motifs['non_msa'][struct_type] * 100
                                   if total_motifs['non_msa'][struct_type] > 0 else 0
                                   for count in non_msa_feature_counts[struct_type].values()]

            # Create grouped bar plot
            plot_creator.get_grouped_barplot(
                table_source=f'motif_msa_distribution_{struct_type}',
                groups=[msa_percentages, non_msa_percentages],
                group_labels=list(msa_feature_counts[struct_type].keys()),
                legend_labels=['MSA', 'No MSA'],
                ylabel='Percentage of Motifs [%]',
                title=f'{struct_type.upper()}_Distribution_of_RNA_Motif_Residues_by_MSA_Status',
                yaxis_range=[0, 100]
            )

            # Save statistics
            summary_text = f"\nRNA Motif Distribution Statistics for {struct_type.upper()}:\n" + "=" * 30 + "\n"
            summary_text += f"Total MSA motifs: {total_motifs['msa'][struct_type]}\n"
            summary_text += f"Total non-MSA motifs: {total_motifs['non_msa'][struct_type]}\n\n"

            for feature in msa_feature_counts[struct_type].keys():
                msa_count = msa_feature_counts[struct_type][feature]
                non_msa_count = non_msa_feature_counts[struct_type][feature]
                msa_pct = msa_count / total_motifs['msa'][struct_type] * 100 if total_motifs['msa'][
                                                                                    struct_type] > 0 else 0
                non_msa_pct = non_msa_count / total_motifs['non_msa'][struct_type] * 100 if total_motifs['non_msa'][
                                                                                                struct_type] > 0 else 0

                summary_text += f"{feature}:\n"
                summary_text += f"  MSA: {msa_count} ({msa_pct:.1f}%)\n"
                summary_text += f"  No MSA: {non_msa_count} ({non_msa_pct:.1f}%)\n\n"

            save_analysis_data(summary_text)

    # Update RNAMotif_isDuplex in both tables using filtered data
    update_motif_duplex_status(db, exp_results, 'exp_protein_rna', exp_df['PDBId'].tolist())
    update_motif_duplex_status(db, pred_results, 'pred_protein_rna', pred_df['exp_db_id'].tolist())

    # Pass filtered DataFrames to plotting functions
    create_duplex_comparison_plots(db, plot_creator, exp_df, pred_df)
    create_duplex_status_distribution(db, plot_creator, exp_df, pred_df)
    create_motif_score_correlation_plot(db, plot_creator, exp_df, pred_df)


def calculate_stem_length(stem):
    """Calculate length of a stem based on consecutive numbers"""
    if not stem:
        return 0
        
    length = 0
    current_seq = [stem[0]]
    
    for num in stem[1:]:
        if num == current_seq[-1] + 1:
            current_seq.append(num)
        else:
            length += len(current_seq)
            current_seq = [num]
    
    # Add length of last sequence
    length += len(current_seq)
    return length // 2  # Divide by 2 since we count both sides of the stem

def calculate_hairpin_internal_length(loop_list):
    """Calculate length for hairpin and internal loops by counting all residues"""
    if not loop_list:
        return 0
    # If it's a list of lists (full feature list)
    if isinstance(loop_list[0], list):
        return sum(len(loop) for loop in loop_list)
    # If it's a single loop (list of residues)
    return len(loop_list)


def calculate_multibranch_length(loop_list):
    """Calculate length for multibranch loops: count residues + 1 for each sublist"""
    if not loop_list:
        return 0
    # If it's a list of lists (full feature list)
    if isinstance(loop_list[0], list):
        junctions = len(loop_list)
        residues = sum(len(sublist) for sublist in loop_list if sublist)
        return junctions + residues
    # If it's a single loop (list of residues)
    return 1 + len(loop_list) if loop_list else 1


def calculate_pseudoknot_length(pseudoknot_list):
    """Calculate length for pseudoknots by counting pairs"""
    if not pseudoknot_list:
        return 0

    # If it's a list of tuples (pairs), count the number of pairs
    if isinstance(pseudoknot_list, list) and pseudoknot_list and isinstance(pseudoknot_list[0], tuple):
        return len(pseudoknot_list)

    # If it's a list of lists of tuples (multiple pseudoknots)
    if isinstance(pseudoknot_list, list) and pseudoknot_list and isinstance(pseudoknot_list[0], list):
        return sum(len(pseudoknot) for pseudoknot in pseudoknot_list)

    return 0

def save_analysis_data(content, filename='secondary_feat_analysis.txt'):
    config = StartConfig()
    dir = config.get_parent_folder()
    filepath = os.path.join(dir, filename)
    # print("filepath", filepath)
    with open(filepath, 'a') as f:
        f.write(content + '\n')

def compare_rna_features(pred_table, exp_table, plot_creator):
    """Compare RNA features between predicted and experimental structures"""
    db = DatabaseMethods()
    
    # Get initial data as DataFrames for filtering
    pred_df = pd.read_sql_query("SELECT * FROM pred_protein_rna", db.connection)
    exp_df = pd.read_sql_query("SELECT * FROM exp_protein_rna", db.connection)
    
    # Apply filters if specified in command line args
    single_chain_only = sys.argv[2] == 'single_chain' if len(sys.argv) > 2 else False
    msa_option = sys.argv[3].lower() if len(sys.argv) > 3 else None  # Convert to lowercase
    
    if single_chain_only:
        exp_df, pred_df = db.filter_by_singleChain(exp_df, pred_df)
        
    if msa_option is not None:
        # Handle both +MSA/-MSA and +msa/-msa formats
        msa_option = '+MSA' if msa_option in ['+msa', '+msa'] else '-MSA' if msa_option in ['-msa', '-msa'] else None
        if msa_option:
            pred_df = db.filter_by_msa(pred_df, msa_option)
            # Filter experimental data to match
            exp_df = exp_df[exp_df['PDBId'].isin(pred_df['exp_db_id'])]

    # Print debug information about filtering
    print(f"Total entries in exp_df: {len(exp_df)}")
    print(f"Total entries in pred_df: {len(pred_df)}")

    # Define features to analyze
    features = {
        'RNA_Stems': calculate_stem_length,
        'RNA_HairpinLoops': calculate_hairpin_internal_length,
        'RNA_InternalLoops': calculate_hairpin_internal_length,
        'RNA_MultibranchLoops': calculate_multibranch_length,
        'RNA_Pseudoknots': calculate_pseudoknot_length
    }

    for feature_name, length_calculator in features.items():
        # Use filtered DataFrames instead of new query
        merged_data = pd.merge(
            pred_df[['exp_db_id', feature_name]],
            exp_df[['PDBId', feature_name]],
            left_on='exp_db_id',
            right_on='PDBId'
        )
        
        # Convert to format similar to previous query results
        results = [(row['PDBId'], row[f'{feature_name}_x'], row[f'{feature_name}_y']) 
                  for _, row in merged_data.iterrows()]

        # Analysis containers
        feat_counts = {'pred': [], 'exp': []}
        feat_lengths = {'pred': [], 'exp': []}
        feat_list_pred = []
        feat_list_exp = []
        total_residues = {'pred': 0, 'exp': 0}
        count_pdbs = {}

        # Process results
        for pdb_id, pred_feat_str, exp_feat_str in results:
            try:
                pred_feat = ast.literal_eval(pred_feat_str) if pred_feat_str else []
                exp_feat = ast.literal_eval(exp_feat_str) if exp_feat_str else []

                feat_list_pred.append(pred_feat)
                feat_list_exp.append(exp_feat)

                # Count features
                pred_count = len(pred_feat)
                exp_count = len(exp_feat)

                # Store PDB IDs for correlation plot
                key = (exp_count, pred_count)
                if key not in count_pdbs:
                    count_pdbs[key] = []
                count_pdbs[key].append(pdb_id)

                feat_counts['pred'].append(pred_count)
                feat_counts['exp'].append(exp_count)

                # Calculate lengths based on feature type
                if feature_name == 'RNA_Pseudoknots':
                    # For pseudoknots, calculate length directly from the list
                    if pred_feat:
                        total_length = length_calculator(pred_feat)  # Will count all pairs
                        feat_lengths['pred'].append(total_length)
                        total_residues['pred'] += total_length

                    if exp_feat:
                        total_length = length_calculator(exp_feat)  # Will count all pairs
                        feat_lengths['exp'].append(total_length)
                        total_residues['exp'] += total_length
                else:
                    # For other features, sum lengths of individual features
                    if pred_feat:
                        total_length = 0
                        for single_feat in pred_feat:
                            total_length += length_calculator(single_feat)
                        feat_lengths['pred'].append(total_length)
                        total_residues['pred'] += total_length

                    if exp_feat:
                        total_length = 0
                        for single_feat in exp_feat:
                            total_length += length_calculator(single_feat)
                        feat_lengths['exp'].append(total_length)
                        total_residues['exp'] += total_length

            except (ValueError, SyntaxError) as e:
                print(f"Error processing {pdb_id}: {e}")
                continue

        summary_text = f"\nSummary Statistics for {feature_name}:\n" + "=" * 30 + "\n"
        summary_text += f"Total structures analyzed: {len(results)}\n"
        summary_text += f"Predicted features total: {total_residues['pred']}\n"
        summary_text += f"Experimental features total: {total_residues['exp']}\n"
        summary_text += f"Average predicted length: {np.mean(feat_lengths['pred']):.2f}\n"
        summary_text += f"Average experimental length: {np.mean(feat_lengths['exp']):.2f}\n"
        save_analysis_data(summary_text)

        # Create plots for this feature
        create_feature_plots(
            feature_name,
            feat_counts,
            feat_lengths,
            total_residues,
            count_pdbs,
            feat_list_pred,
            feat_list_exp,
            plot_creator,
            length_calculator
        )

def create_feature_plots(feature_name, feat_counts, feat_lengths, total_residues, count_pdbs,
                        feat_list_pred, feat_list_exp, plot_creator, length_calculator):
    """Create all plots for a given RNA feature"""

    # Get feature display name for plot labels
    feature_display = {
        'RNA_Stems': 'Stem',
        'RNA_HairpinLoops': 'Hairpin Loop',
        'RNA_InternalLoops': 'Internal Loop',
        'RNA_MultibranchLoops': 'Multibranch Loop',
        'RNA_Pseudoknots': 'Pseudoknot'
    }[feature_name]

    # if feature_name == 'RNA_Pseudoknots':
    #     print("\nDebugging length distribution data:")
    #     print(f"feat_lengths['pred']: {feat_lengths['pred']}")
    #     print(f"feat_lengths['exp']: {feat_lengths['exp']}")

    # 1. Distribution of feature counts
    plot_creator.get_violin_plot(
        table_source='rna_secondary_feature',
        data_values=[feat_counts['pred'], feat_counts['exp']],
        labels=['Predicted', 'Experimental'],
        name=f'{feature_name.lower()}_count_distribution',
        yAxis_label=f'Number of {feature_display}s'
    )

    # 2. Distribution of feature lengths
    plot_creator.get_violin_plot(
        table_source='rna_secondary_feature',
        data_values=[feat_lengths['pred'], feat_lengths['exp']],
        labels=['Predicted', 'Experimental'],
        name=f'{feature_name.lower()}_length_distribution',
        yAxis_label=f'{feature_display} Length'
    )

    # 3. Count correlation - Bubble plot
    if len(feat_counts['exp']) > 1 and len(feat_counts['pred']) > 1:
        plot_creator.get_bubble_scatter(
            table_source='rna_secondary_feature',
            xAxis_score=feat_counts['exp'],
            yAxis_score=feat_counts['pred'],
            xAxis_label=f'Experimental {feature_display} Count',
            yAxis_label=f'Predicted {feature_display} Count',
            name=f'{feature_name.lower()}_count_correlation_bubble'
        )

    # 4. Count correlation with PDB IDs - Bubble plot
    unique_points = list(count_pdbs.keys())
    if len(unique_points) > 1:
        plot_creator.get_bubble_scatter(
            table_source='rna_secondary_feature',
            xAxis_score=[x[0] for x in unique_points],
            yAxis_score=[x[1] for x in unique_points],
            xAxis_label=f'Experimental {feature_display} Count',
            yAxis_label=f'Predicted {feature_display} Count',
            name=f'{feature_name.lower()}_count_correlation_bubble_pdbId',
            pdb_ids=[count_pdbs[point] for point in unique_points]
        )

    # 5. Length correlation - Bubble plot
    if len(feat_lengths['exp']) > 1 and len(feat_lengths['pred']) > 1:
        plot_creator.get_bubble_scatter(
            table_source='rna_secondary_feature',
            xAxis_score=feat_lengths['exp'],
            yAxis_score=feat_lengths['pred'],
            xAxis_label=f'Experimental {feature_display} Length',
            yAxis_label=f'Predicted {feature_display} Length',
            name=f'{feature_name.lower()}_length_correlation_bubble'
        )

    # 6. Count vs Length plots
    exp_counts = []
    exp_lengths = []  # Individual lengths, not averages
    pred_counts = []
    pred_lengths = []  # Individual lengths, not averages

    # Calculate data points for structures that have features
    for pred_feat, exp_feat in zip(feat_list_pred, feat_list_exp):
        # Handle experimental structures
        if exp_feat and len(exp_feat) > 0:
            exp_counts.append(len(exp_feat))
            # Calculate total length based on feature type
            if feature_name == 'RNA_Pseudoknots':
                total_length = length_calculator(exp_feat)  # Calculate directly for pseudoknots
            else:
                total_length = sum(length_calculator(feat) for feat in exp_feat)
            exp_lengths.append(total_length)

        # Handle predicted structures
        if pred_feat and len(pred_feat) > 0:
            pred_counts.append(len(pred_feat))
            # Calculate total length based on feature type
            if feature_name == 'RNA_Pseudoknots':
                total_length = length_calculator(pred_feat)  # Calculate directly for pseudoknots
            else:
                total_length = sum(length_calculator(feat) for feat in pred_feat)
            pred_lengths.append(total_length)

    # Create scatter plots
    plot_creator.get_stem_length_count_plot(
        table_source='rna_secondary_feature',
        exp_counts=exp_counts,
        exp_lengths=exp_lengths,  # Now using total lengths per structure
        pred_counts=pred_counts,
        pred_lengths=pred_lengths,  # Now using total lengths per structure
        name=f'{feature_name.lower()}_count_vs_length_comparison',
        feature_name=feature_display
    )

    # Print summary statistics
    print(f"\n{feature_display} Analysis Summary:")
    print(f"Average {feature_display.lower()}s per structure: "
          f"Predicted={np.mean(feat_counts['pred']):.2f}, "
          f"Experimental={np.mean(feat_counts['exp']):.2f}")
    if feat_lengths['pred'] and feat_lengths['exp']:
        print(f"Average {feature_display.lower()} length: "
              f"Predicted={np.mean(feat_lengths['pred']):.2f}, "
              f"Experimental={np.mean(feat_lengths['exp']):.2f}")
    print(f"Total residues in {feature_display.lower()}s: "
          f"Predicted={total_residues['pred']}, "
          f"Experimental={total_residues['exp']}")

if __name__ == "__main__":
    if len(sys.argv) not in [2, 3, 4]:
        print("Usage: python secondaryFeatPlots.py <table_name> [single_chain|all] [+MSA|-MSA]")
        sys.exit(1)

    table = sys.argv[1]
    single_chain_only = sys.argv[2] == 'single_chain' if len(sys.argv) > 2 else False
    msa_option = sys.argv[3].lower() if len(sys.argv) > 3 else None  # Convert to lowercase
    
    if msa_option is not None:
        # Handle both +MSA/-MSA and +msa/-msa formats
        msa_option = '+MSA' if msa_option in ['+msa', '+msa'] else '-MSA' if msa_option in ['-msa', '-msa'] else None

    plot_creator = PlotCreator('rna_secondary_feature', msa_option, single_chain_only)

    config = StartConfig()
    dir = config.get_parent_folder()
    filepath = os.path.join(dir, 'secondary_feat_analysis.txt')
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(f"RNA Secondary Structure Feature Analysis\n")

    if table == "pred_protein_rna":
        compare_rna_features("pred_protein_rna", "exp_protein_rna", plot_creator)
        analyze_and_plot_motifs()
