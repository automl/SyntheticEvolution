import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import ast
import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig
from plots.plotCreator import PlotCreator

config = StartConfig()
# Change plot directory to be in the same folder as plots
plot_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results', 'plots', 'correlation_analysis')
os.makedirs(plot_dir, exist_ok=True)
plot_creator = PlotCreator("correlation_analysis")


def process_list_column(value, column_name=None):
    """Convert string representations of lists to numeric value
    For RNALength, takes the maximum value if it's a list
    For pTM and ipTM columns, takes average if multiple values
    For af3_rna_MSA, returns 0 if any element is 0
    """
    if pd.isna(value) or value == '':
        return None
    try:
        if isinstance(value, str) and value.startswith('['):
            # Parse string representation of list
            lst = ast.literal_eval(value)
            if isinstance(lst, (list, tuple)):
                if column_name == 'RNALength':
                    return float(max(lst))  # Take maximum for RNALength
                if len(lst) == 1:
                    return float(lst[0])
                return sum(float(x) for x in lst) / len(lst)
            # For af3_rna_MSA, return 0 if any element is 0
            if 0 in lst:
                return 0
            return 1
        return float(value)
    except:
        return None


def process_rna_inf(value):
    """Process RNA_INF list into its four components"""
    if pd.isna(value) or value == '':
        return None, None, None, None
    try:
        if isinstance(value, str) and value.startswith('['):
            lst = ast.literal_eval(value)
            if len(lst) == 4:
                return float(lst[0]), float(lst[1]), float(lst[2]), float(lst[3])
    except:
        pass
    return None, None, None, None


def calculate_conformation_similarity(pred_dict, exp_dict):
    """Calculate similarity score between predicted and experimental conformational distributions"""
    if not pred_dict or not exp_dict:
        return None
    
    try:
        # Convert string representations to dictionaries if needed
        if isinstance(pred_dict, str):
            pred_dict = ast.literal_eval(pred_dict)
        if isinstance(exp_dict, str):
            exp_dict = ast.literal_eval(exp_dict)
        
        # Get all unique keys
        all_keys = set(pred_dict.keys()) | set(exp_dict.keys())
        
        # Get total counts for normalization
        pred_total = sum(pred_dict.values())
        exp_total = sum(exp_dict.values())
        
        if pred_total == 0 or exp_total == 0:
            return None
        
        # Calculate normalized difference for each key
        differences = []
        for key in all_keys:
            pred_freq = pred_dict.get(key, 0) / pred_total
            exp_freq = exp_dict.get(key, 0) / exp_total
            differences.append(abs(pred_freq - exp_freq))
        
        # Return similarity score (1 - average difference)
        # Score ranges from 0 (completely different) to 1 (identical distributions)
        return 1 - (sum(differences) / len(differences))
    except:
        return None


def create_correlation_matrix(df, subset_name):
    """Create and save correlation matrix for a data subset"""
    # Skip correlation analysis if there's only one sample
    if len(df) <= 1:
        print(f"Skipping correlation analysis for {subset_name} (n={len(df)}): insufficient samples")
        return None, None

    # Calculate correlation matrix
    corr_matrix = df.corr(method='spearman')

    # Create heatmap
    plt.figure(figsize=(24, 20))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, fmt='.2f',
                square=True, linewidths=0.5, cbar_kws={"shrink": .5})
    plt.title(f'Correlation Matrix - {subset_name} (n={len(df)})', pad=20)
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()

    # Save plot in plot directory
    plot_creator.save_plot(f'correlation_matrix_{subset_name.lower().replace(" ", "_")}')
    plt.close()

    # Save correlation matrix to CSV in plot directory
    csv_path = os.path.join(plot_dir, f'correlation_matrix_{subset_name.lower().replace(" ", "_")}.csv')
    corr_matrix.to_csv(csv_path)

    # Define excluded correlation pairs
    excluded_pairs = {
        frozenset(['Complex_RMDS', 'Protein_RMSD']),
        frozenset(['RNA_RMSD', 'Protein_RMSD']),
        frozenset(['Complex_RMDS', 'RNA_RMSD']),
        frozenset(['Binding_affinity_kd', 'Free_energy']),
        frozenset(['af3_protein_ipTM', 'af3_rna_ipTM']),
        frozenset(['af3_global_pae_avg', 'af3_chain_pair_pae_min']),
        frozenset(['NumberRNAs', 'ProteinLength'])
    }
    
    # Add RNA metric combinations to excluded pairs
    rna_metrics = ['RNA_f1_score', 'RNA_precision', 'RNA_recall', 'RNA_mcc', 'RNA_wl']
    for i in range(len(rna_metrics)):
        for j in range(i + 1, len(rna_metrics)):
            excluded_pairs.add(frozenset([rna_metrics[i], rna_metrics[j]]))

    # Add af3_ranking_score combinations to excluded pairs
    af3_metrics = [
        'af3_protein_pTM', 'af3_rna_pTM', 'af3_protein_ipTM', 'af3_rna_ipTM',
        'af3_protein_pLDDT_avg', 'af3_rna_pLDDT_avg', 'af3_global_pae_avg',
        'af3_chain_pair_pae_min', 'af3_fraction_disordered', 'af3_has_clash'
    ]
    for metric in af3_metrics:
        excluded_pairs.add(frozenset(['af3_ranking_score', metric]))

    # Find and save strong correlations (using absolute values)
    correlations = []
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            var1 = corr_matrix.columns[i]
            var2 = corr_matrix.index[j]
            corr_value = corr_matrix.iloc[i, j]
            
            # Check if correlation pair is in excluded pairs using absolute value
            if (abs(corr_value) > 0.5 and abs(corr_value) < 1 and 
                frozenset([var1, var2]) not in excluded_pairs):
                correlations.append({
                    'Variable 1': var1,
                    'Variable 2': var2,
                    'Correlation': corr_value,
                    'Abs_Correlation': abs(corr_value)
                })

    corr_df = pd.DataFrame(correlations)
    if not correlations:
        print(f"No strong correlations found for {subset_name}")
        return corr_matrix, None

    # Sort by absolute correlation value
    corr_df = corr_df.sort_values('Abs_Correlation', ascending=False)
    
    # Save correlations to CSV and text file in plot directory
    csv_path = os.path.join(plot_dir, f'correlations_{subset_name.lower().replace(" ", "_")}.csv')
    txt_path = os.path.join(plot_dir, f'correlations_{subset_name.lower().replace(" ", "_")}.txt')
    
    corr_df.to_csv(csv_path, index=False)

    with open(txt_path, 'w') as f:
        f.write(f"Strong correlations (|correlation| > 0.5) for {subset_name}:\n\n")
        for _, row in corr_df.iterrows():
            sign = "+" if row['Correlation'] > 0 else "-"
            f.write(f"{row['Variable 1']} - {row['Variable 2']}: {sign}{abs(row['Correlation']):.3f}\n")

    return corr_matrix, corr_df


def get_correlation_matrices(filter_type='rna_length'):
    """
    Create correlation matrices based on different filters.
    
    Args:
        filter_type (str): Type of filter to apply ('all', 'rna_length', 'rna_family', 'domain_names', 
                          'rna_complexity', 'gc_content', 'rna_msa', 'complex_type')
    """
    # Connect to database
    db = DatabaseMethods()

    # Add domain_names and RNAFamily to the query
    pred_columns = [
        "p.exp_db_id", "p.ProteinLength", "p.NumberProteins", "p.RNALength", "p.NumberRNAs",
        "p.RNAMotifLength", "p.Hbond_proteinRNA", "p.vdWbond_proteinRNA",
        "p.ProteinRNAInterfaceArea", "p.ProteinRNAInterfaceRatio", "p.Free_energy",
        "p.Binding_affinity_kd", "p.af3_protein_pTM", "p.af3_rna_pTM", "p.af3_protein_ipTM",
        "p.af3_rna_ipTM", "p.af3_protein_pLDDT_avg", "p.af3_rna_pLDDT_avg",
        "p.af3_global_pae_avg", "p.af3_chain_pair_pae_min", "p.af3_fraction_disordered",
        "p.af3_has_clash", "p.af3_ranking_score", "p.af3_rna_MSA", "p.Complex_RMSD",
        "p.Protein_RMSD", "p.RNA_RMSD", "p.Protein_LDDT", "p.RNA_LDDT", "p.Protein_TM",
        "p.RNA_TM", "p.RNA_DI", "p.RNA_INF",
        "p.RNA_f1_score", "p.RNA_precision", "p.RNA_recall", "p.RNA_mcc", "p.RNA_wl",
        "p.RNAmotif_score", "p.RNAMotif_isDuplex", "e.RNAMotif_isDuplex as e_RNAMotif_isDuplex",
        "p.RNA_isDuplex", "e.RNA_isDuplex as e_RNA_isDuplex",
        "p.RNA_GlycosidicBond as p_RNA_GlycosidicBond", "e.RNA_GlycosidicBond as e_RNA_GlycosidicBond",
        "p.RNA_SugarPucker as p_RNA_SugarPucker", "e.RNA_SugarPucker as e_RNA_SugarPucker",
        "p.RNA_GammaAngle as p_RNA_GammaAngle", "e.RNA_GammaAngle as e_RNA_GammaAngle",
        "e.domain_counts", "e.domain_names", "e.RNAFamily",
        "p.RNA_ElectrostaticPotential as RNA_ElPotential", 
        "e.RNA_ElectrostaticPotential as e_RNA_ElPotential",
        "p.RNA_GC_Content as RNA_GC_Content",
        "p.RNA_AU_Content as RNA_AU_Content",
        "p.RNA_SequenceComplexity", "e.exp_pH as exp_pH"
    ]

    # Query data with JOIN
    query = f"""
    SELECT {', '.join(pred_columns)}
    FROM pred_protein_rna p
    JOIN exp_protein_rna e ON p.exp_db_id = e.PDBId
    """
    df = pd.read_sql_query(query, db.connection)

    # Process the data as before
    df[['RNA_inf', 'RNA_wc', 'RNA_nonWc', 'RNA_stack']] = pd.DataFrame(
        df['RNA_INF'].apply(process_rna_inf).tolist(),
        index=df.index
    )

    # Calculate RNA_ElPotential_Diff
    def calculate_elpotential_diff(row):
        pred_val = row['RNA_ElPotential']
        exp_val = row['e_RNA_ElPotential']
        if pd.isna(pred_val) or pd.isna(exp_val):
            return None
        # Calculate percentage difference and normalize to 0-1
        diff = abs(pred_val - exp_val) / max(abs(pred_val), abs(exp_val))
        return 1 / (1 + np.exp(-5 * (diff - 0.5)))  # Sigmoid normalization

    df['RNA_ElPotential_Diff'] = df.apply(calculate_elpotential_diff, axis=1)

    # Process all numeric columns
    list_columns = ['Complex_RMSD', 'Protein_RMSD', 'RNA_RMSD', 'Hbond_proteinRNA',
                    'vdWbond_proteinRNA', 'af3_chain_pair_pae_min',
                    'af3_protein_pTM', 'af3_rna_pTM', 'af3_protein_ipTM', 'af3_rna_ipTM',
                    'af3_rna_MSA', 'RNALength']

    for col in list_columns:
        df[col] = df[col].apply(lambda x: process_list_column(x, col))

    # Process other columns as before
    df['domain_counts'] = df['domain_counts'].apply(process_domain_counts)
    df['RNADuplex_Mismatch'] = df.apply(calculate_duplex_mismatch, axis=1)
    df['RNAMotifDuplex_Mismatch'] = df.apply(calculate_motif_duplex_mismatch, axis=1)
    df['GlycosidicBond_Similarity'] = df.apply(lambda row: calculate_conformation_similarity(row['p_RNA_GlycosidicBond'], row['e_RNA_GlycosidicBond']), axis=1)
    df['SugarPucker_Similarity'] = df.apply(lambda row: calculate_conformation_similarity(row['p_RNA_SugarPucker'], row['e_RNA_SugarPucker']), axis=1)
    df['GammaAngle_Similarity'] = df.apply(lambda row: calculate_conformation_similarity(row['p_RNA_GammaAngle'], row['e_RNA_GammaAngle']), axis=1)

    # Drop non-numeric columns
    columns_to_drop = ['exp_db_id', 'e_RNA_isDuplex', 'e_RNAMotif_isDuplex',
                      'p_RNA_GlycosidicBond', 'e_RNA_GlycosidicBond',
                      'p_RNA_SugarPucker', 'e_RNA_SugarPucker',
                      'p_RNA_GammaAngle', 'e_RNA_GammaAngle',
                      'RNA_INF', 'domain_names', 'RNAFamily',
                      'e_RNA_ElPotential']  # Drop the experimental ElPotential after calculating diff
    numeric_df = df.drop(columns_to_drop, axis=1)
    numeric_df = numeric_df.apply(pd.to_numeric, errors='coerce')

    # First, process all data without filtering
    if filter_type == 'all':
        create_correlation_matrix(numeric_df, "All_Data")

    elif filter_type == 'rna_length':
        filters = [
            ('RNA Length 0-10', lambda x: x <= 10),
            ('RNA Length 10-20', lambda x: 10 < x <= 20),
            ('RNA Length 20-60', lambda x: 20 < x <= 60),
            ('RNA Length >60', lambda x: x > 60)
        ]
        for name, condition in filters:
            subset = numeric_df[numeric_df['RNALength'].apply(lambda x: condition(x) if pd.notna(x) else False)]
            if len(subset) > 1:
                create_correlation_matrix(subset, name)

    elif filter_type == 'rna_family':
        # Process RNA families
        def normalize_rna_family(x):
            if pd.isna(x):
                return 'None'
            if isinstance(x, str):
                if x.startswith('['):
                    try:
                        lst = ast.literal_eval(x)
                        if all(item == 'None' for item in lst):
                            return 'None'
                        return next((item for item in lst if item != 'None'), 'None')
                    except:
                        return x
            return str(x)

        # Get unique normalized RNA families and their counts
        df['normalized_rna_family'] = df['RNAFamily'].apply(normalize_rna_family)
        family_counts = df['normalized_rna_family'].value_counts()
        
        # Filter families with 3 or more samples and get top 5
        common_families = family_counts[family_counts >= 3].head(5).index.tolist()
        print(f"\nMost common RNA families (min 3 samples):")
        for family in common_families:
            print(f"  {family}: {family_counts[family]} samples")
        
        for family in common_families:
            subset = numeric_df[df['normalized_rna_family'] == family]
            if len(subset) > 1:  # Only process if we have more than one sample
                create_correlation_matrix(subset, f"RNA_Family_{family}")

    elif filter_type == 'domain_names':
        # Process domains
        def extract_domains(x):
            if pd.isna(x):
                return []
            # Filter out "-" and empty strings, and strip space
            return [d.strip() for d in str(x).split(';') if d.strip() and d.strip() != '-']

        # Get all domains and their counts
        domain_counts = {}
        for domains in df['domain_names'].dropna():
            for domain in extract_domains(domains):
                domain_counts[domain] = domain_counts.get(domain, 0) + 1

        # Filter domains with 3 or more samples and get top 5
        common_domains = sorted([(count, domain) for domain, count in domain_counts.items() if count >= 3], 
                              reverse=True)[:5]
        
        print(f"\nMost common domains (min 3 samples):")
        for count, domain in common_domains:
            print(f"  {domain}: {count} samples")

        # Process each common domain
        for _, domain in common_domains:
            subset = numeric_df[df['domain_names'].apply(
                lambda x: domain in extract_domains(x) if pd.notna(x) else False
            )]
            
            if len(subset) > 1:  # Only process if we have more than one sample
                create_correlation_matrix(subset, f"Domain_{domain}")

    elif filter_type == 'rna_complexity':
        # Process RNA sequence complexity ranges
        complexity_ranges = [
            ('RNA Complexity 0.00-0.39', lambda x: 0.00 <= x <= 0.39),
            ('RNA Complexity 0.40-0.69', lambda x: 0.40 <= x <= 0.69),
            ('RNA Complexity 0.70-1.00', lambda x: 0.70 <= x <= 1.00)
        ]
        
        for name, condition in complexity_ranges:
            subset = numeric_df[numeric_df['RNA_SequenceComplexity'].apply(lambda x: condition(x) if pd.notna(x) else False)]
            if len(subset) > 1:
                create_correlation_matrix(subset, name)

    elif filter_type == 'gc_content':
        # Process GC content ranges
        gc_ranges = [
            ('GC Content 0.00-0.50', lambda x: 0.00 <= x <= 0.50),
            ('GC Content 0.51-1.00', lambda x: 0.51 <= x <= 1.00)
        ]
        
        for name, condition in gc_ranges:
            subset = numeric_df[numeric_df['RNA_GC_Content'].apply(lambda x: condition(x) if pd.notna(x) else False)]
            if len(subset) > 1:
                create_correlation_matrix(subset, name)

    elif filter_type == 'au_content':
        # Process GC content ranges
        au_ranges = [
            ('AU Content 0.00-0.50', lambda x: 0.00 <= x <= 0.50),
            ('AU Content 0.51-1.00', lambda x: 0.51 <= x <= 1.00)
        ]

        for name, condition in au_ranges:
            subset = numeric_df[numeric_df['RNA_AU_Content'].apply(lambda x: condition(x) if pd.notna(x) else False)]
            if len(subset) > 1:
                create_correlation_matrix(subset, name)

    elif filter_type == 'rna_msa':
        # Process based on RNA MSA value
        def get_msa_value(x):
            if pd.isna(x):
                return None
            try:
                if isinstance(x, str) and x.startswith('['):
                    lst = ast.literal_eval(x)
                    return 0 if 0 in lst else 1
                return float(x)
            except:
                return None

        df['msa_value'] = df['af3_rna_MSA'].apply(get_msa_value)
        
        # Create subsets for MSA=0 and MSA=1
        for msa_value in [0, 1]:
            subset = numeric_df[df['msa_value'] == msa_value]
            if len(subset) > 1:
                create_correlation_matrix(subset, f"RNA_MSA_{msa_value}")

    elif filter_type == 'complex_type':
        # Process based on number of proteins and RNAs
        def is_single_complex(row):
            return 1 if pd.notna(row['NumberProteins']) and pd.notna(row['NumberRNAs']) and row['NumberProteins'] == 1 and row['NumberRNAs'] == 1 else 0

        df['is_single_complex'] = df.apply(is_single_complex, axis=1)
        
        # Create subsets for single and multi-complex
        for complex_type in [0, 1]:
            subset = numeric_df[df['is_single_complex'] == complex_type]
            if len(subset) > 1:
                name = "Single_Complex" if complex_type == 1 else "Multi_Complex"
                create_correlation_matrix(subset, name)

    return db  # Return the database connection to be closed in main()


def process_domain_counts(x):
    """Process domain_counts column to return total number of domains"""
    if pd.notna(x):
        try:
            # Split by semicolon and filter out any non-numeric values
            counts = [int(count) for count in x.split(';') if count.isdigit()]
            return sum(counts) if counts else None
        except:
            return None
    return None


def calculate_duplex_mismatch(row):
    """Calculate mismatch between predicted and experimental RNA_isDuplex"""
    pred_duplex = row['RNA_isDuplex']
    exp_duplex = row['e_RNA_isDuplex']
    if pd.isna(pred_duplex) or pd.isna(exp_duplex):
        return None
    return 0 if pred_duplex == exp_duplex else 1


def calculate_motif_duplex_mismatch(row):
    """Calculate mismatch between predicted and experimental RNAMotif_isDuplex"""
    pred_duplex = row['RNAMotif_isDuplex']
    exp_duplex = row['e_RNAMotif_isDuplex']
    if pd.isna(pred_duplex) or pd.isna(exp_duplex):
        return None
    return 0 if pred_duplex == exp_duplex else 1


def main():
    """Main function to create correlation matrices"""
    try:
        db = None
        print("\nCreating correlation matrix for all data...")
        db = get_correlation_matrices('all')
        
        print("\nCreating correlation matrices by RNA Length...")
        db = get_correlation_matrices('rna_length')
        
        print("\nCreating correlation matrices by RNA Family...")
        db = get_correlation_matrices('rna_family')
        
        print("\nCreating correlation matrices by Domain Names...")
        db = get_correlation_matrices('domain_names')
        
        print("\nCreating correlation matrices by RNA Sequence Complexity...")
        db = get_correlation_matrices('rna_complexity')
        
        print("\nCreating correlation matrices by GC Content...")
        db = get_correlation_matrices('gc_content')

        # print("\nCreating correlation matrices by AU Content...")
        # db = get_correlation_matrices('au_content')

        print("\nCreating correlation matrices by RNA MSA...")
        db = get_correlation_matrices('rna_msa')
        
        print("\nCreating correlation matrices by Complex Type...")
        db = get_correlation_matrices('complex_type')
    finally:
        if db:
            db.connection.close()


if __name__ == "__main__":
    main()