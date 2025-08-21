import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import ast
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig
import os
from plots.plotCreator import PlotCreator

config = StartConfig()
plot_creator = PlotCreator("correlation_analysis")


def process_list_column(value):
    """Convert string representations of lists to first numeric value"""
    if pd.isna(value) or value == '':
        return None
    try:
        if isinstance(value, str) and value.startswith('['):
            # Parse string representation of list
            lst = ast.literal_eval(value)
            # For pTM and ipTM columns, take average if multiple values
            if isinstance(lst, (list, tuple)):
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


def get_correlation_matrix():
    # Connect to database
    db = DatabaseMethods()

    # Combined query to get both prediction and experimental data
    query = """
        SELECT p.exp_db_id, p.ProteinLength, p.RNALength, 
               p.RNAMotifLength, p.Hbond_proteinRNA, p.vdWbond_proteinRNA, 
               p.ProteinRNAInterfaceArea, p.ProteinRNAInterfaceRatio, p.Free_energy,
               p.Binding_affinity_kd, p.af3_protein_pTM, p.af3_rna_pTM, p.af3_protein_ipTM,
               p.af3_rna_ipTM, p.af3_protein_pLDDT_avg, p.af3_rna_pLDDT_avg,
               p.af3_global_pae_avg, p.af3_chain_pair_pae_min, p.af3_fraction_disordered,
               p.af3_has_clash, p.af3_ranking_score, 
               p.RNAmotif_score, p.RNAMotif_inExp,
               e.RNAMotif_z_score as exp_RNAMotif_z_score, 
               e.BindingIntensity as exp_BindingIntensity
        FROM pred_protein_rna p
        INNER JOIN exp_protein_rna e ON UPPER(p.exp_db_id) = UPPER(e.pdbid)
    """

    # Get data using pd.read_sql_query
    df = pd.read_sql_query(query, db.connection)

    # Process list columns
    list_columns = ['Hbond_proteinRNA',
                    'vdWbond_proteinRNA', 'af3_chain_pair_pae_min',
                    'af3_protein_pTM', 'af3_rna_pTM', 'af3_protein_ipTM', 'af3_rna_ipTM']

    for col in list_columns:
        df[col] = df[col].apply(process_list_column)

    # Drop non-numeric columns and rows with NaN
    df = df.drop('exp_db_id', axis=1)
    df = df.apply(pd.to_numeric, errors='coerce')

    # print("DataFrame info before correlation:")
    # print(df.info())
    # print("\nNumber of non-null values per column:")
    # print(df.count())

    # Drop rows where all values are NaN
    df = df.dropna(how='all')

    # Calculate correlation matrix
    corr_matrix = df.corr(method='spearman')  # Using Spearman correlation for robustness

    # Print correlation matrix info
    print("\nCorrelation matrix shape:", corr_matrix.shape)
    print("Number of non-null values in correlation matrix:", corr_matrix.count().sum())

    # Create heatmap only if we have valid correlations
    if not corr_matrix.isna().all().all():
        plt.figure(figsize=(20, 16))
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, fmt='.2f',
                    square=True, linewidths=0.5, cbar_kws={"shrink": .5})
        plt.title('Correlation Matrix of Protein-RNA Interface Metrics')
        plt.xticks(rotation=90)
        plt.yticks(rotation=0)
        plt.tight_layout()

        # Replace direct plot saving with plot_creator.save_plot
        plot_creator.save_plot('correlation_matrix')
        plt.close()

        # Save correlation matrix to CSV
        corr_matrix.to_csv(os.path.join(config.parent_folder, 'correlation_matrix.csv'))

    # Find correlations (|correlation| > 0.5)
    correlations = []
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if not pd.isna(corr_matrix.iloc[i, j]) and abs(corr_matrix.iloc[i, j]) > 0.5:
                correlations.append({
                    'Variable_1': corr_matrix.columns[i],
                    'Variable_2': corr_matrix.index[j],
                    'Correlation': corr_matrix.iloc[i, j]
                })

    # Save correlations to CSV and text file
    if correlations:
        corr_df = pd.DataFrame(correlations)
        corr_df = corr_df.sort_values('Correlation', key=abs, ascending=False)
        corr_df.to_csv(os.path.join(config.parent_folder, 'correlations.csv'), index=False)

        # Write to text file
        with open(os.path.join(config.parent_folder, 'correlations.txt'), 'w') as f:
            f.write("Strong correlations (|correlation| > 0.5):\n\n")
            for _, row in corr_df.iterrows():
                f.write(f"{row['Variable_1']} - {row['Variable_2']}: {row['Correlation']:.3f}\n")
    else:
        print("No strong correlations found")
        corr_df = pd.DataFrame(columns=['Variable_1', 'Variable_2', 'Correlation'])

    return corr_matrix, corr_df


if __name__ == "__main__":
    try:
        corr_matrix, correlations = get_correlation_matrix()
        if not correlations.empty:
            print("\nStrong correlations (|correlation| > 0.5):")
            print(correlations)
        else:
            print("\nNo strong correlations found")
    except Exception as e:
        print(f"An error occurred: {str(e)}")