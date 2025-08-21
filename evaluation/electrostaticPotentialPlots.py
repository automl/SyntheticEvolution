import os
import sys
import pandas as pd
import numpy as np
from plots.plotCreator import PlotCreator
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig
import ast

config = StartConfig()
database = DatabaseMethods()

class ElectrostaticPotentialPlotter:
    def __init__(self, database, config, msa_option=None, single_chain_only=False):
        """Initialize the plotter with database connection and configuration"""
        self.database = database
        self.config = config
        self.msa_option = msa_option
        self.single_chain_only = single_chain_only
        # Convert boolean msa_option to string format for PlotCreator
        plot_msa_option = '+MSA' if msa_option else '-MSA' if msa_option is not None else None
        self.plot_creator = PlotCreator('electrostatic_potential', plot_msa_option, single_chain_only)
        
        # Get data from database using the connection
        self.pred_df = pd.read_sql_query("SELECT * FROM pred_protein_rna", database.connection)
        self.exp_df = pd.read_sql_query("SELECT * FROM exp_protein_rna", database.connection)
        
        # Filter data based on single_chain option
        if single_chain_only:
            self.exp_df, self.pred_df = self.database.filter_by_singleChain(self.exp_df, self.pred_df)
        
        # Filter data based on MSA option
        if msa_option is not None:
            print("\nFiltering based on MSA option...")
            msa_str = '+MSA' if msa_option else '-MSA'
            self.pred_df = self.database.filter_by_msa(self.pred_df, msa_str)
            # Filter experimental data to match
            self.exp_df = self.exp_df[
                self.exp_df['PDBId'].isin(self.pred_df['exp_db_id'])
            ]
            print(f"Remaining entries after MSA filtering: {len(self.exp_df)}")
        
        # Print available columns for debugging
        print("\nColumns in pred_df:")
        print(self.pred_df.columns.tolist())
        print("\nColumns in exp_df:")
        print(self.exp_df.columns.tolist())
        
        # Merge dataframes with all necessary columns
        columns_to_merge = [
            'PDBId',
            'RNA_ElectrostaticPotential',
            'ProteinRNAInterfaceArea',
            'Binding_affinity_kd',
            'RNA_GlycosidicBond',
            'RNA_SugarPucker',
            'RNA_GammaAngle'
        ]
        
        # Merge dataframes
        self.merged_df = pd.merge(
            self.pred_df,
            self.exp_df[columns_to_merge],
            left_on='exp_db_id',
            right_on='PDBId',
            how='left',
            suffixes=('_pred', '_exp')
        )
        
        # Print columns after merge for debugging
        print("\nColumns after merge:")
        print(self.merged_df.columns.tolist())
        
        # Calculate similarity metrics
        self.calculate_conformation_similarity()
        self.extract_rna_inf_metrics()
        
    def parse_dict_string(self, dict_string):
        """Parse string representation of dictionary into actual dictionary"""
        try:
            if pd.isna(dict_string):
                return {}
            # Remove curly braces and split by comma
            dict_str = dict_string.strip('{}')
            pairs = dict_str.split(',')
            result = {}
            
            for pair in pairs:
                if ':' in pair:
                    key, value = pair.split(':')
                    # Clean up the key and value
                    key = key.strip().strip("'").strip('"')
                    value = int(value.strip())
                    result[key] = value
            return result
        except:
            return {}
            
    def calculate_conformation_similarity(self):
        """Calculate similarity metrics for RNA conformations"""
        # Calculate GlycosidicBond Similarity
        self.merged_df['GlycosidicBond_Similarity'] = self.merged_df.apply(
            lambda row: self.calculate_similarity(
                self.parse_dict_string(row['RNA_GlycosidicBond_pred']),
                self.parse_dict_string(row['RNA_GlycosidicBond_exp'])
            ), axis=1
        )
        
        # Calculate SugarPucker Similarity
        self.merged_df['SugarPucker_Similarity'] = self.merged_df.apply(
            lambda row: self.calculate_similarity(
                self.parse_dict_string(row['RNA_SugarPucker_pred']),
                self.parse_dict_string(row['RNA_SugarPucker_exp'])
            ), axis=1
        )
        
        # Calculate GammaAngle Similarity
        self.merged_df['GammaAngle_Similarity'] = self.merged_df.apply(
            lambda row: self.calculate_similarity(
                self.parse_dict_string(row['RNA_GammaAngle_pred']),
                self.parse_dict_string(row['RNA_GammaAngle_exp'])
            ), axis=1
        )
        
        # Print sample data for debugging
        print("\nSample data for conformation metrics:")
        print("RNA_GlycosidicBond_pred sample:", self.merged_df['RNA_GlycosidicBond_pred'].head())
        print("RNA_GlycosidicBond_exp sample:", self.merged_df['RNA_GlycosidicBond_exp'].head())
        print("GlycosidicBond_Similarity sample:", self.merged_df['GlycosidicBond_Similarity'].head())
    
    def calculate_similarity(self, pred_dict, exp_dict):
        """Calculate similarity between two dictionaries of counts"""
        if not pred_dict or not exp_dict:
            return None
            
        total_min = 0
        total_max = 0
        
        # Get all unique keys
        all_keys = set(pred_dict.keys()) | set(exp_dict.keys())
        
        for key in all_keys:
            pred_count = pred_dict.get(key, 0)
            exp_count = exp_dict.get(key, 0)
            total_min += min(pred_count, exp_count)
            total_max += max(pred_count, exp_count)
        
        return total_min / total_max if total_max > 0 else 0
    
    def extract_rna_inf_metrics(self):
        """Extract RNA INF metrics from the RNA_INF column"""
        def extract_metrics(inf_str):
            try:
                if pd.isna(inf_str):
                    return None, None, None
                # Parse the string into a list
                inf_list = ast.literal_eval(inf_str)
                if len(inf_list) >= 4:
                    return inf_list[1], inf_list[2], inf_list[3]  # wc, nonWc, stack
                return None, None, None
            except:
                return None, None, None
        
        # Apply extraction to RNA_INF column
        self.merged_df[['RNA_wc', 'RNA_nonWc', 'RNA_stack']] = pd.DataFrame(
            self.merged_df['RNA_INF'].apply(extract_metrics).tolist(),
            index=self.merged_df.index
        )

    def process_rmsd(self, rmsd_str):
        """Process RMSD string to get first value"""
        try:
            if isinstance(rmsd_str, str):
                # Handle list-like string
                if '[' in rmsd_str:
                    values = ast.literal_eval(rmsd_str)
                    return float(values[0]) if values else None
                return float(rmsd_str)
            elif isinstance(rmsd_str, (int, float)):
                return float(rmsd_str)
            return None
        except:
            return None
            
    def process_binding_affinity(self, value):
        """Process binding affinity to log scale"""
        try:
            value = float(value)
            if value > 0:  # Only take log of positive values
                return np.log10(value)
            return None
        except:
            return None

    def calculate_differences(self):
        """Calculate differences between predicted and experimental values"""
        # Calculate RNA_ElPotential_Diff
        self.merged_df['RNA_ElPotential_Diff'] = (
            pd.to_numeric(self.merged_df['RNA_ElectrostaticPotential_pred'], errors='coerce') -
            pd.to_numeric(self.merged_df['RNA_ElectrostaticPotential_exp'], errors='coerce')
        )
        
        # Calculate ProteinRNAInterfaceArea_Diff
        self.merged_df['ProteinRNAInterfaceArea_Diff'] = (
            pd.to_numeric(self.merged_df['ProteinRNAInterfaceArea_pred'], errors='coerce') -
            pd.to_numeric(self.merged_df['ProteinRNAInterfaceArea_exp'], errors='coerce')
        )
        
        # Process binding affinity to log scale and calculate difference
        self.merged_df['Binding_affinity_kd_pred_log'] = self.merged_df['Binding_affinity_kd_pred'].apply(self.process_binding_affinity)
        self.merged_df['Binding_affinity_kd_exp_log'] = self.merged_df['Binding_affinity_kd_exp'].apply(self.process_binding_affinity)
        self.merged_df['Binding_affinity_kd_Diff'] = (
            self.merged_df['Binding_affinity_kd_pred_log'] - self.merged_df['Binding_affinity_kd_exp_log']
        )
        
        # Process RNA_RMSD to get first value
        self.merged_df['RNA_RMSD'] = self.merged_df['RNA_RMSD'].apply(self.process_rmsd)
        
        # Remove rows where RNA_ElectrostaticPotential is None
        self.merged_df = self.merged_df.dropna(subset=['RNA_ElPotential_Diff'])
        
        # Print available columns after processing
        print("\nColumns after processing:")
        print(self.merged_df.columns.tolist())
        
    def create_correlation_plots(self):
        """Create all correlation plots for RNA electrostatic potential differences"""
        # First, create direct comparison plot for RNA_ElectrostaticPotential
        x_values = pd.to_numeric(self.merged_df['RNA_ElectrostaticPotential_exp'], errors='coerce')
        y_values = pd.to_numeric(self.merged_df['RNA_ElectrostaticPotential_pred'], errors='coerce')
        mask = ~(x_values.isna() | y_values.isna())
        
        if mask.any():
            corr = np.corrcoef(x_values[mask], y_values[mask])[0, 1]
            self.plot_creator.get_scatterplot(
                table_source='electrostatic_potential',
                xAxis_score=x_values[mask],
                xAxis_label='Experimental RNA Electrostatic Potential',
                yAxis_label=f'Predicted RNA Electrostatic Potential (r={corr:.2f})',
                name='el_potential_pred_vs_exp',
                yAxis_score=y_values[mask]
            )
        
        # List of metrics to plot against RNA_ElPotential_Diff
        metrics = {
            'GlycosidicBond_Similarity': 'Glycosidic Bond Similarity',
            'SugarPucker_Similarity': 'Sugar Pucker Similarity',
            'GammaAngle_Similarity': 'Gamma Angle Similarity',
            'RNA_RMSD': 'RNA RMSD',
            'RNA_LDDT': 'RNA LDDT',
            'RNA_TM': 'RNA TM-score',
            'RNA_GDT_TS': 'RNA GDT-TS',
            'RNA_wc': 'RNA Watson-Crick',
            'RNA_nonWc': 'RNA non-Watson-Crick',
            'RNA_stack': 'RNA Stacking',
            'Binding_affinity_kd_Diff': 'Binding Affinity Difference (log₁₀)',
            'ProteinRNAInterfaceArea_Diff': 'Interface Area Difference'
        }
        
        print("\nGenerating correlation plots...")
        for metric, label in metrics.items():
            if metric in self.merged_df.columns:
                print(f"Processing {metric}...")
                print(f"Sample data for {metric}:")
                print(self.merged_df[metric].head())
                
                x_values = pd.to_numeric(self.merged_df[metric], errors='coerce')
                y_values = pd.to_numeric(self.merged_df['RNA_ElPotential_Diff'], errors='coerce')
                
                # Drop NaN values
                mask = ~(x_values.isna() | y_values.isna())
                x_values = x_values[mask]
                y_values = y_values[mask]
                
                if len(x_values) > 0:
                    # Calculate correlation coefficient
                    corr = np.corrcoef(x_values, y_values)[0, 1]
                    print(f"  Found {len(x_values)} valid data points, correlation: {corr:.2f}")
                    
                    # Use the new difference scatterplot function
                    self.plot_creator.get_difference_scatterplot(
                        table_source='electrostatic_potential',
                        xAxis_score=x_values,
                        yAxis_score=y_values,
                        xAxis_label=label,
                        yAxis_label='RNA Electrostatic Potential Difference',
                        name=f'el_potential_diff_vs_{metric.lower()}',
                        diff_values=y_values  # The difference values determine the color coding
                    )
                else:
                    print(f"  No valid data points for correlation between RNA_ElPotential_Diff and {metric}")
            else:
                print(f"Warning: {metric} not found in dataframe")

def main():
    """Main function to generate all electrostatic potential plots"""
    if len(sys.argv) not in [2, 3]:
        print("Usage: python electrostaticPotentialPlots.py <all|single_chain> [+MSA|-MSA]")
        print("Options:")
        print("  all          : Generate plots for all chain combinations")
        print("  single_chain : Generate plots for single chain pairs only")
        print("  +MSA/+msa    : Plot only structures with MSA")
        print("  -MSA/-msa    : Plot only structures without MSA")
        sys.exit(1)

    # Parse command line arguments
    chain_option = sys.argv[1].lower()
    if chain_option not in ['all', 'single_chain']:
        print("Error: First argument must be either 'all' or 'single_chain'")
        sys.exit(1)

    # Set MSA option
    msa_option = None
    if len(sys.argv) == 3:
        msa_arg = sys.argv[2].lower()
        if msa_arg in ['+msa', '+msa']:
            msa_option = True
        elif msa_arg in ['-msa', '-msa']:
            msa_option = False
        else:
            print("Error: MSA argument must be either '+MSA'/'+msa' or '-MSA'/'-msa'")
            sys.exit(1)

    single_chain_only = (chain_option == 'single_chain')

    print("Initializing ElectrostaticPotentialPlotter...")
    plotter = ElectrostaticPotentialPlotter(database, config, msa_option, single_chain_only)
    
    print("Calculating differences...")
    plotter.calculate_differences()
    
    print("Creating correlation plots...")
    plotter.create_correlation_plots()
    
    print("Done!")

if __name__ == "__main__":
    main() 