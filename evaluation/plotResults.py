import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ast import literal_eval
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig
from plots.plotCreator import PlotCreator
import json
import numpy as np
import re
import argparse
from collections import OrderedDict
import ast


def safe_eval(val):
    if isinstance(val, str):
        val = val.strip()
        if (val.startswith('[') and val.endswith(']')) or (val.startswith('{') and val.endswith('}')):
            try:
                return ast.literal_eval(val)
            except Exception as e:
                print(f"Error in safe_eval: {e} for value: {val}")
                return None
        return val
    return val

def motif_ipTM_match_score(row):
    """
    For multi-chain complexes (NumberRNAs > 1), check if the chain with the largest experimental motif length
    matches the chain with the highest predicted RNA ipTM, based on sequence identity.
    Returns 1 if they match, otherwise returns (100 - percent difference)/100,
    where percent difference is relative to the highest predicted RNA ipTM.
    Returns np.nan for rows with NumberRNAs <= 1 or on error.
    """

    try:
        # Only process if NumberRNAs > 1
        if 'NumberRNAs' not in row or row['NumberRNAs'] <= 1:
            return np.nan

        # Use experimental columns
        motif_col = None
        for col in row.index:
            if 'RNAMotifLength_exp' in col and not pd.isna(row[col]):
                motif_col = col
                break
        if motif_col is None:
            return np.nan
        motif_lengths = safe_eval(row[motif_col]) if isinstance(row[motif_col], str) else row[motif_col]

        chain_pairs_col = None
        for col in row.index:
            if 'ChainIDpairList_proteinRNA_exp' in col and not pd.isna(row[col]):
                chain_pairs_col = col
                break
        if chain_pairs_col is None:
            return np.nan
        chain_pairs = safe_eval(row[chain_pairs_col]) if isinstance(row[chain_pairs_col], str) else row[chain_pairs_col]

        rna_chain_ids = safe_eval(row['RNAChainIDs']) if isinstance(row['RNAChainIDs'], str) else row['RNAChainIDs']
        exp_rna_seqs = safe_eval(row['RNASequence_exp']) if isinstance(row['RNASequence_exp'], str) else row['RNASequence_exp']
        pred_rna_seqs = safe_eval(row['RNASequence_pred']) if isinstance(row['RNASequence_pred'], str) else row['RNASequence_pred']
        af3_rna_iptm = safe_eval(row['af3_rna_ipTM']) if isinstance(row['af3_rna_ipTM'], str) else row['af3_rna_ipTM']
        protein_chain_ids = safe_eval(row['ProteinChainIDs']) if isinstance(row['ProteinChainIDs'], str) else row['ProteinChainIDs']

        # 1. Find index of max motif length
        motif_max_idx = motif_lengths.index(max(motif_lengths))

        # 2. Get the chain pair at that index
        chain_pair = chain_pairs[motif_max_idx]
        # 3. Find the RNA chain ID (not in protein chains)
        rna_chain_id = [cid for cid in chain_pair if cid not in protein_chain_ids][0]
        # 4. Find the index of this RNA chain in RNAChainIDs
        rna_chain_idx = rna_chain_ids.index(rna_chain_id)
        # 5. Get the corresponding sequence
        exp_seq = exp_rna_seqs[rna_chain_idx]
        # 6. Find the index of this sequence in pred_rna_seqs
        pred_seq_idx = pred_rna_seqs.index(exp_seq)
        # 7. Find the index of max ipTM
        iptm_max_idx = af3_rna_iptm.index(max(af3_rna_iptm))
        # 8. If pred_seq_idx matches iptm_max_idx, return 1
        if pred_seq_idx == iptm_max_idx:
            # print("Raw: 100.0 Normalized: 1.0")
            return 1.0
        # 9. Otherwise, return (100 - percent difference)/100, relative to max ipTM
        iptm_at_motif = af3_rna_iptm[pred_seq_idx]
        iptm_max = af3_rna_iptm[iptm_max_idx]
        if iptm_max == 0 or iptm_at_motif is None or iptm_max is None:
            return 0.0
        percent_diff = abs(iptm_max - iptm_at_motif) / abs(iptm_max) * 100
        normalized = (100 - percent_diff) / 100
        # Clamp to [0, 1]
        normalized = max(0.0, min(1.0, normalized))
        return normalized
    except Exception as e:
        print(f"Error in motif_ipTM_match_score: {e}")
        return np.nan


class ResultPlotter:
    def __init__(self, database, config, msa_option=None, single_chain_only=False):
        """Initialize database connection and configuration"""
        self.database = database
        self.config = config
        self.msa_option = msa_option

        # Set up base results directory - remove MSA folder creation at root level
        self.results_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results', 'plots')

        # Create results directory if it doesn't exist
        os.makedirs(self.results_dir, exist_ok=True)

        # Get data from database
        self.pred_df = self.get_predicted_data()
        self.exp_df = self.get_experimental_data()
        
        # Apply filters if specified
        if single_chain_only:
            # Filter based on ChainIDpairList_proteinRNA instead of ProteinChainIDs
            self.pred_df = self.pred_df[self.pred_df['ChainIDpairList_proteinRNA'].apply(self.is_single_chain)]
            # Filter experimental data to match
            self.exp_df = self.exp_df[self.exp_df['PDBId'].isin(self.pred_df['exp_db_id'])]
            
        if self.msa_option:
            self.pred_df = self.database.filter_by_msa(self.pred_df, self.msa_option)
            # Filter experimental data to match
            self.exp_df = self.exp_df[self.exp_df['PDBId'].isin(self.pred_df['exp_db_id'])]
            
        self.merged_df = self.merge_data()

    def get_predicted_data(self):
        """Get predicted data from database"""

        def get_msa_filter(msa_option):
            if msa_option == '+MSA':
                return "AND (af3_rna_MSA = '[1]' OR af3_rna_MSA LIKE '[1,1%')"
            elif msa_option == '-MSA':
                return "AND (af3_rna_MSA = '[0]' OR af3_rna_MSA LIKE '%0%')"
            return ""

        msa_filter = get_msa_filter(self.msa_option)

        if self.config.pred_table == 'pred_rna_rna':
            query = f"""
                    SELECT exp_db_id, FileName, 
                           af3_rna_pTM, af3_rna_ipTM,af3_rna_pLDDT_avg,
                           af3_global_pae_avg, af3_chain_pair_pae_min,
                           af3_fraction_disordered, af3_rna_MSA,
                           Free_energy, Binding_affinity_kd,
                           RNA_GlycosidicBond, RNA_SugarPucker, RNA_GammaAngle,
                           RNALength, RNA_GC_Content, RNA_SequenceComplexity, RNA_Shannon_Entropy_K2, RNA_MFE_Value,
                           RNA_RMSD, RNA_TM, RNA_LDDT,
                           RNA_f1_score, RNA_precision, RNA_recall, RNA_mcc, RNA_wl
                    FROM {self.config.pred_table}
                    WHERE 1=1 {msa_filter}
                    """
        else:
            query = f"""
            SELECT exp_db_id, FileName, 
                   af3_protein_pTM, af3_rna_pTM, 
                   af3_protein_ipTM, af3_rna_ipTM,
                   af3_protein_pLDDT_avg, af3_rna_pLDDT_avg,
                   af3_global_pae_avg, af3_chain_pair_pae_min,
                   af3_fraction_disordered, af3_rna_MSA,
                   Free_energy, Binding_affinity_kd,
                   RNA_GlycosidicBond, RNA_SugarPucker, RNA_GammaAngle,
                   ChainIDpairList_proteinRNA,
                   ProteinRNAInterfaceRatio as pred_interface_ratio, 
                   RNAmotif_score,
                   RNALength, ProteinLength,
                   RNA_GC_Content, RNA_SequenceComplexity, RNA_Shannon_Entropy_K2, RNA_MFE_Value,
                   RNA_RMSD, RNA_TM, Protein_TM, Complex_TM, Complex_RMSD, RNA_LDDT,
                   RNA_f1_score, RNA_precision, RNA_recall, RNA_mcc, RNA_wl,
                   NumberRNAs, RNA_msa_size, RNAChainIDs as RNAChainIDs_pred,
                   RNASequence as RNASequence_pred, AAmatch_score
            FROM {self.config.pred_table}
            WHERE 1=1 {msa_filter}
            """
        return pd.read_sql_query(query, self.database.connection)

    def get_experimental_data(self):
        """Get experimental data from database"""
        query = f"""
        SELECT PDBId,
               Free_energy, Binding_affinity_kd,
               RNA_GlycosidicBond, RNA_SugarPucker, RNA_GammaAngle,
               ProteinRNAInterfaceRatio as exp_interface_ratio,
               RNAMotifLength as RNAMotifLength_exp, ChainIDpairList_proteinRNA as ChainIDpairList_proteinRNA_exp, RNAChainIDs, 
               RNASequence as RNASequence_exp, ProteinChainIDs
        FROM {self.config.ref_table}
        """
        return pd.read_sql_query(query, self.database.connection)

    def is_single_chain(self, chain_pairs):
        """Check if the structure has a single chain pair"""
        try:
            pairs = literal_eval(chain_pairs)
            return len(pairs) == 1
        except:
            return False

    def filter_chains(self, df, single_chain_only):
        """Filter data based on chain pairs"""
        if single_chain_only:
            return df[df['ChainIDpairList_proteinRNA'].apply(self.is_single_chain)]
        return df

    def merge_data(self):
        """Merge predicted and experimental data based on exp_db_id/PDBId"""
        return pd.merge(
            self.pred_df,
            self.exp_df,
            left_on='exp_db_id',
            right_on='PDBId',
            suffixes=('_pred', '_exp')
        )

    def get_max_value_from_list(self, value):
        """Extract maximum value from a string list or return float value"""
        try:
            if isinstance(value, str):
                # Handle list format - strip brackets and split
                value = value.strip('[]').replace(' ', '')
                if ',' in value:
                    values = [float(x) for x in value.split(',')]
                    return max(values)
                return float(value)
            elif isinstance(value, (int, float)):
                return float(value)
            return None
        except (ValueError, TypeError):
            return None

    def plot_af3_metrics(self, single_chain_only=False):
        """Create violin plots for AF3 metrics (excluding PAE) and scatter plots vs RMSD"""
        import matplotlib.pyplot as plt
        # Get filtered data first
        df = self.filter_data(self.pred_df, single_chain_only)

        # Create violin plot using PlotCreator with af3_metrics type
        plot_creator = PlotCreator('af3_metrics', self.msa_option, single_chain_only)

        # Add new plot for Chain-Pair PAE vs RNAmotif_score
        if not df['af3_chain_pair_pae_min'].isna().all() and not df['RNAmotif_score'].isna().all():
            # Handle empty strings and None values properly
            def process_pae(x):
                if pd.isna(x) or x == '':
                    return None
                try:
                    return float(x) if isinstance(x, str) else float(x)
                except (ValueError, TypeError):
                    return None

            pae_values = df['af3_chain_pair_pae_min'].apply(process_pae).dropna()
            valid_indices = pae_values.index.intersection(df.dropna(subset=['RNAmotif_score']).index)

            if len(valid_indices) > 0:
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    yAxis_score=df.loc[valid_indices, 'RNAmotif_score'],
                    yAxis_label='RNA Motif Score',
                    xAxis_label='Chain-Pair PAE',
                    name='chain_pair_pae_vs_motif_score',
                    xAxis_score=pae_values[valid_indices]
                )

        # Merge with experimental data for interface ratio comparison
        merged_df = pd.merge(
            df,
            self.exp_df[['PDBId', 'exp_interface_ratio']],
            left_on='exp_db_id',
            right_on='PDBId',
            how='left'
        )

        # Convert interface ratios to float
        merged_df['pred_interface_ratio'] = pd.to_numeric(merged_df['pred_interface_ratio'], errors='coerce')
        merged_df['exp_interface_ratio'] = pd.to_numeric(merged_df['exp_interface_ratio'], errors='coerce')

        # Calculate interface ratio difference
        merged_df['r_is_diff'] = ((merged_df['pred_interface_ratio'] - merged_df['exp_interface_ratio']) /
                                  merged_df['exp_interface_ratio'] * 100).round(2)

        # Add new plots for interface ratio difference
        # Plot Chain-Pair PAE vs r_is_diff
        if not merged_df['af3_chain_pair_pae_min'].isna().all() and not merged_df['r_is_diff'].isna().all():
            pae_values = merged_df['af3_chain_pair_pae_min'].apply(process_pae).dropna()
            valid_indices = pae_values.index.intersection(merged_df.dropna(subset=['r_is_diff']).index)

            if len(valid_indices) > 0:
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=merged_df.loc[valid_indices, 'r_is_diff'],
                    xAxis_label='Interface Ratio Difference (%)',
                    yAxis_label='Chain-Pair PAE',
                    name='chain_pair_pae_vs_r_is_diff',
                    yAxis_score=pae_values[valid_indices]
                )

        # Plot RNA ipTM vs r_is_diff
        if not merged_df['af3_rna_ipTM'].isna().all() and not merged_df['r_is_diff'].isna().all():
            def process_tm_values(tm_str, rna_length_str):
                if pd.isna(tm_str) or pd.isna(rna_length_str):
                    return None
                try:
                    # Convert TM string to list
                    if isinstance(tm_str, str):
                        tm_values = eval(tm_str)
                    else:
                        tm_values = [tm_str]
                    
                    # Convert RNA length string to list
                    if isinstance(rna_length_str, str):
                        if '[' in rna_length_str:
                            rna_lengths = eval(rna_length_str)
                        else:
                            rna_lengths = [int(rna_length_str)]
                    else:
                        rna_lengths = [int(rna_length_str)]
                    
                    # Find index of first occurrence of longest RNA
                    if len(rna_lengths) > 1:
                        max_length = max(rna_lengths)
                        max_length_idx = rna_lengths.index(max_length)
                        if max_length_idx < len(tm_values):
                            return tm_values[max_length_idx]
                        return None
                    return tm_values[0]
                except (ValueError, TypeError, SyntaxError, IndexError):
                    return None

            # Process both pTM and ipTM values
            merged_df['processed_pTM'] = merged_df.apply(
                lambda row: process_tm_values(row['af3_rna_pTM'], row['RNALength']), axis=1
            )
            merged_df['processed_ipTM'] = merged_df.apply(
                lambda row: process_tm_values(row['af3_rna_ipTM'], row['RNALength']), axis=1
            )

            # Plot ipTM vs r_is_diff
            iptm_values = merged_df['processed_ipTM'].dropna()
            valid_indices = iptm_values.index.intersection(merged_df.dropna(subset=['r_is_diff']).index)

            if len(valid_indices) > 0:
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=merged_df.loc[valid_indices, 'r_is_diff'],
                    xAxis_label='Interface Ratio Difference (%)',
                    yAxis_label='RNA ipTM',
                    name='rna_iptm_vs_r_is_diff',
                    yAxis_score=iptm_values[valid_indices]
                )

        # Get RMSD data from database with necessary columns for filtering
        rmsd_query = f"""
        SELECT exp_db_id, RNA_RMSD, af3_chain_pair_pae_min, 
               ChainIDpairList_proteinRNA, af3_rna_MSA
        FROM {self.config.pred_table}
        WHERE RNA_RMSD IS NOT NULL
        """
        rmsd_df = pd.read_sql_query(rmsd_query, self.database.connection)

        # Apply the same filtering to RMSD data
        rmsd_df = self.filter_data(rmsd_df, single_chain_only)

        if not rmsd_df.empty:
            # Get first RMSD value for each entry
            def extract_first_rmsd(x):
                if pd.isna(x) or x == '':
                    return None
                try:
                    if isinstance(x, str):
                        if '[' in x:
                            values = eval(x)
                            if isinstance(values, list):
                                return float(values[0])
                        return float(x)
                    return float(x)
                except (ValueError, TypeError, SyntaxError, IndexError):
                    return None

            # Get PAE value for each entry
            def extract_pae(x):
                if pd.isna(x) or x == '':
                    return None
                try:
                    if isinstance(x, str):
                        values = eval(x)
                        return float(values)
                    return float(x)
                except (ValueError, TypeError, SyntaxError):
                    return None

            rmsd_df['RNA_RMSD_first'] = rmsd_df['RNA_RMSD'].apply(extract_first_rmsd)
            rmsd_df['maxMotif_chain_pair_pae'] = rmsd_df['af3_chain_pair_pae_min'].apply(extract_pae)

            # Merge with filtered main dataframe
            df = df.merge(rmsd_df[['exp_db_id', 'RNA_RMSD_first', 'maxMotif_chain_pair_pae']],
                          left_on='exp_db_id',
                          right_on='exp_db_id',
                          how='left')

            # Create scatter plots for AF3 metrics vs RMSD
            metrics_to_plot = {
                'af3_rna_pTM': 'RNA pTM',
                'af3_rna_ipTM': 'RNA ipTM',
                'af3_rna_pLDDT_avg': 'RNA pLDDT',
                'af3_global_pae_avg': 'Global PAE',
                'maxMotif_chain_pair_pae': 'Chain-Pair PAE'
            }

            for column, label in metrics_to_plot.items():
                # Get values, handling lists for pTM and ipTM
                if column in ['af3_rna_pTM', 'af3_rna_ipTM']:
                    values = df.apply(
                        lambda row: process_tm_values(row[column], row['RNALength']), 
                        axis=1
                    ).dropna()
                else:
                    values = df[column].dropna()

                # Create scatter plot only if we have both metrics and RMSD values
                valid_indices = values.index.intersection(df.dropna(subset=['RNA_RMSD_first']).index)
                if len(valid_indices) > 0:
                    # Convert values to float before comparison
                    rmsd_values = pd.to_numeric(df.loc[valid_indices, 'RNA_RMSD_first'], errors='coerce')
                    metric_values = pd.to_numeric(values[valid_indices], errors='coerce')
                    
                    # Add check for zero/negative values before log10
                    valid_mask = (rmsd_values > 0) & (metric_values > 0)

                    if valid_mask.any():
                        plot_creator.get_scatterplot(
                            table_source='af3_metrics',
                            xAxis_score=rmsd_values[valid_mask],
                            xAxis_label='RNA RMSD [Å]',
                            yAxis_label=label,
                            name=f'{column}_vs_RMSD',
                            yAxis_score=metric_values[valid_mask]
                        )

        # Helper to process metrics by length
        def process_metric_by_length(metric_str, length_str):
            if pd.isna(metric_str) or pd.isna(length_str):
                return None
            try:
                # Convert metric string to list
                if isinstance(metric_str, str):
                    metric_values = eval(metric_str)
                else:
                    metric_values = [metric_str]
                # Convert length string to list
                if isinstance(length_str, str):
                    if '[' in length_str:
                        lengths = eval(length_str)
                    else:
                        lengths = [int(length_str)]
                else:
                    lengths = [int(length_str)]
                # Find index of first occurrence of max length
                if len(lengths) > 1:
                    max_length = max(lengths)
                    max_length_idx = lengths.index(max_length)
                    if max_length_idx < len(metric_values):
                        return metric_values[max_length_idx]
                    return None
                return metric_values[0]
            except Exception:
                return None

        def average_metric(metric_str):
            if pd.isna(metric_str):
                return np.nan
            try:
                if isinstance(metric_str, str):
                    if '[' in metric_str:
                        values = eval(metric_str)
                        return float(np.mean(values))
                    return float(metric_str)
                return float(metric_str)
            except Exception:
                return np.nan

        def max_metric(metric_str):
            if pd.isna(metric_str):
                return np.nan
            try:
                if isinstance(metric_str, str):
                    if '[' in metric_str:
                        values = eval(metric_str)
                        return float(np.max(values))
                    return float(metric_str)
                return float(metric_str)
            except Exception:
                return np.nan

        # Use OrderedDict to guarantee order
        metric_keys = [
            'Protein pTM', 'RNA pTM', 'Protein ipTM', 'RNA ipTM',
            'Protein pLDDT', 'RNA pLDDT', 'Fraction Disordered'
        ]
        data_dict = OrderedDict()
        data_dict['Protein pTM'] = df['af3_protein_pTM'].apply(average_metric) if 'af3_protein_pTM' in df.columns else pd.Series(dtype=float)
        data_dict['RNA pTM'] = df['af3_rna_pTM'].apply(average_metric) if 'af3_rna_pTM' in df.columns else pd.Series(dtype=float)
        data_dict['Protein ipTM'] = df['af3_protein_ipTM'].apply(max_metric) if 'af3_protein_ipTM' in df.columns else pd.Series(dtype=float)
        data_dict['RNA ipTM'] = df['af3_rna_ipTM'].apply(max_metric) if 'af3_rna_ipTM' in df.columns else pd.Series(dtype=float)
        data_dict['Protein pLDDT'] = pd.to_numeric(df['af3_protein_pLDDT_avg'], errors='coerce') if 'af3_protein_pLDDT_avg' in df.columns else pd.Series(dtype=float)
        data_dict['RNA pLDDT'] = pd.to_numeric(df['af3_rna_pLDDT_avg'], errors='coerce') if 'af3_rna_pLDDT_avg' in df.columns else pd.Series(dtype=float)
        data_dict['Fraction Disordered'] = pd.to_numeric(df['af3_fraction_disordered'], errors='coerce') if 'af3_fraction_disordered' in df.columns else pd.Series(dtype=float)

        # Ensure all metrics are present, even if empty
        for key in metric_keys:
            if key not in data_dict:
                data_dict[key] = pd.Series(dtype=float)

        data_values = [pd.to_numeric(vals, errors='coerce').dropna().reset_index(drop=True) for vals in data_dict.values()]
        labels = list(data_dict.keys())
        plot_creator.get_violin_plot(
            table_source='af3_metrics',
            data_values=data_values,
            labels=labels,
            name='af3_metrics',
            rotation=90
        )

        # Add scatter plot of RNA pTM vs RNA_TM if available
        rna_ptm = df.apply(lambda row: process_metric_by_length(row['af3_rna_pTM'], row['RNALength']) if ('RNALength' in df.columns and not pd.isna(row['af3_rna_pTM']) and not pd.isna(row['RNALength'])) else None, axis=1)
        rna_tm = pd.to_numeric(df['RNA_TM'], errors='coerce') if 'RNA_TM' in df.columns else pd.Series(dtype=float)
        valid = (~rna_ptm.isna()) & (~rna_tm.isna())
        if valid.any():
            plot_creator.get_scatterplot(
                table_source='af3_metrics',
                xAxis_score=rna_ptm[valid],
                xAxis_label='RNA pTM',
                yAxis_label='RNA TM',
                name='rna_pTM_vs_rna_TM',
                yAxis_score=rna_tm[valid],
                xAxis_limit=(0, 1),
                yAxis_limit=(0, 1)
            )

        # Add scatter plot of Protein pTM vs Protein_TM if available
        if 'af3_protein_pTM' in df.columns and 'Protein_TM' in df.columns and 'ProteinLength' in df.columns:
            print("@£werd")
            protein_ptm = df.apply(lambda row: process_metric_by_length(row['af3_protein_pTM'], row['ProteinLength']) if (not pd.isna(row['af3_protein_pTM']) and not pd.isna(row['ProteinLength'])) else None, axis=1)
            protein_tm = pd.to_numeric(df['Protein_TM'], errors='coerce')
            print("dsf", protein_ptm)
            print("protein_tm", protein_tm)
            valid = (~protein_ptm.isna()) & (~protein_tm.isna())
            if valid.any():
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=protein_ptm[valid],
                    xAxis_label='Protein pTM',
                    yAxis_label='Protein TM',
                    name='protein_pTM_vs_protein_TM',
                    yAxis_score=protein_tm[valid],
                    xAxis_limit=(0, 1),
                    yAxis_limit=(0, 1)
                )

        # Add scatter plot of max af3_rna_ipTM vs RNAmotif_score if available
        if 'af3_rna_ipTM' in df.columns and 'RNAmotif_score' in df.columns:
            rna_iptm_max = df['af3_rna_ipTM'].apply(max_metric)
            motif_score = pd.to_numeric(df['RNAmotif_score'], errors='coerce')
            valid = (~rna_iptm_max.isna()) & (~motif_score.isna())
            if valid.any():
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    yAxis_score=rna_iptm_max[valid],
                    yAxis_label='max RNA ipTM',
                    xAxis_label='RNA motif score',
                    name='rna_ipTMmax_vs_rnamotif_score',
                    xAxis_score=motif_score[valid]
                )

        # Add scatter plot of rna_iptm_max vs Complex_TM if available
        if 'af3_rna_ipTM' in df.columns and 'Complex_TM' in df.columns:
            rna_iptm_max = df['af3_rna_ipTM'].apply(max_metric)
            complex_tm = pd.to_numeric(df['Complex_TM'], errors='coerce')
            valid = (~rna_iptm_max.isna()) & (~complex_tm.isna())
            if valid.any():
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=rna_iptm_max[valid],
                    xAxis_label='max RNA ipTM',
                    yAxis_label='Complex TM',
                    name='rna_ipTMmax_vs_ComplexTM',
                    yAxis_score=complex_tm[valid],
                    xAxis_limit=(0, 1),
                    yAxis_limit=(0, 1)
                )

        # Add scatter plot of chain_pair_pae vs Complex_TM if available
        if 'af3_chain_pair_pae_min' in df.columns and 'Complex_TM' in df.columns:
            chain_pair_pae = pd.to_numeric(df['af3_chain_pair_pae_min'], errors='coerce')
            complex_tm = pd.to_numeric(df['Complex_TM'], errors='coerce')
            valid = (~chain_pair_pae.isna()) & (~complex_tm.isna())
            if valid.any():
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=chain_pair_pae[valid],
                    xAxis_label='Chain Pair PAE',
                    yAxis_label='Complex TM',
                    name='chain_pair_pae_vs_ComplexTM',
                    yAxis_score=complex_tm[valid],
                    yAxis_limit=(0, 1)
                )

        # Add scatter plot of rna_iptm_max vs Complex_RMSD if available
        if 'af3_rna_ipTM' in df.columns and 'Complex_RMSD' in df.columns:
            rna_iptm_max = df['af3_rna_ipTM'].apply(max_metric)
            def first_value(val):
                if pd.isna(val):
                    return np.nan
                try:
                    if isinstance(val, str) and '[' in val:
                        values = eval(val)
                        return float(values[0])
                    return float(val)
                except Exception:
                    return np.nan
            complex_rmsd = df['Complex_RMSD'].apply(first_value)
            valid = (~rna_iptm_max.isna()) & (~complex_rmsd.isna())
            if valid.any():
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=rna_iptm_max[valid],
                    xAxis_label='max RNA ipTM',
                    yAxis_label='Complex RMSD',
                    name='rna_ipTMmax_vs_ComplexRMSD',
                    yAxis_score=complex_rmsd[valid]
                )

        # Add scatter plot of max RNA ipTM vs binding affinity if available
        if 'af3_rna_ipTM' in df.columns and 'Binding_affinity_kd' in df.columns:
            rna_iptm_max = df['af3_rna_ipTM'].apply(max_metric)
            def process_binding_affinity(val):
                if pd.isna(val):
                    return np.nan
                try:
                    if isinstance(val, str):
                        return float(val)
                    return float(val)
                except Exception:
                    return np.nan
            binding_affinity = df['Binding_affinity_kd'].apply(process_binding_affinity)
            valid = (~rna_iptm_max.isna()) & (~binding_affinity.isna())
            if valid.any():
                # Convert to log scale if values are very small
                if binding_affinity[valid].min() < 1e-10:
                    # Handle zero values before log conversion
                    non_zero_mask = binding_affinity > 0
                    binding_affinity_log = pd.Series(index=binding_affinity.index)
                    binding_affinity_log[non_zero_mask] = np.log10(binding_affinity[non_zero_mask])
                    binding_affinity_log[~non_zero_mask] = np.nan
                    valid = valid & non_zero_mask  # Update valid mask to exclude zeros
                    ylabel = 'Log₁₀[Binding Affinity (M)]'
                else:
                    binding_affinity_log = binding_affinity
                    ylabel = 'Binding Affinity (M)'

                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=rna_iptm_max[valid],
                    xAxis_label='max RNA ipTM',
                    yAxis_label=ylabel,
                    name='rna_ipTMmax_vs_binding_affinity',
                    yAxis_score=binding_affinity_log[valid]
                )

        # Add scatter plot of Chain Pair PAE vs binding affinity if available
        if 'af3_chain_pair_pae_min' in df.columns and 'Binding_affinity_kd' in df.columns:
            chain_pair_pae = pd.to_numeric(df['af3_chain_pair_pae_min'], errors='coerce')
            binding_affinity = df['Binding_affinity_kd'].apply(process_binding_affinity)
            valid = (~chain_pair_pae.isna()) & (~binding_affinity.isna())
            if valid.any():
                # Convert to log scale if values are very small
                if binding_affinity[valid].min() < 1e-10:
                    # Handle zero values before log conversion
                    non_zero_mask = binding_affinity > 0
                    binding_affinity_log = pd.Series(index=binding_affinity.index)
                    binding_affinity_log[non_zero_mask] = np.log10(binding_affinity[non_zero_mask])
                    binding_affinity_log[~non_zero_mask] = np.nan
                    valid = valid & non_zero_mask  # Update valid mask to exclude zeros
                    ylabel = 'Log₁₀[Binding Affinity (M)]'
                else:
                    binding_affinity_log = binding_affinity
                    ylabel = 'Binding Affinity (M)'

                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=chain_pair_pae[valid],
                    xAxis_label='Chain Pair PAE',
                    yAxis_label=ylabel,
                    name='chain_pair_pae_vs_binding_affinity',
                    yAxis_score=binding_affinity_log[valid]
                )

        # Add scatter plot of RNA pLDDT vs RNA LDDT if available
        if 'af3_rna_pLDDT_avg' in df.columns and 'RNA_LDDT' in df.columns:
            rna_plddt = pd.to_numeric(df['af3_rna_pLDDT_avg'], errors='coerce')
            rna_lddt = pd.to_numeric(df['RNA_LDDT'], errors='coerce')
            valid = (~rna_plddt.isna()) & (~rna_lddt.isna())
            if valid.any():
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=rna_plddt[valid],
                    xAxis_label='RNA pLDDT',
                    yAxis_label='RNA LDDT',
                    name='rna_plddt_vs_rna_lddt',
                    yAxis_score=rna_lddt[valid],
                    xAxis_limit=(0, 1),
                    yAxis_limit=(0, 1)
                )

        # Add MSA color-coded scatter plot of RNA pLDDT vs RNA LDDT if available
        if all(col in df.columns for col in ['af3_rna_pLDDT_avg', 'RNA_LDDT', 'RNA_msa_size', 'RNAChainIDs_pred', 'RNALength']):
            def classify_msa(row):
                try:
                    msa_dict = row['RNA_msa_size']
                    if isinstance(msa_dict, str):
                        msa_dict = safe_eval(msa_dict)
                    if not isinstance(msa_dict, dict):
                        return np.nan

                    chain_ids = row['RNAChainIDs_pred']
                    if isinstance(chain_ids, str):
                        chain_ids = safe_eval(chain_ids)
                    # Ensure it's a list
                    if not isinstance(chain_ids, (list, tuple)):
                        chain_ids = [chain_ids]

                    rna_lengths = row['RNALength']
                    if isinstance(rna_lengths, str):
                        rna_lengths = safe_eval(rna_lengths)
                    # Ensure it's a list
                    if not isinstance(rna_lengths, (list, tuple)):
                        rna_lengths = [rna_lengths]

                    # If only one chain, just use its MSA value
                    if len(msa_dict) == 1:
                        msa_val = list(msa_dict.values())[0]
                    else:
                        # Find index of longest RNA
                        max_len = max(rna_lengths)
                        max_idx = rna_lengths.index(max_len)
                        chain_id = chain_ids[max_idx]
                        msa_val = msa_dict.get(str(chain_id), 1)

                    return '+MSA' if msa_val > 1 else '-MSA'
                except Exception as e:
                    print(f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                    return np.nan

            msa_status = df.apply(classify_msa, axis=1)
            rna_plddt = pd.to_numeric(df['af3_rna_pLDDT_avg'], errors='coerce')
            rna_lddt = pd.to_numeric(df['RNA_LDDT'], errors='coerce')
            valid = (~rna_plddt.isna()) & (~rna_lddt.isna()) & (~msa_status.isna())
            
            if valid.any():
                color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                plot_creator.get_msa_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=rna_plddt[valid],
                    xAxis_label='RNA pLDDT',
                    name='rna_plddt_vs_rna_lddt_all_chains',
                    yAxis_label='RNA LDDT',
                    yAxis_score=rna_lddt[valid],
                    color_by=msa_status[valid],
                    color_map=color_map
                )

        # Add MSA color-coded scatter plot for Global PAE vs RMSD
        if all(col in df.columns for col in
               ['af3_global_pae_avg', 'Complex_RMSD', 'RNA_msa_size', 'RNAChainIDs_pred', 'RNALength']):
            def classify_msa(row):
                try:
                    msa_dict = row['RNA_msa_size']
                    if isinstance(msa_dict, str):
                        msa_dict = safe_eval(msa_dict)
                    if not isinstance(msa_dict, dict):
                        return np.nan

                    chain_ids = row['RNAChainIDs_pred']
                    if isinstance(chain_ids, str):
                        chain_ids = safe_eval(chain_ids)
                    # Ensure it's a list
                    if not isinstance(chain_ids, (list, tuple)):
                        chain_ids = [chain_ids]

                    rna_lengths = row['RNALength']
                    if isinstance(rna_lengths, str):
                        rna_lengths = safe_eval(rna_lengths)
                    # Ensure it's a list
                    if not isinstance(rna_lengths, (list, tuple)):
                        rna_lengths = [rna_lengths]

                    # If only one chain, just use its MSA value
                    if len(msa_dict) == 1:
                        msa_val = list(msa_dict.values())[0]
                    else:
                        # Find index of longest RNA
                        max_len = max(rna_lengths)
                        max_idx = rna_lengths.index(max_len)
                        chain_id = chain_ids[max_idx]
                        msa_val = msa_dict.get(str(chain_id), 1)

                    return '+MSA' if msa_val > 1 else '-MSA'
                except Exception as e:
                    print(
                        f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                    return np.nan

            def extract_first_rmsd(x):
                if pd.isna(x) or x == '':
                    return None
                try:
                    if isinstance(x, str):
                        if '[' in x:
                            values = eval(x)
                            if isinstance(values, list):
                                return float(values[0])
                        return float(x)
                    return float(x)
                except (ValueError, TypeError, SyntaxError, IndexError):
                    return None

            msa_status = df.apply(classify_msa, axis=1)
            global_pae = pd.to_numeric(df['af3_global_pae_avg'], errors='coerce')
            rmsd_values = df['Complex_RMSD'].apply(extract_first_rmsd)
            valid = (~global_pae.isna()) & (~rmsd_values.isna()) & (~msa_status.isna())

        if valid.any():
            color_map = {'+MSA': 'blue', '-MSA': 'orange'}
            plot_creator.get_msa_scatterplot(
                table_source='af3_metrics',
                xAxis_score=rmsd_values[valid],
                xAxis_label='Complex RMSD [Å]',
                name='af3_global_pae_avg_vs_RMSD_all_chains',
                yAxis_label='Global PAE',
                yAxis_score=global_pae[valid],
                color_by=msa_status[valid],
                color_map=color_map
            )

            # Add MSA color-coded scatter plot for Global PAE vs RNAmotif_score
            if all(col in df.columns for col in
                   ['af3_global_pae_avg', 'RNAmotif_score', 'RNA_msa_size', 'RNAChainIDs_pred', 'RNALength']):
                def classify_msa(row):
                    try:
                        msa_dict = row['RNA_msa_size']
                        if isinstance(msa_dict, str):
                            msa_dict = safe_eval(msa_dict)
                        if not isinstance(msa_dict, dict):
                            return np.nan

                        chain_ids = row['RNAChainIDs_pred']
                        if isinstance(chain_ids, str):
                            chain_ids = safe_eval(chain_ids)
                        # Ensure it's a list
                        if not isinstance(chain_ids, (list, tuple)):
                            chain_ids = [chain_ids]

                        rna_lengths = row['RNALength']
                        if isinstance(rna_lengths, str):
                            rna_lengths = safe_eval(rna_lengths)
                        # Ensure it's a list
                        if not isinstance(rna_lengths, (list, tuple)):
                            rna_lengths = [rna_lengths]

                        # If only one chain, just use its MSA value
                        if len(msa_dict) == 1:
                            msa_val = list(msa_dict.values())[0]
                        else:
                            # Find index of longest RNA
                            max_len = max(rna_lengths)
                            max_idx = rna_lengths.index(max_len)
                            chain_id = chain_ids[max_idx]
                            msa_val = msa_dict.get(str(chain_id), 1)

                        return '+MSA' if msa_val > 1 else '-MSA'
                    except Exception as e:
                        print(
                            f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                        return np.nan

                msa_status = df.apply(classify_msa, axis=1)
                global_pae = pd.to_numeric(df['af3_global_pae_avg'], errors='coerce')
                motif_score = pd.to_numeric(df['RNAmotif_score'], errors='coerce')
                valid = (~global_pae.isna()) & (~motif_score.isna()) & (~msa_status.isna())

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=global_pae[valid],
                        xAxis_label='Global PAE',
                        name='af3_global_pae_avg_vs_RNAmotif_score_all_chains',
                        yAxis_label='RNA Motif Score',
                        yAxis_score=motif_score[valid],
                        color_by=msa_status[valid],
                        color_map=color_map
                    )

            if all(col in df.columns for col in
                   ['RNA_TM', 'RNA_GC_Content', 'RNA_msa_size', 'RNAChainIDs_pred', 'RNALength']):
                def classify_msa(row):
                    try:
                        msa_dict = row['RNA_msa_size']
                        if isinstance(msa_dict, str):
                            msa_dict = safe_eval(msa_dict)
                        if not isinstance(msa_dict, dict):
                            return np.nan

                        chain_ids = row['RNAChainIDs_pred']
                        if isinstance(chain_ids, str):
                            chain_ids = safe_eval(chain_ids)
                        # Ensure it's a list
                        if not isinstance(chain_ids, (list, tuple)):
                            chain_ids = [chain_ids]

                        rna_lengths = row['RNALength']
                        if isinstance(rna_lengths, str):
                            rna_lengths = safe_eval(rna_lengths)
                        # Ensure it's a list
                        if not isinstance(rna_lengths, (list, tuple)):
                            rna_lengths = [rna_lengths]

                        # If only one chain, just use its MSA value
                        if len(msa_dict) == 1:
                            msa_val = list(msa_dict.values())[0]
                        else:
                            # Find index of longest RNA
                            max_len = max(rna_lengths)
                            max_idx = rna_lengths.index(max_len)
                            chain_id = chain_ids[max_idx]
                            msa_val = msa_dict.get(str(chain_id), 1)

                        return '+MSA' if msa_val > 1 else '-MSA'
                    except Exception as e:
                        print(
                            f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                        return np.nan

                msa_status = df.apply(classify_msa, axis=1)
                gc_content = pd.to_numeric(df['RNA_GC_Content'], errors='coerce')
                tm_score = pd.to_numeric(df['RNA_TM'], errors='coerce')
                valid = (~gc_content.isna()) & (~tm_score.isna()) & (~msa_status.isna())

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=gc_content[valid],
                        xAxis_label='RNA GC-content',
                        name='GC_tm',
                        yAxis_label='RNA TM',
                        yAxis_score=tm_score[valid],
                        color_by=msa_status[valid],
                        color_map=color_map
                    )

                lddt_score = pd.to_numeric(df['RNA_LDDT'], errors='coerce')
                valid2 = (~gc_content.isna()) & (~lddt_score.isna()) & (~msa_status.isna())

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=gc_content[valid2],
                        xAxis_label='RNA GC-content',
                        name='GC_lddt',
                        yAxis_label='RNA LDDT',
                        yAxis_score=lddt_score[valid2],
                        color_by=msa_status[valid2],
                        color_map=color_map
                    )

                seq_complex = pd.to_numeric(df['RNA_SequenceComplexity'], errors='coerce')
                valid3 = (~seq_complex.isna()) & (~tm_score.isna()) & (~msa_status.isna())
                valid4 = (~seq_complex.isna()) & (~lddt_score.isna()) & (~msa_status.isna())

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=seq_complex[valid3],
                        xAxis_label='RNA Sequence Complexity',
                        name='seq_complex_tm',
                        yAxis_label='RNA TM',
                        yAxis_score=tm_score[valid3],
                        color_by=msa_status[valid3],
                        color_map=color_map
                    )

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=seq_complex[valid4],
                        xAxis_label='RNA Sequence Complexity',
                        name='seq_complex_lddt',
                        yAxis_label='RNA LDDT',
                        yAxis_score=lddt_score[valid4],
                        color_by=msa_status[valid4],
                        color_map=color_map
                    )

                diversity = pd.to_numeric(df['RNA_Shannon_Entropy_K2'], errors='coerce')
                valid3 = (~diversity.isna()) & (~tm_score.isna()) & (~msa_status.isna())
                valid4 = (~diversity.isna()) & (~lddt_score.isna()) & (~msa_status.isna())

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=diversity[valid3],
                        xAxis_label='RNA Shannon Entropy K2',
                        name='rna_shannon_tm',
                        yAxis_label='RNA TM',
                        yAxis_score=tm_score[valid3],
                        color_by=msa_status[valid3],
                        color_map=color_map
                    )

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=diversity[valid4],
                        xAxis_label='RNA Shannon Entropy K2',
                        name='rna_shannon_lddt',
                        yAxis_label='RNA LDDT',
                        yAxis_score=lddt_score[valid4],
                        color_by=msa_status[valid4],
                        color_map=color_map
                    )

                mfe = pd.to_numeric(df['RNA_MFE_Value'], errors='coerce')
                valid3 = (~diversity.isna()) & (~tm_score.isna()) & (~msa_status.isna())
                valid4 = (~diversity.isna()) & (~lddt_score.isna()) & (~msa_status.isna())

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=mfe[valid3],
                        xAxis_label='RNA MFE',
                        name='rna_mfe_tm',
                        yAxis_label='RNA TM',
                        yAxis_score=tm_score[valid3],
                        color_by=msa_status[valid3],
                        color_map=color_map
                    )

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=mfe[valid4],
                        xAxis_label='RNA MFE',
                        name='rna_mfe_lddt',
                        yAxis_label='RNA LDDT',
                        yAxis_score=lddt_score[valid4],
                        color_by=msa_status[valid4],
                        color_map=color_map
                    )

            # Add MSA color-coded scatter plot for Global PAE vs AAmatch_score
            if all(col in df.columns for col in
                   ['af3_global_pae_avg', 'AAmatch_score', 'RNA_msa_size', 'RNAChainIDs_pred', 'RNALength']):
                def classify_msa(row):
                    try:
                        msa_dict = row['RNA_msa_size']
                        if isinstance(msa_dict, str):
                            msa_dict = safe_eval(msa_dict)
                        if not isinstance(msa_dict, dict):
                            return np.nan

                        chain_ids = row['RNAChainIDs_pred']
                        if isinstance(chain_ids, str):
                            chain_ids = safe_eval(chain_ids)
                        # Ensure it's a list
                        if not isinstance(chain_ids, (list, tuple)):
                            chain_ids = [chain_ids]

                        rna_lengths = row['RNALength']
                        if isinstance(rna_lengths, str):
                            rna_lengths = safe_eval(rna_lengths)
                        # Ensure it's a list
                        if not isinstance(rna_lengths, (list, tuple)):
                            rna_lengths = [rna_lengths]

                        # If only one chain, just use its MSA value
                        if len(msa_dict) == 1:
                            msa_val = list(msa_dict.values())[0]
                        else:
                            # Find index of longest RNA
                            max_len = max(rna_lengths)
                            max_idx = rna_lengths.index(max_len)
                            chain_id = chain_ids[max_idx]
                            msa_val = msa_dict.get(str(chain_id), 1)

                        return '+MSA' if msa_val > 1 else '-MSA'
                    except Exception as e:
                        print(
                            f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                        return np.nan

                msa_status = df.apply(classify_msa, axis=1)
                global_pae = pd.to_numeric(df['af3_global_pae_avg'], errors='coerce')
                aa_match_score = pd.to_numeric(df['AAmatch_score'], errors='coerce')
                valid = (~global_pae.isna()) & (~aa_match_score.isna()) & (~msa_status.isna())

                if valid.any():
                    color_map = {'+MSA': 'blue', '-MSA': 'orange'}
                    plot_creator.get_msa_scatterplot(
                        table_source='af3_metrics',
                        xAxis_score=global_pae[valid],
                        xAxis_label='Global PAE',
                        name='af3_global_pae_avg_vs_AAmatch_score_all_chains',
                        yAxis_label='AA Match Score',
                        yAxis_score=aa_match_score[valid],
                        color_by=msa_status[valid],
                        color_map=color_map
                    )

        # Add scatter plot of chain_pair_pae vs Complex_RMSD if available
        if 'af3_chain_pair_pae_min' in df.columns and 'Complex_RMSD' in df.columns:
            chain_pair_pae = pd.to_numeric(df['af3_chain_pair_pae_min'], errors='coerce')
            def first_value(val):
                if pd.isna(val):
                    return np.nan
                try:
                    if isinstance(val, str) and '[' in val:
                        values = eval(val)
                        return float(values[0])
                    return float(val)
                except Exception:
                    return np.nan
            complex_rmsd = df['Complex_RMSD'].apply(first_value)
            valid = (~chain_pair_pae.isna()) & (~complex_rmsd.isna())
            if valid.any():
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=chain_pair_pae[valid],
                    xAxis_label='Chain Pair PAE',
                    yAxis_label='Complex RMSD',
                    name='chain_pair_pae_vs_ComplexRMSD',
                    yAxis_score=complex_rmsd[valid]
                )

        # Add scatter plot of Fraction Disordered vs RNAmotif_score if available
        if 'af3_fraction_disordered' in df.columns and 'RNAmotif_score' in df.columns:
            fraction_disordered = pd.to_numeric(df['af3_fraction_disordered'], errors='coerce')
            motif_score = pd.to_numeric(df['RNAmotif_score'], errors='coerce')
            valid = (~fraction_disordered.isna()) & (~motif_score.isna())
            if valid.any():
                plot_creator.get_scatterplot(
                    table_source='af3_metrics',
                    xAxis_score=fraction_disordered[valid],
                    xAxis_label='Fraction Disordered',
                    yAxis_label='RNA motif score',
                    name='fraction_disordered_vs_rnamotif_score',
                    yAxis_score=motif_score[valid]
                )

        # Add 3D scatter plot: Chain Pair PAE, max RNA ipTM, RNAmotif_score
        if 'af3_chain_pair_pae_min' in df.columns and 'af3_rna_ipTM' in df.columns and 'RNAmotif_score' in df.columns:
            chain_pair_pae = pd.to_numeric(df['af3_chain_pair_pae_min'], errors='coerce')
            rna_iptm_max = df['af3_rna_ipTM'].apply(max_metric)
            motif_score = pd.to_numeric(df['RNAmotif_score'], errors='coerce')
            valid = (~chain_pair_pae.isna()) & (~rna_iptm_max.isna()) & (~motif_score.isna())
            if valid.any():
                plot_creator.get_3d_scatterplotA(
                    table_source='af3_metrics',
                    xAxis_score=chain_pair_pae[valid],
                    yAxis_score=rna_iptm_max[valid],
                    zAxis_score=motif_score[valid],
                    xAxis_label='Chain Pair PAE (smaller=better)',
                    yAxis_label='max RNA ipTM (bigger=better)',
                    zAxis_label='RNA motif score',
                    name='chain_pair_pae_rna_ipTMmax_vs_rnamotif_score_3d',
                )

        # Add MSA color-coded 3D scatter plot: Chain Pair PAE, max RNA ipTM, RNAmotif_score
        if all(col in df.columns for col in ['af3_chain_pair_pae_min', 'af3_rna_ipTM', 'RNAmotif_score', 'RNA_msa_size', 'RNAChainIDs_pred', 'RNALength']):
            def classify_msa(row):
                try:
                    msa_dict = row['RNA_msa_size']
                    if isinstance(msa_dict, str):
                        msa_dict = safe_eval(msa_dict)
                    if not isinstance(msa_dict, dict):
                        return np.nan

                    chain_ids = row['RNAChainIDs_pred']
                    if isinstance(chain_ids, str):
                        chain_ids = safe_eval(chain_ids)
                    # Ensure it's a list
                    if not isinstance(chain_ids, (list, tuple)):
                        chain_ids = [chain_ids]

                    rna_lengths = row['RNALength']
                    if isinstance(rna_lengths, str):
                        rna_lengths = safe_eval(rna_lengths)
                    # Ensure it's a list
                    if not isinstance(rna_lengths, (list, tuple)):
                        rna_lengths = [rna_lengths]

                    # If only one chain, just use its MSA value
                    if len(msa_dict) == 1:
                        msa_val = list(msa_dict.values())[0]
                    else:
                        # Find index of longest RNA
                        max_len = max(rna_lengths)
                        max_idx = rna_lengths.index(max_len)
                        chain_id = chain_ids[max_idx]
                        msa_val = msa_dict.get(str(chain_id), 1)

                    return '+MSA' if msa_val > 1 else '-MSA'
                except Exception as e:
                    print(f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                    return np.nan

            msa_status = df.apply(classify_msa, axis=1)
            chain_pair_pae = pd.to_numeric(df['af3_chain_pair_pae_min'], errors='coerce')
            rna_iptm_max = df['af3_rna_ipTM'].apply(max_metric)
            motif_score = pd.to_numeric(df['RNAmotif_score'], errors='coerce')
            valid = (~chain_pair_pae.isna()) & (~rna_iptm_max.isna()) & (~motif_score.isna()) & (~msa_status.isna())
            
            if valid.any():
                # Create 3D scatter plot with MSA color coding
                fig = plt.figure(figsize=(12, 8))
                ax = fig.add_subplot(111, projection='3d')
                
                # Separate data by MSA status
                msa_plus_mask = msa_status[valid] == '+MSA'
                msa_minus_mask = msa_status[valid] == '-MSA'
                
                # Plot +MSA points in blue
                if msa_plus_mask.any():
                    ax.scatter(chain_pair_pae[valid][msa_plus_mask], 
                              rna_iptm_max[valid][msa_plus_mask], 
                              motif_score[valid][msa_plus_mask], 
                              c='blue', label='+MSA', alpha=0.7, s=50)
                
                # Plot -MSA points in orange
                if msa_minus_mask.any():
                    ax.scatter(chain_pair_pae[valid][msa_minus_mask], 
                              rna_iptm_max[valid][msa_minus_mask], 
                              motif_score[valid][msa_minus_mask], 
                              c='orange', label='-MSA', alpha=0.7, s=50)
                
                ax.set_xlabel('Chain Pair PAE (smaller=better)')
                ax.set_ylabel('max RNA ipTM (bigger=better)')
                ax.set_zlabel('RNA motif score')
                ax.legend()
                plt.tight_layout()
                plt.savefig(os.path.join(plot_creator.results_dir, 'chain_pair_pae_rna_ipTMmax_vs_rnamotif_score_3d_msa.png'))
                plt.close()

        # Add 3D scatter plot: Chain Pair PAE, max RNA ipTM, Complex_TM
        if 'af3_chain_pair_pae_min' in df.columns and 'af3_rna_ipTM' in df.columns and 'Complex_TM' in df.columns:
            chain_pair_pae = pd.to_numeric(df['af3_chain_pair_pae_min'], errors='coerce')
            rna_iptm_max = df['af3_rna_ipTM'].apply(max_metric)
            complex_tm = pd.to_numeric(df['Complex_TM'], errors='coerce')
            valid = (~chain_pair_pae.isna()) & (~rna_iptm_max.isna()) & (~complex_tm.isna())
            if valid.any():
                plot_creator.get_3d_scatterplotA(
                    table_source='af3_metrics',
                    xAxis_score=chain_pair_pae[valid],
                    yAxis_score=rna_iptm_max[valid],
                    zAxis_score=complex_tm[valid],
                    xAxis_label='Chain Pair PAE (smaller=better)',
                    yAxis_label='max RNA ipTM (bigger=better)',
                    zAxis_label='Complex TM',
                    name='chain_pair_pae_rna_ipTMmax_vs_ComplexTM_3d',
                )

        # Add MSA color-coded 3D scatter plot: Chain Pair PAE, max RNA ipTM, Complex_TM
        if all(col in df.columns for col in ['af3_chain_pair_pae_min', 'af3_rna_ipTM', 'Complex_TM', 'RNA_msa_size', 'RNAChainIDs_pred', 'RNALength']):
            def classify_msa(row):
                try:
                    msa_dict = row['RNA_msa_size']
                    if isinstance(msa_dict, str):
                        msa_dict = safe_eval(msa_dict)
                    if not isinstance(msa_dict, dict):
                        return np.nan

                    chain_ids = row['RNAChainIDs_pred']
                    if isinstance(chain_ids, str):
                        chain_ids = safe_eval(chain_ids)
                    # Ensure it's a list
                    if not isinstance(chain_ids, (list, tuple)):
                        chain_ids = [chain_ids]

                    rna_lengths = row['RNALength']
                    if isinstance(rna_lengths, str):
                        rna_lengths = safe_eval(rna_lengths)
                    # Ensure it's a list
                    if not isinstance(rna_lengths, (list, tuple)):
                        rna_lengths = [rna_lengths]

                    # If only one chain, just use its MSA value
                    if len(msa_dict) == 1:
                        msa_val = list(msa_dict.values())[0]
                    else:
                        # Find index of longest RNA
                        max_len = max(rna_lengths)
                        max_idx = rna_lengths.index(max_len)
                        chain_id = chain_ids[max_idx]
                        msa_val = msa_dict.get(str(chain_id), 1)

                    return '+MSA' if msa_val > 1 else '-MSA'
                except Exception as e:
                    print(f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                    return np.nan

            msa_status = df.apply(classify_msa, axis=1)
            chain_pair_pae = pd.to_numeric(df['af3_chain_pair_pae_min'], errors='coerce')
            rna_iptm_max = df['af3_rna_ipTM'].apply(max_metric)
            complex_tm = pd.to_numeric(df['Complex_TM'], errors='coerce')
            valid = (~chain_pair_pae.isna()) & (~rna_iptm_max.isna()) & (~complex_tm.isna()) & (~msa_status.isna())
            
            if valid.any():
                # Create 3D scatter plot with MSA color coding
                fig = plt.figure(figsize=(12, 8))
                ax = fig.add_subplot(111, projection='3d')
                
                # Separate data by MSA status
                msa_plus_mask = msa_status[valid] == '+MSA'
                msa_minus_mask = msa_status[valid] == '-MSA'
                
                # Plot +MSA points in blue
                if msa_plus_mask.any():
                    ax.scatter(chain_pair_pae[valid][msa_plus_mask], 
                              rna_iptm_max[valid][msa_plus_mask], 
                              complex_tm[valid][msa_plus_mask], 
                              c='blue', label='+MSA', alpha=0.7, s=50)
                
                # Plot -MSA points in orange
                if msa_minus_mask.any():
                    ax.scatter(chain_pair_pae[valid][msa_minus_mask], 
                              rna_iptm_max[valid][msa_minus_mask], 
                              complex_tm[valid][msa_minus_mask], 
                              c='orange', label='-MSA', alpha=0.7, s=50)
                
                ax.set_xlabel('Chain Pair PAE (smaller=better)')
                ax.set_ylabel('max RNA ipTM (bigger=better)')
                ax.set_zlabel('Complex TM')
                ax.legend()
                plt.tight_layout()
                plt.savefig(os.path.join(plot_creator.results_dir, 'chain_pair_pae_rna_ipTMmax_vs_ComplexTM_3d_msa.png'))
                plt.close()

        # Plot Motif ipTM Match Score for multi-chain complexes if possible
        # if all(col in df.columns for col in ['NumberRNAs','RNAMotifLength','ChainIDpairList_proteinRNA_exp','RNAChainIDs','RNASequence_exp','RNASequence_pred','af3_rna_ipTM']):
        multi_chain_df = self.merged_df[self.merged_df['NumberRNAs'] > 1]
        # print(multi_chain_df)
        if not multi_chain_df.empty:
            motif_scores = multi_chain_df.apply(motif_ipTM_match_score, axis=1).dropna()
            n_perfect = (motif_scores == 1.0).sum()
            n_imperfect = (motif_scores != 1.0).sum()
            data = {
                'Perfect Match (1.0)': n_perfect,
                'Imperfect (<1.0)': n_imperfect
            }
            plot_creator.get_barplot(
                table_source='af3_metrics',
                data=data,
                title='Motif ipTM Match Score Distribution',
                xlabel='Motif ipTM Match',
                ylabel='Count',
                bar_color='#FF8F00',
                rotation=0
            )

        import matplotlib.pyplot as plt

        multi_chain_df = self.merged_df[self.merged_df['NumberRNAs'] > 1]
        if not multi_chain_df.empty:
            motif_scores = multi_chain_df.apply(motif_ipTM_match_score, axis=1)

            # Helper for max ipTM
            def get_max_iptm(val):
                if isinstance(val, str):
                    try:
                        vals = ast.literal_eval(val)
                    except Exception:
                        return np.nan
                else:
                    vals = val
                if isinstance(vals, (list, tuple, np.ndarray)):
                    return np.nanmax([float(x) for x in vals if x is not None])
                elif vals is not None:
                    return float(vals)
                return np.nan

            # Helper for chain_pair_pae_min (single value per row)
            def get_chain_pair_pae(val):
                if isinstance(val, str):
                    try:
                        return float(val)
                    except Exception:
                        return np.nan
                elif val is not None:
                    return float(val)
                return np.nan

            # Grouping
            perfect_mask = motif_scores == 1.0
            imperfect_mask = motif_scores != 1.0

            # --- Plot 1: max af3_rna_ipTM ---
            max_iptm = multi_chain_df['af3_rna_ipTM'].apply(get_max_iptm)
            perfect_iptm = max_iptm[perfect_mask].dropna()
            imperfect_iptm = max_iptm[imperfect_mask].dropna()
            data_iptm = [perfect_iptm if not perfect_iptm.empty else pd.Series(dtype=float),
                         imperfect_iptm if not imperfect_iptm.empty else pd.Series(dtype=float)]
            labels = ['Perfect Match (1.0)', 'Imperfect (<1.0)']

            # Percentages
            n_total = len(motif_scores.dropna())
            n_perfect = perfect_mask.sum()
            n_imperfect = imperfect_mask.sum()
            percent_perfect = 100 * n_perfect / n_total if n_total else 0
            percent_imperfect = 100 * n_imperfect / n_total if n_total else 0
            percentages = [percent_perfect, percent_imperfect]

            # Violin plot with annotation below x-ticks
            fig, ax = plt.subplots(figsize=(8, 6))
            parts = ax.violinplot(data_iptm, showmeans=False, showmedians=True)
            ax.set_xticks([1, 2])
            ax.set_xticklabels(labels)
            ax.set_ylabel('max af3_rna_ipTM')
            ax.set_title('max af3_rna_ipTM by Motif Match')
            # Annotate percentages below x-ticks
            for i, pct in enumerate(percentages):
                ax.text(i + 1, -0.08, f'{pct:.1f}%', ha='center', va='top', fontsize=12, color='black',
                        transform=ax.get_xaxis_transform())
            plt.tight_layout()
            plt.savefig(os.path.join(plot_creator.results_dir, 'max_af3_rna_ipTM_by_motif_match.png'))
            plt.close()

            # --- Plot 2: af3_chain_pair_pae_min ---
            chain_pair_pae = multi_chain_df['af3_chain_pair_pae_min'].apply(get_chain_pair_pae)
            perfect_pae = chain_pair_pae[perfect_mask].dropna()
            imperfect_pae = chain_pair_pae[imperfect_mask].dropna()
            data_pae = [perfect_pae if not perfect_pae.empty else pd.Series(dtype=float),
                        imperfect_pae if not imperfect_pae.empty else pd.Series(dtype=float)]

            fig, ax = plt.subplots(figsize=(8, 6))
            parts = ax.violinplot(data_pae, showmeans=False, showmedians=True)
            ax.set_xticks([1, 2])
            ax.set_xticklabels(labels)
            ax.set_ylabel('af3_chain_pair_pae_min')
            ax.set_title('max_af3_rna_ipTM_pae by Motif Match')
            for i, pct in enumerate(percentages):
                ax.text(i + 1, -0.08, f'{pct:.1f}%', ha='center', va='top', fontsize=12, color='black',
                        transform=ax.get_xaxis_transform())
            plt.tight_layout()
            plt.savefig(os.path.join(plot_creator.results_dir, 'af3_chain_pair_pae_min_by_motif_match.png'))
            plt.close()

        from matplotlib.patches import Patch

        multi_chain_df = self.merged_df[self.merged_df['NumberRNAs'] == 1]
        if not multi_chain_df.empty:
            def get_chain_pair_pae(val):
                if isinstance(val, str):
                    try:
                        return float(val)
                    except Exception:
                        return np.nan
                elif val is not None:
                    return float(val)
                return np.nan

            chain_pair_pae = multi_chain_df['af3_chain_pair_pae_min'].apply(get_chain_pair_pae)
            rna_motif_score = pd.to_numeric(multi_chain_df['RNAmotif_score'], errors='coerce')
            complex_tm_score = pd.to_numeric(multi_chain_df['Complex_TM'], errors='coerce')

            # Classification
            pae_le_10_mask = chain_pair_pae <= 10
            pae_gt_10_mask = chain_pair_pae > 10

            motif_score_le_10 = rna_motif_score[pae_le_10_mask].dropna()
            motif_score_gt_10 = rna_motif_score[pae_gt_10_mask].dropna()
            complex_tm_le_10 = complex_tm_score[pae_le_10_mask].dropna()
            complex_tm_gt_10 = complex_tm_score[pae_gt_10_mask].dropna()

            # Percentages
            n_total = len(chain_pair_pae.dropna())
            n_le_10 = pae_le_10_mask.sum()
            n_gt_10 = pae_gt_10_mask.sum()
            percent_le_10 = 100 * n_le_10 / n_total if n_total else 0
            percent_gt_10 = 100 * n_gt_10 / n_total if n_total else 0
            percentages = [percent_le_10, percent_gt_10]

            # Prepare data and positions: for each group, orange (motif), blue (tm)
            data = [
                motif_score_le_10 if not motif_score_le_10.empty else pd.Series(dtype=float),
                complex_tm_le_10 if not complex_tm_le_10.empty else pd.Series(dtype=float),
                motif_score_gt_10 if not motif_score_gt_10.empty else pd.Series(dtype=float),
                complex_tm_gt_10 if not complex_tm_gt_10.empty else pd.Series(dtype=float)
            ]
            positions = [1, 2, 4, 5]
            colors = ['#FF8F00', '#1976D2', '#FF8F00', '#1976D2']

            fig, ax = plt.subplots(figsize=(10, 6))
            parts = ax.violinplot(data, positions=positions, showmeans=False, showmedians=True, widths=0.7)

            # Set correct colors: alternate orange/blue
            for i, pc in enumerate(parts['bodies']):
                pc.set_facecolor(colors[i])
                pc.set_alpha(0.7)

            # X-tick labels and annotation
            ax.set_xticks([1.5, 4.5])
            ax.set_xticklabels(['PAE ≤ 10', 'PAE > 10'])
            ax.set_ylabel('Score')
            ax.set_title('RNAmotif_score and Complex_TM by Chain Pair PAE')
            # Annotate percentages below x-ticks
            for i, pct in enumerate(percentages):
                ax.text(1.5 + 3 * i, -0.08, f'{pct:.1f}%', ha='center', va='top', fontsize=12, color='black',
                        transform=ax.get_xaxis_transform())
            # Add legend
            legend_patches = [Patch(facecolor='#FF8F00', label='RNAmotif_score'),
                              Patch(facecolor='#1976D2', label='Complex_TM')]
            ax.legend(handles=legend_patches, loc='upper right')
            plt.tight_layout()
            plt.savefig(os.path.join(plot_creator.results_dir, 'RNAmotif_and_ComplexTM_by_chain_pair_pae.png'))
            plt.close()

        def classify_msa(row):
            try:
                msa_dict = row['RNA_msa_size']
                if isinstance(msa_dict, str):
                    msa_dict = safe_eval(msa_dict)
                if not isinstance(msa_dict, dict):
                    return np.nan

                chain_ids = row['RNAChainIDs_pred']
                if isinstance(chain_ids, str):
                    chain_ids = safe_eval(chain_ids)
                # Ensure it's a list
                if not isinstance(chain_ids, (list, tuple)):
                    chain_ids = [chain_ids]

                rna_lengths = row['RNALength']
                if isinstance(rna_lengths, str):
                    rna_lengths = safe_eval(rna_lengths)
                # Ensure it's a list
                if not isinstance(rna_lengths, (list, tuple)):
                    rna_lengths = [rna_lengths]

                # If only one chain, just use its MSA value
                if len(msa_dict) == 1:
                    msa_val = list(msa_dict.values())[0]
                else:
                    # Find index of longest RNA
                    max_len = max(rna_lengths)
                    max_idx = rna_lengths.index(max_len)
                    chain_id = chain_ids[max_idx]
                    msa_val = msa_dict.get(str(chain_id), 1)

                return '+MSA' if msa_val > 1 else '-MSA'
            except Exception as e:
                print(
                    f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                return np.nan

        # Prepare data
        df = self.pred_df  # or self.merged_df if you want to include experimental columns
        if all(col in df.columns for col in ['af3_rna_pTM', 'RNA_TM', 'RNA_msa_size', 'RNAChainIDs_pred', 'RNALength']):
            msa_status = df.apply(classify_msa, axis=1)
            rna_ptm = df.apply(lambda row: process_metric_by_length(row['af3_rna_pTM'], row[
                'RNALength']) if 'RNALength' in df.columns and not pd.isna(row['af3_rna_pTM']) and not pd.isna(
                row['RNALength']) else None, axis=1)
            rna_tm = pd.to_numeric(df['RNA_TM'], errors='coerce')
            valid = (~rna_ptm.isna()) & (~rna_tm.isna()) & (~msa_status.isna())

            color_map = {'+MSA': 'blue', '-MSA': 'orange'}

            # Call the color scatterplot
            plot_creator.get_msa_scatterplot(
                table_source='af3_metrics',
                xAxis_score=rna_ptm[valid],
                xAxis_label='RNA pTM',
                name='msa_rna_pTM_vs_rna_TM_all_chains',
                yAxis_label='RNA TM',
                yAxis_score=rna_tm[valid],
                color_by=msa_status[valid],
                color_map=color_map
            )

        # Add MSA color-coded 3D scatter plot: Chain Pair PAE, Global PAE, RNAmotif_score
        if all(col in df.columns for col in
               ['af3_chain_pair_pae_min', 'af3_global_pae_avg', 'RNAmotif_score', 'RNA_msa_size',
                'RNAChainIDs_pred', 'RNALength']):
            def classify_msa(row):
                try:
                    msa_dict = row['RNA_msa_size']
                    if isinstance(msa_dict, str):
                        msa_dict = safe_eval(msa_dict)
                    if not isinstance(msa_dict, dict):
                        return np.nan

                    chain_ids = row['RNAChainIDs_pred']
                    if isinstance(chain_ids, str):
                        chain_ids = safe_eval(chain_ids)
                    if not isinstance(chain_ids, (list, tuple)):
                        chain_ids = [chain_ids]

                    rna_lengths = row['RNALength']
                    if isinstance(rna_lengths, str):
                        rna_lengths = safe_eval(rna_lengths)
                    if not isinstance(rna_lengths, (list, tuple)):
                        rna_lengths = [rna_lengths]

                    if len(msa_dict) == 1:
                        msa_val = list(msa_dict.values())[0]
                    else:
                        max_len = max(rna_lengths)
                        max_idx = rna_lengths.index(max_len)
                        chain_id = chain_ids[max_idx]
                        msa_val = msa_dict.get(str(chain_id), 1)

                    return '+MSA' if msa_val > 1 else '-MSA'
                except Exception as e:
                    print(
                        f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                    return np.nan

            msa_status = df.apply(classify_msa, axis=1)
            chain_pair_pae = pd.to_numeric(df['af3_chain_pair_pae_min'], errors='coerce')
            global_pae = pd.to_numeric(df['af3_global_pae_avg'], errors='coerce')
            motif_score = pd.to_numeric(df['RNAmotif_score'], errors='coerce')
            valid = (~chain_pair_pae.isna()) & (~global_pae.isna()) & (~motif_score.isna()) & (~msa_status.isna())

            if valid.any():
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=(12, 8))
                ax = fig.add_subplot(111, projection='3d')
                msa_plus_mask = msa_status[valid] == '+MSA'
                msa_minus_mask = msa_status[valid] == '-MSA'
                if msa_plus_mask.any():
                    ax.scatter(chain_pair_pae[valid][msa_plus_mask],
                               global_pae[valid][msa_plus_mask],
                               motif_score[valid][msa_plus_mask],
                               c='blue', label='+MSA', alpha=0.7, s=50)
                if msa_minus_mask.any():
                    ax.scatter(chain_pair_pae[valid][msa_minus_mask],
                               global_pae[valid][msa_minus_mask],
                               motif_score[valid][msa_minus_mask],
                               c='orange', label='-MSA', alpha=0.7, s=50)
                ax.set_xlabel('Chain Pair PAE (smaller=better)')
                ax.set_ylabel('Global PAE (smaller=better)')
                ax.set_zlabel('RNA motif score')
                ax.legend()
                plt.tight_layout()
                plt.savefig(os.path.join(plot_creator.results_dir,
                                         'chain_pair_pae_global_pae_vs_rnamotif_score_3d_msa.png'))
                plt.close()

        # Add MSA color-coded 3D scatter plot: Chain Pair PAE, Global PAE, RNA_RMSD
        if all(col in df.columns for col in
               ['af3_chain_pair_pae_min', 'af3_global_pae_avg', 'RNA_RMSD', 'RNA_msa_size', 'RNAChainIDs_pred',
                'RNALength']):

            def parse_rna_rmsd(val):
                if pd.isna(val):
                    return np.nan
                try:
                    if isinstance(val, str):
                        if '[' in val:
                            values = ast.literal_eval(val)
                            if isinstance(values, list) and len(values) > 0:
                                return float(values[0])
                        return float(val)
                    elif isinstance(val, (list, tuple)):
                        return float(val[0]) if len(val) > 0 else np.nan
                    return float(val)
                except Exception:
                    return np.nan

            def classify_msa(row):
                try:
                    msa_dict = row['RNA_msa_size']
                    if isinstance(msa_dict, str):
                        msa_dict = safe_eval(msa_dict)
                    if not isinstance(msa_dict, dict):
                        return np.nan

                    chain_ids = row['RNAChainIDs_pred']
                    if isinstance(chain_ids, str):
                        chain_ids = safe_eval(chain_ids)
                    if not isinstance(chain_ids, (list, tuple)):
                        chain_ids = [chain_ids]

                    rna_lengths = row['RNALength']
                    if isinstance(rna_lengths, str):
                        rna_lengths = safe_eval(rna_lengths)
                    if not isinstance(rna_lengths, (list, tuple)):
                        rna_lengths = [rna_lengths]

                    if len(msa_dict) == 1:
                        msa_val = list(msa_dict.values())[0]
                    else:
                        max_len = max(rna_lengths)
                        max_idx = rna_lengths.index(max_len)
                        chain_id = chain_ids[max_idx]
                        msa_val = msa_dict.get(str(chain_id), 1)

                    return '+MSA' if msa_val > 1 else '-MSA'
                except Exception as e:
                    print(
                        f"Error in classify_msa: {e}. Row data: RNALength={row.get('RNALength')}, RNAChainIDs_pred={row.get('RNAChainIDs_pred')}")
                    return np.nan

            msa_status = df.apply(classify_msa, axis=1)
            chain_pair_pae = pd.to_numeric(df['af3_chain_pair_pae_min'], errors='coerce')
            global_pae = pd.to_numeric(df['af3_global_pae_avg'], errors='coerce')
            rna_rmsd = df['RNA_RMSD'].apply(parse_rna_rmsd)
            valid = (~chain_pair_pae.isna()) & (~global_pae.isna()) & (~rna_rmsd.isna()) & (~msa_status.isna())

            if valid.any():
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=(12, 8))
                ax = fig.add_subplot(111, projection='3d')
                msa_plus_mask = msa_status[valid] == '+MSA'
                msa_minus_mask = msa_status[valid] == '-MSA'
                if msa_plus_mask.any():
                    ax.scatter(chain_pair_pae[valid][msa_plus_mask],
                               global_pae[valid][msa_plus_mask],
                               rna_rmsd[valid][msa_plus_mask],
                               c='blue', label='+MSA', alpha=0.7, s=50)
                if msa_minus_mask.any():
                    ax.scatter(chain_pair_pae[valid][msa_minus_mask],
                               global_pae[valid][msa_minus_mask],
                               rna_rmsd[valid][msa_minus_mask],
                               c='orange', label='-MSA', alpha=0.7, s=50)
                ax.set_xlabel('Chain Pair PAE (smaller=better)')
                ax.set_ylabel('Global PAE (smaller=better)')
                ax.set_zlabel('RNA RMSD')
                ax.legend()
                plt.tight_layout()
                plt.savefig(
                    os.path.join(plot_creator.results_dir, 'chain_pair_pae_global_pae_vs_rna_rmsd_3d_msa.png'))
                plt.close()

    def plot_global_pae_metrics(self, single_chain_only=False):
        """Plot global PAE metrics separately"""
        df = self.filter_data(self.pred_df, single_chain_only)

        # Define metrics to plot (excluding PAE)
        metrics = {
            'Global PAE': 'af3_global_pae_avg'
        }

        # Process data for each metric
        data_dict = {}
        for label, column in metrics.items():
            values = df[column].dropna()
            if len(values) > 0:
                data_dict[label] = values

        # Create violin plot using PlotCreator with af3_metrics type
        plot_creator = PlotCreator('af3_metrics', self.msa_option, single_chain_only)
        plot_creator.set_plot_style()

        # Convert data for violin plot
        data_values = [data_dict[key] for key in data_dict.keys()]
        labels = list(data_dict.keys())

        # Clean data_values to ensure all are numeric
        cleaned_data_values = []
        for vals in data_values:
            numeric_vals = pd.to_numeric(vals, errors='coerce').dropna()
            cleaned_data_values.append(numeric_vals)

        plot_creator.get_violin_plot(
            table_source='pae_metrics',
            data_values=cleaned_data_values,
            labels=labels,
            name='global_pae',
            yAxis_label='PAE Score'
        )

    def plot_bond_metrics(self, single_chain_only=False):
        """Plot RNA bond metrics (glycosidic, sugar pucker, gamma)"""
        plot_creator = PlotCreator('bond_metrics', self.msa_option, single_chain_only)
        plot_creator.set_plot_style()

        # Change filter_chains to filter_data
        df = self.filter_data(self.merged_df, single_chain_only)
        # print("df", self.merged_df)

        # Plot glycosidic bond
        self._plot_glycosidic_bond(plot_creator, df, single_chain_only)

        # Plot sugar pucker
        self._plot_sugar_pucker(plot_creator, df, single_chain_only)

        # Plot gamma angle
        self._plot_gamma_angle(plot_creator, df, single_chain_only)

    def _plot_glycosidic_bond(self, plot_creator, df, single_chain_only):
        """Plot glycosidic bond distributions"""
        categories = ['anti', 'syn']

        # Process glycosidic bond data
        chi_counts_pred = {cat: 0 for cat in categories}
        chi_counts_exp = {cat: 0 for cat in categories}

        for chi_dict_str in df['RNA_GlycosidicBond_pred'].dropna():
            counts = self.parse_dict_string(chi_dict_str)
            for cat in categories:
                chi_counts_pred[cat] += counts.get(cat, 0)

        for chi_dict_str in df['RNA_GlycosidicBond_exp'].dropna():
            counts = self.parse_dict_string(chi_dict_str)
            for cat in categories:
                chi_counts_exp[cat] += counts.get(cat, 0)

        # Calculate percentages
        total_pred = sum(chi_counts_pred.values())
        total_exp = sum(chi_counts_exp.values())

        chi_percent_pred = {cat: (count / total_pred * 100 if total_pred > 0 else 0)
                            for cat, count in chi_counts_pred.items()}
        chi_percent_exp = {cat: (count / total_exp * 100 if total_exp > 0 else 0)
                           for cat, count in chi_counts_exp.items()}

        # Create plot
        plt.figure(figsize=(8, 6))
        x = range(len(categories))
        plt.bar(x, [chi_percent_exp[cat] for cat in categories], color='orange', label='Experimental', width=0.35)
        plt.bar([i + 0.35 for i in x], [chi_percent_pred[cat] for cat in categories],
                color='green', alpha=0.7, label='Predicted', width=0.35)
        plt.xticks([i + 0.175 for i in x], categories)
        plt.xlabel("Chi Orientation")
        plt.ylabel("Percentage (%)")
        plt.legend()
        plt.tight_layout()

        suffix = '_single_chain' if single_chain_only else ''
        plot_creator.save_plot(f'glycosidic_bond{suffix}')

    def _plot_sugar_pucker(self, plot_creator, df, single_chain_only):
        """Plot sugar pucker distributions"""
        categories = ["C3'-endo", "C2'-endo", "C3'-exo", "C4'-exo", "C2'-exo", "O4'-endo", "O4'-exo"]

        # Process data
        pucker_counts_pred = {cat: 0 for cat in categories}
        pucker_counts_exp = {cat: 0 for cat in categories}

        for pucker_dict_str in df['RNA_SugarPucker_pred'].dropna():
            counts = self.parse_dict_string(pucker_dict_str)
            for cat in categories:
                pucker_counts_pred[cat] += counts.get(cat, 0)

        for pucker_dict_str in df['RNA_SugarPucker_exp'].dropna():
            counts = self.parse_dict_string(pucker_dict_str)
            for cat in categories:
                pucker_counts_exp[cat] += counts.get(cat, 0)

        # Calculate percentages and create plot
        total_pred = sum(pucker_counts_pred.values())
        total_exp = sum(pucker_counts_exp.values())

        pucker_percent_pred = {cat: (count / total_pred * 100 if total_pred > 0 else 0)
                               for cat, count in pucker_counts_pred.items()}
        pucker_percent_exp = {cat: (count / total_exp * 100 if total_exp > 0 else 0)
                              for cat, count in pucker_counts_exp.items()}

        plt.figure(figsize=(12, 6))
        x = range(len(categories))
        plt.bar(x, [pucker_percent_exp[cat] for cat in categories], color='orange', label='Experimental', width=0.35)
        plt.bar([i + 0.35 for i in x], [pucker_percent_pred[cat] for cat in categories],
                color='green', alpha=0.7, label='Predicted', width=0.35)
        plt.xticks([i + 0.175 for i in x], categories, rotation=45, ha='right')
        plt.xlabel("Sugar Pucker Type")
        plt.ylabel("Percentage (%)")
        plt.legend()
        plt.tight_layout()

        suffix = '_single_chain' if single_chain_only else ''
        plot_creator.save_plot(f'sugar_pucker{suffix}')

    def _plot_gamma_angle(self, plot_creator, df, single_chain_only):
        """Plot gamma angle distributions"""
        categories = ['g+', 'g-', 't']

        # Process data
        gamma_counts_pred = {cat: 0 for cat in categories}
        gamma_counts_exp = {cat: 0 for cat in categories}

        for gamma_dict_str in df['RNA_GammaAngle_pred'].dropna():
            counts = self.parse_dict_string(gamma_dict_str)
            for cat in categories:
                gamma_counts_pred[cat] += counts.get(cat, 0)

        for gamma_dict_str in df['RNA_GammaAngle_exp'].dropna():
            counts = self.parse_dict_string(gamma_dict_str)
            for cat in categories:
                gamma_counts_exp[cat] += counts.get(cat, 0)

        # Calculate percentages and create plot
        total_pred = sum(gamma_counts_pred.values())
        total_exp = sum(gamma_counts_exp.values())

        gamma_percent_pred = {cat: (count / total_pred * 100 if total_pred > 0 else 0)
                              for cat, count in gamma_counts_pred.items()}
        gamma_percent_exp = {cat: (count / total_exp * 100 if total_exp > 0 else 0)
                             for cat, count in gamma_counts_exp.items()}

        plt.figure(figsize=(8, 6))
        x = range(len(categories))
        plt.bar(x, [gamma_percent_exp[cat] for cat in categories], color='orange', label='Experimental', width=0.35)
        plt.bar([i + 0.35 for i in x], [gamma_percent_pred[cat] for cat in categories],
                color='green', alpha=0.7, label='Predicted', width=0.35)
        plt.xticks([i + 0.175 for i in x], categories)
        plt.xlabel("Gamma Angle Conformation")
        plt.ylabel("Percentage (%)")
        plt.legend()
        plt.tight_layout()

        suffix = '_single_chain' if single_chain_only else ''
        plot_creator.save_plot(f'gamma_angle{suffix}')

    def plot_rna_sugar_pucker(self, single_chain_only=False):
        """Plot sugar pucker distributions as percentages"""
        df = self.filter_data(self.merged_df, single_chain_only)
        categories = ["C3'-endo", "C2'-endo", "C3'-exo", "C4'-exo", "C2'-exo", "O4'-endo", "O4'-exo"]

        # Parse dictionary strings and aggregate counts
        pucker_counts_pred = {cat: 0 for cat in categories}
        pucker_counts_exp = {cat: 0 for cat in categories}

        for pucker_dict_str in df['RNA_SugarPucker_pred'].dropna():
            counts = self.parse_dict_string(pucker_dict_str)
            for cat in categories:
                pucker_counts_pred[cat] += counts.get(cat, 0)

        for pucker_dict_str in df['RNA_SugarPucker_exp'].dropna():
            counts = self.parse_dict_string(pucker_dict_str)
            for cat in categories:
                pucker_counts_exp[cat] += counts.get(cat, 0)

        # Calculate percentages
        total_pred = sum(pucker_counts_pred.values())
        total_exp = sum(pucker_counts_exp.values())

        pucker_percent_pred = {cat: (count / total_pred * 100 if total_pred > 0 else 0)
                               for cat, count in pucker_counts_pred.items()}
        pucker_percent_exp = {cat: (count / total_exp * 100 if total_exp > 0 else 0)
                              for cat, count in pucker_counts_exp.items()}

        # Create plot
        plt.figure(figsize=(12, 6))
        x = range(len(categories))
        plt.bar(x, [pucker_percent_exp[cat] for cat in categories], color='yellow', label='Experimental', width=0.35)
        plt.bar([i + 0.35 for i in x], [pucker_percent_pred[cat] for cat in categories],
                color='green', alpha=0.7, label='Predicted', width=0.35)
        plt.xticks([i + 0.175 for i in x], categories, rotation=45, ha='right')
        plt.xlabel("Sugar Pucker Type")
        plt.ylabel("Percentage (%)")
        plt.legend()
        plt.tight_layout(pad=2.0)

        suffix = '_single_chain' if single_chain_only else ''
        plt.savefig(os.path.join(self.results_dir, f'sugar_pucker{suffix}.png'))
        plt.close()

    def plot_rna_gamma_angle(self, single_chain_only=False):
        """Plot gamma angle distributions as percentages"""
        df = self.filter_data(self.merged_df, single_chain_only)
        categories = ['g+', 'g-', 't']

        # Parse dictionary strings and aggregate counts
        gamma_counts_pred = {cat: 0 for cat in categories}
        gamma_counts_exp = {cat: 0 for cat in categories}

        for gamma_dict_str in df['RNA_GammaAngle_pred'].dropna():
            counts = self.parse_dict_string(gamma_dict_str)
            for cat in categories:
                gamma_counts_pred[cat] += counts.get(cat, 0)

        for gamma_dict_str in df['RNA_GammaAngle_exp'].dropna():
            counts = self.parse_dict_string(gamma_dict_str)
            for cat in categories:
                gamma_counts_exp[cat] += counts.get(cat, 0)

        # Calculate percentages
        total_pred = sum(gamma_counts_pred.values())
        total_exp = sum(gamma_counts_exp.values())

        gamma_percent_pred = {cat: (count / total_pred * 100 if total_pred > 0 else 0)
                              for cat, count in gamma_counts_pred.items()}
        gamma_percent_exp = {cat: (count / total_exp * 100 if total_exp > 0 else 0)
                             for cat, count in gamma_counts_exp.items()}

        # Create plot
        plt.figure(figsize=(8, 6))
        x = range(len(categories))
        plt.bar(x, [gamma_percent_exp[cat] for cat in categories], color='yellow', label='Experimental', width=0.35)
        plt.bar([i + 0.35 for i in x], [gamma_percent_pred[cat] for cat in categories],
                color='green', alpha=0.7, label='Predicted', width=0.35)
        plt.xticks([i + 0.175 for i in x], categories)
        plt.xlabel("Gamma Angle Conformation")
        plt.ylabel("Percentage (%)")
        plt.legend()
        plt.tight_layout(pad=2.0)

        suffix = '_single_chain' if single_chain_only else ''
        plt.savefig(os.path.join(self.results_dir, f'gamma_angle{suffix}.png'))
        plt.close()

    def plot_interface_metrics(self, single_chain_only=False):
        """Plot interface metrics (free energy and binding affinity)"""
        plot_creator = PlotCreator('binding_affinity', self.msa_option, single_chain_only)

        # Plot free energy
        df = self.filter_chains(self.merged_df, True)  # Always use single_chain=True

        plt.figure()  # Create new figure for free energy plot
        # Plot free energy if data exists
        if df['Free_energy_exp'].notna().sum() > 0 and df['Free_energy_pred'].notna().sum() > 0:
            plot_creator.get_scatterplot(
                table_source='pred_vs_exp',
                xAxis_score=df['Free_energy_exp'].values,
                xAxis_label='Experimental Free Energy (kcal/mol)',
                yAxis_label='Predicted Free Energy (kcal/mol)',
                name='free_energy',
                yAxis_score=df['Free_energy_pred'].values
            )
        else:
            print("No valid free energy data to plot")
            plt.close()  # Close empty figure if no data

        # Plot binding affinity separately
        plt.figure()  # Create new figure for binding affinity plot
        kd_exp = df['Binding_affinity_kd_exp'].apply(self.parse_binding_affinity)
        kd_pred = df['Binding_affinity_kd_pred'].apply(self.parse_binding_affinity)

        # Remove None values and sort
        mask = kd_exp.notna() & kd_pred.notna()
        kd_exp = kd_exp[mask]
        kd_pred = kd_pred[mask]

        if len(kd_exp) > 0:
            # Convert to log scale if values are very small
            if kd_exp.min() < 1e-10 or kd_pred.min() < 1e-10:
                kd_exp = np.log10(kd_exp)
                kd_pred = np.log10(kd_pred)
                xlabel = 'Log₁₀[Experimental Kd (M)]'
                ylabel = 'Log₁₀[Predicted Kd (M)]'
            else:
                xlabel = 'Experimental Kd (M)'
                ylabel = 'Predicted Kd (M)'

            plot_creator.get_scatterplot(
                table_source='pred_vs_exp',
                xAxis_score=kd_exp.values,
                xAxis_label=xlabel,
                yAxis_label=ylabel,
                name='binding_affinity',
                yAxis_score=kd_pred.values
            )
        else:
            print("No valid binding affinity data to plot")
            plt.close()  # Close empty figure if no data

    def parse_dict_string(self, dict_string):
        """Parse string representation of dictionary into actual dictionary"""
        try:
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

    def parse_pae_list(self, pae_string):
        """Parse PAE list string and return maximum value"""
        try:
            pae_list = literal_eval(pae_string)
            return max(pae_list) if pae_list else 0
        except:
            return 0

    def parse_binding_affinity(self, value):
        """Parse binding affinity values, handling various formats"""
        try:
            if pd.isna(value):
                return None
            if isinstance(value, (int, float)):
                return float(value)
            # Handle scientific notation strings
            if 'e' in value.lower():
                return float(value)
            # Handle other string formats
            return float(value)
        except (ValueError, TypeError):
            return None

    def filter_by_msa(self, df, msa_option):
        """Filter dataframe based on MSA option"""
        if msa_option not in ['+MSA', '-MSA']:
            return df

        def get_msa_value(msa_str):
            try:
                if pd.isna(msa_str):
                    return None
                if isinstance(msa_str, str):
                    values = eval(msa_str)
                    return 1 if any(v == 1 for v in values) else 0
                return msa_str
            except:
                return None

        df['msa_value'] = df['af3_rna_MSA'].apply(get_msa_value)
        msa_filter = df['msa_value'] == (1 if msa_option == '+MSA' else 0)
        return df[msa_filter]

    def plot_domain_rmsd(self, single_chain_only=False):
        """Create bar chart of average RMSD by protein domain"""
        try:
            # Check if domain_names column exists
            check_query = """
            SELECT COUNT(*) 
            FROM pragma_table_info('exp_protein_rna') 
            WHERE name='domain_names'
            """
            cursor = self.database.connection.cursor()
            cursor.execute(check_query)
            has_domain = bool(cursor.fetchone()[0])

            if not has_domain:
                print("domain_names column not found in exp_protein_rna table")
                return

            # Get domain and RMSD data from SQLite
            query = """
            SELECT e.domain_names, p.RNA_RMSD
            FROM exp_protein_rna e
            JOIN pred_protein_rna p ON e.PDBId = p.exp_db_id
            WHERE e.domain_names IS NOT NULL 
            AND p.RNA_RMSD IS NOT NULL
            """
            cursor.execute(query)
            rows = cursor.fetchall()

            if not rows:
                print("No domain and RMSD data found")
                return

            # Process the data
            domain_rmsd_data = []
            for domain_names, rmsd in rows:
                if domain_names and rmsd:
                    # Split domain names
                    domains = domain_names.split(';')
                    # Get RMSD value using the extract_rmsd function
                    rmsd_value = self.extract_rmsd(rmsd)
                    if rmsd_value is not None:
                        # Add each domain with its RMSD
                        for domain in domains:
                            if domain.strip():  # Skip empty domains
                                domain_rmsd_data.append((domain.strip(), rmsd_value))

            if not domain_rmsd_data:
                print("No valid domain-RMSD pairs found")
                return

            # Convert to DataFrame for easier processing
            df = pd.DataFrame(domain_rmsd_data, columns=['domain', 'rmsd'])

            # Calculate statistics for each domain
            domain_stats = df.groupby('domain').agg({
                'rmsd': ['mean', 'std', 'count']
            }).reset_index()
            domain_stats.columns = ['domain', 'mean', 'std', 'count']

            # Sort by mean RMSD
            domain_stats = domain_stats.sort_values('mean')

            # Create plot using PlotCreator
            plot_creator = PlotCreator('domain_metrics', self.msa_option, single_chain_only)

            # Create larger bar plot with adjusted dimensions
            plt.figure(figsize=(20, 10))  # Increased figure size
            bars = plt.bar(
                range(len(domain_stats)),
                domain_stats['mean'],
                yerr=domain_stats['std'],
                capsize=5
            )

            # Customize plot with larger text and more space
            plt.xticks(
                range(len(domain_stats)),
                domain_stats['domain'],
                rotation=90,
                ha='right',
                fontsize=12  # Increased font size
            )
            plt.xlabel('Protein Domains', fontsize=14, labelpad=10)  # Added padding
            plt.ylabel('Average RNA RMSD [Å]', fontsize=14, labelpad=10)

            # Remove count labels from bars
            # Adjust margins to prevent label cutoff
            plt.margins(x=0.01)  # Reduce horizontal margins
            plt.tight_layout()

            # Save plot
            plot_creator.save_plot('domain_rmsd_distribution')
            plt.close()

        except Exception as e:
            print(f"Error creating domain RMSD plot: {e}")

    def extract_rmsd(self, rmsd_str):
        """Extract first RMSD value from string or float"""
        if pd.isna(rmsd_str):
            return None
        try:
            if isinstance(rmsd_str, (int, float)):
                return float(rmsd_str)
            # Handle list format
            if isinstance(rmsd_str, str) and '[' in rmsd_str:
                values = [float(x.strip()) for x in rmsd_str.strip('[]').split(',')]
                return values[0]  # Take first value for RMSD
            return float(rmsd_str)
        except:
            print(f"Could not parse RMSD: {rmsd_str}")
            return None

    def plot_domain_coverage(self, single_chain_only=False):
        """Plot how well predicted interaction residues align with domain positions"""
        try:
            # Get domain mappings first
            mapping_query = """
            SELECT DISTINCT domain_ids, domain_names 
            FROM exp_protein_rna 
            WHERE domain_ids IS NOT NULL 
            AND domain_names IS NOT NULL
            """
            cursor = self.database.connection.cursor()
            cursor.execute(mapping_query)
            mapping_rows = cursor.fetchall()

            # Create domain ID to name mapping
            domain_map = {}
            for ids, names in mapping_rows:
                id_list = ids.split(';')
                name_list = names.split(';')
                if len(id_list) == len(name_list):
                    domain_map.update(dict(zip(id_list, name_list)))

            # Get domain and AAMotif data with MSA column
            query = """
            SELECT e.PDBId, e.domain_ids, e.domain_pos, p.AAMotif, p.ChainIDpairList_proteinRNA, p.af3_rna_MSA
            FROM exp_protein_rna e
            JOIN pred_protein_rna p ON e.PDBId = p.exp_db_id
            WHERE e.domain_pos IS NOT NULL 
            AND p.AAMotif IS NOT NULL
            """
            cursor.execute(query)
            rows = cursor.fetchall()

            # Convert to DataFrame for filtering
            df = pd.DataFrame(rows, columns=['PDBId', 'domain_ids', 'domain_pos', 'AAMotif',
                                             'ChainIDpairList_proteinRNA', 'af3_rna_MSA'])

            # Apply filters
            df = self.filter_data(df, single_chain_only)

            if df.empty:
                print("No data after filtering")
                return

            # Process filtered data
            rows = df.values.tolist()

            # Process the data
            domain_coverage_data = []
            domain_coverage_types = {}
            pdb_coverage_data = {}
            pdb_categories = {
                'all_in_all': 0,  # All residues in all domains (multiple domains only)
                'partial_in_all': 0,  # Some residues in all domains (multiple domains only)
                'all_in_one': 0,  # All residues in one domain (any number of domains)
                'partial_in_one': 0,  # Some residues in one domain (any number of domains)
                'outside_all': 0  # No residues in any domain
            }

            for pdb_id, domain_ids, domain_pos, aa_motif, chain_pairs, msa in rows:
                try:
                    # Parse domain positions
                    domain_ranges = {}
                    domain_id_list = domain_ids.split(';')
                    domain_pos_list = domain_pos.split(';')

                    # Skip entries with invalid domain positions
                    if all(pos.strip() == '-' for pos in domain_pos_list):
                        continue

                    for d_id, d_pos in zip(domain_id_list, domain_pos_list):
                        d_id = d_id.strip()
                        if d_pos.strip() == '-':
                            continue

                        try:
                            # Extract all ranges for this domain
                            ranges = []
                            if ':' in d_pos:
                                range_str = d_pos.split(':')[1]
                            else:
                                range_str = d_pos

                            # Handle multiple ranges separated by comma
                            for r in range_str.split(','):
                                start, end = map(int, r.strip().split('-'))
                                ranges.append((start, end))

                            domain_ranges[d_id] = ranges
                            # print(f"Processing {pdb_id}:")
                            # print(f"  Domain {d_id}: ranges {ranges}")

                        except ValueError:
                            continue

                    if not domain_ranges:  # Skip if no valid domain ranges found
                        continue

                    # Parse AAMotif indices
                    # try:
                    # print(f"\nProcessing {pdb_id}:")
                    # print(f"Raw AAMotif data: {aa_motif}")

                    if isinstance(aa_motif, str):
                        # Handle string representation of list
                        if aa_motif.startswith('[') and aa_motif.endswith(']'):
                            aa_list = literal_eval(aa_motif)
                        else:
                            aa_list = [aa_motif]
                    else:
                        aa_list = [str(aa_motif)]

                    if not isinstance(aa_list, list):
                        aa_list = [aa_list]

                    # Process each motif list separately
                    for motif_list in aa_list:
                        aa_indices = []
                        if isinstance(motif_list, str):
                            # Extract numbers between parentheses
                            numbers = re.findall(r'\((\d+)\)', motif_list)
                            if numbers:
                                aa_indices = [int(num) for num in numbers]
                                # print(f"  AAMotif residue indices: {aa_indices}")
                            else:
                                # print(f"  No residue indices found in: {motif_list}")
                                continue

                        if aa_indices:
                            # Track coverage for this PDB
                            domains_with_all = 0  # Domains containing all residues
                            domains_with_partial = 0  # Domains containing some residues
                            total_domains = len(domain_ranges)
                            max_coverage = 0  # Initialize max coverage for domain plots

                            for domain_id, ranges in domain_ranges.items():
                                # Initialize coverage type counters if not exists
                                if domain_id not in domain_coverage_types:
                                    domain_coverage_types[domain_id] = {
                                        'within': 0,
                                        'partial': 0,
                                        'outside': 0
                                    }

                                # Count residues within this domain
                                indices_within = sum(1 for idx in aa_indices
                                                     if any(start <= idx <= end
                                                            for start, end in ranges))

                                # Update domain coverage types (for existing plots)
                                if indices_within == len(aa_indices):
                                    domain_coverage_types[domain_id]['within'] += 1
                                    domains_with_all += 1
                                elif indices_within > 0:
                                    domain_coverage_types[domain_id]['partial'] += 1
                                    domains_with_partial += 1
                                else:
                                    domain_coverage_types[domain_id]['outside'] += 1

                                # Calculate coverage for domain plots
                                coverage = indices_within / len(aa_indices)
                                domain_coverage_data.append((domain_id, coverage))
                                max_coverage = max(max_coverage, coverage)

                            # Store maximum coverage for PDB scatter plot
                            pdb_coverage_data[pdb_id] = max_coverage

                            # Categorize this PDB sample
                            if total_domains > 1:  # Multiple domains
                                if domains_with_all == total_domains:
                                    pdb_categories['all_in_all'] += 1
                                elif domains_with_partial == total_domains:
                                    pdb_categories['partial_in_all'] += 1
                                elif domains_with_all > 0:
                                    pdb_categories['all_in_one'] += 1
                                elif domains_with_partial > 0:
                                    pdb_categories['partial_in_one'] += 1
                                else:
                                    pdb_categories['outside_all'] += 1
                            else:  # Single domain
                                if domains_with_all == 1:
                                    pdb_categories['all_in_one'] += 1
                                elif domains_with_partial == 1:
                                    pdb_categories['partial_in_one'] += 1
                                else:
                                    pdb_categories['outside_all'] += 1

                except Exception as e:
                    print(f"Error processing {pdb_id}: {e}")
                    continue

            if not domain_coverage_data:
                print("No valid domain coverage data found")
                return

            # Create plots
            plot_creator = PlotCreator('domain_metrics', self.msa_option, single_chain_only)

            # Create coverage distribution plot
            plt.figure(figsize=(20, 10))
            domains = list(domain_coverage_types.keys())
            # Map domain IDs to names
            domain_names = [domain_map.get(d, d) for d in domains]

            within_counts = [domain_coverage_types[d]['within'] for d in domains]
            partial_counts = [domain_coverage_types[d]['partial'] for d in domains]
            outside_counts = [domain_coverage_types[d]['outside'] for d in domains]

            # Calculate percentages
            totals = np.array(within_counts) + np.array(partial_counts) + np.array(outside_counts)
            within_pct = 100 * np.array(within_counts) / totals
            partial_pct = 100 * np.array(partial_counts) / totals
            outside_pct = 100 * np.array(outside_counts) / totals

            x = range(len(domains))
            plt.bar(x, within_pct, label='Within Domain', color='green')
            plt.bar(x, partial_pct, bottom=within_pct, label='Partially Within', color='yellow')
            plt.bar(x, outside_pct, bottom=within_pct + partial_pct, label='Outside Domain', color='red')

            # Use domain names instead of IDs
            plt.xticks(x, domain_names, rotation=90, ha='right', fontsize=12)
            plt.xlabel('Protein Domains', fontsize=14, labelpad=10)
            plt.ylabel('Percentage of Samples (%)', fontsize=14, labelpad=10)
            plt.legend()
            plt.margins(x=0.01)
            plt.tight_layout()
            plot_creator.save_plot('domain_coverage_distribution')
            plt.close()

            # Create domain interaction coverage plot
            plt.figure(figsize=(20, 10))
            df = pd.DataFrame(domain_coverage_data, columns=['domain_id', 'coverage'])

            # Map domain IDs to names in DataFrame - fixed inplace warning
            df = df.assign(
                domain_name=df['domain_id'].map(domain_map).fillna(df['domain_id'])
            )

            domain_stats = df.groupby('domain_name').agg({
                'coverage': ['mean', 'std', 'count']
            }).reset_index()
            domain_stats.columns = ['domain_id', 'mean', 'std', 'count']
            domain_stats = domain_stats.sort_values('mean', ascending=False)

            plt.bar(
                range(len(domain_stats)),
                domain_stats['mean'],
                yerr=domain_stats['std'],
                capsize=5
            )

            plt.xticks(
                range(len(domain_stats)),
                domain_stats['domain_id'],
                rotation=90,
                ha='right',
                fontsize=12
            )
            plt.xlabel('Protein Domains', fontsize=14, labelpad=10)
            plt.ylabel('Average Coverage of Predicted Interactions', fontsize=14, labelpad=10)
            plt.ylim(0, 1)
            plt.margins(x=0.01)
            plt.tight_layout()
            plot_creator.save_plot('domain_interaction_coverage')
            plt.close()

            # Create PDB coverage scatter plot
            plt.figure(figsize=(15, 8))
            pdb_ids = list(pdb_coverage_data.keys())
            coverages = list(pdb_coverage_data.values())

            plt.scatter(range(len(pdb_ids)), coverages, alpha=0.6)
            plt.xticks(
                range(len(pdb_ids)),
                pdb_ids,
                rotation=90,
                ha='right',
                fontsize=10
            )
            plt.xlabel('PDB ID', fontsize=12)
            plt.ylabel('Maximum Domain Coverage', fontsize=12)
            plt.ylim(0, 1)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plot_creator.save_plot('pdb_domain_coverage')
            plt.close()

            # Add new plot for PDB categories
            plt.figure(figsize=(12, 8))
            categories = [
                'All residues\nin all domains\n(multiple domains)',
                'Partial residues\nin all domains\n(multiple domains)',
                'All residues in\none domain',
                'Partial residues in\none domain',
                'Outside all\ndomains'
            ]
            counts = [
                pdb_categories['all_in_all'],
                pdb_categories['partial_in_all'],
                pdb_categories['all_in_one'],
                pdb_categories['partial_in_one'],
                pdb_categories['outside_all']
            ]

            # Calculate percentages
            total_samples = sum(counts)
            percentages = [count / total_samples * 100 for count in counts]

            # Create bar plot
            bars = plt.bar(range(len(categories)), percentages)

            # Add count labels on top of bars
            for i, (bar, count) in enumerate(zip(bars, counts)):
                plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                         f'n={count}',
                         ha='center', va='bottom')

            plt.xticks(range(len(categories)), categories, rotation=45, ha='right')
            # plt.xlabel('Coverage Category', fontsize=12)
            plt.ylabel('Percentage of Samples (%)', fontsize=12)
            # plt.title('Distribution of Domain Coverage Patterns', fontsize=14)
            plt.tight_layout()
            plot_creator.save_plot('pdb_coverage_categories')
            plt.close()

        except Exception as e:
            print(f"Error creating domain coverage plot: {e}")

    def filter_data(self, df, single_chain_only):
        """Filter data based on chain pairs and MSA option"""
        filtered_df = df.copy()

        # Apply single chain filter if requested
        if single_chain_only:
            filtered_df = filtered_df[filtered_df['ChainIDpairList_proteinRNA'].apply(self.is_single_chain)]

        # Apply MSA filter if specified
        if self.msa_option == '+MSA':
            filtered_df = filtered_df[
                (filtered_df['af3_rna_MSA'] == '[1]') |
                (filtered_df['af3_rna_MSA'].str.startswith('[1,1', na=False))
                ]
        elif self.msa_option == '-MSA':
            filtered_df = filtered_df[
                (filtered_df['af3_rna_MSA'] == '[0]') |
                (filtered_df['af3_rna_MSA'].str.contains('0', na=False))
                ]

        return filtered_df

    def plot_rna_metrics(self, single_chain_only=False):
        """Create scatter plots for RNA metrics vs length and RMSD"""
        # Get filtered data first
        df = self.filter_data(self.pred_df, single_chain_only)

        # Create plot creator
        plot_creator = PlotCreator('rna_metrics', self.msa_option, single_chain_only)

        # Process RNA length values
        def extract_length(length_str):
            if pd.isna(length_str):
                return None
            try:
                if isinstance(length_str, (int, float)):
                    return int(length_str)
                if isinstance(length_str, str) and '[' in length_str:
                    values = [int(x.strip()) for x in length_str.strip('[]').split(',')]
                    return max(values)
                return int(length_str)
            except:
                return None

        # Process RMSD values
        def extract_rmsd(rmsd_str):
            if pd.isna(rmsd_str):
                return None
            try:
                if isinstance(rmsd_str, (int, float)):
                    return float(rmsd_str)
                if isinstance(rmsd_str, str) and '[' in rmsd_str:
                    values = [float(x.strip()) for x in rmsd_str.strip('[]').split(',')]
                    return values[0]
                return float(rmsd_str)
            except:
                return None

        # Apply transformations
        df['RNALength'] = df['RNALength'].apply(extract_length)
        df['RNA_RMSD'] = df['RNA_RMSD'].apply(extract_rmsd)

        # Drop rows with None values and convert to numeric
        df = df.dropna(subset=['RNALength', 'RNA_RMSD'])
        df['RNALength'] = pd.to_numeric(df['RNALength'])
        df['RNA_RMSD'] = pd.to_numeric(df['RNA_RMSD'])

        # Define metrics to plot
        metrics = {
            'F1-Score': 'RNA_f1_score',
            'Precision': 'RNA_precision',
            'Recall': 'RNA_recall',
            'MCC': 'RNA_mcc',
            'WL': 'RNA_wl'
        }

        # # Plot metrics vs length
        # for metric_name, metric_col in metrics.items():
        #     if not df[metric_col].isna().all():
        #         plot_creator.get_scatterplot(
        #             table_source='rna_metrics',
        #             xAxis_score=df['RNALength'],
        #             xAxis_label='RNA Length',
        #             yAxis_label=metric_name,
        #             name=f'{metric_col}_vs_length',
        #             yAxis_score=df[metric_col]
        #         )

        fig_length, axes_length = plt.subplots(2, 3, figsize=(15, 10))
        axes_length = axes_length.flatten()

        for idx, (metric_name, metric_col) in enumerate(metrics.items()):
            ax = axes_length[idx]
            plot_creator.scatter_plot(
                df['RNALength'],
                df[metric_col],
                ax=ax,
                xlabel='RNA Length',
                ylabel=metric_name
            )

        # Remove extra subplot
        if len(axes_length) > len(metrics):
            fig_length.delaxes(axes_length[-1])
        fig_length.tight_layout()

        # # Plot metrics vs RMSD
        # df_rmsd = df.dropna(subset=['RNA_RMSD'])
        # for metric_name, metric_col in metrics.items():
        #     if not df_rmsd[metric_col].isna().all():
        #         plot_creator.get_scatterplot(
        #             table_source='rna_metrics',
        #             xAxis_score=df_rmsd['RNA_RMSD'],
        #             xAxis_label='RNA RMSD',
        #             yAxis_label=metric_name,
        #             name=f'{metric_col}_vs_rmsd',
        #             yAxis_score=df_rmsd[metric_col]
        #         )

        fig_rmsd, axes_rmsd = plt.subplots(2, 3, figsize=(15, 10))
        axes_rmsd = axes_rmsd.flatten()

        # Filter out None RMSD values
        df_rmsd = df.dropna(subset=['RNA_RMSD'])

        for idx, (metric_name, metric_col) in enumerate(metrics.items()):
            ax = axes_rmsd[idx]
            plot_creator.scatter_plot(
                df_rmsd['RNA_RMSD'],
                df_rmsd[metric_col],
                ax=ax,
                xlabel='RNA RMSD',
                ylabel=metric_name
            )

        # Remove extra subplot
        if len(axes_rmsd) > len(metrics):
            fig_rmsd.delaxes(axes_rmsd[-1])

        fig_rmsd.tight_layout()

        # Save plots using plt.savefig directly
        plt.figure(fig_length.number)
        plt.savefig(os.path.join(plot_creator.results_dir, 'rna_metrics_vs_length.png'))

        plt.figure(fig_rmsd.number)
        plt.savefig(os.path.join(plot_creator.results_dir, 'rna_metrics_vs_rmsd.png'))

        plt.close('all')

    def analyze_domain_relationships(self, single_chain_only=False):
        """Analyze relationships between domain counts and various metrics"""
        plot_creator = PlotCreator('domain_count', self.msa_option, single_chain_only)

        # Query to get required data
        query = """
        SELECT p.Binding_affinity_kd, p.af3_rna_MSA, p.RNAmotif_score,
               p.ProteinRNAInterfaceArea, p.RNAMotifLength, e.domain_counts
        FROM pred_protein_rna p
        JOIN exp_protein_rna e ON p.exp_db_id = e.PDBId
        WHERE e.domain_counts IS NOT NULL
        """
        results = self.database.execute_query(query)

        # Process functions remain the same...
        def process_domain_counts(x):
            if pd.notna(x):
                try:
                    counts = [int(count) for count in x.split(';') if count.isdigit()]
                    return sum(counts) if counts else None
                except:
                    return None
            return None

        def process_motif_length(x):
            try:
                if '[' in x:
                    lengths = ast.literal_eval(x)
                    return sum(lengths)
                else:
                    return int(x)
            except:
                return None

        def process_binding_affinity(x):
            try:
                if x and x.lower() != 'none':
                    value = float(x)
                    return np.log10(value) if value > 0 else None
                return None
            except:
                return None

        def process_msa_status(x):
            try:
                if isinstance(x, str):
                    msa_list = ast.literal_eval(x)
                    if len(msa_list) == 1:
                        return 2 if msa_list[0] == 1 else 0
                    elif all(v == 1 for v in msa_list):
                        return 2
                    elif all(v == 0 for v in msa_list):
                        return 0
                    elif 1 in msa_list:
                        return 1
                    return 0
                return 0
            except:
                return 0

        # Create separate paired data for each plot
        motif_score_data = {'domain_counts': [], 'motif_score': []}
        interface_area_data = {'domain_counts': [], 'interface_area': []}
        binding_affinity_data = {'domain_counts': [], 'binding_affinity': []}
        
        # Collect data for MSA analysis
        msa_domain_data = {'domain_counts': [], 'msa_status': []}

        # Process results and create paired data
        for row in results:
            domain_count = process_domain_counts(row[5])
            
            if domain_count is not None:
                # For motif score plot
                motif_score = row[2]
                if motif_score is not None:
                    motif_score_data['domain_counts'].append(domain_count)
                    motif_score_data['motif_score'].append(motif_score)
                
                # For interface area plot
                interface_area = row[3]
                if interface_area is not None:
                    interface_area_data['domain_counts'].append(domain_count)
                    interface_area_data['interface_area'].append(interface_area)
                
                # For binding affinity plot
                binding_aff = process_binding_affinity(row[0])
                if binding_aff is not None:
                    binding_affinity_data['domain_counts'].append(domain_count)
                    binding_affinity_data['binding_affinity'].append(binding_aff)
                
                # For MSA analysis
                msa_status = process_msa_status(row[1])
                if msa_status is not None:
                    msa_domain_data['domain_counts'].append(domain_count)
                    msa_domain_data['msa_status'].append(msa_status)

        # Create scatter plots with paired data
        fig_scatter, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))

        # RNAmotif score vs domain count
        if motif_score_data['domain_counts']:
            plot_creator.scatter_plot(
                x=motif_score_data['domain_counts'],
                y=motif_score_data['motif_score'],
                ax=ax1,
                xlabel='Number of Domains',
                ylabel='RNA Motif Score'
            )

        # Interface area vs domain count
        if interface_area_data['domain_counts']:
            plot_creator.scatter_plot(
                x=interface_area_data['domain_counts'],
                y=interface_area_data['interface_area'],
                ax=ax2,
                xlabel='Number of Domains',
                ylabel='Protein-RNA Interface Area'
            )

        # Binding affinity vs domain count
        if binding_affinity_data['domain_counts']:
            min_aff = min(binding_affinity_data['binding_affinity'])
            ylabel = 'Log₁₀[Binding Affinity (M)]' if min_aff < 1e-10 else 'Binding Affinity (M)'
            plot_creator.scatter_plot(
                x=binding_affinity_data['domain_counts'],
                y=binding_affinity_data['binding_affinity'],
                ax=ax3,
                xlabel='Number of Domains',
                ylabel=ylabel
            )

        fig_scatter.tight_layout()
        plt.savefig(os.path.join(plot_creator.results_dir, 'domain_analysis_scatter_plots.png'))
        plt.close(fig_scatter)

        # Create boxplot for binding affinity vs domain count
        if binding_affinity_data['domain_counts']:
            unique_counts = sorted(set(binding_affinity_data['domain_counts']))
            affinity_by_count = {count: [] for count in unique_counts}
            for count, aff in zip(binding_affinity_data['domain_counts'], binding_affinity_data['binding_affinity']):
                affinity_by_count[count].append(aff)

            min_aff = min(min(vals) for vals in affinity_by_count.values() if vals)
            yAxis_label = 'Log₁₀[Binding Affinity (M)]' if min_aff < 1e-10 else 'Binding Affinity (M)'

            plot_creator.get_boxplot(
                table_source='domain_count',
                data_values=[affinity_by_count[count] for count in unique_counts],
                labels=[f'Domains: {count}' for count in unique_counts],
                name='binding_affinity_vs_domains',
                yAxis_label=yAxis_label,
                rotation=45
            )

        # Create segmented bar plot for MSA status vs domain count distribution
        if msa_domain_data['domain_counts']:
            msa_domain_counts = {
                'No RNA MSA': [],
                'Mixed RNA MSA': [],
                'Full RNA MSA': []
            }

            msa_labels = {0: 'No RNA MSA', 1: 'Mixed RNA MSA', 2: 'Full RNA MSA'}

            for msa, count in zip(msa_domain_data['msa_status'], msa_domain_data['domain_counts']):
                msa_domain_counts[msa_labels[msa]].append(count)

            segment_labels = sorted(set(msa_domain_data['domain_counts']))

            plot_creator.get_segmented_barplot(
                table_source='domain_count',
                data_dict=msa_domain_counts,
                category_labels=['No RNA MSA', 'Mixed RNA MSA', 'Full RNA MSA'],
                segment_labels=segment_labels,
                xlabel='MSA Status',
                ylabel='Number of Structures',
                name='domain_count_distribution_by_msa',
                rotation=0
            )


def _standardize_msa_option(msa_option):
    """Standardize MSA option to uppercase and ensure proper format"""
    if not msa_option:
        return None
    msa_option = msa_option.upper()
    if msa_option in ['+MSA', '-MSA', '+msa', '-msa']:
        # Check for minus sign explicitly
        return '-MSA' if msa_option.startswith('-') else '+MSA'
    return None


def main():
    if len(sys.argv) not in [2, 3]:
        print("Usage: python plotResults.py <all|single_chain> [+MSA|-MSA]")
        print("Options:")
        print("  all          : Generate plots for all chain combinations")
        print("  single_chain : Generate plots for single chain pairs only")
        print("  +MSA/+msa    : Plot only structures with MSA")
        print("  -MSA/-msa    : Plot only structures without MSA")
        sys.exit(1)

    plot_type = sys.argv[1]
    msa_option = _standardize_msa_option(sys.argv[2]) if len(sys.argv) == 3 else None
    single_chain_only = (plot_type == 'single_chain')

    # Initialize database and config
    database = DatabaseMethods()
    config = StartConfig()
    plotter = ResultPlotter(database, config, msa_option, single_chain_only)

    # Create plots
    plotter.plot_af3_metrics(single_chain_only)
    # plotter.plot_global_pae_metrics(single_chain_only)
    # plotter.plot_bond_metrics(single_chain_only)
    # plotter.plot_rna_metrics(single_chain_only)
    #
    # # if single_chain_only:
    # plotter.plot_interface_metrics(single_chain_only)
    #
    # domain = True
    # if domain:
    #     plotter.plot_domain_rmsd(single_chain_only)
    #     plotter.plot_domain_coverage(single_chain_only)
    #     plotter.analyze_domain_relationships(single_chain_only)

    print(f"\nPlots generated for {'single chain pairs' if single_chain_only else 'all chains'}")
    if msa_option:
        print(f"MSA filter: {msa_option}")
    print(f"Results saved in: {os.path.dirname(plotter.results_dir)}")

    # Close database connection
    database.close_connection()


if __name__ == "__main__":
    main()