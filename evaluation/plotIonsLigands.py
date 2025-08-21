import pandas as pd
import ast
import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from database.startConfig import StartConfig
from database.databaseMethods import DatabaseMethods
from plots.plotCreator import PlotCreator
import numpy as np
from collections import Counter, defaultdict

def extract_items(df, col):
    """Extract all ion/ligand names from a DataFrame column."""
    items = []
    for val in df[col]:
        try:
            group = ast.literal_eval(val) if isinstance(val, str) else val
            for x in group:
                if isinstance(x, list) and len(x) > 1:
                    items.append(x[1])
        except Exception:
            continue
    return items

def get_samples_with_item(df, col, item):
    """Return a boolean mask for rows where the item is present in the given column."""
    mask = []
    for val in df[col]:
        found = False
        try:
            group = ast.literal_eval(val) if isinstance(val, str) else val
            for x in group:
                if isinstance(x, list) and len(x) > 1 and x[1] == item:
                    found = True
                    break
        except Exception:
            pass
        mask.append(found)
    return pd.Series(mask, index=df.index)

def classify_ions_ligands(ions, ligands):
    import ast

    try:
        ions = ast.literal_eval(ions) if isinstance(ions, str) else ions
        ligands = ast.literal_eval(ligands) if isinstance(ligands, str) else ligands
    except Exception:
        return 'unknown'

    # Case: both empty
    if not ions and not ligands:
        return 'none'

    # Check if any sublist starts with 0 → disallowed
    def is_disallowed(group):
        return any(isinstance(x, list) and x and x[0] == 0 for x in group)

    if is_disallowed(ions) or is_disallowed(ligands):
        return 'unsupported'

    # Check if all non-empty sublists start with 1 → allowed
    def is_allowed(group):
        return all(isinstance(x, list) and x and x[0] == 1 for x in group)

    if is_allowed(ions) and is_allowed(ligands) and (ions or ligands):
        return 'supported'

    return 'other'

def classify_presence(ions, ligands):
    import ast
    try:
        ions_list = ast.literal_eval(ions) if isinstance(ions, str) else ions
        ligands_list = ast.literal_eval(ligands) if isinstance(ligands, str) else ligands
    except Exception:
        return 'unknown'

    has_ions = bool(ions_list)
    has_ligands = bool(ligands_list)

    if not has_ions and not has_ligands:
        return 'none'
    elif has_ions and not has_ligands:
        return 'ions'
    elif not has_ions and has_ligands:
        return 'ligand'
    elif has_ions and has_ligands:
        return 'ions/ligands'
    
    return 'unknown'

def is_empty_list(val):
    try:
        group = ast.literal_eval(val) if isinstance(val, str) else val
        return isinstance(group, list) and len(group) == 0
    except Exception:
        return False

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

def main():
    config = StartConfig()
    db = DatabaseMethods()
    pred_table = config.pred_table
    query = f"""
        SELECT 
            p.exp_db_id,
            p.ContactList as pred_contacts,
            p.NumberRNAs,
            p.RNALength,
            p.RNAChainIDs,
            p.ChainIDpairList_proteinRNA,
            p.RNAMotifLength,
            e.ContactList as exp_contacts,
            e.Ions,
            e.Ligands,
            p.af3_global_pae_avg,
            p.af3_chain_pair_pae_min,
            p.Complex_RMSD,
            p.RNA_RMSD,
            p.RNA_TM,
            p.RNA_LDDT,
            p.RNAmotif_score,
            p.Complex_TM,
            p.AAmatch_score,
            p.RNA_DI,
            p.Complex_LDDT
        FROM {pred_table} p
        LEFT JOIN exp_protein_rna e ON p.exp_db_id = e.PDBId
        WHERE p.ContactList IS NOT NULL AND e.ContactList IS NOT NULL
    """
    df = pd.read_sql_query(query, db.connection)
    db.close_connection()

    df['category'] = df.apply(lambda row: classify_ions_ligands(row['Ions'], row['Ligands']), axis=1)

    metrics = [
        ('af3_global_pae_avg', 'Global PAE'),
        ('af3_chain_pair_pae_min', 'Chain Pair PAE'),
        ('Complex_RMSD', 'Complex RMSD'),
        ('RNA_RMSD', 'RNA RMSD'),
        ('RNA_TM', 'RNA TM'),
        ('RNA_LDDT', 'RNA LDDT'),
        ('RNAmotif_score', 'RNA Motif Score'),
        ('Complex_TM', 'Complex TM'),
        ('Complex_LDDT', 'Complex LDDT'),
        ('AAmatch_score', 'AA Match Score'),
        ('RNA_DI', 'RNA_DI')
    ]

    plot_creator = PlotCreator('ions_ligands_analysis', None, False)

    # ions = "[[0, 'SO4'], [1, 'CL']]"
    # ligands = "[]"
    # print(classify_ions_ligands(ions, ligands))

    for metric, label in metrics:
        data = []
        labels = []
        for cat in ['none', 'supported', 'unsupported']:
            print(f"Category: {cat}, count: {len(df[df['category'] == cat])}")

        for cat in ['none', 'supported', 'unsupported']:
            if metric in ['Complex_RMSD', 'RNA_RMSD']:
                vals = df[df['category'] == cat][metric].apply(parse_rna_rmsd).dropna()
                vals = pd.to_numeric(vals, errors='coerce').dropna()
            else:
                vals = pd.to_numeric(df[df['category'] == cat][metric], errors='coerce').dropna()
            if len(vals) > 0:
                data.append(vals)
                labels.append(cat)
        if data:
            plot_creator.get_grouped_violin_plot(
                table_source='ions_ligands',
                data_values=data,
                labels=labels,
                name=f'{metric}_by_ions_ligands',
                yAxis_label=label
            )

    # --- NEW PLOTS BY PRESENCE CATEGORY ---
    df['presence_category'] = df.apply(lambda row: classify_presence(row['Ions'], row['Ligands']), axis=1)

    for metric, label in metrics:
        data = []
        labels = []
        for cat in ['none', 'ions', 'ligand', 'ions/ligands']:
            category_df = df[df['presence_category'] == cat]
            if category_df.empty:
                continue

            if metric in ['Complex_RMSD', 'RNA_RMSD']:
                vals = category_df[metric].apply(parse_rna_rmsd).dropna()
            else:
                vals = pd.to_numeric(category_df[metric], errors='coerce').dropna()
            
            if not vals.empty:
                data.append(vals)
                labels.append(cat)
        
        if data:
            plot_creator.get_grouped_violin_plot(
                table_source='ions_ligands_presence',
                data_values=data,
                labels=labels,
                name=f'{metric}_by_presence',
                yAxis_label=label
            )

    # --- IONS ---
    ion_names = extract_items(df, 'Ions')
    ion_counts = Counter(ion_names)
    frequent_ions = [ion for ion, count in ion_counts.items() if count > 2]

    for metric, label in metrics:
        data = []
        labels_ = []
        # Frequent ions
        for ion in frequent_ions:
            mask = get_samples_with_item(df, 'Ions', ion)
            vals = df.loc[mask, metric]
            if metric in ['Complex_RMSD', 'RNA_RMSD']:
                vals = vals.apply(parse_rna_rmsd)
            vals = pd.to_numeric(vals, errors='coerce').dropna()
            if len(vals) > 0:
                data.append(vals)
                labels_.append(ion)
        # Add group for samples with no ions
        mask_no_ions = df['Ions'].apply(is_empty_list)
        vals = df.loc[mask_no_ions, metric]
        if metric in ['Complex_RMSD', 'RNA_RMSD']:
            vals = vals.apply(parse_rna_rmsd)
        vals = pd.to_numeric(vals, errors='coerce').dropna()
        if len(vals) > 0:
            data.append(vals)
            labels_.append('none')
        if data:
            plot_creator.get_grouped_violin_plot(
                table_source='ions_ligands',
                data_values=data,
                labels=labels_,
                name=f'{metric}_by_frequent_ions',
                yAxis_label=label
            )

    # --- LIGANDS ---
    ligand_names = extract_items(df, 'Ligands')
    ligand_counts = Counter(ligand_names)
    frequent_ligands = [lig for lig, count in ligand_counts.items() if count > 4]

    for metric, label in metrics:
        data = []
        labels_ = []
        # Frequent ligands
        for lig in frequent_ligands:
            mask = get_samples_with_item(df, 'Ligands', lig)
            vals = df.loc[mask, metric]
            if metric in ['Complex_RMSD', 'RNA_RMSD']:
                vals = vals.apply(parse_rna_rmsd)
            vals = pd.to_numeric(vals, errors='coerce').dropna()
            if len(vals) > 0:
                data.append(vals)
                labels_.append(lig)
        # Add group for samples with no ligands
        mask_no_ligands = df['Ligands'].apply(is_empty_list)
        vals = df.loc[mask_no_ligands, metric]
        if metric in ['Complex_RMSD', 'RNA_RMSD']:
            vals = vals.apply(parse_rna_rmsd)
        vals = pd.to_numeric(vals, errors='coerce').dropna()
        if len(vals) > 0:
            data.append(vals)
            labels_.append('none')
        if data:
            plot_creator.get_grouped_violin_plot(
                table_source='ions_ligands',
                data_values=data,
                labels=labels_,
                name=f'{metric}_by_frequent_ligands',
                yAxis_label=label,
                rotation=90
            )

if __name__ == '__main__':
    main() 