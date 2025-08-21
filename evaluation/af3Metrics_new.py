import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import glob
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from structures.proteinRNAdnaComplex import RNAProteinDNAcomplex
import numpy as np
import ast
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig
from structures.rna import RNA
from structures.dna import DNA
from structures.protein import Protein

database_methods = DatabaseMethods()
startConfig = StartConfig()
sourceTable = startConfig.pred_table
iLDDTs = []
rna_iLDDTs = []
dna_iLDDTs = []


# Based on: https://colab.research.google.com/drive/1YB2378jqsZMNkVPonF1Mxrp-jn188w16?usp=sharing&pli=1#scrollTo=YZAVOqx6dig5
def extract_data(json_file_path):
    try:
        with open(json_file_path, 'r') as file:
            return json.load(file)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
        return None


# Convert a number to Excel-style column label: A, B, C, ..., Z, AA, AB, ..., AZ, BA, etc
def number_to_column_label(num):
    label = ''
    while num > 0:
        num, remainder = divmod(num - 1, 26)
        label = chr(65 + remainder) + label
    return label


def generate_and_save_heatmap(data, title, xlabel, ylabel, cbar_label, output_file_path):
    plt.figure(figsize=(12, 8))
    sns.heatmap(data, cmap='coolwarm', cbar_kws={'label': cbar_label}, square=True)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(output_file_path)
    plt.close()


# pLDDTs for all entities in one figure
def plot_and_save_plddts_combined(df, output_file_path):
    # Group the DataFrame by the 'atom_chain_ids' column into chains
    grouped = df.groupby('atom_chain_ids')
    plt.figure(figsize=(18, 10))
    for chain_id, chain_df in grouped:
        plt.plot(chain_df.index, chain_df['atom_plddts'], label=f'Chain {chain_id}')
    plt.title('pLDDTs per atom for all entities')
    plt.xlabel('Atom id')
    plt.ylabel('Atom pLDDT')
    plt.legend()
    plt.savefig(output_file_path)
    plt.close()


def insert_data_to_db(file_name, iptm_metrics, pae_metrics):
    if sourceTable == 'pred_protein_rna':
        protein_plddt_avg = round(np.mean(pae_metrics['protein_plddt']) / 100.0, 2) if pae_metrics['protein_plddt'] else 0
        rna_plddt_avg = round(np.mean(pae_metrics['rna_plddt']) / 100.0, 2) if pae_metrics['rna_plddt'] else 0
        
        database_methods.update_or_insert('pred_protein_rna',
            ['af3_protein_pTM', 'af3_rna_pTM', 
             'af3_protein_ipTM', 'af3_rna_ipTM',
             'af3_protein_pLDDT_avg', 'af3_rna_pLDDT_avg',
             'af3_global_pae_avg', 'af3_chain_pair_pae_min',
             'af3_fraction_disordered', 'af3_has_clash', 'af3_ranking_score'],
            (json.dumps(iptm_metrics['protein_ptm']),
             json.dumps(iptm_metrics['rna_ptm']),
             json.dumps(iptm_metrics['protein_iptm']),
             json.dumps(iptm_metrics['rna_iptm']),
             protein_plddt_avg,
             rna_plddt_avg,
             pae_metrics['global_pae_avg'],
             json.dumps(pae_metrics['chain_pair_pae_min']),
             iptm_metrics['fraction_disordered'],
             iptm_metrics['has_clash'],
             iptm_metrics['af3_ranking_score']),
                                          condition=f"FileName = '{file_name}'")

    elif sourceTable == 'pred_protein_dna':
        protein_plddt_avg = round(np.mean(pae_metrics['protein_plddt']) / 100.0, 2) if pae_metrics['protein_plddt'] else 0
        dna_plddt_avg = round(np.mean(pae_metrics['dna_plddt']) / 100.0, 2) if pae_metrics['dna_plddt'] else 0
        
        database_methods.update_or_insert('pred_protein_dna',
            ['af3_protein_pTM', 'af3_dna_pTM', 
             'af3_protein_ipTM', 'af3_dna_ipTM',
             'af3_protein_pLDDT_avg', 'af3_dna_pLDDT_avg',
             'af3_global_pae_avg', 'af3_chain_pair_pae_min',
             'af3_fraction_disordered', 'af3_has_clash', 'af3_ranking_score'],
            (json.dumps(iptm_metrics['protein_ptm']),
             json.dumps(iptm_metrics['dna_ptm']),
             json.dumps(iptm_metrics['protein_iptm']),
             json.dumps(iptm_metrics['dna_iptm']),
             protein_plddt_avg,
             dna_plddt_avg,
             pae_metrics['global_pae_avg'],
             json.dumps(pae_metrics['chain_pair_pae_min']),
             iptm_metrics['fraction_disordered'],
             iptm_metrics['has_clash'],
             iptm_metrics['af3_ranking_score']),
            condition=f"FileName = '{file_name}'")

    elif sourceTable == 'pred_rna_rna':
        rna_plddt_avg = round(np.mean(pae_metrics['rna_plddt']) / 100.0, 2) if pae_metrics['rna_plddt'] else 0
        
        database_methods.update_or_insert('pred_rna_rna',
            ['af3_rna_pTM', 'af3_rna_ipTM', 'af3_rna_pLDDT_avg',
             'af3_global_pae_avg', 'af3_chain_pair_pae_min',
             'af3_fraction_disordered', 'af3_has_clash', 'af3_ranking_score'],
            (json.dumps(iptm_metrics['rna_ptm']),
             json.dumps(iptm_metrics['rna_iptm']),
             rna_plddt_avg,
             pae_metrics['global_pae_avg'],
             json.dumps(pae_metrics['chain_pair_pae_min']),
             iptm_metrics['fraction_disordered'],
             iptm_metrics['has_clash'],
             iptm_metrics['af3_ranking_score']),
                                          condition=f"FileName = '{file_name}'")

    elif sourceTable == 'pred_protein_rna_dna':
        protein_plddt_avg = round(np.mean(pae_metrics['protein_plddt']) / 100.0, 2) if pae_metrics['protein_plddt'] else 0
        rna_plddt_avg = round(np.mean(pae_metrics['rna_plddt']) / 100.0, 2) if pae_metrics['rna_plddt'] else 0
        dna_plddt_avg = round(np.mean(pae_metrics['dna_plddt']) / 100.0, 2) if pae_metrics['dna_plddt'] else 0
        
        database_methods.update_or_insert('pred_protein_rna_dna',
            ['af3_protein_pTM', 'af3_rna_pTM', 'af3_dna_pTM',
             'af3_protein_ipTM', 'af3_rna_ipTM', 'af3_dna_ipTM',
             'af3_protein_pLDDT_avg', 'af3_rna_pLDDT_avg', 'af3_dna_pLDDT_avg',
             'af3_global_pae_avg', 'af3_chain_pair_pae_min',
             'af3_fraction_disordered', 'af3_has_clash', 'af3_ranking_score'],
            (json.dumps(iptm_metrics['protein_ptm']),
             json.dumps(iptm_metrics['rna_ptm']),
             json.dumps(iptm_metrics['dna_ptm']),
             json.dumps(iptm_metrics['protein_iptm']),
             json.dumps(iptm_metrics['rna_iptm']),
             json.dumps(iptm_metrics['dna_iptm']),
             protein_plddt_avg,
             rna_plddt_avg,
             dna_plddt_avg,
             pae_metrics['global_pae_avg'],
             json.dumps(pae_metrics['chain_pair_pae_min']),
             iptm_metrics['fraction_disordered'],
             iptm_metrics['has_clash'],
             iptm_metrics['af3_ranking_score']),
                                          condition=f"FileName = '{file_name}'")


def get_interfaceAtomID_from_db(file_path):
    def query_rna_sequence(file_name):
        parts = file_name.split('_')
        proteinRNA = RNAProteinDNAcomplex.get_proteinRNAdnaComplex_from_db_predTable(id=parts[1].upper(),
                                                                                     file_name=file_name)
        all_interface_id = proteinRNA.get_interfaceAtom_ids()
        rna_interface_id = proteinRNA.get_RNA_interfaceAtom_ids()
        if sourceTable == 'pred_protein_rna_dna':
            protDNA_interface_id = proteinRNA.get_dnaProt_interfaceAtom_ids()  # check if protDNA_interface_id necessary
            dna_interface_id = proteinRNA.get_DNA_interfaceAtom_ids()
        else:
            protDNA_interface_id = None
            dna_interface_id = None

        return all_interface_id, rna_interface_id, protDNA_interface_id, dna_interface_id

    file_name = os.path.basename(file_path)
    cif_file_name1 = file_name.replace('_full_data_', '_model_').replace('.json', '.cif')
    file = os.path.splitext(file_name)[0]
    all_interfaceAtomID, rna_interfaceAtomID, protDNA_interface_id, dna_interface_id = query_rna_sequence(
        cif_file_name1)
    # print("interfaceAtomID", interfaceAtomID)

    if not all_interfaceAtomID or not rna_interfaceAtomID:  # TO-DO: use cifParser.py
        print("Parse .cif file first before extracting interactions.")

    return all_interfaceAtomID, rna_interfaceAtomID, protDNA_interface_id, dna_interface_id


def derive_cif_file_name(json_file_name):
    if '_summary_confidences_' in json_file_name:
        return json_file_name.replace('_summary_confidences_', '_model_').replace('.json', '.cif')
    elif '_full_data_' in json_file_name:
        return json_file_name.replace('_full_data_', '_model_').replace('.json', '.cif')
    else:
        raise ValueError(
            f"{json_file_name} does not contain a recognized pattern: _summary_confidences_ or _full_data_.")


def plot_violin(iLDDTs, folder_path, title):
    plt.figure()
    sns.violinplot(data=iLDDTs, inner="point")
    plt.ylabel(title, fontsize=18)
    title = title.replace(' ', '_').replace('/', '_')
    save_path = os.path.join(folder_path, f"{title}.png")
    plt.savefig(save_path)
    plt.close()


def identify_chain_type(chain_id, file_name):
    """
    Identify if a chain is RNA, DNA, or protein based on the structures package.

    Args:
        chain_id: The chain identifier
        file_name: The CIF file name

    Returns:
        str: 'rna', 'dna', or 'protein'
    """
    parts = file_name.split('_')
    pdb_id = parts[1].upper()

    # Check RNA first
    rna = RNA.get_rna_from_db(id=pdb_id, file_name=file_name)
    if rna and chain_id in rna.get_rna_chain_IDs():
        return 'rna'

    # Check DNA
    dna = DNA.get_dna_from_db(id=pdb_id, file_name=file_name)
    if dna and chain_id in dna.get_dna_chain_IDs():
        return 'dna'

    # Check Protein
    protein = Protein.get_protein_from_db(id=pdb_id, file_name=file_name)
    if protein and chain_id in protein.get_protein_chain_IDs():
        return 'protein'

    return 'unknown'


def get_af3Metrics_paths_iptm(folder_path, json_file_path):
    data = extract_data(json_file_path)
    if data:
        # Extract basic metrics
        chain_ptm = data.get('chain_ptm', [])
        chain_iptm = data.get('chain_iptm', [])
        chain_pair_iptm = data.get('chain_pair_iptm', [])
        chain_pair_pae_min_matrix = data.get('chain_pair_pae_min', [])
        # print("chain_pair_pae_min_matrix", chain_pair_pae_min_matrix)

        # Basic metrics
        fraction_disordered = data.get('fraction_disordered', 0)
        has_clash = data.get('has_clash', False)
        iptm = data.get('iptm', 0)
        af3_ranking_score = data.get('ranking_score', 0)

        # Create matrices
        df_iptm_matrix = pd.DataFrame(chain_pair_iptm)
        df_pae_min_matrix = pd.DataFrame(chain_pair_pae_min_matrix) #for visualisation only
        # print("df_pae_min_matrix", df_pae_min_matrix)

        # Set labels for both matrices
        chain_labels = [number_to_column_label(i) for i in range(1, len(df_iptm_matrix.index) + 1)]
        df_iptm_matrix.index = df_iptm_matrix.columns = chain_labels
        df_pae_min_matrix.index = df_pae_min_matrix.columns = chain_labels

        # Get CIF file name for chain type identification
        cif_file_name = derive_cif_file_name(os.path.basename(json_file_path))

        # Identify chain types and separate metrics
        protein_ptm = []
        rna_ptm = []
        dna_ptm = []
        chain_types = {}

        for i, chain in enumerate(chain_labels):
            chain_type = identify_chain_type(chain, cif_file_name)
            chain_types[chain] = chain_type

            if chain_type == 'protein':
                protein_ptm.append(chain_ptm[i])
            elif chain_type == 'rna':
                rna_ptm.append(chain_ptm[i])
            elif chain_type == 'dna':
                dna_ptm.append(chain_ptm[i])

        # Calculate average PTM
        ptm_avg = np.mean(chain_ptm) if chain_ptm else 0

        # Calculate average ipTM per chain type
        protein_iptm = []
        rna_iptm = []
        dna_iptm = []

        for i, chain in enumerate(chain_labels):
            chain_type = identify_chain_type(chain, cif_file_name)
            if chain_type == 'protein':
                protein_iptm.append(chain_iptm[i])
            elif chain_type == 'rna':
                rna_iptm.append(chain_iptm[i])
            elif chain_type == 'dna':
                dna_iptm.append(chain_iptm[i])

        metrics = {
            # Single values
            "ptm_avg": round(ptm_avg, 2),
            "iptm": round(iptm, 2),
            "protein_iptm_avg": round(np.mean(protein_iptm), 2) if protein_iptm else 0,
            "rna_iptm_avg": round(np.mean(rna_iptm), 2) if rna_iptm else 0,
            "dna_iptm_avg": round(np.mean(dna_iptm), 2) if dna_iptm else 0,

            # Lists
            "protein_ptm": [round(x, 2) for x in protein_ptm],
            "rna_ptm": [round(x, 2) for x in rna_ptm],
            "dna_ptm": [round(x, 2) for x in dna_ptm],
            "protein_iptm": [round(x, 2) for x in protein_iptm],
            "rna_iptm": [round(x, 2) for x in rna_iptm],
            "dna_iptm": [round(x, 2) for x in dna_iptm],
            "chain_ptm": [round(x, 2) for x in chain_ptm],
            "chain_iptm": [round(x, 2) for x in chain_iptm],
            "chain_pair_pae_min": chain_pair_pae_min_matrix,

            # Other metrics
            "fraction_disordered": round(fraction_disordered, 2),
            "has_clash": has_clash,
            "af3_ranking_score": round(af3_ranking_score, 2),

            # For visualization
            "df_iptm_matrix": df_iptm_matrix,
            "df_pae_min_matrix": df_pae_min_matrix,
            "chain_labels": chain_labels,
            "chain_types": chain_types
        }

        return metrics

    return None


def get_pae_value_for_chain_pair(chain_pair, pae_matrix):
    """Get PAE value from matrix for specific chain pair"""
    try:
        if not chain_pair or len(chain_pair) != 2:
            return None
            
        chain1, chain2 = chain_pair
        
        # Convert matrix to numpy array if it isn't already
        pae_array = np.array(pae_matrix)
        
        # If array is 1D, reshape it to 2D square matrix
        if len(pae_array.shape) == 1:
            n = int(np.sqrt(len(pae_array)))  # Calculate matrix size
            pae_array = pae_array.reshape(n, n)
        
        # Get indices based on alphabetical order (A=0, B=1, etc.)
        chain1_idx = ord(chain1) - ord('A')
        chain2_idx = ord(chain2) - ord('A')
        
        # Check if indices are within matrix bounds
        if chain1_idx < 0 or chain2_idx < 0 or chain1_idx >= len(pae_array) or chain2_idx >= len(pae_array[0]):
            print(f"Chain indices out of bounds: {chain1}({chain1_idx}), {chain2}({chain2_idx})")
            return None
        
        # Get PAE value from matrix
        return float(pae_array[chain1_idx, chain2_idx])
        
    except (IndexError, TypeError, ValueError) as e:
        print(f"Error getting PAE value for chain pair {chain_pair}: {e}")
        return None


def get_af3Metrics_paths_pae(cif_name, json_file_path, iptm_metrics):
    """Process PAE metrics from JSON file"""
    # try:
    data = extract_data(json_file_path)
    if not data:
        return None

    atom_chain_ids = data.get('atom_chain_ids', [])
    atom_plddts = data.get('atom_plddts', [])
    contact_probs = data.get('contact_probs', [])
    pae_data = data.get('pae', [])

    # Get chain_pair_pae_min from iptm_metrics instead of full_data JSON
    chain_pair_pae_min = iptm_metrics.get('chain_pair_pae_min')

    # Initialize metrics dictionary
    metrics = {}

    # Calculate per-chain pLDDT statistics using chain types from iptm_metrics
    df_plddts = pd.DataFrame({'atom_chain_ids': atom_chain_ids, 'atom_plddts': atom_plddts})

    protein_plddts = []
    rna_plddts = []
    dna_plddts = []
    chain_plddt_stats = {}

    # Calculate global PAE average with rounding
    global_pae_avg = 0
    if pae_data and len(pae_data) > 0:
        all_pae_values = [val for row in pae_data for val in row]
        global_pae_avg = round(np.mean(all_pae_values), 2) if all_pae_values else 0

    chain_types = iptm_metrics['chain_types']

    for chain, chain_df in df_plddts.groupby('atom_chain_ids'):
        plddts = chain_df['atom_plddts'].tolist()
        chain_type = chain_types.get(chain, 'unknown')

        # Calculate pLDDT stats for this chain
        chain_plddt_stats[chain] = {
            'avg_plddt': np.mean(plddts) if plddts else 0,
            'min_plddt': min(plddts) if plddts else 0,
            'max_plddt': max(plddts) if plddts else 0
        }

        if chain_type == 'protein':
            protein_plddts.extend(plddts)
        elif chain_type == 'rna':
            rna_plddts.extend(plddts)
        elif chain_type == 'dna':
            dna_plddts.extend(plddts)

    # Get chain pair with longest motif
    pdb_id = os.path.basename(json_file_path).split('_')[1]
    rna = RNA.get_rna_from_db(id=pdb_id, file_name=cif_name) #!don't forget the filename for pred_protein_rna
    protein = Protein.get_protein_from_db(id=pdb_id, file_name=cif_name)
    if rna.get_rna_chain_number() == 1 and protein.get_protein_chain_number() == 1:
        # print("protein.protein_chain_IDs", protein.protein_chain_IDs)
        # print("get_rna_chain_IDs", rna.get_rna_chain_IDs())
        chain_pair = [f"{protein.protein_chain_IDs}", f"{rna.get_rna_chain_IDs()}"]
        # print("%%%%chain_pair", chain_pair)
        chain_pair_pae = get_pae_value_for_chain_pair(chain_pair, chain_pair_pae_min)
    else:
        chain_pair = rna.get_chain_pair_with_longest_motif()
        # print("%%%%chain_pair", chain_pair)
        chain_pair_pae = get_pae_value_for_chain_pair(chain_pair, chain_pair_pae_min)
        # print("£££££££££chain_pair_pae", chain_pair_pae)

    final_metrics = {
        # Lists
        "protein_plddt": [round(x, 2) for x in protein_plddts],
        "rna_plddt": [round(x, 2) for x in rna_plddts],
        "dna_plddt": [round(x, 2) for x in dna_plddts],

        # Chain-specific stats
        "chain_plddt_stats": {chain: {
            'avg_plddt': round(stats['avg_plddt'], 2),
            'min_plddt': round(stats['min_plddt'], 2),
            'max_plddt': round(stats['max_plddt'], 2)
        } for chain, stats in chain_plddt_stats.items()},
        "global_pae_avg": global_pae_avg,
        "chain_pair_pae_min": chain_pair_pae,

        # For visualization
        "contact_probs": contact_probs,
        "df_plddts": df_plddts,
        "pae_data": pae_data
    }

    return final_metrics

    # except Exception as e:
    #     print(f"Error processing PAE metrics: {e}")
    #     return None


def process_interface_metrics(interface_atom_ids_str, rna_interfaceAtomID_str,
                              protDNA_interface_id_str, dna_interface_id_str, atom_plddts):
    metrics = {
        "interface_plddt_avg": 0,
        "rna_interface_plddt_avg": 0,
        "dna_interface_plddt_avg": 0
    }

    if interface_atom_ids_str:
        # Process interface atoms
        interface_atom_ids = process_interface_ids(interface_atom_ids_str)
        metrics["interface_plddt_avg"] = calculate_plddt_avg(interface_atom_ids, atom_plddts)

        # Process RNA interface atoms
        if rna_interfaceAtomID_str:
            rna_interface_ids = process_interface_ids(rna_interfaceAtomID_str)
            metrics["rna_interface_plddt_avg"] = calculate_plddt_avg(rna_interface_ids, atom_plddts)

        # Process DNA interface atoms if available
        if protDNA_interface_id_str and dna_interface_id_str:
            dna_interface_ids = process_interface_ids(dna_interface_id_str)
            metrics["dna_interface_plddt_avg"] = calculate_plddt_avg(dna_interface_ids, atom_plddts)

    return metrics


def process_interface_ids(interface_ids_str):
    ids = ast.literal_eval(interface_ids_str)
    if any(isinstance(i, list) for i in ids):
        return [item for sublist in ids for item in sublist]
    return ids


def calculate_plddt_avg(atom_ids, atom_plddts):
    plddts = [atom_plddts[i] for i in atom_ids if 0 <= i < len(atom_plddts)]
    return np.mean(plddts) if plddts else 0


def create_af3_visualizations(json_file_path, iptm_metrics, pae_metrics, folder_path):
    fig = plt.figure(figsize=(30, 20))
    gs = GridSpec(3, 2, figure=fig)

    # Commented out plots
    # create_iptm_heatmap(fig, gs[0, 0], iptm_metrics['df_iptm_matrix'])
    # create_contact_heatmap(fig, gs[0, 1], pae_metrics['contact_probs'])

    # Active plots
    create_pae_min_heatmap(fig, gs[1, 0], iptm_metrics['df_pae_min_matrix'])
    create_pae_heatmap(fig, gs[1, 1], pae_metrics['pae_data'])
    create_plddt_plot(fig, gs[2, 0], pae_metrics['df_plddts'])
    create_stats_plot(fig, gs[2, 1], iptm_metrics, pae_metrics)

    # Save the visualization
    save_visualization(fig, folder_path, json_file_path)


def create_iptm_heatmap(fig, position, df_iptm_matrix):
    ax = fig.add_subplot(position)
    sns.heatmap(df_iptm_matrix, cmap='coolwarm', ax=ax)
    ax.set_title('ipTM Matrix Heatmap', fontsize=28)


def create_pae_min_heatmap(fig, position, df_pae_min_matrix):
    ax = fig.add_subplot(position)
    sns.heatmap(df_pae_min_matrix, cmap='viridis', ax=ax)
    ax.set_title('Chain-pair Minimum PAE Heatmap', fontsize=28)


def create_pae_heatmap(fig, position, pae_data):
    ax = fig.add_subplot(position)
    sns.heatmap(pae_data, cmap='viridis', ax=ax)
    ax.set_title('Full PAE Matrix', fontsize=28)
    ax.set_xlabel('Residue', fontsize=18)
    ax.set_ylabel('Residue', fontsize=18)


def create_plddt_plot(fig, position, df_plddts):
    ax = fig.add_subplot(position)
    grouped = df_plddts.groupby('atom_chain_ids')
    for chain_id, chain_df in grouped:
        ax.plot(chain_df.index, chain_df['atom_plddts'], label=f'Chain {chain_id}')
    ax.set_title('pLDDTs per Chain', fontsize=28)
    ax.set_xlabel('Atom Index', fontsize=18)
    ax.set_ylabel('pLDDT', fontsize=18)
    ax.legend()


def create_contact_heatmap(fig, position, contact_probs):
    ax = fig.add_subplot(position)
    sns.heatmap(contact_probs, cmap='viridis', ax=ax)
    ax.set_title('Contact Probability Heatmap', fontsize=28)
    ax.set_xlabel('Residue', fontsize=18)
    ax.set_ylabel('Residue', fontsize=18)


def create_stats_plot(fig, position, iptm_metrics, pae_metrics):
    ax = fig.add_subplot(position)

    # Prepare data for plotting
    chains = iptm_metrics['chain_labels']
    data = []

    for i, chain in enumerate(chains):
        plddt = pae_metrics['chain_plddt_stats'][chain]['avg_plddt'] / 100.0
        data.append({
            'Chain': chain,
            'pTM': iptm_metrics['chain_ptm'][i],
            'ipTM': iptm_metrics['chain_iptm'][i],
            'pLDDT': plddt
        })

    # Create DataFrame
    df = pd.DataFrame(data)

    # Plot chain-specific metrics
    df.plot(x='Chain', kind='bar', ax=ax, width=0.8)

    # Add ranking score as text - moved to upper left to avoid legend overlap
    af3_ranking_score = iptm_metrics['af3_ranking_score']
    ax.text(0.05, 0.95, f'Ranking Score: {af3_ranking_score:.2f}',
            transform=ax.transAxes, ha='left', va='top',
            bbox=dict(facecolor='white', alpha=0.8))

    ax.set_title('Chain-wise Metrics', fontsize=28)
    ax.set_xlabel('Chain', fontsize=18)
    ax.set_ylabel('Score', fontsize=18)
    ax.legend(fontsize=12, loc='upper right')  # Explicitly set legend to upper right
    ax.grid(True, alpha=0.3)

    # Adjust y-axis to show full range
    ax.set_ylim(0, 1.0)

    # Rotate x-axis labels for better readability
    plt.setp(ax.get_xticklabels(), rotation=45)


def save_visualization(fig, folder_path, json_file_path):
    pdb_id = os.path.basename(json_file_path).split('_')[1]
    file_base_graphAF3Metrics = f"fold_{pdb_id}_graphAF3Metrics.png"
    save_path = os.path.join(folder_path, file_base_graphAF3Metrics)
    plt.subplots_adjust(hspace=0.3, wspace=0.3)  # Added wspace adjustment
    plt.savefig(save_path, bbox_inches='tight', dpi=300)  # Added bbox_inches='tight' for better spacing
    plt.close(fig)


def print_metrics(iptm_metrics, pae_metrics):
    print("\n=== AF3 Metrics Summary ===")
    print(f"Global PAE Average: {pae_metrics['global_pae_avg']:.2f}")
    # print("Chain-pair PAE min (row with highest sum):", [f"{x:.2f}" for x in iptm_metrics['chain_pair_pae_min']])
    print(f"Chain-pair PAE min: {pae_metrics['chain_pair_pae_min']}")
    print(f"AF3 Ranking Score: {iptm_metrics['af3_ranking_score']:.2f}")
    print(f"Fraction Disordered: {iptm_metrics['fraction_disordered']:.2f}")
    print(f"Has Clash: {iptm_metrics['has_clash']}")

    print("\n--- Protein Metrics ---")
    if iptm_metrics['protein_ptm']:
        print("Protein pTM values:", [f"{x:.2f}" for x in iptm_metrics['protein_ptm']])
        print("Protein ipTM values:", [f"{x:.2f}" for x in iptm_metrics['protein_iptm']])
        if pae_metrics['protein_plddt']:
            protein_plddt_avg = np.mean(pae_metrics['protein_plddt']) / 100.0
            print(f"Protein pLDDT average: {protein_plddt_avg:.2f}")
    else:
        print("No protein chains found")

    print("\n--- RNA Metrics ---")
    if iptm_metrics['rna_ptm']:
        print("RNA pTM values:", [f"{x:.2f}" for x in iptm_metrics['rna_ptm']])
        print("RNA ipTM values:", [f"{x:.2f}" for x in iptm_metrics['rna_iptm']])
        if pae_metrics['rna_plddt']:
            rna_plddt_avg = np.mean(pae_metrics['rna_plddt']) / 100.0
            print(f"RNA pLDDT average: {rna_plddt_avg:.2f}")
    else:
        print("No RNA chains found")

    print("\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python af3Metrics_new.py <af_folder_path> [--no-plots]")
        print("Options:")
        print("  --no-plots : Skip generating plots")
        sys.exit(1)

    folder_path = sys.argv[1]
    json_file_paths_iptm = glob.glob(os.path.join(folder_path, '*_summary_confidences_[0-4].json'))
    json_file_paths_pae = glob.glob(os.path.join(folder_path, '*_full_data_[0-4].json'))
    skip_plots = "--no-plots" in sys.argv

    for json_file_path in json_file_paths_iptm:
        # Process iPTM metrics
        iptm_metrics = get_af3Metrics_paths_iptm(folder_path, json_file_path)
        if not iptm_metrics:
            continue

        # Find matching PAE file
        cif_name = derive_cif_file_name(os.path.basename(json_file_path))
        print("!£$@$cif_name", cif_name)
        matching_pae_files = [f for f in json_file_paths_pae
                              if derive_cif_file_name(os.path.basename(f)) == cif_name]

        if not matching_pae_files:
            continue

        # Process PAE metrics
        pae_metrics = get_af3Metrics_paths_pae(cif_name, matching_pae_files[0], iptm_metrics)
        if not pae_metrics:
            continue

        # Print metrics
        # print(f"\nProcessing: {os.path.basename(json_file_path)}")
        # print_metrics(iptm_metrics, pae_metrics)

        # Create visualizations
        if not skip_plots:
            create_af3_visualizations(json_file_path, iptm_metrics, pae_metrics, folder_path)

        # Update database
        insert_data_to_db(cif_name, iptm_metrics, pae_metrics)






































