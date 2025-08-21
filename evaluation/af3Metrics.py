import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import glob
import sys
from structures.proteinRNAdnaComplex import RNAProteinDNAcomplex
import numpy as np
import ast
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig

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
    heatmap = sns.heatmap(data, cmap='coolwarm', cbar_kws={'label': cbar_label}, square=True)
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

def insert_data_to_db(file_name, ptm, iptm, fraction_disordered, iLDDT, rna_iLDDT, dna_iLDDT):
    if sourceTable == 'pred_protein_rna':
        database_methods.update_or_insert('pred_protein_rna',
                                          ['ptm', 'iptm', 'fraction_disordered', 'iLDDT', 'RNA_iLDDT'],
                                          (ptm, iptm, fraction_disordered, iLDDT, rna_iLDDT),
                                          condition=f"FileName = '{file_name}'")

    elif sourceTable == 'pred_rna_rna':
        database_methods.update_or_insert('pred_rna_rna',
                                          ['ptm', 'iptm', 'fraction_disordered', 'iLDDT', 'RNA_iLDDT'],
                                          (ptm, iptm, fraction_disordered, iLDDT, rna_iLDDT),
                                          condition=f"FileName = '{file_name}'")

    elif sourceTable == 'pred_protein_dna':
        database_methods.update_or_insert('pred_protein_dna',
                                          ['ptm', 'iptm', 'fraction_disordered', 'iLDDT', 'DNA_iLDDT'],
                                          (ptm, iptm, fraction_disordered, iLDDT, rna_iLDDT),
                                          condition=f"FileName = '{file_name}'")

    elif sourceTable == 'pred_protein_rna_dna':
        database_methods.update_or_insert('pred_rna_rna',
                                          ['ptm', 'iptm', 'fraction_disordered', 'iLDDT', 'RNA_iLDDT', 'DNA_iLDDT'],
                                          (ptm, iptm, fraction_disordered, iLDDT, rna_iLDDT, dna_iLDDT),
                                          condition=f"FileName = '{file_name}'")

def get_interfaceAtomID_from_db(file_path):
    def query_rna_sequence(file_name):
        parts = file_name.split('_')
        proteinRNA = RNAProteinDNAcomplex.get_proteinRNAdnaComplex_from_db_predTable(id=parts[1].upper(), file_name=file_name)
        all_interface_id = proteinRNA.get_interfaceAtom_ids()
        rna_interface_id = proteinRNA.get_RNA_interfaceAtom_ids()
        if sourceTable == 'pred_protein_rna_dna':
            protDNA_interface_id = proteinRNA.get_dnaProt_interfaceAtom_ids() # check if protDNA_interface_id necessary
            dna_interface_id = proteinRNA.get_DNA_interfaceAtom_ids()
        else:
            protDNA_interface_id = None
            dna_interface_id = None

        return all_interface_id, rna_interface_id, protDNA_interface_id, dna_interface_id

    file_name = os.path.basename(file_path)
    cif_file_name1 = file_name.replace('_full_data_', '_model_').replace('.json', '.cif')
    file = os.path.splitext(file_name)[0]
    all_interfaceAtomID, rna_interfaceAtomID, protDNA_interface_id, dna_interface_id = query_rna_sequence(cif_file_name1)
    # print("interfaceAtomID", interfaceAtomID)

    if not all_interfaceAtomID or not rna_interfaceAtomID: # TO-DO: use cifParser.py
        print("Parse .cif file first before extracting interactions.")

    return all_interfaceAtomID, rna_interfaceAtomID, protDNA_interface_id, dna_interface_id

def derive_cif_file_name(json_file_name):
    if '_summary_confidences_' in json_file_name:
        return json_file_name.replace('_summary_confidences_', '_model_').replace('.json', '.cif')
    elif '_full_data_' in json_file_name:
        return json_file_name.replace('_full_data_', '_model_').replace('.json', '.cif')
    else:
        raise ValueError(f"{json_file_name} does not contain a recognized pattern: _summary_confidences_ or _full_data_.")

def plot_violin(iLDDTs, folder_path, title):
    plt.figure()
    sns.violinplot(data=iLDDTs, inner="point")
    # plt.title(f'Distribution of {title}', fontsize=18)
    plt.ylabel(title, fontsize=18)
    title = title.replace(' ', '_').replace('/', '_')
    save_path = os.path.join(folder_path, f"{title}.png")
    plt.savefig(save_path)
    plt.close()

# def get_rnaMotif_from_db(cif_file_name):
#     def query_rna_motif(cif_file_name):
#         parts = cif_file_name.split('_')
#         rna = RNA.get_rna_from_db(id=parts[1].upper(), is_pred=True, file_name=cif_file_name)
#         result = rna.get_rna_motif()
#         return result
#
#     # Query the RNA sequence
#     rna_motif = query_rna_motif(cif_file_name)
#
#     if not rna_motif: # TO-DO: use cifParser.py
#         print("Extract interactions with extractInteractions.py.")
#
#     return rna_motif
#
# def create_rna_motif_logo(cif_file_name, rna_interface_plddts):
#     rna_motif_data_str = get_rnaMotif_from_db(cif_file_name)
#     motif_data = []
#     elements = rna_motif_data_str.split(', ') # rna_motif_data_str to list of tuples (nucleotide, position)
#     for element in elements:
#         nucleotide = element[0]
#         position = int(element[2:-1])
#         motif_data.append((nucleotide, position))
#     print(motif_data)
#
#     positions = [pos for _, pos in motif_data] # Initialize the rna_motif_data dictionary with zeroes
#     motif_data = {
#         'Position': positions,
#         'A': [0] * len(positions),
#         'C': [0] * len(positions),
#         'G': [0] * len(positions),
#         'U': [0] * len(positions)
#     }
#     print(motif_data)
#     print(rna_interface_plddts)
#     # Set frequencies based on the motif data and corresponding iLDDT scores
#     for (nucleotide, pos), score in zip(motif_data, rna_interface_plddts):
#         idx = positions.index(pos)
#         rna_motif_data[nucleotide][idx] = score
#
#     print(rna_motif_data)
#
#     # Create a DataFrame
#     logomaker_df = pd.DataFrame(rna_motif_data)
#     print(logomaker_df)
#     # Set the 'Position' column as the index
#     logomaker_df.set_index('Position', inplace=True)
#
#     # Create a Logo object
#     logo = logomaker.Logo(logomaker_df)
#     plt.title("RNA Motif Sequence Logo")
#     plt.xlabel("Position")
#     plt.ylabel("Frequency")
#     title = title.replace(' ', '_').replace('/', '_')
#     save_path = os.path.join(folder_path, f"{title}.png")
#     plt.savefig(save_path)

def get_af3Metrics_paths_iptm(folder_path):
    # json_file_paths_iptm = glob.glob(os.path.join(folder_path, '*_summary_confidences_[0-4].json'))
    #
    # for json_file_path in json_file_paths_iptm:
    data = extract_data(json_file_path)
    if data:
        chain_iptm = data['chain_iptm']
        chain_pair_iptm = data['chain_pair_iptm']
        fraction_disordered = data['fraction_disordered']
        # has_clash = data['has_clash']
        iptm = data['iptm']
        # num_recycles = data['num_recycles']
        ptm = data['ptm']
        # ranking_score = data['ranking_score']
        chain_ptm = data['chain_ptm']
        df_chain_stats = pd.DataFrame({'chain_ptm': chain_ptm, 'chain_iptm': chain_iptm})
        df_chain_stats.index = [number_to_column_label(i) for i in range(1, len(df_chain_stats.index) + 1)]
        df_iptm_matrix = pd.DataFrame(chain_pair_iptm)
        df_iptm_matrix.index = [number_to_column_label(i) for i in range(1, len(df_iptm_matrix.index) + 1)]
        df_iptm_matrix.columns = [number_to_column_label(i) for i in range(1, len(df_iptm_matrix.columns) + 1)]
        # # file_base = os.path.basename(json_file_path).replace('_summary_confidences_', '_iptm_heatmap_')
        # # output_file_path = os.path.join(folder_path, file_base.replace('.json', '.png'))
        # # generate_and_save_heatmap(df_iptm_matrix, 'ipTM Matrix Heatmap', 'Entity id', 'Entity id', 'ipTM', output_file_path)

    return ptm, iptm, fraction_disordered, df_iptm_matrix

def get_af3Metrics_paths_pae(folder_path, df_iptm_matrix):
    # json_file_paths_pae = glob.glob(os.path.join(folder_path, '*_full_data_[0-4].json'))
    # df_iptm_matrix = None
    iLDDT = None

    # for json_file_path in json_file_paths_pae:
    data = extract_data(json_file_path)
    if data:
        pae_data = data['pae']
        # token_res_ids = data['token_res_ids']
        atom_chain_ids = data.get('atom_chain_ids', [])
        atom_plddts = data.get('atom_plddts', [])
        contact_data = data.get('contact_probs', [])
        # file_base = os.path.basename(json_file_path).replace('_full_data_', '_pae_heatmap_')
        # output_file_path = os.path.join(folder_path, file_base.replace('.json', '.png'))
        # generate_and_save_heatmap(pae_data, 'Predicted Alignment Error (PAE) Heatmap', 'Scored residue', 'Aligned residue', 'Expected Position Error (Ångströms)', output_file_path)
        # file_base_plddts = os.path.basename(json_file_path).replace('_full_data_', '_plddts_')
        # if atom_chain_ids and atom_plddts:
        #     df_plddts = pd.DataFrame({'atom_chain_ids': atom_chain_ids, 'atom_plddts': atom_plddts})
        #     output_file_path_plddts = os.path.join(folder_path, file_base_plddts.replace('.json', '.png'))
        #     plot_and_save_plddts_combined(df_plddts, output_file_path_plddts)
        # file_base_contact = os.path.basename(json_file_path).replace('_full_data_', '_contact_probs_heatmap_')
        # if contact_data:
        #     output_file_path_contact = os.path.join(folder_path, file_base_contact.replace('.json', '.png'))
        #     generate_and_save_heatmap(contact_data, 'Contact Probability Heatmap', 'Residue id', 'Residue id', 'Contact Probability', output_file_path_contact)

        fig = plt.figure(figsize=(30, 20))
        gs = GridSpec(3, 2, figure=fig)

        ax1 = fig.add_subplot(gs[0, 0])
        sns.heatmap(df_iptm_matrix, cmap='coolwarm', ax=ax1)
        ax1.set_title('ipTM Matrix Heatmap', fontsize=28)
        ax1.set_xlabel('Entity id', fontsize=18)
        ax1.set_ylabel('Entity id', fontsize=18)

        ax2 = fig.add_subplot(gs[0, 1])
        sns.heatmap(pae_data, cmap='viridis', cbar_kws={'label': 'Expected Position Error (Ångströms)'}, square=True,
                    ax=ax2)
        ax2.set_title('Predicted Alignment Error (PAE) Heatmap', fontsize=28)
        ax2.set_xlabel('Scored residue', fontsize=18)
        ax2.set_ylabel('Aligned residue', fontsize=18)

        ax3 = fig.add_subplot(gs[1, 1])
        sns.heatmap(contact_data, cmap='viridis', cbar_kws={'label': 'Contact Probability'}, square=True, ax=ax3)
        ax3.set_title('Contact Probability Heatmap', fontsize=28)
        ax3.set_xlabel('Residue id', fontsize=18)
        ax3.set_ylabel('Residue id', fontsize=18)

        df_plddts = pd.DataFrame({'atom_chain_ids': atom_chain_ids, 'atom_plddts': atom_plddts})
        grouped = df_plddts.groupby('atom_chain_ids')
        subdfs = {}
        for i, (chain, chain_df) in enumerate(grouped):
            letters = ''
            j = i
            while j >= 0:
                letters = chr(65 + j % 26) + letters
                j = j // 26 - 1
            subdfs[letters] = chain_df
        sorted_keys = sorted(subdfs.keys(), key=lambda x: (len(x), x))
        ax4 = fig.add_subplot(gs[1, 0])
        for key in sorted_keys:
            chain_df = subdfs[key]
            ax4.plot(chain_df.index, chain_df['atom_plddts'], label=f'Entity {key}')
        ax4.set_title('pLDDTs per atom', fontsize=28)
        ax4.set_xlabel('Atom id', fontsize=18)
        ax4.set_ylabel('Atom pLDDT', fontsize=18)
        ax4.legend()

        # Extract and plot interface atom pLDDTs
        cif_file_name = derive_cif_file_name(os.path.basename(json_file_path))
        interface_atom_ids_str, rna_interfaceAtomID_str, protDNA_interface_id_str, dna_interface_id_str = get_interfaceAtomID_from_db(cif_file_name) # rna or dna
        # if sourceTable == 'pred_protein_rna' or sourceTable == 'pred_rna_rna' or sourceTable == 'pred_protein_dna':
        # print("sourceTable", sourceTable)
        # print("rna_interfaceAtomID_str", rna_interfaceAtomID_str)
        # print("interface_atom_ids_str", interface_atom_ids_str)
        if interface_atom_ids_str == "":
            interface_atom_ids = [0]
            interface_plddts = [0]
            iLDDT = 0
            rna_interface_plddts = [0]
            rna_iLDDT = 0
            print("No interface atoms found.")
        else:
            # Transform from string to list
            interface_atom_ids = ast.literal_eval(interface_atom_ids_str)
            rna_interfaceAtomIDs = ast.literal_eval(rna_interfaceAtomID_str)

            # If more than one lists: join lists into a single list
            interface_atom_ids = [item for sublist in interface_atom_ids for item in sublist] if any(
                isinstance(i, list) for i in interface_atom_ids) else interface_atom_ids
            rna_interface_atom_ids = [item for sublist in rna_interfaceAtomIDs for item in sublist] if any(
                isinstance(i, list) for i in rna_interfaceAtomIDs) else rna_interfaceAtomIDs
            # print("rna_interface_atom_ids", rna_interface_atom_ids)

            # Calculate iLDDT and rna_iLDDT
            interface_plddts = [atom_plddts[i] for i in interface_atom_ids if 0 <= i < len(atom_plddts)]
            print("interface_plddts", interface_plddts)
            iLDDT = round(np.mean(interface_plddts), 2)
            print("iLDDT", iLDDT)
            rna_interface_plddts = [atom_plddts[i] for i in rna_interface_atom_ids if 0 <= i < len(atom_plddts)]
            rna_iLDDT = round(np.mean(rna_interface_plddts), 2)
            print("rna_iLDDT", rna_iLDDT)

        if sourceTable == 'pred_protein_rna' or sourceTable == 'pred_rna_rna' or sourceTable == 'pred_protein_dna':
            dna_iLDDT = 0

        # else:
        #     interface_atom_ids_str = interface_atom_ids_str.strip('[]')
        #     print("interface_atom_ids_str", interface_atom_ids_str)
        #     rna_interfaceAtomID_str = rna_interfaceAtomID_str.strip('[]') # rna or dna
        #     print("rna_interfaceAtomID_str", rna_interfaceAtomID_str)

        # if sourceTable == 'pred_protein_rna_dna':
        #     print("ERRORERRORERRORERRORERROR")
        # # Filter out any empty strings before converting them to integers
        #     interface_atom_ids = [int(x) for x in interface_atom_ids_str.split(',') if x.strip()]
        #     rna_interfaceAtomID_str = [int(x) for x in rna_interfaceAtomID_str.split(',') if x.strip()]
        #     if sourceTable == 'pred_protein_rna_dna':
        #         print("protDNA_interface_id_str", protDNA_interface_id_str)
        #         protDNA_interface_id_str = protDNA_interface_id_str.strip('[]')
        #         dna_interface_id_str = dna_interface_id_str.strip('[]')  # rna or dna
        #         print("dna_interface_id_str", dna_interface_id_str)
        #         # Filter out any empty strings before converting them to integers
        #         protDNA_interface_atom_ids = [int(x) for x in protDNA_interface_id_str.split(',') if x.strip()]
        #         dna_interface_id_str = [int(x) for x in dna_interface_id_str.split(',') if x.strip()]
        #     # valid_interface_atom_ids = [i for i in interface_atom_ids if 0 <= i < len(atom_plddts)] # Ensure indices are within the valid range
        #     # print(f"Atoms in total: {len(atom_plddts)}")
        #     # print(f"Interface_atom_ids: {interface_atom_ids}")
        #     if len(interface_atom_ids) > 2:
        #         interface_plddts = [atom_plddts[i] for i in interface_atom_ids if 0 <= i < len(atom_plddts)]
        #         iLDDT = round(np.mean(interface_plddts), 2)
        #         rna_interface_plddts = [atom_plddts[i] for i in rna_interfaceAtomID_str if 0 <= i < len(atom_plddts)]
        #         # print("rna_interfaceAtomID_str", rna_interfaceAtomID_str)
        #         rna_iLDDT = round(np.mean(rna_interface_plddts), 2)
        #         if sourceTable == 'pred_protein_rna' or sourceTable == 'pred_rna_rna' or sourceTable == 'pred_protein_dna':
        #             dna_iLDDT = 0
        elif sourceTable == 'pred_protein_rna_dna':
            # Transform from string to list
            protDNA_interface_atom_ids = ast.literal_eval(protDNA_interface_id_str)
            dna_interface_id_str = ast.literal_eval(dna_interface_id_str)

            # If more than one lists: join lists into a single list
            protDNA_interface_atom_id = [item for sublist in protDNA_interface_atom_ids for item in sublist] if any(
                isinstance(i, list) for i in protDNA_interface_atom_ids) else protDNA_interface_atom_ids
            dna_interface_id = [item for sublist in dna_interface_id_str for item in sublist] if any(
                isinstance(i, list) for i in dna_interface_id_str) else dna_interface_id_str
            print("dna_interface_id_str", dna_interface_id_str)

            # Calculate iLDDT and dna_iLDDT
            interface_plddts = [atom_plddts[i] for i in protDNA_interface_atom_id if 0 <= i < len(atom_plddts)]
            dna_interface_plddts = [atom_plddts[i] for i in dna_interface_id if 0 <= i < len(atom_plddts)]
            print("dna_interface_plddts", dna_interface_id)
            dna_iLDDT = round(np.mean(dna_interface_plddts), 2)
            print("dna_iLDDT", dna_iLDDT)
            dna_filtered_plddts = [plddt for plddt in dna_interface_plddts if plddt != 0]
            iLDDTav = np.average(dna_filtered_plddts) if dna_filtered_plddts else None
            if iLDDTav is not None:
                dna_iLDDTs.append(dna_iLDDT)
            print("dna_iLDDTs", dna_iLDDTs)
            # print(f"Average pLDDT for interface atoms: {iLDDT}, "
            #       f"Average rna pLDDT for interface atoms: {rna_iLDDT}, "
            #           f"min pLDDT for interface atoms: {min(interface_plddts)}, "
            #           f"max pLDDT for interface atoms: {max(interface_plddts)}, "
            #           f"median pLDDT for interface atoms: {np.median(interface_plddts)}")
        else:
            interface_atom_ids = [0]
            interface_plddts = [0]
            iLDDT = 0
            rna_interface_plddts = [0]
            rna_iLDDT = 0
            dna_iLDDT = 0
            print("No interface atoms found.")
        if len(interface_atom_ids) != len(interface_plddts):
            print(f"Error: length of plddts ({len(interface_plddts)}) unequal to length of atom ids ({len(interface_atom_ids)})")
            min_length = min(len(interface_atom_ids), len(interface_plddts))
            interface_atom_ids = interface_atom_ids[:min_length]
            interface_plddts = interface_plddts[:min_length]

        filtered_plddts = [plddt for plddt in interface_plddts if plddt != 0]
        iLDDTm = np.average(filtered_plddts) if filtered_plddts else None
        if iLDDTm is not None:
            iLDDTs.append(iLDDT)
        # print("iLDDTs", iLDDTs)
        # interface_plddts_list = list(zip(valid_interface_atom_ids, interface_plddts))
        # print(interface_plddts_list)

        rna_filtered_plddts = [plddt for plddt in rna_interface_plddts if plddt != 0]
        iLDDTav = np.average(rna_filtered_plddts) if rna_filtered_plddts else None
        if iLDDTav is not None:
            rna_iLDDTs.append(rna_iLDDT)

        ax5 = fig.add_subplot(gs[2, 0])
        ax5.plot(interface_atom_ids, interface_plddts, marker='o')
        ax5.set_title('Interface Atom pLDDTs', fontsize=28)
        ax5.set_xlabel('Atom id', fontsize=18)
        ax5.set_ylabel('pLDDT', fontsize=18)

        interface_plddt_avg = np.mean(interface_plddts)
        ax5.text(0.95, 0.95, f'Avg pLDDT: {interface_plddt_avg:.2f}',
                 horizontalalignment='right', verticalalignment='top',
                 transform=ax5.transAxes, fontsize=14, bbox=dict(facecolor='white', alpha=0.5))

        file_base_graphAF3Metrics = os.path.basename(json_file_path).replace('_full_data_', '_graphAF3Metrics_')
        plt.subplots_adjust(hspace=0.3) # adjust vertical space between graphs
        save_path = os.path.join(folder_path, os.path.join(folder_path, file_base_graphAF3Metrics.replace('.json', '.png')))
        plt.savefig(save_path)
        plt.close(fig)
    else:
        iLDDT = 0
        rna_interface_plddts = [0]
        rna_iLDDT = 0
        dna_iLDDT = 0

    return iLDDT, iLDDTs, rna_iLDDT, rna_iLDDTs, rna_interface_plddts, dna_iLDDT

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python af3Metrics.py <af_folder_path>")
        sys.exit(1)
    folder_path = sys.argv[1]
    json_file_paths_iptm = glob.glob(os.path.join(folder_path, '*_summary_confidences_[0-4].json'))
    json_file_paths_pae = glob.glob(os.path.join(folder_path, '*_full_data_[0-4].json'))

    for json_file_path in json_file_paths_iptm:
        # base_name = os.path.splitext(os.path.basename(json_file_path))[0]
        # file_prefix = '_'.join(base_name.split('_')[:-3]) # fold_q3u8w9_neg6_s1
        ptm, iptm, fraction_disordered, df_iptm_matrix = get_af3Metrics_paths_iptm(folder_path)
        cif_file_name_iptm = derive_cif_file_name(os.path.basename(json_file_path))
        json_file_paths_pae = glob.glob(os.path.join(folder_path, '*_full_data_[0-4].json'))
        for json_file_path in json_file_paths_pae:
            cif_file_name = derive_cif_file_name(os.path.basename(json_file_path))
            # base_name = os.path.splitext(os.path.basename(json_file_path))[0]
            # file_prefix_pae = '_'.join(base_name.split('_')[:-3]) # fold_q3u8w9_neg6_s1
            if cif_file_name == cif_file_name_iptm:
                print(cif_file_name)
                iLDDT, iLDDTs, rna_iLDDT, rna_iLDDTs, rna_interface_plddts, dna_iLDDT = get_af3Metrics_paths_pae(folder_path, df_iptm_matrix)
                insert_data_to_db(cif_file_name, ptm, iptm, fraction_disordered, iLDDT, rna_iLDDT, dna_iLDDT)
                # print(cif_file_name, ptm, iptm, fraction_disordered, iLDDT, rna_iLDDT)
                break
                # create_rna_motif_logo(cif_file_name, rna_interface_plddts)
    # print("iLDDTs", iLDDTs)
    plot_violin(iLDDTs, folder_path, "Average iLDDTs")
    if sourceTable == 'pred_protein_rna' or sourceTable == 'pred_rna_rna' or sourceTable == 'pred_protein_rna_dna':
        plot_violin(rna_iLDDTs, folder_path, "Average RNA iLDDTs")
    elif sourceTable == 'pred_protein_dna':
        plot_violin(dna_iLDDTs, folder_path, "Average DNA iLDDTs")
    if sourceTable == 'pred_protein_rna_dna':
        plot_violin(dna_iLDDTs, folder_path, "Average DNA iLDDTs")

    database_methods.close_connection()