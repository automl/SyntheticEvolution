import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from evaluationMetrics import RNAMetrics
import ast
from structures.rna import RNA
from structures.dna import DNA
from structures.protein import Protein
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig

database_methods = DatabaseMethods()
config = StartConfig()

def pair_af_pdb_files(dssr_files): #Pair AF and PDB files based on common prefix before '_af_' or '_pdb_'.
    af_files = {}
    pdb_files = {}

    for f in dssr_files:
        if '_af_' in f:
            prefix = f.split('_af_')[0]
            af_files[prefix] = f
        elif '_pdb_' in f:
            prefix = f.split('_pdb_')[0]
            pdb_files[prefix] = f

    common_prefixes = set(af_files.keys()).intersection(set(pdb_files.keys())) # Identify common prefixes between AF and PDB files
    paired_files = [(af_files[prefix], pdb_files[prefix]) for prefix in common_prefixes]
    return paired_files

def load_dssr_data(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

# Extract all/wc/nonWC base pairs
# def extract_type_bp(dssr_data):
#     all_base_pairs = dssr_data['pairs']
#     # print("before", all_base_pairs)
#     wc_base_pairs = extract_dssr_bp([pair for pair in all_base_pairs if pair['name'] == 'WC'])
#     # print("wc_base_pairs", wc_base_pairs)
#     non_wc_base_pairs = extract_dssr_bp([pair for pair in all_base_pairs if pair['name'] not in ['WC']]) # nonWC or pair['name'] != 'WC'???
#     # print("non_wc_base_pairs", non_wc_base_pairs)
#     stacking_base_pairs = [
#         ','.join([f"{stack['nts_short'][i]}({nt.split('/')[-1]})" for i, nt in enumerate(stack['nts_long'].split(','))])
#         for stack in dssr_data['stacks']
#     ]
#     # print("stacking_base_pairs", stacking_base_pairs)
#     all_base_pairs = extract_dssr_bp(all_base_pairs)
#     # print("AFTER", all_base_pairs)
#     return all_base_pairs, wc_base_pairs, non_wc_base_pairs, stacking_base_pairs

def extract_type_bp(dssr_data):
    all_base_pairs = []
    wc_base_pairs = []
    non_wc_base_pairs = []
    stacking_base_pairs = []

    # Check if 'pairs' is directly in dssr_data or if it's nested within other objects
    if isinstance(dssr_data, dict):
        if 'pairs' in dssr_data:
            all_base_pairs = dssr_data['pairs']
        else:
            # Check if 'pairs' is nested within other items in dssr_data
            for item in dssr_data.values():
                if isinstance(item, dict) and 'pairs' in item:
                    all_base_pairs.extend(item['pairs'])
                elif isinstance(item, list):
                    for sub_item in item:
                        if isinstance(sub_item, dict) and 'pairs' in sub_item:
                            all_base_pairs.extend(sub_item['pairs'])

    # Categorise base pairs into WC and non-WC
    wc_base_pairs = extract_dssr_bp([pair for pair in all_base_pairs if pair['name'] == 'WC'])
    non_wc_base_pairs = extract_dssr_bp([pair for pair in all_base_pairs if pair['name'] != 'WC'])

    # Extract stacking pairs if 'stacks' exists in the data
    if 'stacks' in dssr_data:
        stacking_base_pairs = [
            ','.join([f"{stack['nts_short'][i]}({nt.split('/')[-1]})" for i, nt in enumerate(stack['nts_long'].split(','))])
            for stack in dssr_data['stacks']
        ]

    all_base_pairs = extract_dssr_bp(all_base_pairs)

    return all_base_pairs, wc_base_pairs, non_wc_base_pairs, stacking_base_pairs



def extract_dssr_bp(dssr_base_pairs):
    base_pairs = []
    for pair in dssr_base_pairs:
        if '-' in pair['bp']:
            nt1, nt2 = pair['bp'].split('-')
        elif '+' in pair['bp']:
            nt1, nt2 = pair['bp'].split('+')
        else:
            nt1, nt2 = pair['bp'], None
        base_pairs.append({'nt1': nt1, 'nt2': nt2})
    return base_pairs

# based on RNA-Assessment def INF(self, src_struct, trg_struct, type):
def calculate_inf(ref_base_pairs, pred_base_pairs, stacking = None): # INF comparing base pairs from a source and target structure.

    def pair_set(base_pairs):
        return set((pair['nt1'], pair['nt2']) for pair in base_pairs)

    if stacking:
        src_pair_set = set(ref_base_pairs)
        trg_pair_set = set(pred_base_pairs)
        # print("src_pair_set", src_pair_set)
        # print("trg_pair_set", trg_pair_set)

        if not (src_pair_set & trg_pair_set):

            def convert_stacking_base_pairs(stacking_base_pairs): # Without indices
                converted_pairs = []
                for pair in stacking_base_pairs:
                    nucleotides = ''.join([nt.split('(')[0] for nt in pair.split(',')])
                    converted_pairs.append(nucleotides)
                return converted_pairs

            src_pair_set = set(convert_stacking_base_pairs(ref_base_pairs))
            trg_pair_set = set(convert_stacking_base_pairs(pred_base_pairs))
            # print("src_pair_set", src_pair_set)
            # print("trg_pair_set", trg_pair_set)
    else:
        src_pair_set = pair_set(ref_base_pairs)
        trg_pair_set = pair_set(pred_base_pairs)

    TP = len(src_pair_set & trg_pair_set)  # Intersection
    FN = len(src_pair_set - trg_pair_set)  # Elements in src but not in trg
    FP = len(trg_pair_set - src_pair_set)  # Elements in trg but not in src

    if TP == 0 and (FP == 0 or FN == 0):
        INF = -1.0
    else: # Precision and sensitivity
        PPV = float(TP) / (float(TP) + float(FP)) if (TP + FP) > 0 else 0
        STY = float(TP) / (float(TP) + float(FN)) if (TP + FN) > 0 else 0
        INF = (PPV * STY) ** 0.5

    return INF

def get_RMSD_from_DB():
    pass

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <dssr_directory_path>")
        sys.exit(1)
    directory = sys.argv[1]

    dssr_files = [f for f in os.listdir(directory) if f.endswith('.json')]
    paired_files = pair_af_pdb_files(dssr_files)
    # print(paired_files)

    for af_file, pdb_file in paired_files:
        print(af_file, " and ", pdb_file)
        inf = []
        pdb_id = os.path.basename(pdb_file).split('.')[0][:4].upper()
        pdb_dssr_data = load_dssr_data(os.path.join(directory, pdb_file))
        af_dssr_data = load_dssr_data(os.path.join(directory, af_file))
        pdb_all_bp, pdb_wc_bp, pdb_nonWC_bp, pdb_stacking_bp = extract_type_bp(pdb_dssr_data)
        af_all_bp, af_wc_bp, af_nonWC_bp, af_stacking_bp = extract_type_bp(af_dssr_data)
        # print("All Base Pairs:", extract_type_bp(af_dssr_data))
        all_inf_value = round(calculate_inf(pdb_all_bp, af_all_bp),2)
        inf.append(all_inf_value)
        # if all_inf_value != -1: #would it be better to filter them out at a later stage? Are we updating the list instead of single value?
        #     inf.append(all_inf_value)
        print(f"All INF: {all_inf_value}")
        wc_inf_value = round(calculate_inf(pdb_wc_bp, af_wc_bp),2)
        inf.append(wc_inf_value)
        print(f"WC INF: {wc_inf_value}")
        nonWC_inf_value = round(calculate_inf(pdb_nonWC_bp, af_nonWC_bp),2)
        inf.append(nonWC_inf_value)
        print(f"nonWC INF: {nonWC_inf_value}")
        stacking_inf_value = round(calculate_inf(pdb_stacking_bp, af_stacking_bp, True),2)
        inf.append(stacking_inf_value)
        print(f"Stacking INF: {stacking_inf_value}")
        rna_metrics = RNAMetrics(pdb_id)
        rna_metrics.load_from_db(config.pred_table)
        rna_rmsd = rna_metrics.get_rna_rmsd() #rna_rna solved with get_rna_rmsd_values and  if isinstance(rna_rmsd, list): rna_rmsd_value = float(rna_rmsd[0]) else rna_rmsd_value = float(rna_rmsd)
        # if rna_rmsd.startswith("["):
        #     rna_rmsd_list = ast.literal_eval(rna_rmsd)
        #     if isinstance(rna_rmsd_list, list) and len(rna_rmsd_list) > 0:
        #         rna_rmsd_value = float(rna_rmsd_list[0])
        # elif rna_rmsd.strip():
        #     rna_rmsd_value = float(rna_rmsd)
        print(rna_metrics.rna_rmsd, rna_metrics.id)
        print("rna_rmsd",rna_rmsd, rna_rmsd.strip() )
        
        # Handle empty or None rna_rmsd values
        if not rna_rmsd or rna_rmsd.strip() == '':
            print(f"Warning: No RMSD value found for {pdb_id}")
            rna_rmsd_value = None
        elif rna_rmsd.startswith("["):
            rna_rmsd_list = ast.literal_eval(rna_rmsd)
            if isinstance(rna_rmsd_list, list) and len(rna_rmsd_list) > 0:
                rna_rmsd_value = float(rna_rmsd_list[0])
            else:
                rna_rmsd_value = None
        else:
            try:
                rna_rmsd_value = float(rna_rmsd)
            except ValueError:
                print(f"Warning: Invalid RMSD value for {pdb_id}: {rna_rmsd}")
                rna_rmsd_value = None

        # print(pdb_id)
        if all_inf_value == 0 or rna_rmsd_value is None:
            di_all = None
        else:
            di_all = round((rna_rmsd_value / all_inf_value), 2)
            print(f"DI all: {di_all}")

        file_name_cif_af = RNA.get_filename_from_id(config.pred_table, pdb_id)
        rna = RNA.get_rna_from_db(id=pdb_id, file_name=file_name_cif_af)
        protein = Protein.get_protein_from_db(id=pdb_id, file_name=file_name_cif_af)
        dna = DNA.get_dna_from_db(id=pdb_id, file_name=file_name_cif_af)
        inf_values = str(inf)
        # print("££££££, prot, rna, dna", protein.is_protein, rna.is_rna, dna.is_dna)
        if protein.is_protein and rna.is_rna and not dna.is_dna:
            database_methods.update_or_insert('pred_protein_rna',
                                              ['RNA_INF', 'RNA_DI'],
                                              (inf_values, di_all),
                                              condition=f"exp_db_id = '{pdb_id}'")

        if protein.is_protein and rna.is_rna and dna.is_dna:
            database_methods.update_or_insert('pred_protein_rna_dna',
                                              ['RNA_INF', 'RNA_DI'],
                                              (inf_values, di_all),
                                              condition=f"exp_db_id = '{pdb_id}'")

        if not protein.is_protein and rna.is_rna and not dna.is_dna:
            database_methods.update_or_insert('pred_rna_rna',
                                              ['RNA_INF', 'RNA_DI'],
                                              (inf_values, di_all),
                                              condition=f"exp_db_id = '{pdb_id}'")

    database_methods.close_connection()
