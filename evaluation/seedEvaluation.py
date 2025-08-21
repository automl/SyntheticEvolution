import sqlite3
import json
import statistics
import os
import ast

database_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'database', 'rbpDatabase.db')
# print("database_path", database_path)

def calculate_mean_std(values):
    # filtered_values = [v for v in values if v is not None]  # Filter out None values
    filtered_values = [v for v in values if v is not None and isinstance(v, (int, float))]
    if not filtered_values:  # If the filtered list is empty, return None for both mean and std
        return None, None
    mean_val = statistics.mean(filtered_values)
    std_val = statistics.stdev(filtered_values) if len(filtered_values) > 1 else 0.0  # Standard deviation requires at least two values
    return round(mean_val, 2), round(std_val, 2)


def calculate_mean_std_str(values):
    # Filter out None, non-numeric, and invalid values
    filtered_values = []
    for val in values:
        try:
            if val is not None and val != '':
                filtered_values.append(float(val))
        except ValueError:
            continue

    if not filtered_values:
        return (None, None)

    mean_val = round(statistics.mean(filtered_values), 2)
    std_val = round(statistics.stdev(filtered_values), 2) if len(filtered_values) > 1 else 0.0

    return (mean_val, std_val)

def process_rmsd_column(rmsd_column):
    # RMSD column = a list of two values
    # rmsd_1, rmsd_2 = [], []
    # for rmsd_entry in rmsd_column:
    #     if rmsd_entry:
    #         values = json.loads(rmsd_entry)  # Convert the string into a list
    #         rmsd_1.append(values[0])
    #         rmsd_2.append(values[1])
    # return calculate_mean_std(rmsd_1), calculate_mean_std(rmsd_2)
    rmsd_1, rmsd_2 = [], []
    for rmsd_entry in rmsd_column:
        if rmsd_entry:
            if isinstance(rmsd_entry, str):
                try:
                    values = json.loads(rmsd_entry)  # Convert the string into a list
                    if isinstance(values, list) and len(values) >= 2:
                        rmsd_1.append(values[0])
                        rmsd_2.append(values[1])
                except json.JSONDecodeError:
                    # Handle case where rmsd_entry is a string but not a JSON format
                    continue
            elif isinstance(rmsd_entry, (list, tuple)) and len(rmsd_entry) >= 2:
                rmsd_1.append(rmsd_entry[0])
                rmsd_2.append(rmsd_entry[1])
            elif isinstance(rmsd_entry, float):  # Handle float values (if applicable)
                # If only a single float value is provided, you can decide how to handle it
                rmsd_1.append(rmsd_entry)
                rmsd_2.append(0)  # Or some other default value for the second RMSD value
            else:
                continue  # Skip if it's not in the expected format
    return calculate_mean_std(rmsd_1), calculate_mean_std(rmsd_2)

def process_rna_inf_column(rna_inf_column):
    # Mean and std for each of the 4 RNA_INF values, disregarding any -1 values.
    inf_1, inf_2, inf_3, inf_4 = [], [], [], []
    for inf_entry in rna_inf_column:
        if inf_entry:
            values = json.loads(inf_entry)
            # Append values only if they are not -1
            if values[0] != -1:
                inf_1.append(values[0])
            if values[1] != -1:
                inf_2.append(values[1])
            if values[2] != -1:
                inf_3.append(values[2])
            if values[3] != -1:
                inf_4.append(values[3])
    return (calculate_mean_std(inf_1), calculate_mean_std(inf_2),
            calculate_mean_std(inf_3), calculate_mean_std(inf_4))



def process_chainidpair_column(chainid_column):
    item_count = []
    for pair in chainid_column:
        if pair:
            item_count.append(len(json.loads(pair)))  # Count the number of items in the list
    return calculate_mean_std(item_count)


def process_bond_column(bond_column):
    # H- and vdW-bonds
    bond_values = []
    for bond_entry in bond_column:
        if bond_entry:
            bond_values.extend(json.loads(bond_entry))  # Extend list with all values in the list
    return calculate_mean_std(bond_values)


conn = sqlite3.connect(database_path)
cursor = conn.cursor()

cursor.execute("""
    SELECT ChainIDpairList_proteinRNA, Hbond_proteinRNA, vdWbond_proteinRNA, 
           ProteinRNAInterfaceArea, ProteinRNAInterfaceRatio, ptm, iptm, iLDDT, RNA_iLDDT, 
           All_RMSD, Protein_RMSD, RNA_RMSD, Protein_LDDT, RNA_LDDT, Protein_TM, RNA_TM, 
           Protein_GDT_TS, RNA_GDT_TS, RNA_INF, RNA_DI, RNALength
    FROM pred_protein_rna
    WHERE LENGTH(ProteinChainIDs) = 1 AND LENGTH(RNAChainIDs) = 1
""")

rows = cursor.fetchall()

chainid_column = []
hbond_column = []
vdwbond_column = []
protein_rna_interface_area = []
protein_rna_interface_ratio = []
ptm_column = []
iptm_column = []
ilddt_column = []
rna_ilddt_column = []
all_rmsd_column = []
protein_rmsd_column = []
rna_rmsd_column = []
protein_lddt_column = []
rna_lddt_column = []
protein_tm_column = []
rna_tm_column = []
protein_gdt_ts_column = []
rna_gdt_ts_column = []
rna_inf_column = []
rna_di_column = []
rna_length_column = []

for row in rows:
    (chainidpair, hbond, vdwbond, pr_area, pr_ratio, ptm, iptm, ilddt, rna_ilddt,
     all_rmsd, protein_rmsd, rna_rmsd, protein_lddt, rna_lddt, protein_tm, rna_tm,
     protein_gdt_ts, rna_gdt_ts, rna_inf, rna_di, rna_length) = row

    chainid_column.append(chainidpair)
    hbond_column.append(hbond)
    vdwbond_column.append(vdwbond)
    protein_rna_interface_area.append(pr_area)
    protein_rna_interface_ratio.append(pr_ratio)
    ptm_column.append(ptm)
    iptm_column.append(iptm)
    ilddt_column.append(ilddt)
    rna_ilddt_column.append(rna_ilddt)
    all_rmsd_column.append(all_rmsd)
    protein_rmsd_column.append(protein_rmsd)
    rna_rmsd_column.append(rna_rmsd)
    protein_lddt_column.append(protein_lddt)
    rna_lddt_column.append(rna_lddt)
    protein_tm_column.append(protein_tm)
    rna_tm_column.append(rna_tm)
    protein_gdt_ts_column.append(protein_gdt_ts)
    rna_gdt_ts_column.append(rna_gdt_ts)
    rna_inf_column.append(rna_inf)
    rna_di_column.append(rna_di)
    if not isinstance(rna_length, int):
        rna_length = ast.literal_eval(rna_length)
        rna_lengths = [int(length) for length in rna_length]
        max_length = max(rna_lengths)
        rna_length_column.append(max_length)
    else:
        rna_length_column.append(rna_length)

chainid_mean, chainid_std = process_chainidpair_column(chainid_column)
hbond_mean, hbond_std = process_bond_column(hbond_column)
vdwbond_mean, vdwbond_std = process_bond_column(vdwbond_column)

protein_rna_interface_area_mean, protein_rna_interface_area_std = calculate_mean_std(protein_rna_interface_area)
protein_rna_interface_ratio_mean, protein_rna_interface_ratio_std = calculate_mean_std(protein_rna_interface_ratio)
ptm_mean, ptm_std = calculate_mean_std(ptm_column)
iptm_mean, iptm_std = calculate_mean_std(iptm_column)
ilddt_mean, ilddt_std = calculate_mean_std(ilddt_column)
rna_ilddt_mean, rna_ilddt_std = calculate_mean_std(rna_ilddt_column)
rna_length_mean, rna_length_std = calculate_mean_std(rna_length_column)

all_rmsd_mean_1, all_rmsd_mean_2 = process_rmsd_column(all_rmsd_column)
protein_rmsd_mean_1, protein_rmsd_mean_2 = process_rmsd_column(protein_rmsd_column)
rna_rmsd_mean_1, rna_rmsd_mean_2 = process_rmsd_column(rna_rmsd_column)

protein_lddt_mean, protein_lddt_std = calculate_mean_std_str(protein_lddt_column)
rna_lddt_mean, rna_lddt_std = calculate_mean_std_str(rna_lddt_column)
protein_tm_mean, protein_tm_std = calculate_mean_std_str(protein_tm_column)
rna_tm_mean, rna_tm_std = calculate_mean_std_str(rna_tm_column)
protein_gdt_ts_mean, protein_gdt_ts_std = calculate_mean_std_str(protein_gdt_ts_column)
rna_gdt_ts_mean, rna_gdt_ts_std = calculate_mean_std_str(rna_gdt_ts_column)
rna_di_mean, rna_di_std = calculate_mean_std_str(rna_di_column)

rna_inf_means = process_rna_inf_column(rna_inf_column)

print("rna_length Mean:", rna_length_mean, "STD:", rna_length_std)
print("ChainID Pair List Mean:", chainid_mean, "STD:", chainid_std)
print("Hbond Mean:", hbond_mean, "STD:", hbond_std)
print("vdW Bond Mean:", vdwbond_mean, "STD:", vdwbond_std)
print("Protein-RNA Interface Area Mean:", protein_rna_interface_area_mean, "STD:", protein_rna_interface_area_std)
print("Protein-RNA Interface Ratio Mean:", protein_rna_interface_ratio_mean, "STD:", protein_rna_interface_ratio_std)
print("PTM Mean:", ptm_mean, "STD:", ptm_std)
print("IPTM Mean:", iptm_mean, "STD:", iptm_std)
print("iLDDT Mean:", ilddt_mean, "STD:", ilddt_std)
print("RNA iLDDT Mean:", rna_ilddt_mean, "STD:", rna_ilddt_std)
print("All RMSD Mean (1st value):", all_rmsd_mean_1)
print("All RMSD Mean (2nd value):", all_rmsd_mean_2)
print("Protein RMSD Mean (1st value):", protein_rmsd_mean_1)
print("Protein RMSD Mean (2nd value):", protein_rmsd_mean_2)
print("RNA RMSD Mean (1st value):", rna_rmsd_mean_1)
print("RNA RMSD Mean (2nd value):", rna_rmsd_mean_2)
print("Protein LDDT Mean:", protein_lddt_mean, "STD:", protein_lddt_std)
print("RNA LDDT Mean:", rna_lddt_mean, "STD:", rna_lddt_std)
print("Protein TM Mean:", protein_tm_mean, "STD:", protein_tm_std)
print("RNA TM Mean:", rna_tm_mean, "STD:", rna_tm_std)
print("Protein GDT-TS Mean:", protein_gdt_ts_mean, "STD:", protein_gdt_ts_std)
print("RNA GDT-TS Mean:", rna_gdt_ts_mean, "STD:", rna_gdt_ts_std)

print("RNA INF Means:", rna_inf_means)
print("RNA DI Mean:", rna_di_mean, "STD:", rna_di_std)

conn.close()
