import json
import numpy as np
import glob
import os
import sys
import re

def calculate_avg_pae_from_pae_array(json_file_path):
    with open(json_file_path, 'r') as f:
        data = json.load(f)

    # Extract the PAE array
    pae_array = data.get('pae', [])

    # Flatten the PAE array if it's 2D or nested
    flattened_pae = [pae for sublist in pae_array for pae in sublist] if isinstance(pae_array[0], list) else pae_array

    # Calculate the average PAE
    avg_pae = np.mean(flattened_pae) if flattened_pae else None
    return round(avg_pae, 2)


# Function to retrieve metrics from the summary_confidences JSON file (fold_xxxx_s1_summary_confidences_0.json)
def get_af3Metrics_paths_iptm(summary_confidences_file):
    # Load the JSON data
    with open(summary_confidences_file, 'r') as f:
        data = json.load(f)

    # Extract relevant metrics from the JSON structure
    ptm = data.get('ptm', 0.0)
    iptm = data.get('iptm', 0.0)
    fraction_disordered = data.get('fraction_disordered', 0.0)
    ranking_score = data.get('ranking_score', 0.0)
    num_recycles = data.get('num_recycles', 0)
    has_clash = data.get('has_clash', 0.0)

    return ptm, iptm, fraction_disordered, ranking_score, num_recycles, has_clash


# Function to retrieve iLDDT and other metrics from the PAE file
def get_af3Metrics_paths_pae(folder_path):
    # Placeholder function to get iLDDT and other metrics from your files
    iLDDT = 0.75
    rna_iLDDT = 0.80
    iLDDTs = [0.75, 0.80]  # Example list
    rna_iLDDTs = [0.76, 0.81]  # Example list
    rna_interface_plddts = [0.78, 0.79]  # Example list
    dna_iLDDT = 0.74
    return iLDDT, iLDDTs, rna_iLDDT, rna_iLDDTs, rna_interface_plddts, dna_iLDDT


# Function to insert data into the database
def insert_data_to_db(file_name, ptm, iptm, fraction_disordered, iLDDT, rna_iLDDT, ranking_score, num_recycles,
                      has_clash, avg_pae):
    # Assuming you have a `database_methods` module or method to perform database insertions
    print(f"Inserting data for file {file_name}:")
    print(f"ptm={ptm}, iptm={iptm}, fraction_disordered={fraction_disordered}, iLDDT={iLDDT}, rna_iLDDT={rna_iLDDT}, "
          f"ranking_score={ranking_score}, num_recycles={num_recycles}, has_clash={has_clash}, avg_pae={avg_pae}")
    # Example database insertion
    # database_methods.update_or_insert(...) - Uncomment when ready to insert into your database


# Main function to process JSON files
def process_json_files(folder_path):
    # Get the JSON file paths for PAE and other metrics
    json_file_paths_iptm = glob.glob(os.path.join(folder_path, '*_summary_confidences_[0-4].json'))
    json_file_paths_pae = glob.glob(os.path.join(folder_path, '*_full_data_[0-4].json'))

    for json_file_path in json_file_paths_iptm:
        # Retrieve metrics from the summary_confidences file
        ptm, iptm, fraction_disordered, ranking_score, num_recycles, has_clash = get_af3Metrics_paths_iptm(
            json_file_path)

        # Derive the corresponding CIF file name from the JSON file
        cif_file_name_iptm = derive_cif_file_name(os.path.basename(json_file_path))
        print("cif_file_name_iptm", cif_file_name_iptm)
        for json_file_path_pae in json_file_paths_pae:
            cif_file_name = derive_cif_file_name(os.path.basename(json_file_path_pae))
            print(cif_file_name)
            if cif_file_name == cif_file_name_iptm:
                print(f"Processing file: {cif_file_name}")

                # Calculate avg_pae from the PAE array in the full data JSON
                avg_pae = calculate_avg_pae_from_pae_array(json_file_path_pae)

                # Get other metrics (like iLDDT, RNA iLDDT, etc.) from the PAE file
                iLDDT, iLDDTs, rna_iLDDT, rna_iLDDTs, rna_interface_plddts, dna_iLDDT = get_af3Metrics_paths_pae(
                    folder_path)

                # Insert the calculated values into the database
                insert_data_to_db(cif_file_name, ptm, iptm, fraction_disordered, iLDDT, rna_iLDDT,
                                  ranking_score, num_recycles, has_clash, avg_pae)
                break

def derive_cif_file_name(filename):
    # Use regular expression to match the protein code pattern (letters and digits)
    match = re.match(r"^(?:fold_)?([a-zA-Z0-9]+)_s1_(?:summary_confidences|full_data)_0.json$", filename)
    if match:
        # Extract the protein code
        return match.group(1)
    else:
        raise ValueError(f"Filename {filename} does not match expected pattern")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python af3Metrics.py <af_folder_path>")
        sys.exit(1)
    folder_path = sys.argv[1]
    process_json_files(folder_path)
