import os
import json
import matplotlib.pyplot as plt
from collections import Counter

dssr_folder = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/dssr'

chi_counts_pred_total = Counter()
chi_counts_ref_total = Counter()
pucker_counts_pred_total = Counter()
pucker_counts_ref_total = Counter()

def classify_chi_angle(chi_angle):
    """Classify chi angle as 'anti' or 'syn'."""
    if (90 <= chi_angle <= 180) or (-180 <= chi_angle <= -90):
        return 'anti'
    elif -90 <= chi_angle <= 90:
        return 'syn'
    return 'unknown'


def classify_sugar_pucker(pucker):
    """Classify sugar pucker conformation."""
    if pucker == "C3'-endo":
        return 'C3\'-endo'
    elif pucker == "C2'-endo":
        return 'C2\'-endo'
    elif pucker == "C3'-exo":
        return 'C3\'-exo'
    elif pucker == "C4'-exo":
        return 'C4\'-exo'
    return 'other'


# Process each pair of files in the dssr folder
for filename in os.listdir(dssr_folder):
    # Only proceed if it's a predicted file
    if not filename.endswith('_af_dssr.json'):
        continue

    # Extract PDB ID and create paths for predicted and reference files
    pdb_id = filename.split('_')[0]
    predicted_file = os.path.join(dssr_folder, f"{pdb_id}_af_dssr.json")
    reference_file = os.path.join(dssr_folder, f"{pdb_id}_pdb_dssr.json")

    # Load and parse the predicted file
    with open(predicted_file) as f:
        predicted_data = json.load(f)
    for nucleotide in predicted_data.get('nts', []):
        # Classify and count chi orientation
        chi_angle = nucleotide.get('chi', None)
        if chi_angle is not None:
            orientation = classify_chi_angle(chi_angle)
            chi_counts_pred_total[orientation] += 1

        # Classify and count sugar pucker type
        sugar_pucker = nucleotide.get('puckering', 'N/A')
        pucker_class = classify_sugar_pucker(sugar_pucker)
        pucker_counts_pred_total[pucker_class] += 1

    # Load and parse the reference file
    with open(reference_file) as f:
        reference_data = json.load(f)
    for nucleotide in reference_data.get('nts', []):
        # Classify and count chi orientation
        chi_angle = nucleotide.get('chi', None)
        if chi_angle is not None:
            orientation = classify_chi_angle(chi_angle)
            chi_counts_ref_total[orientation] += 1

        # Classify and count sugar pucker type
        sugar_pucker = nucleotide.get('puckering', 'N/A')
        pucker_class = classify_sugar_pucker(sugar_pucker)
        pucker_counts_ref_total[pucker_class] += 1

# Calculate percentages for chi orientation and sugar pucker for all samples in the table
total_chi_pred = sum(chi_counts_pred_total.values())
total_chi_ref = sum(chi_counts_ref_total.values())
total_pucker_pred = sum(pucker_counts_pred_total.values())
total_pucker_ref = sum(pucker_counts_ref_total.values())

# Convert counts to percentages
chi_percent_pred = {k: (v / total_chi_pred * 100) for k, v in chi_counts_pred_total.items()}
chi_percent_ref = {k: (v / total_chi_ref * 100) for k, v in chi_counts_ref_total.items()}
pucker_percent_pred = {k: (v / total_pucker_pred * 100) for k, v in pucker_counts_pred_total.items()}
pucker_percent_ref = {k: (v / total_pucker_ref * 100) for k, v in pucker_counts_ref_total.items()}

# Plotting aggregated percentage results for all pdb_ids
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Chi orientation plot (aggregated percentage)
ax1.bar(chi_percent_ref.keys(), chi_percent_ref.values(), color='red', label='Reference')
ax1.bar(chi_percent_pred.keys(), chi_percent_pred.values(), color='yellow', alpha=0.7, label='Predicted')
# ax1.set_title("Cumulative Chi Orientation Percentages (Syn vs Anti) for All PDB IDs")
ax1.set_xlabel("Chi Orientation")
ax1.set_ylabel("Percentage (%)")
ax1.legend()

# Sugar pucker plot (aggregated percentage)
ax2.bar(pucker_percent_ref.keys(), pucker_percent_ref.values(), color='red', label='Reference')
ax2.bar(pucker_percent_pred.keys(), pucker_percent_pred.values(), color='yellow', alpha=0.7, label='Predicted')
# ax2.set_title("Cumulative Sugar Pucker Type Percentages (C3'-endo vs C2'-endo) for All PDB IDs")
ax2.set_xlabel("Sugar Pucker Type")
ax2.set_ylabel("Percentage (%)")
ax2.legend()

# Display the plots
# plt.suptitle("Aggregate Chi and Pucker Classifications as Percentages for All PDB IDs")
plt.tight_layout()
plt.show()



