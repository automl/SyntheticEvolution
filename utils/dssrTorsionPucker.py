import json

# Load DSSR JSON file
with open('/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/done/a2021_s1/dssr/7r9g_pdb_dssr.json') as f:
    dssr_data = json.load(f)

# Lists to store classified nucleotides based on glycosidic bond and sugar pucker
syn_nucleotides = []
anti_nucleotides = []
c3_endo_nucleotides = []
c2_endo_nucleotides = []
other_pucker_nucleotides = []

# Function to classify glycosidic bond orientation based on chi angle
# https://x3dna.org/highlights/the-chi-x-torsion-angle-characterizes-base-sugar-relative-orientation
def classify_chi_angle(chi_angle):
    if (90 <= chi_angle <= 180) or (-180 <= chi_angle <= -90):
        return 'anti'
    elif -90 <= chi_angle <= 90:
        return 'syn'
    else:
        return 'other'


# Function to classify sugar pucker conformation
def classify_sugar_pucker(pucker):
    if pucker == "C3'-endo":
        return 'C3\'-endo'
    elif pucker == "C2'-endo":
        return 'C2\'-endo'
    else:
        return 'other'


# Parse nucleotides in the JSON data
for nucleotide in dssr_data.get('nts', []):
    # Extract nucleotide ID
    nt_id = nucleotide.get('nt_id')

    # Extract chi angle and classify glycosidic bond orientation
    chi_angle = nucleotide.get('chi', None)
    orientation = classify_chi_angle(chi_angle) if chi_angle is not None else 'unknown'

    # Append to appropriate list based on glycosidic bond orientation
    if orientation == 'anti':
        anti_nucleotides.append({'nt_id': nt_id, 'chi': chi_angle})
    elif orientation == 'syn':
        syn_nucleotides.append({'nt_id': nt_id, 'chi': chi_angle})

    # Extract sugar pucker conformation and classify
    sugar_pucker = nucleotide.get('puckering', 'N/A')
    pucker_class = classify_sugar_pucker(sugar_pucker)

    # Append to appropriate list based on sugar pucker conformation
    if pucker_class == 'C3\'-endo':
        c3_endo_nucleotides.append({'nt_id': nt_id, 'pucker': sugar_pucker})
    elif pucker_class == 'C2\'-endo':
        c2_endo_nucleotides.append({'nt_id': nt_id, 'pucker': sugar_pucker})
    elif pucker_class == 'C2\'-endo':
        c2_endo_nucleotides.append({'nt_id': nt_id, 'pucker': sugar_pucker})
    else:
        other_pucker_nucleotides.append({'nt_id': nt_id, 'pucker': sugar_pucker})

    # Extract backbone torsion angles
    backbone_angles = {key: nucleotide.get(key, 'N/A') for key in
                       ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta']}

    # Print details for each nucleotide
    print(f"Nucleotide ID: {nt_id}")
    print(f"  Chi angle: {chi_angle} ({orientation})")
    print(f"  Sugar pucker: {sugar_pucker} ({pucker_class})")
    print("  Backbone angles:")
    for angle_name, angle_value in backbone_angles.items():
        print(f"    {angle_name}: {angle_value}")
    print("\n")

# Summary output of glycosidic bond and sugar pucker classifications
print("Summary of Glycosidic Bond Orientations:")
print(f"Anti conformation nucleotides ({len(anti_nucleotides)}): {anti_nucleotides}")
print(f"Syn conformation nucleotides ({len(syn_nucleotides)}): {syn_nucleotides}")

print("\nSummary of Sugar Pucker Conformations:")
print(f"C3'-endo nucleotides ({len(c3_endo_nucleotides)}): {c3_endo_nucleotides}")
print(f"C2'-endo nucleotides ({len(c2_endo_nucleotides)}): {c2_endo_nucleotides}")
print(f"Other sugar puckers ({len(other_pucker_nucleotides)}): {other_pucker_nucleotides}")
