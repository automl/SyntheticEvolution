import sqlite3
from collections import Counter

# Path to your database
db_path = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/database/rbpDatabase.db'

# Connect to the database
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Get all PDB IDs from both tables
cursor.execute("SELECT PDBId FROM exp_protein_rna")
exp_all = [row[0].strip().lower() for row in cursor.fetchall()]
exp_ids = set(exp_all)
exp_counter = Counter(exp_all)

cursor.execute("SELECT exp_db_id FROM pred_protein_rna")
pred_all = [row[0].strip().lower() for row in cursor.fetchall()]
pred_ids = set(pred_all)
pred_counter = Counter(pred_all)

# Compute differences
missing_in_exp = pred_ids - exp_ids
missing_in_pred = exp_ids - pred_ids

# Report missing
print("❌ In pred_protein_rna but missing from exp_protein_rna:")
for item in sorted(missing_in_exp):
    print("  -", item)
print("Total missing in exp:", len(missing_in_exp))

print("\n❌ In exp_protein_rna but missing from pred_protein_rna:")
for item in sorted(missing_in_pred):
    print("  -", item)
print("Total missing in pred:", len(missing_in_pred))

# Report duplicates
print("\n⚠️ Duplicates in exp_protein_rna:")
for item, count in exp_counter.items():
    if count > 1:
        print(f"  - {item}: {count} entries")

print("\n⚠️ Duplicates in pred_protein_rna:")
for item, count in pred_counter.items():
    if count > 1:
        print(f"  - {item}: {count} entries")

# Close connection
conn.close()