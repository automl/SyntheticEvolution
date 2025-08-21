import sqlite3
import pandas as pd
import os

# Database and TSV file paths
DB_PATH = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/database/rbpDatabase.db"  # Change to your actual database path
TSV_PATH = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/a2021_nonX/a2021nonX_exp_protein_rna.added_domain_infos.tsv"  # Change to the actual path of your TSV file

# Connect to SQLite database
conn = sqlite3.connect(DB_PATH)
cursor = conn.cursor()

# Step 1: Add new columns to the table if they don't exist
columns_to_add = [
    "ALTER TABLE exp_protein_rna ADD COLUMN domain_ids TEXT",
    "ALTER TABLE exp_protein_rna ADD COLUMN domain_names TEXT",
    "ALTER TABLE exp_protein_rna ADD COLUMN domain_counts TEXT",
    "ALTER TABLE exp_protein_rna ADD COLUMN domain_pos TEXT"
]

for query in columns_to_add:
    try:
        cursor.execute(query)
    except sqlite3.OperationalError:
        pass  # Ignore if column already exists

# Step 2: Load TSV file into a Pandas DataFrame
df = pd.read_csv(TSV_PATH, sep="\t")

# Step 3: Update exp_protein table with domain information
for _, row in df.iterrows():
    pdb_id = row["PDBId"]  # Adjust with ProteinID for rnacompete data
    domain_ids = row["domain_ids"]
    domain_names = row["domain_names"]
    domain_counts = row["domain_counts"]
    domain_pos = row["domain_pos"]

    cursor.execute("""
        UPDATE exp_protein_rna
        SET domain_ids = ?, domain_names = ?, domain_counts = ?, domain_pos = ?
        WHERE PDBId = ?
    """, (domain_ids, domain_names, domain_counts, domain_pos, pdb_id))

# Commit and close the database connection
conn.commit()
conn.close()

print("Database updated successfully!")

# # Connect to SQLite database
# conn = sqlite3.connect(DB_PATH)
# cursor = conn.cursor()
#
# # Load the new TSV file into a Pandas DataFrame
# df = pd.read_csv(TSV_PATH, sep="\t")
#
# # Step 1: Ensure that only existing PDBId rows are updated
# for _, row in df.iterrows():
#     pdb_id = row["PDBId"]  # Adjust if the column name differs
#     domain_ids = row["domain_ids"]
#     domain_names = row["domain_names"]
#     domain_counts = row["domain_counts"]
#     domain_pos = row["domain_pos"]
#
#     # Check if the PDBId exists in the database
#     cursor.execute("SELECT COUNT(*) FROM exp_protein_rna WHERE PDBId = ?", (pdb_id,))
#     exists = cursor.fetchone()[0]
#
#     if exists:  # Only update if the PDBId exists in the table
#         cursor.execute("""
#             UPDATE exp_protein_rna
#             SET domain_ids = ?, domain_names = ?, domain_counts = ?, domain_pos = ?
#             WHERE PDBId = ?
#         """, (domain_ids, domain_names, domain_counts, domain_pos, pdb_id))
#
# # Commit and close the database connection
# conn.commit()
# conn.close()
#
# print("Database updated successfully!")