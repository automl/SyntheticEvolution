import sqlite3

# Connect to your SQLite database
conn = sqlite3.connect('rbpDatabase.db')
cursor = conn.cursor()

tables_to_drop = ['exp_protein_rna', 'pred_protein_rna', 'exp_protein_dna', 'pred_protein_dna',
                  'exp_protein_rna_dna', 'pred_protein_rna_dna', 'exp_rna_rna', 'pred_rna_rna',
                  'unified_exp', 'unified_pred']

# Generate the DROP TABLE statements for each table
drop_table_statements = '; '.join([f"DROP TABLE IF EXISTS {table}" for table in tables_to_drop])

try:
    cursor.executescript(drop_table_statements)
    print("Tables deleted successfully.")
except sqlite3.Error as e:
    print(f"An error occurred: {e}")
finally:
    # Commit the changes and close the connection
    conn.commit()
    conn.close()
