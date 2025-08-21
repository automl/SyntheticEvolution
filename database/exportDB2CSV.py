import sqlite3
import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from database.databaseMethods import DatabaseMethods

database = DatabaseMethods()

def convert_table(table):

    conn = sqlite3.connect(database.db_path)
    cursor = conn.cursor()
    cursor.execute(f'SELECT * FROM {table}')
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    results_path = os.path.join(project_root, "results")
    os.makedirs(results_path, exist_ok=True) # Make directory if it doesn't exist
    save_path = os.path.join(results_path, f"{table}.csv")

    rows = cursor.fetchall()
    print('Saving csv to', save_path)
    with open(save_path, 'w', newline='', encoding='utf-8') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow([i[0] for i in cursor.description])

        for row in rows:
            formatted_row = [
                f"{item:.2f}" if isinstance(item, float) else item
                for item in row
            ]
            csv_writer.writerow(formatted_row)

    print('done.')

    conn.close()

def convert_all_dbs(source_folder, table_name):
    project_root = os.path.abspath(source_folder)
    results_path = os.path.join(project_root, "results")
    os.makedirs(results_path, exist_ok=True)

    for filename in os.listdir(source_folder):
        if filename.endswith('.db'):
            db_path = os.path.join(source_folder, filename)
            csv_name = os.path.splitext(filename)[0] + '.csv'
            save_path = os.path.join(results_path, csv_name)

            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()

            try:
                cursor.execute(f"SELECT * FROM {table_name}")
                rows = cursor.fetchall()
                print(save_path)
                with open(save_path, 'w', newline='', encoding='utf-8') as file:
                    csv_writer = csv.writer(file)
                    csv_writer.writerow([i[0] for i in cursor.description])
                    for row in rows:
                        formatted_row = [
                            f"{item:.2f}" if isinstance(item, float) else item
                            for item in row
                        ]
                        csv_writer.writerow(formatted_row)

                print(f"Saved: {csv_name}")
            except sqlite3.Error as e:
                print(f"Failed to process {filename}: {e}")
            finally:
                conn.close()

# convert_all_dbs("/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/9May2025_config/synthetic80-81/combined_DB", "pred_protein_rna")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <desired_pred_table> in quotation marks. ")
        # f"To select all tables, just leave out the <desired_pred_table> parameter.")
        sys.exit(1)
    table = sys.argv[1]
    convert_table(table)
    # print(f'Data exported to {csv_file}')