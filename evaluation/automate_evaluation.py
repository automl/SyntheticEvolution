import os
import shutil
import subprocess

# 1. change source_folder to your data/sMSA folder in fold_nonserver_config.py in dataPrep package
# 2. in AF3InterfaceEval/data/sMSA you will now have configX folders with the fold_...cif files in af folder
# 3. AF3InterfaceEval/data/protein_rna and AF3InterfaceEval/data/rna_rna folders are already prepared
# 4. change 3. if you have different pdbs and run python parsePDBCIF2DB.py .../data/protein_rna/pdb
# and then copy the new rbpDatabase.db into the /data/protein_rna/ or /rna_rna folders
# 5. optionally run "CONFIGS: save csv from RESULTS" from fold_nonserver_config.py to have everything in 1 folder - change parent_folder
# - sorry fold_nonserver_config is still a mess will functionalise it after the initial draft

# Example:
# /DATA/CONFIGs folder: /Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/2207_config0-10/config0/protein_rna/af and config0/rna_rna/af
# /DATA/protein_rna and /rna_rna folders /Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data/protein_rna/pdb and config.json
# and rbpDatabase.db run only with python parsePDBCIF2DB.py and save in /DATA/protein_rna and /rna_rna folders

def setup_paths():
    base_path = "/home/fred/current_projects/github/af3_evaluation/AF3InterfaceEval"
    return {
        'base': base_path,
        'data': os.path.join(base_path, "data"),
        'database': os.path.join(base_path, "database"),
        'evaluation': os.path.join(base_path, "evaluation"),
        'results': os.path.join(base_path, "results"),
        'configs': os.path.join(base_path, "data/large_scale_hpo_config_predictions")
    }

def process_config_subfolder(config_num, subfolder, paths):
    print(f"\nProcessing config{config_num}/{subfolder}...")

    # Always copy both main dbs to the database directory before processing
    main_db_paths = [
        "/home/fred/current_projects/github/af3_evaluation/AF3InterfaceEval/data/protein_rna/rbpDatabase.db",
        "/home/fred/current_projects/github/af3_evaluation/AF3InterfaceEval/data/rna_rna/rbpDatabase.db"
    ]
    for src_db in main_db_paths:
        if os.path.exists(src_db):
            shutil.copy2(src_db, os.path.join(paths['database'], "rbpDatabase.db"))

    config_folder = os.path.join(paths['configs'], f"config{config_num}")
    af_folder = os.path.join(config_folder, subfolder, "af")
    pdb_folder = os.path.join(paths['data'], f"{subfolder}/pdb")

    # Copy database file
    db_src = os.path.join(paths['data'], f"{subfolder}/rbpDatabase.db")
    db_dst = os.path.join(paths['database'], "rbpDatabase.db")
    shutil.copy2(db_src, db_dst)

    # Copy config file
    config_src = os.path.join(paths['data'], f"{subfolder}/config.json")
    config_dst = os.path.join(paths['database'], "config.json")
    shutil.copy2(config_src, config_dst)

    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = paths['base']
        print("Running parseAFcif2DB.py...")
        subprocess.run([
            "python", 
            os.path.join(paths['database'], "parseAFcif2DB.py"),
            af_folder
        ], env=env, check=True)

        print("Running extractEvalMetrics.py...")
        subprocess.run([
            "python",
            os.path.join(paths['evaluation'], "extractEvalMetrics.py"),
            pdb_folder,
            af_folder
        ], env=env, check=True)

        print("Running exportDB2CSV.py...")
        table_name = f"pred_{subfolder}"
        subprocess.run([
            "python",
            os.path.join(paths['database'], "exportDB2CSV.py"),
            table_name
        ], env=env, check=True)

        # Create results directory in config folder if it doesn't exist
        config_results_dir = os.path.join(config_folder, "results")
        os.makedirs(config_results_dir, exist_ok=True)

        # Move and rename CSV
        csv_src = os.path.join(paths['results'], f"{table_name}.csv")
        csv_dst = os.path.join(config_results_dir, f"config{config_num}_{table_name}.csv")
        if os.path.exists(csv_src):
            shutil.move(csv_src, csv_dst)
            print(f"Saved: {csv_dst}")

        # Move DB file
        db_dst_final = os.path.join(config_results_dir, f"config{config_num}_{table_name}_rbpDatabase.db")
        if os.path.exists(db_dst):
            shutil.move(db_dst, db_dst_final)

        print(f"Successfully processed config{config_num}/{subfolder}")
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error processing config{config_num}/{subfolder}: {str(e)}")
        return False
    except Exception as e:
        print(f"Unexpected error processing config{config_num}/{subfolder}: {str(e)}")
        return False

def main():
    paths = setup_paths()
    # Desired range:
    start = 71
    end = 91

    for config_num in range(start, end):
        for subfolder in ["protein_rna", "rna_rna"]:
            success = process_config_subfolder(config_num, subfolder, paths)
            if not success:
                print(f"Failed to process config{config_num}/{subfolder}. Stopping.")
                break

if __name__ == "__main__":
    main() 