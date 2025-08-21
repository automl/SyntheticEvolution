import os
import sys
import subprocess
from database.startConfig import StartConfig

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # One level above utils package.

def execute_command(script_path, *args):
    command = f"{script_path} {' '.join(args)}"
    # if '.py' in command:
        # print('Found a .py in', command)
        # print(' removing -m...')
        # command.replace(' -m ', ' ')
    print(f"Executing: {command}")
    subprocess.run(command, shell=True, check=True, cwd=project_root)

def process_folder(folder, pdb_command, af_command,  pdb_path, af_path):
    print('folder in process folder:', folder, 'pdb' in folder, 'af' in folder)
    if "pdb" in folder.split('/')[-1]:
        if pdb_command:
            execute_command(f"python {pdb_command} {pdb_path}")

    if "af" in folder.split('/')[-1]:
        if af_command:
            execute_command(f"python {af_command} {af_path}")

def main():
    if len(sys.argv) != 3:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <parent_folder> <+dssr|-dssr>")
        print("Example: python mainExecution.py ./data +dssr")
        print("Options:")
        print("  +dssr: Process DSSR evaluation")
        print("  -dssr: Skip DSSR evaluation")
        sys.exit(1)

    parent_folder = sys.argv[1]
    dssr_flag = sys.argv[2].lower()

    database_path = "database"
    eval_path = "evaluation"

    init_command = os.path.join(project_root, database_path, "startConfig.py")
    execute_command(f"python {init_command} {parent_folder}")

    config = StartConfig()
    pdb_path = config.ref_folder
    af_path = config.pred_folder
    # parent_folder = config.parent_folder
    subfolders = [os.path.join(parent_folder, subfolder) for subfolder in ["pdb", "af"]]
    dssr_path = os.path.join(parent_folder, 'dssr')

    try:
        # addDomain from database and csv_addMSA from utils FIRST!!!
        if dssr_flag == '+dssr' and os.path.exists(dssr_path):
            print("DSSR evaluation required...")
        else:
            for folder in subfolders:
                print("folder", folder, ("pdb" in folder))
                pdb_command = None
                af_command = None
                if "pdb" in folder.split('/')[-1]:
                    pdb_command = os.path.join(project_root, database_path, "parsePDBCIF2DB.py")
                else:
                    af_command = os.path.join(project_root, database_path, "parseAFcif2DB.py")

                process_folder(folder, pdb_command, af_command, pdb_path, af_path)

            for folder in subfolders:
                if "pdb" in folder.split('/')[-1]:
                    pass
                if "af" in folder.split('/')[-1]:
                    execute_command(f"python {os.path.join(project_root, eval_path, 'extractEvalMetrics.py')} {pdb_path} {af_path}")
    except FileNotFoundError as e:
            print(e)

    csv_command = os.path.join(project_root, database_path, "exportDB2CSV.py")
    
    try:
        execute_command(f"python {csv_command} {config.pred_table}")
        execute_command(f"python {csv_command} {config.ref_table}")
    except Exception as e:
        print(e)

if __name__ == "__main__":
    main()