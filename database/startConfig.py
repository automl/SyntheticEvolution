import os
import sys
import json

class StartConfig:
    _instance = None  # Singleton instance for global access
    config_file_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'database', 'config.json')

    def __new__(cls, parent_folder=None):
        # Singleton pattern, ensures only one instance is created
        if cls._instance is None:
            cls._instance = super(StartConfig, cls).__new__(cls)
            cls._instance.load_config(parent_folder)
        return cls._instance

    def set_config(self, parent_folder):
        self.parent_folder = parent_folder
        self.pred_folder = os.path.join(parent_folder, 'af')
        self.ref_folder = os.path.join(parent_folder, 'pdb')
        self.pred_table = self.get_predicted_table_name()
        self.ref_table = self.get_reference_table_name()
        self.db_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'database', 'rbpDatabase.db')

    def clear_config(self):
        self.parent_folder = None
        self.pred_folder = None
        self.ref_folder = None
        self.pred_table = None
        self.ref_table = None
    def get_reference_path(self):
        if not os.path.exists(self.ref_folder):
            print(f"Error: 'pdb' subfolder does not exist in '{parent_folder}'")
            return None
        return self.ref_folder
    def get_predicted_path(self):
        if not os.path.exists(self.pred_folder):
            print(f"Error: 'af' subfolder does not exist in '{parent_folder}'")
            return None
        return self.pred_folder
    def get_predicted_table_name(self):
        folder_name = os.path.basename(self.parent_folder)  # Get the last part of the parent folder path
        return f"pred_{folder_name}"
    def get_reference_table_name(self):
        folder_name = os.path.basename(self.parent_folder)
        return f"exp_{folder_name}"
    def get_database_path(self):
        return self.db_path
    def get_parent_folder(self):
        return self.parent_folder
    def load_config(self, parent_folder=None):
        if os.path.exists(self.config_file_path):
            with open(self.config_file_path, 'r') as file:
                config_data = json.load(file)
                self.set_config(config_data['parent_folder'])
        elif parent_folder:
            self.set_config(parent_folder)
            self.save_config()  # Save initial configuration to JSON

    def save_config(self):
        config_data = {
            'parent_folder': self.parent_folder,
            'pred_folder': self.pred_folder,
            'ref_folder': self.ref_folder,
            'pred_table': self.pred_table,
            'ref_table': self.ref_table,
            'db_path': self.db_path
        }
        with open(self.config_file_path, 'w') as file:
            json.dump(config_data, file, indent=4)

# from startConfig import StartConfig
# table_config = StartConfig()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <parent_folder>")
        sys.exit(1)

    parent_folder = sys.argv[1]
    start_config = StartConfig(parent_folder=parent_folder)
    start_config.save_config()

    # Example output
    print("Configuration Loaded:")
    print(f"Parent Folder: {start_config.parent_folder}")
    print(f"Prediction Folder: {start_config.pred_folder}")
    print(f"Reference Folder: {start_config.ref_folder}")
    print(f"Prediction Table: {start_config.pred_table}")
    print(f"Reference Table: {start_config.ref_table}")
    print(f"Database Path: {start_config.db_path}")