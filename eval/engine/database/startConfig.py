import os
import sys
import json

class StartConfig:
    _instance = None  # Singleton instance for global access
    # config filename is env-overridable so concurrent run_eval processes (SLURM array
    # shards) never share/clobber a single config.json. Default preserves old behavior.
    config_file_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'database',
                                    os.environ.get('EVAL_CONFIG_NAME', 'config.json'))

    def __new__(cls, parent_folder=None):
        # Singleton pattern, ensures only one instance is created
        if cls._instance is None:
            cls._instance = super(StartConfig, cls).__new__(cls)
            cls._instance.load_config(parent_folder)
        return cls._instance

    def set_config(self, classification=None, base_data_path=None, parent_folder=None):
        """Set configuration from classification. Assumes classification is always correct and available.
        
        Can be called with:
        - classification and base_data_path (new method)
        - parent_folder (backward compatibility - extracts classification from folder name)
        - single positional argument treated as parent_folder if it's a path
        """
        # Handle backward compatibility: if first arg is a path (contains / or \), treat as parent_folder
        if classification and (os.sep in str(classification) or '/' in str(classification) or '\\' in str(classification)):
            parent_folder = classification
            classification = None
        
        if parent_folder:
            # Backward compatibility: extract classification from parent_folder
            folder_name = os.path.basename(parent_folder)
            base_data_path = os.path.dirname(parent_folder)
            classification = folder_name
        
        if not classification or not base_data_path:
            raise ValueError("Either (classification and base_data_path) or parent_folder must be provided")
        
        self.classification = classification
        self.base_data_path = base_data_path
        self.parent_folder = os.path.join(base_data_path, classification)
        self.pred_folder = os.path.join(self.parent_folder, 'af')
        self.ref_folder = os.path.join(self.parent_folder, 'pdb')
        self.pred_table = f"pred_{classification}"
        self.ref_table = f"exp_{classification}"
        
        # Set db_path (can be overridden in config)
        if not hasattr(self, 'db_path') or self.db_path is None:
            self.db_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'database', 'rbpDatabase.db')

    def clear_config(self):
        self.parent_folder = None
        self.pred_folder = None
        self.ref_folder = None
        self.pred_table = None
        self.ref_table = None
        self.classification = None
        self.base_data_path = None
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
        """Get predicted table name from classification."""
        return f"pred_{self.classification}"
    
    def get_reference_table_name(self):
        """Get reference table name from classification."""
        return f"exp_{self.classification}"
    def get_database_path(self):
        return self.db_path
    def get_parent_folder(self):
        return self.parent_folder
    def load_config(self, parent_folder=None):
        """Load configuration from config.json. Assumes classification is always present and correct."""
        if os.path.exists(self.config_file_path):
            with open(self.config_file_path, 'r') as file:
                config_data = json.load(file)
                classification = config_data['classification']
                base_data_path = config_data.get('base_data_path', 
                    os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data'))
                db_path = config_data.get('db_path', 
                    os.path.join(os.path.dirname(os.path.dirname(__file__)), 'database', 'rbpDatabase.db')
                )
                self.set_config(classification, base_data_path)
                self.db_path = db_path
        elif parent_folder:
            # Extract classification from parent_folder for backward compatibility
            folder_name = os.path.basename(parent_folder)
            base_data_path = os.path.dirname(parent_folder)
            self.set_config(folder_name, base_data_path)
            self.save_config()  # Save initial configuration to JSON

    def save_config(self):
        """Save configuration in classification-based format."""
        config_data = {
            'classification': self.classification,
            'base_data_path': self.base_data_path,
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