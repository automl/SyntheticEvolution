import numpy as np
import sys
import os

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from database.databaseMethods import DatabaseMethods
from structures.proteinRNAdnaComplex import RNAProteinDNAcomplex

class BindingAffinityEstimator:
    def __init__(self):
        self.database = DatabaseMethods()
        
        # Energy parameters with specific references
        self.energy_params = {
            'protein_rna': {
                # Core parameters
                'HBOND_ENERGY': -0.8,     # Bahadur et al. (2008) JMB 384:1195-1206
                'VDW_ENERGY': -0.15,      # Jones et al. (2001) NAR 29:943-954
                'SURFACE_ENERGY': -0.025,  # Bahadur & Zacharias (2008) Biophys J 94:3066-3074
                
                # Electrostatic factors
                'ELECTROSTATIC_FACTOR': 1.2,  # Ellis et al. (2007) JMB 366:1639-1653
                'IONIC_STRENGTH_FACTOR': 0.9,  # Draper (2004) RNA 10:335-343
                
                # RNA-specific recognition
                'BASE_SPECIFIC_FACTOR': 1.1,   # Allers & Shamoo (2001) JMB 311:75-86
                'SHAPE_COMPLEMENTARITY': 0.9,  # Chen & Varani (2013) FEBS Lett 587:1094-1102
                'ACCESSIBILITY_FACTOR': 1.05,  # Ellis & Jones (2008) NAR 36:5134-5146
                'SEQUENCE_SPECIFICITY': 1.1    # Kligun & Mandel-Gutfreund (2013) NAR 41:8434-8443
            },
            'rna_rna': {
                # Parameters based on:
                # - Turner et al. (2010) Nucleic Acids Res. 38, D280-D282
                # - Mathews et al. (2004) PNAS 101, 7287-7292
                'HBOND_ENERGY': -1.0,     # Strong H-bonds in RNA base pairing
                'VDW_ENERGY': -0.12,      # Base stacking contributions
                'SURFACE_ENERGY': -0.02,   # RNA-RNA interface burial
                'ELECTROSTATIC_FACTOR': 1.4,  # Strong phosphate-phosphate repulsion
                'STACK_FACTOR': 1.3,          # Base stacking contribution
                'HELIX_FACTOR': 1.1           # RNA helical structure contribution
            },
            'protein_dna': {
                # Parameters based on:
                # - Rohs et al. (2010) Nature 461, 1248-1253
                # - Luscombe et al. (2001) Nucleic Acids Res. 29, 2860-2874
                'HBOND_ENERGY': -0.7,     # DNA groove-specific H-bonds
                'VDW_ENERGY': -0.14,      # Shape complementarity
                'SURFACE_ENERGY': -0.022,  # DNA major/minor groove contacts
                'ELECTROSTATIC_FACTOR': 1.1,  # DNA phosphate backbone effect
                'GROOVE_FACTOR': 1.2,         # Major/minor groove specificity
                'CONFORMATIONAL_FACTOR': 0.95  # DNA conformational changes
            }
        }
        
        self.R = 0.001987  # Gas constant in kcal/(mol·K)
        self.T = 298       # Temperature in Kelvin

    def get_complex_data(self, table_name, identifier):
        """Get interface data from database based on table name"""
        column_mappings = {
            'protein_rna': ['Hbond_proteinRNA', 'vdWbond_proteinRNA', 'ProteinRNAInterfaceArea'],
            'rna_rna': ['Hbond_rnaRNA', 'vdWbond_rnaRNA', 'RNARNAInterfaceArea'],
            'protein_dna': ['Hbond_proteinDNA', 'vdWbond_proteinDNA', 'ProteinDNAInterfaceArea']
        }

        prefix = 'pred_' if table_name.startswith('pred_') else 'exp_'
        complex_type = table_name.replace(prefix, '')
        
        if complex_type not in column_mappings:
            raise ValueError(f"Unsupported table type: {table_name}")

        columns = column_mappings[complex_type]
        id_column = 'exp_db_id' if prefix == 'pred_' else 'PDBId'
        
        # Use get_row instead of fetch_one
        query_result = self.database.get_row(table_name, 
            columns,
            f"{id_column} = '{identifier}'")

        if not query_result:
            return None

        return {
            'type': complex_type,
            'hbonds': query_result[columns[0]],
            'vdw_bonds': query_result[columns[1]],
            'interface_area': float(query_result[columns[2]])
        }

    def calculate_free_energy(self, complex_data):
        """Calculate binding free energy based on complex type with specific adjustments"""
        complex_type = complex_data['type']
        params = self.energy_params[complex_type]
        
        # Convert string lists to actual lists if needed and sum up the values
        if isinstance(complex_data['hbonds'], str):
            hbonds_list = eval(complex_data['hbonds'])
        else:
            hbonds_list = complex_data['hbonds']
        num_hbonds = sum(hbonds_list)

        if isinstance(complex_data['vdw_bonds'], str):
            vdw_list = eval(complex_data['vdw_bonds'])
        else:
            vdw_list = complex_data['vdw_bonds']
        num_vdw = sum(vdw_list)

        interface_area = complex_data['interface_area']

        # Base energy calculations
        hbond_energy = num_hbonds * params['HBOND_ENERGY']
        vdw_energy = num_vdw * params['VDW_ENERGY']
        surface_energy = interface_area * params['SURFACE_ENERGY']

        # Apply complex-specific adjustments
        if complex_type == 'protein_rna':
            hbond_energy *= (params['ELECTROSTATIC_FACTOR'] * params['IONIC_STRENGTH_FACTOR'])
            hbond_energy *= (params['BASE_SPECIFIC_FACTOR'] * params['SEQUENCE_SPECIFICITY'])
            surface_energy *= (params['SHAPE_COMPLEMENTARITY'] * params['ACCESSIBILITY_FACTOR'])

        total_energy = hbond_energy + vdw_energy + surface_energy
        return total_energy, {
            'complex_type': complex_type,
            'hbond_contribution': hbond_energy,
            'vdw_contribution': vdw_energy,
            'surface_contribution': surface_energy
        }

    def calculate_binding_affinity(self, free_energy):
        """Calculate binding affinity (Kd) from free energy"""
        kd = np.exp(free_energy / (self.R * self.T))
        return kd

    def store_binding_data(self, identifier, results, is_predicted=True):
        """
        Store binding affinity results in database
        Args:
            identifier: exp_db_id or PDBId
            results: dictionary containing binding parameters
            is_predicted: boolean to determine table prefix (pred_ or exp_)
        """
        table_prefix = 'pred_' if is_predicted else 'exp_'
        id_column = 'exp_db_id' if is_predicted else 'PDBId'
        
        # Prepare data for storage
        data_columns = [
            'Free_energy',
            'Binding_affinity_kd'
        ]
        
        data_values = (
            results['free_energy'],
            float(results['binding_affinity_kd'].replace('e', 'E'))
        )

        try:
            self.database.update_or_insert(
                f"{table_prefix}{results['complex_type']}",
                data_columns,
                data_values,
                condition=f"{id_column} = '{identifier}'"
            )
        except Exception as e:
            print(f"Error storing data for {identifier}: {e}")

    def estimate_binding_parameters(self, table_name, identifier):
        """Main function to estimate binding parameters"""
        complex_data = self.get_complex_data(table_name, identifier)
        if not complex_data:
            raise ValueError(f"No data found for {identifier} in {table_name}")

        free_energy, energy_components = self.calculate_free_energy(complex_data)
        kd = self.calculate_binding_affinity(free_energy)

        results = {
            'free_energy': round(free_energy, 2),
            'binding_affinity_kd': f"{kd:.2e}",
            'complex_type': energy_components['complex_type'],
            'num_hbonds': sum(eval(complex_data['hbonds']) if isinstance(complex_data['hbonds'], str) 
                             else complex_data['hbonds']),
            'num_vdw': sum(eval(complex_data['vdw_bonds']) if isinstance(complex_data['vdw_bonds'], str) 
                          else complex_data['vdw_bonds']),
            'interface_area': complex_data['interface_area']
        }
        
        return results

def main():
    if len(sys.argv) != 2:
        print("Usage: python estimateBindingAffinity.py <table_name>")
        print("Example: python estimateBindingAffinity.py pred_protein_rna")
        print("Supported tables: pred_protein_rna, pred_rna_rna, pred_protein_dna,")
        print("                 exp_protein_rna, exp_rna_rna, exp_protein_dna")
        sys.exit(1)

    table_name = sys.argv[1]
    valid_tables = ['pred_protein_rna', 'pred_rna_rna', 'pred_protein_dna',
                   'exp_protein_rna', 'exp_rna_rna', 'exp_protein_dna']
    
    if table_name not in valid_tables:
        print(f"Error: Invalid table name. Must be one of: {', '.join(valid_tables)}")
        sys.exit(1)

    estimator = BindingAffinityEstimator()
    is_predicted = table_name.startswith('pred_')
    
    # Get all identifiers from the table
    id_column = 'exp_db_id' if is_predicted else 'PDBId'
    rows = estimator.database.select(table_name, [id_column])
    
    if not rows:
        print(f"Error: No entries found in {table_name}")
        sys.exit(1)

    print(f"\nProcessing {len(rows)} entries from {table_name}...")
    
    for row in rows:
        identifier = row[0]  # First column contains the ID
        try:
            results = estimator.estimate_binding_parameters(table_name, identifier)
            estimator.store_binding_data(identifier, results, is_predicted)

            print(f"\nID: {identifier}")
            print(f"Free Energy (ΔG): {results['free_energy']} kcal/mol")
            print(f"Binding Affinity (Kd): {results['binding_affinity_kd']} M")
            print("-" * 50)
        except Exception as e:
            print(f"Error processing {identifier}: {e}")
            continue

    # print(f"\nProcessing complete. Processed {len(rows)} entries.")

if __name__ == "__main__":
    main()
