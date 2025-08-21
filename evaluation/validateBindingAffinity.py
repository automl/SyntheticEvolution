import numpy as np
from database.databaseMethods import DatabaseMethods
import os
import sys
from scipy import stats

class BindingAffinityValidator:
    def __init__(self):
        self.database = DatabaseMethods()

    def get_binding_data(self, file_name):
        """Get both predicted and experimental data for comparison"""
        tables = ['protein_rna', 'rna_rna', 'protein_dna']
        
        for table in tables:
            # Try to get predicted data
            pred_data = self.database.get_row(f'pred_{table}',
                ['Free_energy', 'Binding_affinity_kd'],
                f"FileName = '{file_name}'")
            
            # Try to get experimental data
            exp_data = self.database.get_row(f'exp_{table}',
                ['Free_energy', 'Binding_affinity_kd'],
                f"FileName = '{file_name}'")
            
            if pred_data and exp_data:
                return {
                    'complex_type': table,
                    'predicted': pred_data,
                    'experimental': exp_data
                }
        
        return None

    def calculate_metrics(self, pred_data, exp_data):
        """Calculate validation metrics"""
        metrics = {
            'absolute_error_energy': abs(pred_data['Free_energy'] - exp_data['Free_energy']),
            'relative_error_energy': abs((pred_data['Free_energy'] - exp_data['Free_energy']) / exp_data['Free_energy']) * 100,
            'absolute_error_kd': abs(pred_data['Binding_affinity_kd'] - exp_data['Binding_affinity_kd']),
            'relative_error_kd': abs((pred_data['Binding_affinity_kd'] - exp_data['Binding_affinity_kd']) / exp_data['Binding_affinity_kd']) * 100
        }
        return metrics

    def validate_binding_parameters(self, file_name):
        """Main validation function"""
        data = self.get_binding_data(file_name)
        if not data:
            raise ValueError(f"No matching predicted and experimental data found for {file_name}")

        metrics = self.calculate_metrics(data['predicted'], data['experimental'])
        
        return {
            'complex_type': data['complex_type'],
            'predicted_energy': data['predicted']['Free_energy'],
            'experimental_energy': data['experimental']['Free_energy'],
            'predicted_kd': data['predicted']['Binding_affinity_kd'],
            'experimental_kd': data['experimental']['Binding_affinity_kd'],
            'metrics': metrics
        }

def main():
    if len(sys.argv) != 2:
        print("Usage: python validateBindingAffinity.py <cif_file>")
        sys.exit(1)

    cif_file = os.path.basename(sys.argv[1])
    validator = BindingAffinityValidator()
    
    try:
        results = validator.validate_binding_parameters(cif_file)
        
        print(f"\nComplex Type: {results['complex_type'].replace('_', '-')}")
        print("\nFree Energy Comparison:")
        print(f"Predicted: {results['predicted_energy']:.2f} kcal/mol")
        print(f"Experimental: {results['experimental_energy']:.2f} kcal/mol")
        print(f"Absolute Error: {results['metrics']['absolute_error_energy']:.2f} kcal/mol")
        print(f"Relative Error: {results['metrics']['relative_error_energy']:.1f}%")
        
        print("\nBinding Affinity Comparison:")
        print(f"Predicted Kd: {results['predicted_kd']:.2e} M")
        print(f"Experimental Kd: {results['experimental_kd']:.2e} M")
        print(f"Absolute Error: {results['metrics']['absolute_error_kd']:.2e} M")
        print(f"Relative Error: {results['metrics']['relative_error_kd']:.1f}%")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main() 