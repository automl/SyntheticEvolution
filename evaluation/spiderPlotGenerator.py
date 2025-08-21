import pandas as pd
import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from plots.plotCreator import PlotCreator

# Change plot directory to be in the same folder as plots
plot_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results', 'plots', 'correlation_analysis')
os.makedirs(plot_dir, exist_ok=True)

# Define metric groups
METRIC_GROUPS = {
    'sequence': ["RNA_GC_Content", "RNA_SequenceComplexity", "domain_counts"],
    'length': ["ProteinLength", "NumberProteins", "RNALength", "NumberRNAs"],
    'motif': ["RNAMotifLength", "RNAmotif_score"],
    'interaction': ["Hbond_proteinRNA", "vdWbond_proteinRNA", "ProteinRNAInterfaceArea", 
                   "ProteinRNAInterfaceRatio", "Free_energy", "Binding_affinity_kd", 
                   "RNA_ElPotential", "RNA_ElPotential_Diff"],
    'bond': ["GlycosidicBond_Similarity", "SugarPucker_Similarity", "GammaAngle_Similarity"],
    'af3': ["af3_protein_pTM", "af3_rna_pTM", "af3_protein_ipTM", "af3_rna_ipTM", 
            "af3_protein_pLDDT_avg", "af3_rna_pLDDT_avg", "af3_global_pae_avg", 
            "af3_chain_pair_pae_min", "af3_fraction_disordered", "af3_has_clash", 
            "af3_ranking_score", "af3_rna_MSA"],
    'structure': ["All_RMSD", "Protein_RMSD", "RNA_RMSD", "Protein_LDDT", "RNA_LDDT", 
                 "Protein_TM", "RNA_TM", "Protein_GDT_TS", "RNA_GDT_TS"],
    'rna_secondary': ["RNA_DI", "RNA_inf", "RNA_wc", "RNA_nonWc", "RNA_stack", 
                     "RNA_f1_score", "RNA_precision", "RNA_recall", "RNA_mcc", "RNA_wl", 
                     "RNAMotif_isDuplex", "RNA_isDuplex", "RNADuplex_Mismatch", 
                     "RNAMotifDuplex_Mismatch"]
}

def create_spider_plot_for_subset(correlations_dict, title, filter_type):
    """Create spider plot for a subset of correlations"""
    print(f"\nCreating spider plot for {title}")
    print(f"Input correlations: {correlations_dict}")
    
    if not correlations_dict:
        print(f"No correlations found for {title}")
        return

    # First, collect all correlation pairs from all subsets
    all_pairs = set()
    for subset_name, subset_data in correlations_dict.items():
        if subset_data is not None:  # Changed condition to check for None
            all_pairs.update(pair for pair, _ in subset_data)
    
    if not all_pairs:
        print(f"No correlation pairs found for {title}")
        return
    
    print(f"Found {len(all_pairs)} unique correlation pairs")
    
    # For each subset, get correlation values for all pairs
    enhanced_correlations = {}
    for subset_name, subset_data in correlations_dict.items():
        print(f"\nProcessing subset: {subset_name}")
        # Try different possible file name formats
        possible_filenames = [
            f'correlation_matrix_{subset_name.lower().replace(" ", "_")}.csv',
            f'correlation_matrix_domain_{subset_name.lower()}.csv',
            f'correlation_matrix_rna_family_{subset_name.lower()}.csv'
        ]
        
        csv_path = None
        for filename in possible_filenames:
            path = os.path.join(plot_dir, filename)
            if os.path.exists(path):
                csv_path = path
                print(f"Found correlation matrix file: {filename}")
                break
        
        if csv_path:
            corr_matrix = pd.read_csv(csv_path, index_col=0)
            print(f"Matrix shape: {corr_matrix.shape}")
            
            # Get correlation values for all pairs from the correlation matrix
            subset_correlations = []
            for pair in all_pairs:
                var1, var2 = pair
                corr_value = 0  # Default value if correlation not found
                
                # Try to find the correlation value in the matrix
                if var1 in corr_matrix.columns and var2 in corr_matrix.index:
                    corr_value = corr_matrix.loc[var2, var1]
                elif var2 in corr_matrix.columns and var1 in corr_matrix.index:
                    corr_value = corr_matrix.loc[var1, var2]
                
                # Store both the correlation value and its sign
                subset_correlations.append((pair, corr_value))
            
            # Always add correlations for this subset
            enhanced_correlations[subset_name] = subset_correlations
            print(f"Added {len(subset_correlations)} correlations for {subset_name}")
        else:
            print(f"Warning: Could not find correlation matrix file for {subset_name}")
            print(f"Tried files: {possible_filenames}")
    
    if enhanced_correlations:  # Only create plot if we have data
        print(f"\nCreating spider plot with data for {len(enhanced_correlations)} subsets:")
        for name, corrs in enhanced_correlations.items():
            print(f"- {name}: {len(corrs)} correlations")
        plot_creator = PlotCreator("correlation")
        plot_creator.create_spider_plot(enhanced_correlations, title, filter_type)
    else:
        print(f"No valid correlations found for {title} after processing")


def get_top_correlations(subset_name):
    """Get top 10 correlations from a correlation matrix CSV file"""
    # Handle both formats of filenames
    base_name = subset_name.lower().replace(" ", "_")
    possible_paths = [
        os.path.join(plot_dir, f'correlation_matrix_{base_name}.csv'),
        os.path.join(plot_dir, f'correlation_matrix_all_data.csv') if base_name == "all_data" else None,
        os.path.join(plot_dir, base_name + '.csv')
    ]
    
    csv_path = None
    for path in possible_paths:
        if path and os.path.exists(path):
            csv_path = path
            break
            
    if csv_path and os.path.exists(csv_path):
        corr_matrix = pd.read_csv(csv_path, index_col=0)
        correlations = []
        
        # Try different thresholds in order
        thresholds = [0.8, 0.5, 0.3]
        
        for threshold in thresholds:
            if not correlations:  # Only try lower threshold if no correlations found yet
                for i in range(len(corr_matrix.columns)):
                    for j in range(i):
                        var1 = corr_matrix.columns[i]
                        var2 = corr_matrix.index[j]
                        corr_value = corr_matrix.iloc[j, i]
                        if abs(corr_value) > threshold and abs(corr_value) < 1:
                            correlations.append({
                                'Variable 1': var1,
                                'Variable 2': var2,
                                'Correlation': corr_value,
                                'Abs_Correlation': abs(corr_value)
                            })
                
                if correlations:
                    print(f"Found correlations above {threshold} for {subset_name}")
        
        if correlations:
            corr_df = pd.DataFrame(correlations)
            corr_df = corr_df.sort_values('Abs_Correlation', ascending=False)
            top_10 = corr_df.head(10)
            return [((row['Variable 1'], row['Variable 2']), row['Correlation']) 
                   for _, row in top_10.iterrows()]
        else:
            print(f"No correlations found above 0.3 for {subset_name}")
    return None


def create_custom_group_spider_plots(correlations_dict, title, filter_type):
    """Create spider plots showing correlations between specific metric groups"""
    if not correlations_dict:
        print(f"No correlations found for {title}")
        return

    # Define custom group combinations
    custom_groups = [
        {
            'name': 'interaction_bond_structure',
            'groups': ['interaction', 'bond', 'structure'],
            'title': 'Interaction-Bond-Structure Correlations'
        }
    ]

    for custom_group in custom_groups:
        print(f"\nProcessing custom group: {custom_group['name']}...")
        group_correlations = {}
        
        # Get all metrics from the specified groups
        target_metrics = []
        for group_name in custom_group['groups']:
            target_metrics.extend(METRIC_GROUPS[group_name])
        
        for subset_name, subset_data in correlations_dict.items():
            if not subset_data:
                continue
                
            # Read the correlation matrix for this subset
            csv_filename = f'correlation_matrix_{subset_name.lower().replace(" ", "_")}.csv'
            csv_path = os.path.join(plot_dir, csv_filename)
            
            if not os.path.exists(csv_path):
                continue
                
            corr_matrix = pd.read_csv(csv_path, index_col=0)
            strongest_correlations = []
            
            # Find correlations between metrics from different groups
            for i, metric1 in enumerate(target_metrics):
                if metric1 not in corr_matrix.index:
                    continue
                    
                for j, metric2 in enumerate(target_metrics[i+1:], i+1):
                    if metric2 not in corr_matrix.columns:
                        continue
                        
                    # Check if metrics are from different groups
                    group1 = next(g for g in custom_group['groups'] if metric1 in METRIC_GROUPS[g])
                    group2 = next(g for g in custom_group['groups'] if metric2 in METRIC_GROUPS[g])
                    
                    if group1 != group2:
                        # Try both orientations of the correlation
                        corr_value = 0
                        if metric1 in corr_matrix.index and metric2 in corr_matrix.columns:
                            corr_value = corr_matrix.loc[metric1, metric2]
                        elif metric2 in corr_matrix.index and metric1 in corr_matrix.columns:
                            corr_value = corr_matrix.loc[metric2, metric1]
                            
                        if abs(corr_value) > 0.3:  # Only include strong correlations
                            # Always store pairs in a consistent order (alphabetically)
                            pair = tuple(sorted([metric1, metric2]))
                            strongest_correlations.append((pair, corr_value))
            
            if strongest_correlations:
                group_correlations[subset_name] = strongest_correlations
        
        if group_correlations:
            group_title = f"{title} - {custom_group['title']}"
            print(f"Creating spider plot for {group_title}")
            
            # Get all unique correlation pairs
            all_pairs = set()
            for subset_data in group_correlations.values():
                all_pairs.update(pair for pair, _ in subset_data)
            
            # Create enhanced correlations including all subsets
            enhanced_correlations = {}
            for subset_name in correlations_dict.keys():
                csv_filename = f'correlation_matrix_{subset_name.lower().replace(" ", "_")}.csv'
                csv_path = os.path.join(plot_dir, csv_filename)
                
                if os.path.exists(csv_path):
                    corr_matrix = pd.read_csv(csv_path, index_col=0)
                    subset_correlations = []
                    
                    for pair in all_pairs:
                        var1, var2 = pair
                        corr_value = 0  # Default value
                        
                        # Try both orientations
                        if var1 in corr_matrix.index and var2 in corr_matrix.columns:
                            corr_value = corr_matrix.loc[var1, var2]
                        elif var2 in corr_matrix.index and var1 in corr_matrix.columns:
                            corr_value = corr_matrix.loc[var2, var1]
                        
                        subset_correlations.append((pair, corr_value))
                    
                    enhanced_correlations[subset_name] = subset_correlations
            
            if enhanced_correlations:
                plot_creator = PlotCreator("correlation")
                plot_creator.create_spider_plot(enhanced_correlations, group_title, f"{filter_type}_{custom_group['name']}")


def create_custom_correlation_spider_plots(correlations_dict, title, filter_type):
    """Create spider plots for specific correlation pairs between RNA structure metrics and bond metrics"""
    if not correlations_dict:
        print(f"No correlations found for {title}")
        return

    # Define the specific correlation pairs we want to analyze
    custom_pairs = [
        ('GlycosidicBond_Similarity', 'RNA_stack'),
        ('SugarPucker_Similarity', 'RNA_stack'),
        ('GammaAngle_Similarity', 'RNA_stack'),
        ('GlycosidicBond_Similarity', 'RNA_nonWc'),
        ('SugarPucker_Similarity', 'RNA_nonWc'),
        ('GammaAngle_Similarity', 'RNA_nonWc'),
        ('GlycosidicBond_Similarity', 'RNA_wc'),
        ('SugarPucker_Similarity', 'RNA_wc'),
        ('GammaAngle_Similarity', 'RNA_wc')
    ]

    # Initialize dictionary to store correlations for each subset
    custom_correlations = {}
    
    for subset_name, subset_data in correlations_dict.items():
        if not subset_data:
            continue
            
        # Read the correlation matrix for this subset
        csv_filename = f'correlation_matrix_{subset_name.lower().replace(" ", "_")}.csv'
        csv_path = os.path.join(plot_dir, csv_filename)
        
        if not os.path.exists(csv_path):
            continue
            
        corr_matrix = pd.read_csv(csv_path, index_col=0)
        subset_correlations = []
        
        # Get correlations for each custom pair
        for pair in custom_pairs:
            var1, var2 = pair
            corr_value = 0  # Default value
            
            # Try both orientations of the correlation
            if var1 in corr_matrix.index and var2 in corr_matrix.columns:
                corr_value = corr_matrix.loc[var1, var2]
            elif var2 in corr_matrix.index and var1 in corr_matrix.columns:
                corr_value = corr_matrix.loc[var2, var1]
            
            # Store the pair in a consistent order (alphabetically)
            sorted_pair = tuple(sorted([var1, var2]))
            subset_correlations.append((sorted_pair, corr_value))
        
        if subset_correlations:
            custom_correlations[subset_name] = subset_correlations
    
    if custom_correlations:
        plot_title = f"{title} - RNA Structure-Bond Correlations"
        print(f"Creating spider plot for {plot_title}")
        plot_creator = PlotCreator("correlation")
        plot_creator.create_spider_plot(custom_correlations, plot_title, f"{filter_type}_rna_structure_bond")


def generate_spider_plots(filter_type='rna_length'):
    """Generate spider plots from existing correlation matrices"""
    spider_plot_data = {}
    
    if filter_type == 'all':
        # For unfiltered data, just get correlations from the all data matrix
        correlations = get_top_correlations("All_Data")
        if correlations:
            spider_plot_data["All_Data"] = correlations
            print(f"Added All_Data to spider plot data")

    elif filter_type == 'rna_length':
        length_ranges = ['0-10', '10-20', '20-60', '>60']
        for length_range in length_ranges:
            name = f'RNA Length {length_range}'
            correlations = get_top_correlations(name)
            if correlations:
                spider_plot_data[name] = correlations

    elif filter_type == 'rna_family':
        # Read all CSV files in the directory to find RNA family correlations
        for filename in os.listdir(plot_dir):
            if filename.startswith('correlation_matrix_rna_family_'):
                family_name = filename.replace('correlation_matrix_rna_family_', '').replace('.csv', '')
                correlations = get_top_correlations(f"RNA_Family_{family_name}")
                if correlations:
                    spider_plot_data[family_name] = correlations
                    print(f"Added {family_name} to spider plot data")

    elif filter_type == 'domain_names':
        # Read all CSV files in the directory to find domain correlations
        for filename in os.listdir(plot_dir):
            if filename.startswith('correlation_matrix_domain_'):
                domain_name = filename.replace('correlation_matrix_domain_', '').replace('.csv', '')
                correlations = get_top_correlations(f"Domain_{domain_name}")
                if correlations:
                    spider_plot_data[domain_name] = correlations
                    print(f"Added {domain_name} to spider plot data")

    elif filter_type == 'rna_complexity':
        complexity_ranges = ['0.00-0.39', '0.40-0.69', '0.70-1.00']
        for complexity_range in complexity_ranges:
            name = f'RNA Complexity {complexity_range}'
            correlations = get_top_correlations(name)
            if correlations:
                spider_plot_data[name] = correlations

    elif filter_type == 'gc_content':
        gc_ranges = ['0.00-0.50', '0.51-1.00']
        for gc_range in gc_ranges:
            name = f'GC Content {gc_range}'
            correlations = get_top_correlations(name)
            if correlations:
                spider_plot_data[name] = correlations

    elif filter_type == 'rna_msa':
        # Process RNA MSA values (0 and 1)
        for msa_value in [0, 1]:
            name = f'RNA_MSA_{msa_value}'
            correlations = get_top_correlations(name)
            if correlations:
                spider_plot_data[name] = correlations
                print(f"Added {name} to spider plot data")

    elif filter_type == 'complex_type':
        # Process complex types (Single and Multi)
        complex_types = ['Single_Complex', 'Multi_Complex']
        for complex_type in complex_types:
            correlations = get_top_correlations(complex_type)
            if correlations:
                spider_plot_data[complex_type] = correlations
                print(f"Added {complex_type} to spider plot data")

    # Create regular spider plot if we have any data
    if spider_plot_data:
        title = f"{filter_type.replace('_', ' ').title()} Correlations"
        print(f"\nCreating regular spider plot for {title} with {len(spider_plot_data)} subsets")
        create_spider_plot_for_subset(spider_plot_data, title, filter_type)
        
        # Create custom correlation spider plots
        print(f"\nCreating custom correlation spider plots for {title}")
        create_custom_correlation_spider_plots(spider_plot_data, title, filter_type)
        
        # Create group-based spider plots
        print(f"\nCreating group-based spider plots for {title}")
        create_group_based_spider_plots(spider_plot_data, title, filter_type)
        
        # Create custom group spider plots
        print(f"\nCreating custom group spider plots for {title}")
        create_custom_group_spider_plots(spider_plot_data, title, filter_type)
    else:
        print(f"No correlations found for {filter_type}")


def create_group_based_spider_plots(correlations_dict, title, filter_type):
    """Create spider plots based on metric groups for each subset"""
    if not correlations_dict:
        print(f"No correlations found for {title}")
        return

    # For each metric group
    for group_name, group_metrics in METRIC_GROUPS.items():
        print(f"\nProcessing {group_name} metrics...")
        group_correlations = {}
        
        for subset_name, subset_data in correlations_dict.items():
            if not subset_data:
                continue
                
            # Read the correlation matrix for this subset
            csv_filename = f'correlation_matrix_{subset_name.lower().replace(" ", "_")}.csv'
            csv_path = os.path.join(plot_dir, csv_filename)
            
            if not os.path.exists(csv_path):
                continue
                
            corr_matrix = pd.read_csv(csv_path, index_col=0)
            strongest_correlations = []
            used_metrics = set()
            
            for metric in group_metrics:
                if metric not in corr_matrix.index:
                    continue
                    
                # Find strongest correlation with any metric not in the same group
                best_corr = 0
                best_pair = None
                
                for other_metric in corr_matrix.columns:
                    # Skip if other metric is in the same group or already used
                    if other_metric in group_metrics or other_metric in used_metrics:
                        continue
                        
                    corr_value = corr_matrix.loc[metric, other_metric]
                    if abs(corr_value) > 0.3 and abs(corr_value) > abs(best_corr):
                        best_corr = corr_value
                        best_pair = (metric, other_metric)  # Store as tuple
                
                if best_pair:
                    used_metrics.add(best_pair[1])
                    strongest_correlations.append((best_pair, best_corr))
            
            if strongest_correlations:
                group_correlations[subset_name] = strongest_correlations
        
        if group_correlations:
            group_title = f"{title} - {group_name.replace('_', ' ').title()} Metrics"
            print(f"Creating spider plot for {group_title}")
            
            # Get all unique correlation pairs
            all_pairs = set()
            for subset_data in group_correlations.values():
                all_pairs.update(pair for pair, _ in subset_data)
            
            # Create enhanced correlations including all subsets
            enhanced_correlations = {}
            for subset_name in correlations_dict.keys():
                csv_filename = f'correlation_matrix_{subset_name.lower().replace(" ", "_")}.csv'
                csv_path = os.path.join(plot_dir, csv_filename)
                
                if os.path.exists(csv_path):
                    corr_matrix = pd.read_csv(csv_path, index_col=0)
                    subset_correlations = []
                    
                    for pair in all_pairs:
                        var1, var2 = pair
                        corr_value = 0  # Default value
                        
                        if var1 in corr_matrix.columns and var2 in corr_matrix.index:
                            corr_value = corr_matrix.loc[var2, var1]
                        elif var2 in corr_matrix.columns and var1 in corr_matrix.index:
                            corr_value = corr_matrix.loc[var1, var2]
                        
                        subset_correlations.append((pair, corr_value))
                    
                    enhanced_correlations[subset_name] = subset_correlations
            
            if enhanced_correlations:
                plot_creator = PlotCreator("correlation")
                plot_creator.create_spider_plot(enhanced_correlations, group_title, f"{filter_type}_{group_name}")

def generate_custom_spider_plots(plot_title, file_title, selected_classes=None, custom_pairs=None):
    """
    Generate a combined spider plot for selected classes and correlation pairs.
    Args:
        selected_classes (list): List of class types to include. If None, includes all classes.
        custom_pairs (list): List of tuples containing pairs of metrics to correlate.
                           If None, uses default pairs.
    """
    print("\nGenerating combined custom correlation spider plot...")
    
    # Define the default correlation pairs if none provided
    if custom_pairs is None:
        custom_pairs = [
            ('GlycosidicBond_Similarity', 'RNA_stack'),
            ('SugarPucker_Similarity', 'RNA_stack'),
            ('GammaAngle_Similarity', 'RNA_stack'),
            ('GlycosidicBond_Similarity', 'RNA_nonWc'),
            ('SugarPucker_Similarity', 'RNA_nonWc'),
            ('GammaAngle_Similarity', 'RNA_nonWc'),
            ('GlycosidicBond_Similarity', 'RNA_wc'),
            ('SugarPucker_Similarity', 'RNA_wc'),
            ('GammaAngle_Similarity', 'RNA_wc')
        ]
    
    # Initialize dictionary to store all correlations
    all_correlations = {}
    
    # Define all subclasses with their colors
    subclass_definitions = {
        'all': {
            'ranges': ['All_Data'],
            'prefix': '',
            'colors': ['#000000']  # Black color for all data
        },
        'rna_length': {
            'ranges': ['0-10', '10-20', '20-60', '>60'],
            'prefix': 'RNA Length ',
            'colors': ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728']  # Blue, Green, Orange, Red
        },
        'rna_complexity': {
            'ranges': ['0.00-0.39', '0.40-0.69', '0.70-1.00'],
            'prefix': 'RNA Complexity ',
            'colors': ['#9467bd', '#8c564b', '#e377c2']  # Purple, Brown, Pink
        },
        'gc_content': {
            'ranges': ['0.00-0.50', '0.51-1.00'],
            'prefix': 'GC Content ',
            'colors': ['#7f7f7f', '#bcbd22']  # Gray, Yellow-Green
        },
        'rna_msa': {
            'ranges': ['0', '1'],
            'prefix': 'RNA MSA ',
            'colors': ['#17becf', '#9edae5']  # Light Blue, Pale Blue
        },
        'complex_type': {
            'ranges': ['Single_Complex', 'Multi_Complex'],
            'prefix': '',
            'colors': ['#393b79', '#637939']  # Dark Blue, Dark Green
        },
        'rna_family': {
            'ranges': None,  # Will be populated from files
            'prefix': 'RNA Family ',
            'colors': ['#8c6d31', '#843c39', '#7b4173', '#5254a3', '#31a354']  # Various colors
        },
        'domain_names': {
            'ranges': None,  # Will be populated from files
            'prefix': 'Domain ',
            'colors': ['#e7ba52', '#ad494a', '#d6616b', '#7b4173', '#a55194']  # Various colors
        }
    }
    
    # If no specific classes selected, use all
    if selected_classes is None:
        selected_classes = list(subclass_definitions.keys())
    
    # Process each selected subclass
    color_index = 0  # Keep track of color assignments
    for class_type in selected_classes:
        if class_type not in subclass_definitions:
            print(f"Warning: Unknown class type {class_type}")
            continue
            
        definition = subclass_definitions[class_type]
        
        if class_type in ['rna_family', 'domain_names']:
            # Handle special cases that need to read from files
            file_prefix = 'correlation_matrix_rna_family_' if class_type == 'rna_family' else 'correlation_matrix_domain_'
            file_count = 0
            for filename in os.listdir(plot_dir):
                if filename.startswith(file_prefix):
                    name = filename.replace(file_prefix, '').replace('.csv', '')
                    csv_path = os.path.join(plot_dir, filename)
                    
                    if os.path.exists(csv_path):
                        corr_matrix = pd.read_csv(csv_path, index_col=0)
                        subset_correlations = []
                        
                        for pair in custom_pairs:
                            var1, var2 = pair
                            corr_value = 0
                            
                            if var1 in corr_matrix.index and var2 in corr_matrix.columns:
                                corr_value = corr_matrix.loc[var1, var2]
                            elif var2 in corr_matrix.index and var1 in corr_matrix.columns:
                                corr_value = corr_matrix.loc[var2, var1]
                            
                            sorted_pair = tuple(sorted([var1, var2]))
                            subset_correlations.append((sorted_pair, corr_value))
                        
                        if subset_correlations:
                            subset_name = f"{definition['prefix']}{name}"
                            color = definition['colors'][file_count % len(definition['colors'])]
                            all_correlations[subset_name] = {
                                'correlations': subset_correlations,
                                'color': color
                            }
                            file_count += 1
        else:
            # Handle regular cases with predefined ranges
            for i, range_value in enumerate(definition['ranges']):
                subset_name = f"{definition['prefix']}{range_value}"
                csv_filename = f'correlation_matrix_{subset_name.lower().replace(" ", "_")}.csv'
                csv_path = os.path.join(plot_dir, csv_filename)
                
                if os.path.exists(csv_path):
                    corr_matrix = pd.read_csv(csv_path, index_col=0)
                    subset_correlations = []
                    
                    for pair in custom_pairs:
                        var1, var2 = pair
                        corr_value = 0
                        
                        if var1 in corr_matrix.index and var2 in corr_matrix.columns:
                            corr_value = corr_matrix.loc[var1, var2]
                        elif var2 in corr_matrix.index and var1 in corr_matrix.columns:
                            corr_value = corr_matrix.loc[var2, var1]
                        
                        sorted_pair = tuple(sorted([var1, var2]))
                        subset_correlations.append((sorted_pair, corr_value))
                    
                    if subset_correlations:
                        color = definition['colors'][i % len(definition['colors'])]
                        all_correlations[subset_name] = {
                            'correlations': subset_correlations,
                            'color': color
                        }
    
    if all_correlations:
        print(f"Creating combined spider plot with {len(all_correlations)} subsets")
        plot_creator = PlotCreator("correlation")
        # Pass the correlations and colors to the spider plot creation
        plot_creator.create_spider_plot(all_correlations, plot_title, file_title, use_custom_colors=True)
    else:
        print("No correlations found for any class")

def main():
    """Main function to generate all spider plots"""
    # Comment out the regular spider plot generation if you only want the custom plots

    print("\nGenerating spider plots for all data...")
    generate_spider_plots('all')
    
    print("\nGenerating spider plots for RNA Length correlations...")
    generate_spider_plots('rna_length')
    
    print("\nGenerating spider plots for RNA Family correlations...")
    generate_spider_plots('rna_family')
    
    print("\nGenerating spider plots for Domain correlations...")
    generate_spider_plots('domain_names')
    
    print("\nGenerating spider plots for RNA Sequence Complexity correlations...")
    generate_spider_plots('rna_complexity')
    
    print("\nGenerating spider plots for GC Content correlations...")
    generate_spider_plots('gc_content')
    
    print("\nGenerating spider plots for RNA MSA correlations...")
    generate_spider_plots('rna_msa')
    
    print("\nGenerating spider plots for Complex Type correlations...")
    generate_spider_plots('complex_type')

    
    # custom_pairs = [
    #     ('GlycosidicBond_Similarity', 'RNA_stack'),
    #     ('SugarPucker_Similarity', 'RNA_stack'),
    #     ('GammaAngle_Similarity', 'RNA_stack'),
    #     ('GlycosidicBond_Similarity', 'RNA_nonWc'),
    #     ('SugarPucker_Similarity', 'RNA_nonWc'),
    #     ('GammaAngle_Similarity', 'RNA_nonWc'),
    #     ('GlycosidicBond_Similarity', 'RNA_wc'),
    #     ('SugarPucker_Similarity', 'RNA_wc'),
    #     ('GammaAngle_Similarity', 'RNA_wc')
    # ]

    custom_pairs = [
        ('RNA_ElPotential_Diff', 'af3_rna_MSA'),
        ('RNA_ElPotential_Diff', 'GlycosidicBond_Similarity'),
        ('RNA_ElPotential_Diff', 'SugarPucker_Similarity'),
        ('RNA_ElPotential_Diff', 'GammaAngle_Similarity'),
        ('RNA_ElPotential_Diff', 'RNA_nonWc'),
        ('RNA_ElPotential_Diff', 'RNA_stack'),
        ('RNA_ElPotential_Diff', 'RNA_wc'),
        ('RNA_ElPotential_Diff', 'RNA_wl'),
        ('RNA_ElPotential_Diff', 'RNA_RMSD'),
        ('RNA_ElPotential_Diff', 'RNA_LDDT'),
        ('RNA_ElPotential_Diff', 'RNA_TM'),
        ('RNA_ElPotential_Diff', 'RNA_GDT_TS'),
        ('RNA_ElPotential', 'Binding_affinity_kd'),
        ('RNA_ElPotential', 'ProteinRNAInterfaceRatio')
    ]
    
    # Generate the combined custom spider plot with selected classes and custom pairs
    selected_classes = ['all','rna_length', 'rna_family', 'rna_complexity', 'gc_content', 'rna_msa']
    plot_title = "El. Potential Correlations"
    file_title = "el_potential_correlations"
    generate_custom_spider_plots(plot_title, file_title, selected_classes, custom_pairs)


if __name__ == "__main__":
    main() 