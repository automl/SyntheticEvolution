import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import os
from database.startConfig import StartConfig
from plots.plotCreator import PlotCreator

config = StartConfig()
plot_creator = PlotCreator("feature_importance")

def load_correlation_data(csv_file):
    """Load data from correlation matrix CSV file and prepare it for analysis"""
    # Read the correlation matrix
    df = pd.read_csv(csv_file, index_col=0)
    
    # Convert correlation matrix to long format
    df_long = df.reset_index().melt(id_vars='index', var_name='Variable 2', value_name='Correlation')
    df_long = df_long.rename(columns={'index': 'Variable 1'})
    
    # Remove self-correlations and duplicate pairs
    df_long = df_long[df_long['Variable 1'] != df_long['Variable 2']]
    df_long = df_long[~df_long[['Variable 1', 'Variable 2']].apply(frozenset, axis=1).duplicated()]
    
    # Add all variables as columns with their values
    for col in df.columns:
        df_long[col] = df_long.apply(
            lambda row: df.loc[row['Variable 1'], col] if row['Variable 1'] != col else df.loc[row['Variable 2'], col],
            axis=1
        )
    
    return df_long

def prepare_features_and_target(df, target_variable):
    """Prepare feature matrix X and target variable y"""
    # Remove the target variable and non-feature columns from features
    non_feature_columns = ['Variable 1', 'Variable 2', 'Correlation']
    feature_columns = [col for col in df.columns if col not in non_feature_columns and col != target_variable]
    
    X = df[feature_columns]
    y = df[target_variable]
    
    return X, y

def train_random_forest(X, y):
    """Train a Random Forest model"""
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X, y)
    return model

def calculate_shap_values(model, X):
    """Calculate SHAP values"""
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X)
    return shap_values

def plot_shap_summary(shap_values, X, target_variable, subset_name):
    """Create and save SHAP summary plot"""
    plt.figure(figsize=(12, 8))
    shap.summary_plot(shap_values, X, show=False)
    plt.title(f'SHAP Feature Importance for {target_variable}\n{subset_name}')
    
    # Save plot
    plot_name = f'shap_summary_{target_variable}_{subset_name.lower().replace(" ", "_")}'
    plot_creator.save_plot(plot_name)
    plt.close()

def analyze_feature_importance(csv_file):
    """Analyze feature importance for all parameters"""
    # Extract subset name from filename
    subset_name = os.path.basename(csv_file).replace('correlation_matrix_', '').replace('.csv', '')
    
    # Load data
    df = load_correlation_data(csv_file)
    
    # Get all numeric columns as potential targets
    numeric_columns = df.select_dtypes(include=[np.number]).columns.tolist()
    
    results = []
    for target in numeric_columns:
        print(f"  Analyzing target: {target}")
        
        # Prepare data
        X, y = prepare_features_and_target(df, target)
        
        if len(X.columns) == 0 or len(y) == 0:
            print(f"  Skipping {target} - insufficient data")
            continue
            
        # Handle missing values
        X = X.fillna(X.mean())
        y = y.fillna(y.mean())
        
        if X.empty or y.empty or y.isna().all():
            print(f"  Skipping {target} - all values are NaN")
            continue
        
        try:
            # Train model
            model = train_random_forest(X, y)
            
            # Calculate SHAP values
            shap_values = calculate_shap_values(model, X)
            
            # Plot SHAP summary
            plot_shap_summary(shap_values, X, target, subset_name)
            
            # Calculate feature importance
            feature_importance = pd.DataFrame({
                'Feature': X.columns,
                'Importance': np.abs(shap_values).mean(0)
            })
            feature_importance = feature_importance.sort_values('Importance', ascending=False)
            
            # Save feature importance to CSV
            output_file = os.path.join(
                config.parent_folder,
                f'feature_importance_{target}_{subset_name.lower().replace(" ", "_")}.csv'
            )
            feature_importance.to_csv(output_file, index=False)
            
            # Store top 5 important features
            top_features = feature_importance.head(5)
            results.append({
                'Target': target,
                'Subset': subset_name,
                'Top_Features': top_features['Feature'].tolist(),
                'Importance_Scores': top_features['Importance'].tolist()
            })
        except Exception as e:
            print(f"  Error analyzing {target}: {str(e)}")
    
    return results

def main():
    """Main function to run feature importance analysis"""
    # Find all correlation matrix CSV files
    correlation_files = []
    for root, _, files in os.walk(config.parent_folder):
        for file in files:
            if file.startswith('correlation_matrix_') and file.endswith('.csv'):
                correlation_files.append(os.path.join(root, file))
    
    # Analyze each file
    all_results = []
    for csv_file in correlation_files:
        print(f"\nAnalyzing feature importance for {os.path.basename(csv_file)}...")
        results = analyze_feature_importance(csv_file)
        all_results.extend(results)
    
    # Save summary of all results
    summary_df = pd.DataFrame(all_results)
    summary_file = os.path.join(config.parent_folder, 'feature_importance_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    
    # Create summary text file with grouped results
    with open(os.path.join(config.parent_folder, 'feature_importance_summary.txt'), 'w') as f:
        f.write("Feature Importance Analysis Summary\n\n")
        
        # Group results by subset
        by_subset = {}
        for result in all_results:
            subset = result['Subset']
            if subset not in by_subset:
                by_subset[subset] = []
            by_subset[subset].append(result)
        
        # Write grouped results
        for subset, results in by_subset.items():
            f.write(f"\n{'='*50}\n")
            f.write(f"Subset: {subset}\n")
            f.write(f"{'='*50}\n\n")
            
            for result in results:
                f.write(f"\nTarget: {result['Target']}\n")
                f.write("Top 5 Important Features:\n")
                for feature, score in zip(result['Top_Features'], result['Importance_Scores']):
                    f.write(f"  - {feature}: {score:.4f}\n")
                f.write("\n" + "-"*30 + "\n")

if __name__ == "__main__":
    main() 