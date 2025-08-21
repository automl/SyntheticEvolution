import os
import re
import statistics
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def set_plot_style():
    plt.rcParams.update({
        'figure.figsize': (8, 6),  # Adjust the base figure size here (proportional to 10:6)
        'font.size': 20,            # Set a standard font size for all plots
        'axes.titlesize': 20,       # Title size
        'axes.labelsize': 18,       # X and Y label size
        'xtick.labelsize': 16,      # X-axis tick size
        'ytick.labelsize': 16,      # Y-axis tick size
        'legend.fontsize': 16,      # Legend font size
        'figure.autolayout': True   # Automatically adjust layout for better fitting
    })

def search_and_calculate_per_directory(directory):
    # Define the patterns to search for
    pattern = re.compile(r"\(base\)|\(sugar\)|\(phosphate\)|\(side_chain\)|\(backbone\)")

    # Initialize a dictionary to store statistics for each file
    file_statistics = {}

    # Iterate through all .txt files in the specified directory
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as file:
                content = file.read()

                # Initialize counters for the current file
                counts = {
                    "(base)": 0,
                    "(sugar)": 0,
                    "(phosphate)": 0,
                    "(side_chain)": 0,
                    "(backbone)": 0
                }

                # Find all occurrences in the file
                matches = pattern.findall(content)

                # Update the counts based on matches found
                for match in matches:
                    if match == "(base)":
                        counts["(base)"] += 1
                    elif match == "(sugar)":
                        counts["(sugar)"] += 1
                    elif match == "(phosphate)":
                        counts["(phosphate)"] += 1
                    elif match == "(side_chain)":
                        counts["(side_chain)"] += 1
                    elif match == "(backbone)":
                        counts["(backbone)"] += 1

                # Calculate the backbone (sugar + phosphate) and protein backbone count
                rna_backbone_count = counts["(sugar)"] + counts["(phosphate)"]
                protein_backbone_count = counts["(backbone)"]

                # Prepare a dictionary to store statistics for this file
                file_stats = {
                    "RNA base": counts["(base)"],
                    "RNA sugar": counts["(sugar)"],
                    "RNA phosphate": counts["(phosphate)"],
                    # "RNA backbone": rna_backbone_count,
                    "Protein side chain": counts["(side_chain)"],
                    "Protein backbone": protein_backbone_count
                }

                # Store values for statistical calculations
                for key in file_stats:
                    if key not in file_statistics:
                        file_statistics[key] = []
                    file_statistics[key].append(file_stats[key])

    # Calculate overall statistics per interaction type across all files
    overall_statistics = {}
    for key, values in file_statistics.items():
        mean = statistics.mean(values) if values else 0
        std_dev = statistics.stdev(values) if len(values) > 1 else 0
        median = statistics.median(values) if values else 0

        # Calculate lower and upper quartiles
        q1 = statistics.quantiles(values, n=4)[0] if values else 0
        q3 = statistics.quantiles(values, n=4)[2] if values else 0
        iqr = q3 - q1

        # Determine outliers
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        outliers = [v for v in values if v < lower_bound or v > upper_bound]

        overall_statistics[key] = {
            "mean": mean,
            "std_dev": std_dev,
            "median": median,
            "lower_quartile": q1,
            "upper_quartile": q3,
            "outliers": outliers
        }

    return overall_statistics

def plot_statistics(statistics):
    set_plot_style()
    # Prepare data for plotting
    data = {
        'Interaction Type': [],
        'Lower Quartile': [],
        'Upper Quartile': []
    }

    for interaction_type, stats in statistics.items():
        data['Interaction Type'].append(interaction_type)
        data['Lower Quartile'].append(stats['lower_quartile'])
        data['Upper Quartile'].append(stats['upper_quartile'])

    df = pd.DataFrame(data)

    # Create a box plot using seaborn
    # plt.figure(figsize=(12, 8))
    box_plot = sns.boxplot(data=df.melt(id_vars='Interaction Type', value_vars=['Lower Quartile', 'Upper Quartile']),
                           x='Interaction Type', y='value', color='#FF8F00', fliersize=0) ###48A860 FF8F00

    # plt.title('Box Plot of Interaction Types with Quartiles')
    plt.ylabel(' ')
    plt.xlabel(' ')
    plt.xticks(rotation=45)
    plt.tight_layout()
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    plot_path = os.path.join(project_root, "results")
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    title = f"ref_backbone_interactions"
    save_path = os.path.join(plot_path, f"{title}.png")
    plt.savefig(save_path)
    plt.close()

    # Print overall statistics
    for interaction_type, stats in statistics.items():
        print(f"{interaction_type.capitalize()} Statistics:")
        print(f"Mean: {stats['mean']:.2f}")
        print(f"Standard Deviation: {stats['std_dev']:.2f}")
        print(f"Median: {stats['median']:.2f}")
        print(f"Lower Quartile: {stats['lower_quartile']:.2f}")
        print(f"Upper Quartile: {stats['upper_quartile']:.2f}")
        print(f"Outliers: {stats['outliers']}")
        print()


# Example usage
directory_path = "/data/protein_rna/done/a2021_united/pdb"
overall_stats = search_and_calculate_per_directory(directory_path)
plot_statistics(overall_stats)

# # Print the overall statistics
# for interaction_type, stats in overall_stats.items():
#     print(f"{interaction_type.capitalize()} Statistics:")
#     print(f"Mean: {stats['mean']:.2f}")
#     print(f"Standard Deviation: {stats['std_dev']:.2f}")
#     print(f"Median: {stats['median']:.2f}")
#     print()


