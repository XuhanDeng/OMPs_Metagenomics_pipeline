import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.stats import zscore
import numpy as np
import sys
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate clustered heatmaps for drug-associated genes')
parser.add_argument('--input', type=str, required=True, help='Input CSV file with all genes')
parser.add_argument('--output-dir', type=str, required=True, help='Output directory for heatmap PDFs')
parser.add_argument('--top-n', type=int, required=True, help='Number of top genes to keep per drug')
args = parser.parse_args()

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
# Set working directory to parent (snakemake root)
os.chdir(os.path.dirname(script_dir))

# Load the data
all_data = pd.read_csv(args.input)
all_data["Ko_EC"] = all_data["KO"].astype(str) + "_" + all_data["EC_number"].astype(str)

# Get the current column order
cols = list(all_data.columns)

# Find the index of 'Ko_EC' and 'EC_number'
ko_ec_index = cols.index('Ko_EC')
ec_number_index = cols.index('EC_number')

# Remove 'Ko_EC' from its original position
cols.pop(ko_ec_index)

# Insert 'Ko_EC' after 'EC_number'
cols.insert(ec_number_index + 1, 'Ko_EC')

# Reassign the DataFrame with the new column order
all_data = all_data[cols]

# Filter to top N genes per drug, sorted by sum
top_n_data = (
    all_data.sort_values(by="sum", ascending=False)
    .groupby("Drug", group_keys=False)
    .apply(lambda x: x.drop_duplicates(subset="EC_number").head(args.top_n), include_groups=False)
    .reset_index(drop=True)
)

print(f"Filtered from {len(all_data)} total genes to {len(top_n_data)} top {args.top_n} genes per drug")


def drug_heatmap(df, drug_name, output_dir):
    """Generate clustered heatmap for a specific drug"""

    # Make a copy to avoid modifying original
    df = df.copy()

    # Set 'Ko_EC' as the index
    df.set_index('Ko_EC', inplace=True)

    # Find sample columns (numeric columns after metadata columns)
    # Exclude: Drug, EC_third, btpathway, KO, gene_count, KO_definition, EC_number, Ko_EC, sum
    metadata_cols = ['Drug', 'EC_third', 'btpathway', 'KO', 'gene_count', 'KO_definition', 'EC_number', 'sum']
    sample_cols = [col for col in df.columns if col not in metadata_cols]

    # Extract sample columns
    samples = df[sample_cols]

    # Check if we have enough data
    if samples.shape[0] < 2 or samples.shape[1] < 2:
        print(f"Warning: Not enough data for {drug_name}, skipping heatmap generation")
        return None

    # Perform z-score normalization (row-wise, across samples)
    samples_normalized = samples.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x, axis=1)

    # Remove rows with all NaN values
    samples_normalized = samples_normalized.dropna(how='all')

    if samples_normalized.shape[0] < 2:
        print(f"Warning: Not enough valid data for {drug_name} after normalization, skipping")
        return None

    # Create a clustered heatmap with dendrograms
    clustermap = sns.clustermap(
        samples_normalized.T,
        cmap='viridis',  # Color map
        annot=False,  # Annotate cells with values
        fmt=".2f",  # Format for annotations
        linewidths=.5,  # Width of lines between cells
        cbar_pos=(1, 0.5, .05, 0.3),
        figsize=(15, 5),  # Size of the figure
        row_cluster=True,  # Cluster rows (genes)
        col_cluster=True,  # Cluster columns (samples)
        dendrogram_ratio=0.2,  # Adjust the size of dendrograms
    )

    # Add a title
    plt.suptitle(drug_name + " associated genes", y=1.02)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Save the clustered heatmap as a PDF
    output_path = os.path.join(output_dir, drug_name + '_clustered_heatmap_with_tree.pdf')
    clustermap.savefig(output_path, format='pdf', bbox_inches='tight')
    plt.close()

    print(f"Generated heatmap for {drug_name}: {output_path}")
    return output_path


# Remove 'sum' column if it exists
if 'sum' in top_n_data.columns:
    top_n_data.drop("sum", axis=1, inplace=True)

# Get unique drug names
drug_name_all = list(set(top_n_data.loc[:, "Drug"]))

print(f"Generating heatmaps for {len(drug_name_all)} drugs (top {args.top_n} genes per drug)...")

# Generate heatmap for each drug
generated_files = []
for i, drug_name in enumerate(drug_name_all):
    print(f"Processing {i+1}/{len(drug_name_all)}: {drug_name}")
    drug_choose = top_n_data.loc[top_n_data["Drug"] == drug_name, :]
    result = drug_heatmap(drug_choose, drug_name, args.output_dir)
    if result:
        generated_files.append(result)

print(f"\nCompleted! Generated {len(generated_files)} heatmaps in {args.output_dir}")
