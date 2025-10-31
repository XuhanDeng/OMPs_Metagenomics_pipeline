import pandas as pd
import os
import sys
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Envipath analysis - Match EC numbers with drug biodegradation pathways')
parser.add_argument('--input', type=str, required=True, help='Input KO abundance TSV file')
parser.add_argument('--bt-rules', type=str, required=True, help='BT rules CSV or Excel file')
parser.add_argument('--output-dir', type=str, required=True, help='Output directory')
args = parser.parse_args()

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
# Set working directory to parent (snakemake root)
os.chdir(os.path.dirname(script_dir))

# Load the KO abundance file (TSV format)
ko_abundance = pd.read_csv(
    args.input,
    sep="\t",
    header=0
)

# Load the bt rules file (contains Drug, btpathway, EC information)
# Support both Excel and CSV formats
if args.bt_rules.endswith('.xlsx') or args.bt_rules.endswith('.xls'):
    drug_envipath = pd.read_excel(args.bt_rules, header=None)
else:
    drug_envipath = pd.read_csv(
        args.bt_rules,
        header=None,
        encoding='utf-8-sig'  # Handle BOM if present
    )

# Rename the columns of the drug_envipath DataFrame for clarity
drug_envipath.columns = ["Drug", "btpathway", "EC"]

# Create a new column "EC_number" by prefixing the "EC" column with "EC:"
drug_envipath["EC_number"] = "EC:" + drug_envipath["EC"].astype(str)

# Initialize an empty list to store processed DataFrames
processed_data = []

# Iterate over each row in the drug_envipath DataFrame
for i in range(drug_envipath.shape[0]):
    # Extract the EC number for the current row
    ec_number = drug_envipath.loc[i, "EC_number"]

    # Filter the ko_abundance DataFrame to find rows where the "KO_definition" column contains the current EC number
    tmp = ko_abundance[ko_abundance["KO_definition"].str.contains(ec_number, na=False, regex=False)].copy()

    if tmp.empty:
        continue

    # Add the current EC number, btpathway, and Drug as new columns to the filtered DataFrame
    tmp["EC_third"] = ec_number
    tmp["btpathway"] = drug_envipath.loc[i, "btpathway"]
    tmp["Drug"] = drug_envipath.loc[i, "Drug"]

    # Append the processed DataFrame to the list
    processed_data.append(tmp)

# Concatenate all processed DataFrames into a single DataFrame
if processed_data:
    combined_pathway = pd.concat(processed_data, ignore_index=True)

    # Extract ALL EC numbers from the "KO_definition" column
    # The EC numbers are in the format [EC:x.x.x.x EC:y.y.y.y ...]
    def extract_all_ec_numbers(text):
        import re
        # Find all EC numbers in the text
        ec_matches = re.findall(r'EC:[\d\.]+', str(text))
        return ' '.join(ec_matches) if ec_matches else None

    combined_pathway["EC_number"] = combined_pathway["KO_definition"].apply(extract_all_ec_numbers)

    # Calculate sum of abundance across all samples (skip first 3 columns: KO, gene_count, KO_definition)
    # Also skip the last 3 columns we just added: EC_third, btpathway, Drug
    sample_columns = combined_pathway.columns[3:-3]
    # Convert to numeric, coercing errors to NaN
    combined_pathway[sample_columns] = combined_pathway[sample_columns].apply(pd.to_numeric, errors='coerce')
    combined_pathway["sum"] = combined_pathway[sample_columns].sum(axis=1)

    # Sort by Drug, EC_third, and sum
    envipath_ko = combined_pathway.sort_values(by=["Drug", "EC_third", "sum"], ascending=[True, True, False])

    # Write output to CSV file
    os.makedirs(args.output_dir, exist_ok=True)
    envipath_ko_path = f"{args.output_dir}/envipath_ko.csv"

    envipath_ko.to_csv(envipath_ko_path, index=False)

    print(f"Processing complete! Output saved to:")
    print(f"  - {envipath_ko_path}")
    print(f"\nTotal rows in envipath_ko: {len(envipath_ko)}")
    print(f"\nDrugs processed: {', '.join(envipath_ko['Drug'].unique())}")
else:
    print("No matching data found between bt rules and KO abundance file.")
    sys.exit(1)




