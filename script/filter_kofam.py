#!/usr/bin/env python3
"""
Filter KOfam annotations to keep only the best hit per gene.
Filters based on:
1. Score > threshold
2. E-value <= 1e-5
3. Best hit = highest (score - threshold) difference
"""

import pandas as pd
import sys

def filter_kofam_annotations(input_file, output_file, log_file):
    """Filter KOfam annotations to best hits per gene."""

    try:
        # Read the file, skip first 2 lines (headers and separator line)
        df = pd.read_csv(input_file, sep='\t', skiprows=2, header=None)

        # Assign column names based on KOfam output format
        df.columns = ['marker', 'gene', 'KO', 'threshold', 'score', 'evalue', 'definition']

        # Remove rows where marker is '*' (no KO assignment)
        df = df[df['marker'] != '*'].copy()

        # Remove leading/trailing whitespace from all string columns
        for col in df.columns:
            if df[col].dtype == 'object':
                df[col] = df[col].str.strip()

        # Convert numeric columns
        df['threshold'] = pd.to_numeric(df['threshold'], errors='coerce')
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
        df['evalue'] = pd.to_numeric(df['evalue'], errors='coerce')

        # Remove rows with missing values in critical columns
        df = df.dropna(subset=['gene', 'threshold', 'score', 'evalue'])

        # Calculate score difference
        df['score_diff'] = df['score'] - df['threshold']

        # Filter: score > threshold AND evalue <= 1e-5
        df_filtered = df[(df['score_diff'] > 0) & (df['evalue'] <= 1e-5)].copy()

        # Keep only the best hit per gene (highest score_diff)
        df_best = df_filtered.loc[df_filtered.groupby('gene')['score_diff'].idxmax()]

        # Sort by gene name
        df_best = df_best.sort_values('gene')

        # Drop the score_diff column before writing
        df_best = df_best.drop(columns=['score_diff'])

        # Write header
        with open(output_file, 'w') as f:
            f.write('#\tgene name\tKO\tthrshld\tscore\tE-value\t"KO definition"\n')

        # Write filtered results (without index, with tab separator)
        df_best.to_csv(output_file, sep='\t', index=False, header=False, mode='a')

        # Log summary
        with open(log_file, 'w') as f:
            f.write(f"Total input lines: {len(df)}\n")
            f.write(f"Lines passing filters: {len(df_filtered)}\n")
            f.write(f"Best hits per gene: {len(df_best)}\n")

        return 0

    except Exception as e:
        with open(log_file, 'w') as f:
            f.write(f"Error: {str(e)}\n")
        return 1

if __name__ == "__main__":
    # Snakemake provides these variables
    input_file = snakemake.input.annotations
    output_file = snakemake.output.filtered
    log_file = snakemake.log.out

    sys.exit(filter_kofam_annotations(input_file, output_file, log_file))
