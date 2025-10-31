#!/usr/bin/env python3
"""
Aggregate gene abundance by KO number
Genes with the same KO annotation will be summed together
"""

import pandas as pd
import sys
import os

def aggregate_by_ko(annotated_file, output_file):
    """
    Aggregate abundance by KO number (sum genes with same KO)

    Args:
        annotated_file: TSV file with annotations + abundance
        output_file: Output file with KO-level aggregated abundance
    """
    print(f"Reading annotated abundance data from {annotated_file}...", flush=True)
    df = pd.read_csv(annotated_file, sep='\t', index_col=0)
    print(f"  Loaded {len(df)} genes", flush=True)

    # Identify annotation columns and sample columns
    annotation_cols = ['KO', 'threshold', 'score', 'evalue', 'KO_definition']
    sample_cols = [col for col in df.columns if col not in annotation_cols]
    print(f"  Found {len(sample_cols)} samples", flush=True)

    # Remove genes without annotation
    df_annotated = df[df['KO'] != 'no_annotation'].copy()
    unannotated_count = len(df) - len(df_annotated)
    print(f"\n  Removing {unannotated_count} unannotated genes", flush=True)
    print(f"  Keeping {len(df_annotated)} annotated genes", flush=True)

    # Group by KO and sum abundance values
    print(f"\nAggregating by KO number...", flush=True)

    # For abundance columns: sum
    ko_abundance = df_annotated.groupby('KO')[sample_cols].sum()

    # For annotation columns: keep first occurrence
    ko_annotation = df_annotated.groupby('KO')[['KO_definition']].first()

    # Also count how many genes per KO
    ko_gene_count = df_annotated.groupby('KO').size()
    ko_gene_count.name = 'gene_count'

    # Merge everything together
    result = pd.concat([ko_gene_count, ko_annotation, ko_abundance], axis=1)

    print(f"\nResults:", flush=True)
    print(f"  Unique KO numbers: {len(result)}", flush=True)
    print(f"  Genes per KO (mean): {ko_gene_count.mean():.1f}", flush=True)
    print(f"  Genes per KO (median): {ko_gene_count.median():.0f}", flush=True)
    print(f"  Genes per KO (max): {ko_gene_count.max()}", flush=True)

    # Save to output
    print(f"\nSaving KO-aggregated abundance matrix...", flush=True)
    result.to_csv(output_file, sep='\t')

    print(f"✓ Output saved to: {output_file}\n", flush=True)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python aggregate_by_ko.py <annotated_abundance_file> <output_file>")
        print("  annotated_abundance_file: File with annotations + abundance")
        print("  output_file: Output file with KO-level aggregated abundance")
        sys.exit(1)

    annotated_file = sys.argv[1]
    output_file = sys.argv[2]

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    aggregate_by_ko(annotated_file, output_file)

    print("✓ KO aggregation completed successfully")
