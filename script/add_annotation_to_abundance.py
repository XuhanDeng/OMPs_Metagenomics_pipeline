#!/usr/bin/env python3
"""
Add KOfam annotations to abundance matrix and aggregate by KO
"""

import pandas as pd
import sys
import os

def add_annotations_and_aggregate(abundance_file, annotation_file, output_annotated, output_ko_aggregated):
    """
    Merge KOfam annotations with abundance data, then aggregate by KO

    Args:
        abundance_file: TSV file with gene abundance (contigName as index)
        annotation_file: KOfam annotation file
        output_annotated: Output file with gene-level annotations + abundance
        output_ko_aggregated: Output file with KO-level aggregated abundance
    """
    print(f"Reading abundance data from {abundance_file}...", flush=True)
    abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
    print(f"  Loaded {len(abundance_df)} genes × {len(abundance_df.columns)} samples", flush=True)

    print(f"\nReading KOfam annotations from {annotation_file}...", flush=True)
    # Read annotation file, skip comment lines
    anno_df = pd.read_csv(
        annotation_file,
        sep='\t',
        comment='#',
        header=None,
        names=['gene_name', 'KO', 'threshold', 'score', 'evalue', 'KO_definition']
    )
    print(f"  Loaded {len(anno_df)} annotations", flush=True)

    # Set gene_name as index for merging
    anno_df.set_index('gene_name', inplace=True)

    # Merge annotations with abundance (right join to keep all genes)
    print(f"\nMerging annotations with abundance data...", flush=True)
    merged_df = anno_df.join(abundance_df, how='right')

    # Reorder columns: annotations first, then abundance values
    annotation_cols = ['KO', 'threshold', 'score', 'evalue', 'KO_definition']
    sample_cols = [col for col in merged_df.columns if col not in annotation_cols]
    merged_df = merged_df[annotation_cols + sample_cols]

    # Fill missing annotation values
    merged_df['KO'].fillna('no_annotation', inplace=True)
    merged_df['KO_definition'].fillna('no_annotation', inplace=True)

    # Statistics for gene-level data
    annotated = (merged_df['KO'] != 'no_annotation').sum()
    print(f"\n✓ Gene-level results:", flush=True)
    print(f"  Total genes: {len(merged_df)}", flush=True)
    print(f"  Annotated genes: {annotated} ({annotated/len(merged_df)*100:.1f}%)", flush=True)
    print(f"  Unannotated genes: {len(merged_df) - annotated} ({(len(merged_df)-annotated)/len(merged_df)*100:.1f}%)", flush=True)

    # Save gene-level annotated file
    print(f"\nSaving gene-level annotated abundance matrix...", flush=True)
    merged_df.to_csv(output_annotated, sep='\t')
    print(f"✓ Saved to: {output_annotated}", flush=True)

    # ===== Aggregate by KO =====
    print(f"\n{'='*60}", flush=True)
    print(f"Aggregating by KO number...", flush=True)

    # Remove genes without annotation
    df_annotated = merged_df[merged_df['KO'] != 'no_annotation'].copy()
    print(f"  Removing {len(merged_df) - len(df_annotated)} unannotated genes", flush=True)
    print(f"  Keeping {len(df_annotated)} annotated genes for aggregation", flush=True)

    # Group by KO and sum abundance values
    ko_abundance = df_annotated.groupby('KO')[sample_cols].sum()

    # For annotation columns: keep first occurrence
    ko_annotation = df_annotated.groupby('KO')[['KO_definition']].first()

    # Also count how many genes per KO
    ko_gene_count = df_annotated.groupby('KO').size()
    ko_gene_count.name = 'gene_count'

    # Merge everything together
    result = pd.concat([ko_gene_count, ko_annotation, ko_abundance], axis=1)

    print(f"\n✓ KO-level results:", flush=True)
    print(f"  Unique KO numbers: {len(result)}", flush=True)
    print(f"  Genes per KO (mean): {ko_gene_count.mean():.1f}", flush=True)
    print(f"  Genes per KO (median): {ko_gene_count.median():.0f}", flush=True)
    print(f"  Genes per KO (max): {ko_gene_count.max()}", flush=True)

    # Save KO-level aggregated file
    print(f"\nSaving KO-aggregated abundance matrix...", flush=True)
    result.to_csv(output_ko_aggregated, sep='\t')
    print(f"✓ Saved to: {output_ko_aggregated}", flush=True)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python add_annotation_to_abundance.py <abundance_file> <annotation_file> <output_annotated> <output_ko_aggregated>")
        print("  abundance_file: Merged abundance TSV (TPM/Mean/Count)")
        print("  annotation_file: KOfam annotation file")
        print("  output_annotated: Output file with gene-level annotations")
        print("  output_ko_aggregated: Output file with KO-level aggregated abundance")
        sys.exit(1)

    abundance_file = sys.argv[1]
    annotation_file = sys.argv[2]
    output_annotated = sys.argv[3]
    output_ko_aggregated = sys.argv[4]

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_annotated), exist_ok=True)
    os.makedirs(os.path.dirname(output_ko_aggregated), exist_ok=True)

    add_annotations_and_aggregate(abundance_file, annotation_file, output_annotated, output_ko_aggregated)

    print("\n✓ Annotation and aggregation completed successfully")
