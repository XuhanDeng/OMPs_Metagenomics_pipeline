#!/usr/bin/env python3
"""
Merge abundance files (TPM, mean coverage, read counts) from multiple samples
Using pairwise merging to minimize memory usage
"""

import pandas as pd
import sys
import os
import tempfile

def merge_two_files(file1, file2, sample1, sample2):
    """
    Merge two abundanc  e files into a single dataframe

    Args:
        file1, file2: File paths
        sample1, sample2: Sample names

    Returns:
        Merged dataframe
    """
    # Read both files
    df1 = pd.read_csv(file1, sep='\t')
    df1.columns = ['contigName', sample1]
    df1.set_index('contigName', inplace=True)

    df2 = pd.read_csv(file2, sep='\t')
    df2.columns = ['contigName', sample2]
    df2.set_index('contigName', inplace=True)

    # Merge with outer join
    merged = df1.join(df2, how='outer')
    merged.fillna(0, inplace=True)

    # Clean up
    del df1, df2

    return merged

def merge_abundance_files_pairwise(file_list, output_file, metric_name):
    """
    Merge multiple abundance files using pairwise merging (like merge sort)
    Only keeps 2 files in memory at a time

    Args:
        file_list: List of file paths to merge
        output_file: Output file path
        metric_name: Name of the metric (TPM, Mean, Read Count)
    """
    print(f"Merging {len(file_list)} {metric_name} files using pairwise approach...", flush=True)
    print(f"This will use temporary files to minimize memory usage\n", flush=True)

    # Create temp directory for intermediate files
    temp_dir = tempfile.mkdtemp(prefix=f'merge_{metric_name}_')
    print(f"Temporary directory: {temp_dir}\n", flush=True)

    # Initialize working list with original files
    current_files = []
    for filepath in file_list:
        sample_name = os.path.basename(os.path.dirname(filepath))
        current_files.append((filepath, sample_name, False))  # (path, name, is_temp)

    iteration = 1

    # Iteratively merge pairs until only one file remains
    while len(current_files) > 1:
        print(f"=== Iteration {iteration} ===", flush=True)
        print(f"Files to merge: {len(current_files)}", flush=True)

        next_files = []

        # Process pairs
        i = 0
        while i < len(current_files):
            if i + 1 < len(current_files):
                # Merge two files
                file1, sample1, is_temp1 = current_files[i]
                file2, sample2, is_temp2 = current_files[i + 1]

                print(f"  [{i//2 + 1}] Merging {sample1} + {sample2}...", end='', flush=True)

                # Read file1
                df = pd.read_csv(file1, sep='\t')
                if is_temp1:
                    # Temp file already has index as first column
                    df1 = df.set_index(df.columns[0])
                else:
                    # Original file has 2 columns: contigName and values
                    df.columns = ['contigName', sample1]
                    df1 = df.set_index('contigName')
                del df

                # Read file2
                df = pd.read_csv(file2, sep='\t')
                if is_temp2:
                    # Temp file already has index as first column
                    df2 = df.set_index(df.columns[0])
                else:
                    # Original file has 2 columns: contigName and values
                    df.columns = ['contigName', sample2]
                    df2 = df.set_index('contigName')
                del df

                # Merge
                merged = df1.join(df2, how='outer')
                merged.fillna(0, inplace=True)

                print(f" {len(merged)} genes × {len(merged.columns)} samples", flush=True)

                # Clean up
                del df1, df2

                # Delete temp files
                if is_temp1:
                    os.remove(file1)
                if is_temp2:
                    os.remove(file2)

                # Save merged result to temp file
                temp_file = os.path.join(temp_dir, f'merged_{iteration}_{i//2}.tsv')
                merged.to_csv(temp_file, sep='\t')
                del merged

                next_files.append((temp_file, f"merged_{iteration}_{i//2}", True))
                i += 2

            else:
                # Odd file out, carry forward to next iteration
                print(f"  Carrying forward: {current_files[i][1]}", flush=True)
                next_files.append(current_files[i])
                i += 1

        current_files = next_files
        iteration += 1
        print(flush=True)

    # Final file is the merged result
    final_file, _, _ = current_files[0]
    print(f"Moving final result to {output_file}...", flush=True)

    # Read final file and save to output
    final_df = pd.read_csv(final_file, sep='\t')
    final_df.to_csv(output_file, sep='\t', index=False)

    print(f"✓ Merged matrix: {len(final_df)} genes × {len(final_df.columns) - 1} samples", flush=True)
    print(f"✓ Output saved to: {output_file}", flush=True)

    # Clean up temp directory
    import shutil
    shutil.rmtree(temp_dir)
    print(f"✓ Cleaned up temporary directory\n", flush=True)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python merge_abundance.py <metric_type> <output_file> <input_files...>")
        print("  metric_type: TPM, Mean, or Count")
        sys.exit(1)

    metric_type = sys.argv[1]
    output_file = sys.argv[2]
    input_files = sys.argv[3:]

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Merge files using pairwise approach
    merge_abundance_files_pairwise(input_files, output_file, metric_type)

    print(f"✓ {metric_type} merge completed successfully")
