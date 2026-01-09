#!/usr/bin/env python3
"""
Script to calculate weighted averages from multiple xlsx files.

Each xlsx file should have columns: index, R1, R1err, R2, R2err, Rnoe, Rnoeerr
The script calculates weighted averages using inverse variance weighting (1/err^2).
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys


def calculate_weighted_average(values, errors):
    """
    Calculate weighted average using inverse variance weighting.

    weighted_avg = sum(value / error^2) / sum(1 / error^2)
    weighted_err = sqrt(1 / sum(1 / error^2))

    Args:
        values: array of values
        errors: array of errors (weights)

    Returns:
        tuple: (weighted_average, weighted_error)
    """
    # Remove NaN values
    mask = ~(np.isnan(values) | np.isnan(errors) | (errors == 0))
    values_clean = values[mask]
    errors_clean = errors[mask]

    if len(values_clean) == 0:
        return np.nan, np.nan

    # Inverse variance weights: w = 1/sigma^2
    weights = 1.0 / (errors_clean ** 2)

    # Weighted average
    weighted_avg = np.sum(values_clean * weights) / np.sum(weights)

    # Error on weighted average
    weighted_err = np.sqrt(1.0 / np.sum(weights))

    return weighted_avg, weighted_err


def process_xlsx_files(directory='.', output_file='weighted_averages.xlsx'):
    """
    Process all xlsx files in directory and calculate weighted averages.

    Args:
        directory: directory containing xlsx files
        output_file: output file name for results
    """
    dir_path = Path(directory)

    # Find all xlsx files
    xlsx_files = list(dir_path.glob('*.xlsx'))

    # Filter out the output file if it exists
    xlsx_files = [f for f in xlsx_files if f.name != output_file]

    if not xlsx_files:
        print(f"No xlsx files found in {directory}")
        return

    print(f"Found {len(xlsx_files)} xlsx file(s):")
    for f in xlsx_files:
        print(f"  - {f.name}")

    # Expected columns
    expected_cols = ['index', 'R1', 'R1err', 'R2', 'R2err', 'Rnoe', 'Rnoeerr']

    # Load all files
    all_data = []
    for xlsx_file in xlsx_files:
        print(f"\nLoading {xlsx_file.name}...")
        try:
            df = pd.read_excel(xlsx_file)

            # Check if columns match
            if not all(col in df.columns for col in expected_cols):
                print(f"  Warning: {xlsx_file.name} missing expected columns")
                print(f"  Expected: {expected_cols}")
                print(f"  Found: {list(df.columns)}")
                continue

            # Select only expected columns
            df = df[expected_cols]
            all_data.append(df)
            print(f"  Loaded {len(df)} rows")

        except Exception as e:
            print(f"  Error loading {xlsx_file.name}: {e}")
            continue

    if not all_data:
        print("\nNo valid data loaded. Exiting.")
        return

    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal rows combined: {len(combined_df)}")

    # Get unique indices
    unique_indices = combined_df['index'].unique()
    print(f"Unique indices found: {len(unique_indices)}")

    # Calculate weighted averages for each index
    results = []

    for idx in unique_indices:
        idx_data = combined_df[combined_df['index'] == idx]

        # Calculate weighted averages for R1, R2, Rnoe
        r1_avg, r1_err = calculate_weighted_average(
            idx_data['R1'].values, idx_data['R1err'].values
        )
        r2_avg, r2_err = calculate_weighted_average(
            idx_data['R2'].values, idx_data['R2err'].values
        )
        rnoe_avg, rnoe_err = calculate_weighted_average(
            idx_data['Rnoe'].values, idx_data['Rnoeerr'].values
        )

        results.append({
            'index': idx,
            'R1': r1_avg,
            'R1err': r1_err,
            'R2': r2_avg,
            'R2err': r2_err,
            'Rnoe': rnoe_avg,
            'Rnoeerr': rnoe_err,
            'n_measurements': len(idx_data)
        })

    # Create results dataframe
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('index')

    # Save to xlsx
    output_path = dir_path / output_file
    results_df.to_excel(output_path, index=False)
    print(f"\nWeighted averages saved to: {output_path}")

    # Display summary
    print("\nSummary of weighted averages:")
    print(results_df.to_string(index=False))

    return results_df


if __name__ == '__main__':
    # Get directory from command line argument or use current directory
    directory = sys.argv[1] if len(sys.argv) > 1 else '.'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'weighted_averages.xlsx'

    print(f"Processing xlsx files in: {directory}")
    print(f"Output file: {output_file}\n")

    process_xlsx_files(directory, output_file)
