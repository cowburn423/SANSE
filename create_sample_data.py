#!/usr/bin/env python3
"""
Create sample xlsx files for testing the weighted average calculation.
"""

import pandas as pd
import numpy as np
from pathlib import Path


def create_sample_xlsx_files(output_dir='sample_data', n_files=3):
    """
    Create sample xlsx files with the expected format.

    Args:
        output_dir: directory to save sample files
        n_files: number of sample files to create
    """
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    # Set random seed for reproducibility
    np.random.seed(42)

    # Common indices for all files
    indices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    for file_num in range(1, n_files + 1):
        # Generate realistic NMR relaxation data with noise
        data = {
            'index': indices,
            'R1': [],
            'R1err': [],
            'R2': [],
            'R2err': [],
            'Rnoe': [],
            'Rnoeerr': []
        }

        for idx in indices:
            # Base values (typical for protein NMR)
            r1_base = 1.5 + 0.3 * np.random.randn()  # ~1.5 s^-1
            r2_base = 10.0 + 2.0 * np.random.randn()  # ~10 s^-1
            rnoe_base = 0.7 + 0.1 * np.random.randn()  # ~0.7

            # Add measurement noise
            r1 = r1_base + 0.05 * np.random.randn()
            r2 = r2_base + 0.5 * np.random.randn()
            rnoe = rnoe_base + 0.02 * np.random.randn()

            # Errors (smaller for more precise measurements)
            r1_err = 0.02 + 0.01 * np.random.rand()
            r2_err = 0.2 + 0.1 * np.random.rand()
            rnoe_err = 0.01 + 0.005 * np.random.rand()

            data['R1'].append(r1)
            data['R1err'].append(r1_err)
            data['R2'].append(r2)
            data['R2err'].append(r2_err)
            data['Rnoe'].append(rnoe)
            data['Rnoeerr'].append(rnoe_err)

        # Create DataFrame
        df = pd.DataFrame(data)

        # Save to xlsx
        output_file = output_path / f'sample_data_{file_num}.xlsx'
        df.to_excel(output_file, index=False)
        print(f"Created: {output_file}")

    print(f"\nCreated {n_files} sample xlsx files in '{output_dir}/' directory")


if __name__ == '__main__':
    create_sample_xlsx_files()
