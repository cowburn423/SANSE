# Weighted Average XLSX Processor

This tool processes multiple Excel (.xlsx) files containing measurement data and calculates weighted averages using inverse variance weighting.

**Available in both Python and MATLAB versions.**

## Overview

The script identifies all `.xlsx` files in a directory, loads them as tables with specific columns, and produces weighted averages for each measurement using the error columns as weights.

## Expected Data Format

Each Excel file should contain the following columns:
- `index`: Identifier for each measurement point
- `R1`: R1 relaxation rate measurement
- `R1err`: Error/uncertainty in R1 measurement
- `R2`: R2 relaxation rate measurement
- `R2err`: Error/uncertainty in R2 measurement
- `Rnoe`: NOE (Nuclear Overhauser Effect) measurement
- `Rnoeerr`: Error/uncertainty in Rnoe measurement

## Weighted Average Calculation

The script uses **inverse variance weighting**, which is the standard method for combining measurements with different uncertainties:

```
weighted_avg = Σ(value_i / error_i²) / Σ(1 / error_i²)
weighted_err = √(1 / Σ(1 / error_i²))
```

This method gives more weight to measurements with smaller errors (higher precision).

## Installation

### Python

Install required Python packages:

```bash
pip3 install pandas openpyxl
```

### MATLAB

No additional packages required. The MATLAB scripts use built-in functions.

## Usage

### Python Usage

### Basic Usage

Process all `.xlsx` files in the current directory:

```bash
python3 weighted_average_xlsx.py
```

### Specify Directory

Process files in a specific directory:

```bash
python3 weighted_average_xlsx.py /path/to/data/directory
```

### Custom Output File

Specify both input directory and output filename:

```bash
python3 weighted_average_xlsx.py /path/to/data output_results.xlsx
```

## Example

The repository includes a sample data generator:

```bash
# Create sample data
python3 create_sample_data.py

# Process the sample data
python3 weighted_average_xlsx.py sample_data
```

This will:
1. Create 3 sample `.xlsx` files in the `sample_data/` directory
2. Process all files and calculate weighted averages
3. Save results to `sample_data/weighted_averages.xlsx`

### MATLAB Usage

#### Basic Usage

Process all `.xlsx` files in the current directory:

```matlab
results = weighted_average_xlsx();
```

#### Specify Directory

Process files in a specific directory:

```matlab
results = weighted_average_xlsx('path/to/data/directory');
```

#### Custom Output File

Specify both input directory and output filename:

```matlab
results = weighted_average_xlsx('path/to/data', 'output_results.xlsx');
```

#### Example

The repository includes a sample data generator and example script:

```matlab
% Create sample data
create_sample_data('sample_data_matlab', 3);

% Process the sample data
results = weighted_average_xlsx('sample_data_matlab', 'weighted_averages_matlab.xlsx');

% Or simply run the complete example
example_weighted_average;
```

This will:
1. Create 3 sample `.xlsx` files in the `sample_data_matlab/` directory
2. Process all files and calculate weighted averages
3. Save results to `sample_data_matlab/weighted_averages_matlab.xlsx`
4. Display results and create plots

## Output

The script generates an Excel file containing:
- `index`: The measurement point identifier
- `R1`, `R1err`: Weighted average and error for R1
- `R2`, `R2err`: Weighted average and error for R2
- `Rnoe`, `Rnoeerr`: Weighted average and error for Rnoe
- `n_measurements`: Number of measurements combined for each index

The output also displays a summary table in the terminal showing the calculated weighted averages.

## Features

- **Automatic file discovery**: Finds all `.xlsx` files in the specified directory
- **Data validation**: Checks for expected columns and reports issues
- **Robust error handling**: Handles missing values and zero errors gracefully
- **Combines multiple files**: Merges data from all files before calculating averages
- **Grouped by index**: Calculates separate weighted averages for each unique index value
- **Comprehensive output**: Includes error estimates and measurement counts

## Notes

- The script automatically excludes the output file when searching for input files
- NaN values and zero errors are automatically filtered out before calculation
- Each index can have different numbers of measurements across files
- The weighted error provides a statistically valid estimate of uncertainty

## Files in Repository

### Python Files
- `weighted_average_xlsx.py` - Main Python script for processing xlsx files
- `create_sample_data.py` - Python script to generate sample data
- `sample_data/` - Directory containing Python-generated sample files

### MATLAB Files
- `weighted_average_xlsx.m` - Main MATLAB function for processing xlsx files
- `create_sample_data.m` - MATLAB function to generate sample data
- `example_weighted_average.m` - Complete MATLAB example with plots
- `sample_data_matlab/` - Directory for MATLAB-generated sample files (created on first run)
