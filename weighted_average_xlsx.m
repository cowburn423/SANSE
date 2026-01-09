function results = weighted_average_xlsx(directory, output_file)
% WEIGHTED_AVERAGE_XLSX Calculate weighted averages from multiple xlsx files
%
% Usage:
%   results = weighted_average_xlsx()
%   results = weighted_average_xlsx(directory)
%   results = weighted_average_xlsx(directory, output_file)
%
% Inputs:
%   directory   - Directory containing xlsx files (default: current directory)
%   output_file - Output filename (default: 'weighted_averages.xlsx')
%
% Output:
%   results     - Table containing weighted averages
%
% Description:
%   Processes all xlsx files in the specified directory with columns:
%   index, R1, R1err, R2, R2err, Rnoe, Rnoeerr
%
%   Calculates weighted averages using inverse variance weighting:
%   weighted_avg = sum(value / error^2) / sum(1 / error^2)
%   weighted_err = sqrt(1 / sum(1 / error^2))

    % Default parameters
    if nargin < 1 || isempty(directory)
        directory = '.';
    end

    if nargin < 2 || isempty(output_file)
        output_file = 'weighted_averages.xlsx';
    end

    fprintf('Processing xlsx files in: %s\n', directory);
    fprintf('Output file: %s\n\n', output_file);

    % Find all xlsx files
    xlsx_files = dir(fullfile(directory, '*.xlsx'));

    % Filter out the output file if it exists
    xlsx_files = xlsx_files(~strcmp({xlsx_files.name}, output_file));

    if isempty(xlsx_files)
        fprintf('No xlsx files found in %s\n', directory);
        results = [];
        return;
    end

    fprintf('Found %d xlsx file(s):\n', length(xlsx_files));
    for i = 1:length(xlsx_files)
        fprintf('  - %s\n', xlsx_files(i).name);
    end
    fprintf('\n');

    % Expected columns
    expected_cols = {'index', 'R1', 'R1err', 'R2', 'R2err', 'Rnoe', 'Rnoeerr'};

    % Load all files
    all_data = [];

    for i = 1:length(xlsx_files)
        filepath = fullfile(directory, xlsx_files(i).name);
        fprintf('Loading %s...\n', xlsx_files(i).name);

        try
            % Read the xlsx file
            data = readtable(filepath);

            % Check if columns match
            if ~all(ismember(expected_cols, data.Properties.VariableNames))
                fprintf('  Warning: %s missing expected columns\n', xlsx_files(i).name);
                fprintf('  Expected: ');
                fprintf('%s ', expected_cols{:});
                fprintf('\n  Found: ');
                fprintf('%s ', data.Properties.VariableNames{:});
                fprintf('\n');
                continue;
            end

            % Select only expected columns
            data = data(:, expected_cols);

            % Append to all_data
            all_data = [all_data; data];

            fprintf('  Loaded %d rows\n', height(data));

        catch ME
            fprintf('  Error loading %s: %s\n', xlsx_files(i).name, ME.message);
            continue;
        end
    end

    if isempty(all_data)
        fprintf('\nNo valid data loaded. Exiting.\n');
        results = [];
        return;
    end

    fprintf('\nTotal rows combined: %d\n', height(all_data));

    % Get unique indices
    unique_indices = unique(all_data.index);
    fprintf('Unique indices found: %d\n\n', length(unique_indices));

    % Initialize results arrays
    n_indices = length(unique_indices);
    result_index = zeros(n_indices, 1);
    result_R1 = zeros(n_indices, 1);
    result_R1err = zeros(n_indices, 1);
    result_R2 = zeros(n_indices, 1);
    result_R2err = zeros(n_indices, 1);
    result_Rnoe = zeros(n_indices, 1);
    result_Rnoeerr = zeros(n_indices, 1);
    result_n_measurements = zeros(n_indices, 1);

    % Calculate weighted averages for each index
    for i = 1:n_indices
        idx = unique_indices(i);
        idx_mask = all_data.index == idx;
        idx_data = all_data(idx_mask, :);

        % Calculate weighted averages
        [r1_avg, r1_err] = calculate_weighted_average(...
            idx_data.R1, idx_data.R1err);
        [r2_avg, r2_err] = calculate_weighted_average(...
            idx_data.R2, idx_data.R2err);
        [rnoe_avg, rnoe_err] = calculate_weighted_average(...
            idx_data.Rnoe, idx_data.Rnoeerr);

        % Store results
        result_index(i) = idx;
        result_R1(i) = r1_avg;
        result_R1err(i) = r1_err;
        result_R2(i) = r2_avg;
        result_R2err(i) = r2_err;
        result_Rnoe(i) = rnoe_avg;
        result_Rnoeerr(i) = rnoe_err;
        result_n_measurements(i) = height(idx_data);
    end

    % Create results table
    results = table(result_index, result_R1, result_R1err, ...
                    result_R2, result_R2err, result_Rnoe, result_Rnoeerr, ...
                    result_n_measurements, ...
                    'VariableNames', {'index', 'R1', 'R1err', 'R2', 'R2err', ...
                                     'Rnoe', 'Rnoeerr', 'n_measurements'});

    % Sort by index
    results = sortrows(results, 'index');

    % Save to xlsx
    output_path = fullfile(directory, output_file);
    writetable(results, output_path);
    fprintf('Weighted averages saved to: %s\n\n', output_path);

    % Display summary
    fprintf('Summary of weighted averages:\n');
    disp(results);
end


function [weighted_avg, weighted_err] = calculate_weighted_average(values, errors)
% CALCULATE_WEIGHTED_AVERAGE Calculate weighted average using inverse variance
%
% Inputs:
%   values - Array of measurement values
%   errors - Array of measurement errors (weights)
%
% Outputs:
%   weighted_avg - Weighted average
%   weighted_err - Error on weighted average
%
% Method:
%   weighted_avg = sum(value / error^2) / sum(1 / error^2)
%   weighted_err = sqrt(1 / sum(1 / error^2))

    % Remove NaN values and zero errors
    mask = ~(isnan(values) | isnan(errors) | errors == 0);
    values_clean = values(mask);
    errors_clean = errors(mask);

    if isempty(values_clean)
        weighted_avg = NaN;
        weighted_err = NaN;
        return;
    end

    % Inverse variance weights: w = 1/sigma^2
    weights = 1.0 ./ (errors_clean .^ 2);

    % Weighted average
    weighted_avg = sum(values_clean .* weights) / sum(weights);

    % Error on weighted average
    weighted_err = sqrt(1.0 / sum(weights));
end
