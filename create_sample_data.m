function create_sample_data(output_dir, n_files)
% CREATE_SAMPLE_DATA Create sample xlsx files for testing weighted average
%
% Usage:
%   create_sample_data()
%   create_sample_data(output_dir)
%   create_sample_data(output_dir, n_files)
%
% Inputs:
%   output_dir - Directory to save sample files (default: 'sample_data_matlab')
%   n_files    - Number of sample files to create (default: 3)
%
% Description:
%   Creates sample xlsx files with the expected format for testing
%   the weighted_average_xlsx function. Files contain realistic NMR
%   relaxation data with noise.

    % Default parameters
    if nargin < 1 || isempty(output_dir)
        output_dir = 'sample_data_matlab';
    end

    if nargin < 2 || isempty(n_files)
        n_files = 3;
    end

    % Create output directory
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Set random seed for reproducibility
    rng(42);

    % Common indices for all files
    indices = (1:10)';

    for file_num = 1:n_files
        % Initialize data arrays
        n_rows = length(indices);
        R1_data = zeros(n_rows, 1);
        R1err_data = zeros(n_rows, 1);
        R2_data = zeros(n_rows, 1);
        R2err_data = zeros(n_rows, 1);
        Rnoe_data = zeros(n_rows, 1);
        Rnoeerr_data = zeros(n_rows, 1);

        for i = 1:n_rows
            % Base values (typical for protein NMR)
            r1_base = 1.5 + 0.3 * randn();  % ~1.5 s^-1
            r2_base = 10.0 + 2.0 * randn();  % ~10 s^-1
            rnoe_base = 0.7 + 0.1 * randn();  % ~0.7

            % Add measurement noise
            R1_data(i) = r1_base + 0.05 * randn();
            R2_data(i) = r2_base + 0.5 * randn();
            Rnoe_data(i) = rnoe_base + 0.02 * randn();

            % Errors (smaller for more precise measurements)
            R1err_data(i) = 0.02 + 0.01 * rand();
            R2err_data(i) = 0.2 + 0.1 * rand();
            Rnoeerr_data(i) = 0.01 + 0.005 * rand();
        end

        % Create table
        data = table(indices, R1_data, R1err_data, R2_data, R2err_data, ...
                     Rnoe_data, Rnoeerr_data, ...
                     'VariableNames', {'index', 'R1', 'R1err', 'R2', ...
                                      'R2err', 'Rnoe', 'Rnoeerr'});

        % Save to xlsx
        output_file = fullfile(output_dir, sprintf('sample_data_%d.xlsx', file_num));
        writetable(data, output_file);
        fprintf('Created: %s\n', output_file);
    end

    fprintf('\nCreated %d sample xlsx files in ''%s/'' directory\n', ...
            n_files, output_dir);
end
