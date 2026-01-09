% EXAMPLE_WEIGHTED_AVERAGE - Example script for weighted average calculation
%
% This script demonstrates how to use the weighted_average_xlsx function
% to process multiple xlsx files and calculate weighted averages.
%
% The script will:
% 1. Create sample data files
% 2. Process the files to calculate weighted averages
% 3. Display the results

%% Clear workspace
clear all;
close all;
clc;

%% Step 1: Create sample data
fprintf('========================================\n');
fprintf('Step 1: Creating sample data\n');
fprintf('========================================\n\n');

% Create sample data in a MATLAB-specific directory
output_dir = 'sample_data_matlab';
n_files = 3;

create_sample_data(output_dir, n_files);

%% Step 2: Calculate weighted averages
fprintf('\n========================================\n');
fprintf('Step 2: Calculating weighted averages\n');
fprintf('========================================\n\n');

% Process the sample data
results = weighted_average_xlsx(output_dir, 'weighted_averages_matlab.xlsx');

%% Step 3: Display and analyze results
fprintf('\n========================================\n');
fprintf('Step 3: Analysis\n');
fprintf('========================================\n\n');

if ~isempty(results)
    fprintf('Number of unique indices: %d\n', height(results));
    fprintf('Total measurements per index: %d\n', results.n_measurements(1));

    fprintf('\nMean values across all indices:\n');
    fprintf('  R1:   %.4f ± %.4f\n', mean(results.R1), std(results.R1));
    fprintf('  R2:   %.4f ± %.4f\n', mean(results.R2), std(results.R2));
    fprintf('  Rnoe: %.4f ± %.4f\n', mean(results.Rnoe), std(results.Rnoe));

    fprintf('\nMean weighted errors:\n');
    fprintf('  R1err:   %.4f\n', mean(results.R1err));
    fprintf('  R2err:   %.4f\n', mean(results.R2err));
    fprintf('  Rnoeerr: %.4f\n', mean(results.Rnoeerr));

    %% Optional: Create plots
    if exist('results', 'var') && height(results) > 0
        figure('Position', [100, 100, 1200, 400]);

        % Plot R1
        subplot(1, 3, 1);
        errorbar(results.index, results.R1, results.R1err, 'o-', 'LineWidth', 1.5);
        xlabel('Index');
        ylabel('R1 (s^{-1})');
        title('R1 Relaxation Rate');
        grid on;

        % Plot R2
        subplot(1, 3, 2);
        errorbar(results.index, results.R2, results.R2err, 'o-', 'LineWidth', 1.5);
        xlabel('Index');
        ylabel('R2 (s^{-1})');
        title('R2 Relaxation Rate');
        grid on;

        % Plot Rnoe
        subplot(1, 3, 3);
        errorbar(results.index, results.Rnoe, results.Rnoeerr, 'o-', 'LineWidth', 1.5);
        xlabel('Index');
        ylabel('R_{NOE}');
        title('NOE');
        grid on;

        sgtitle('Weighted Average Results');

        fprintf('\nPlots created.\n');
    end
else
    fprintf('No results to display.\n');
end

fprintf('\n========================================\n');
fprintf('Example completed!\n');
fprintf('========================================\n');
