% problem2_solver.m
% This script runs the experiments for Problem 2 of PS3.
% It calls a modified version of Main_GSSA_1.m (named run_gssa_experiment.m)
% to test various numerical stabilization methods.

clc;
clear all;
close all;

% Add paths to the original code and our new modified method
% This path order ensures our 'methods' folder is checked first,
% though under the new assumption, it is not critical for casing.
addpath('methods/');
addpath('gssa_1_agent_capital/');

disp('Starting Problem 2 GSSA experiments...');
total_tic = tic;

% Define all experimental cases based on PS3 and Num_Stab_Approx.m
% RM maps to Num_Stab_Approx.m:
%   2: LS-SVD
%   3: LAD-PP
%   4: LAD-DP
%   5: RLS-Tikhonov
%   6: RLS-TSVD
%   7: RLAD-PP
%   8: RLAD-DP
%
% 'penalty' is interpreted as:
%   - For Tikhonov/RLAD (RM 5, 7, 8): The exponent for 'eta' (e.g., -7 for 10^-7)
%   - For TSVD (RM 6): The exponent for 'kappa' (e.g., 7 for 10^7)
%   - For others (RM 2, 3, 4): Not used (can set to 0 or NaN)
%
% 'normalize' is 1 (on) or 0 (off).

experiments = {
    % Case 1: Compare LAD Methods
    'LAD-PP (Base)',              3, 0,    1;
    'LAD-DP (Base)',              4, 0,    1;
    'RLAD-PP (eta=10^-7)',        7, -7,   1;
    'RLAD-DP (eta=10^-7)',        8, -7,   1;
    'RLAD-PP (eta=10^-4)',        7, -4,   1;
    'RLAD-DP (eta=10^-4)',        8, -4,   1;

    % Case 2: Compare Tikhonov and SVD
    'RLS-Tikhonov (eta=10^-7)',   5, -7,   1;
    'RLS-Tikhonov (eta=10^-4)',   5, -4,   1;
    'LS-SVD (Base)',              2, 0,    1;

    % Case 3: Compare SVD and Truncated SVD (TSVD)
    'RLS-TSVD (kappa=10^6)',      6, 6,    1;
    'RLS-TSVD (kappa=10^7)',      6, 7,    1; % Default from main_gssa_1
    'RLS-TSVD (kappa=10^8)',      6, 8,    1;

    % Case 4: Role of Normalization
    'LAD-DP (No Norm)',           4, 0,    0;
    'RLS-Tikhonov (No Norm)',     5, -7,   0;
    'RLS-TSVD (No Norm)',         6, 7,    0
};

% Initialize results table
num_experiments = size(experiments, 1);
results = cell(num_experiments + 1, 5);
results(1,:) = {'Case Description', 'CPU Time (s)', 'Mean Error (log10)', 'Max Error (log10)', 'Parameters (RM, penalty, norm)'};

% Run the experiment loop
for i = 1:num_experiments
    desc = experiments{i, 1};
    rm_val = experiments{i, 2};
    penalty_val = experiments{i, 3};
    norm_val = experiments{i, 4};
    
    param_str = sprintf('RM=%d, penalty=%d, norm=%d', rm_val, penalty_val, norm_val);
    fprintf('Running case %d/%d: %s...\n', i, num_experiments, desc);
    
    try
        case_tic = tic;
        % Call our modified function
        [time, e_mean, e_max] = run_gssa_experiment(rm_val, penalty_val, norm_val);
        case_time = toc(case_tic);
        
        % Store results
        results{i+1, 1} = desc;
        results{i+1, 2} = case_time;
        results{i+1, 3} = e_mean;
        results{i+1, 4} = e_max;
        results{i+1, 5} = param_str;
        
    catch ME
        % Handle potential errors during a run
        fprintf('--- ERROR in case %s: %s ---\n', desc, ME.message);
        results{i+1, 1} = [desc ' (FAILED)'];
        results{i+1, 2} = toc(case_tic);
        results{i+1, 3} = NaN;
        results{i+1, 4} = NaN;
        results{i+1, 5} = param_str;
    end
end

total_time = toc(total_tic);
fprintf('\n--- All experiments complete. Total time: %.2f seconds. ---\n\n', total_time);

% Display the final results table
results_table = cell2table(results(2:end,:), 'VariableNames', results(1,:));
disp(results_table);

% Save results to a file
% Ensure the 'results/problem2' directory exists
if ~exist('results/problem2', 'dir')
   mkdir('results/problem2');
end
writetable(results_table, 'results/problem2/p2_stabilization_results.txt', 'Delimiter', '\t');
writetable(results_table, 'results/problem2/p2_stabilization_results.csv');

disp('Results saved to homeworks/hw3/results/problem2/');