%% ========================================================================
% Master Script for ECON-81360 - Problem Set 1
% Author: Ebrahim Ebrami
% Date: Fall 2025
%
% This script runs all parts of the problem set 1 sequentially and saves
% the Command Window output to a log file using the 'diary' command.
% =========================================================================

clear; clc; close all;

% Create a directory for the final results if it doesn't exist
if ~exist('results', 'dir')
    mkdir('results');
end

% Define the log file path
log_filepath = 'results/problem_set_1_log.txt';

% If an old log file exists, delete it to ensure a clean run
if exist(log_filepath, 'file')
    delete(log_filepath);
end

% Start logging to the text file.
diary(log_filepath);

fprintf('============================================================\n');
fprintf('RUNNING ECON-81360 - PROBLEM SET 1\n');
fprintf('Author: Ebrahim Ebrami\n');
fprintf('Execution started at: %s\n', string(datetime('now')));
fprintf('============================================================\n');

% --- Run Problem 1 ---
fprintf('\n--- Starting Problem 1: Function Approximation ---\n\n');
tic;
main_problem1_approximation;
% --- BUG FIX: Added 'toc' variable to the fprintf command ---
fprintf('\n--- Problem 1 Finished. Elapsed time: %.2f seconds ---\n', toc);

% --- Run Problem 2 ---
fprintf('\n\n--- Starting Problem 2: Stochastic Growth Model ---\n\n');
tic;
main_problem2_growth_model;
% --- BUG FIX: Added 'toc' variable to the fprintf command ---
fprintf('\n--- Problem 2 Finished. Elapsed time: %.2f seconds ---\n', toc);

fprintf('\n============================================================\n');
fprintf('EXECUTION COMPLETE\n');
fprintf('============================================================\n');

% Stop logging
diary off;

% Print the final confirmation message
fprintf('\nLog file successfully generated: results/problem_set_1_log.txt\n');

% Close all open figure windows for a clean exit.
close all;