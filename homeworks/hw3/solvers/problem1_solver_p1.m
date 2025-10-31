% problem1_solver_p1.m
% This script is saved in /hw3/solvers/ and run from /hw3/
% It solves Problem 1, Part 1 of PS3

clear; clc;

% --- Setup ---
fprintf('Running Problem 1, Part 1...\n');
n = 2000;
iterations = 1000;

% Create directory for results (path relative to /hw3)
results_dir = 'results/problem1';
output_file = fullfile(results_dir, 'p1_part1_timing.txt');
output_plot = fullfile(results_dir, 'p1_part1_timing_comparison.png');

% Generate standard random matrix and vector
A = randn(n, n);
b = randn(n, 1);

% --- (a) Gaussian Elimination (A\b) ---
fprintf('  Benchmarking (a) Gaussian Elimination (A\\b)...\n');
tic;
for i = 1:iterations
    x_gauss = A\b;
end
time_gauss = toc;

% --- (b) LU Decomposition ---
fprintf('  Benchmarking (b) LU Decomposition...\n');
tic;
[L, U, P] = lu(A); % Pre-compute factorization
for i = 1:iterations
    y_lu = L\(P*b);
    x_lu = U\y_lu;
end
time_lu = toc;

% --- (c) Matrix Inverse (inv(A)) ---
fprintf('  Benchmarking (c) Matrix Inverse...\n');
tic;
A_inv = inv(A); % Pre-compute inverse
for i = 1:iterations
    x_inv = A_inv*b;
end
time_inv = toc;

% --- (d) Cholesky Decomposition ---
% Cholesky requires a symmetric positive definite (SPD) matrix.
fprintf('  Generating SPD matrix for Cholesky...\n');
A_spd = A'*A; 
b_spd = randn(n, 1); % Create a new b for this system

fprintf('  Benchmarking (d) Cholesky Decomposition...\n');
tic;
U_chol = chol(A_spd); % Pre-compute factorization
for i = 1:iterations
    y_chol = U_chol'\b_spd; 
    x_chol = U_chol\y_chol;
end
time_chol = toc;

% --- Save Text Results ---
fid = fopen(output_file, 'w');
fprintf(fid, '--- Problem 1, Part 1: Direct Solver Benchmark ---\n');
fprintf(fid, 'Parameters: n = %d, Iterations = %d\n\n', n, iterations);
fprintf(fid, '(a) Gaussian (A\\b):    %.6f seconds\n', time_gauss);
fprintf(fid, '(b) LU Decomposition:   %.6f seconds\n', time_lu);
fprintf(fid, '(c) Matrix Inverse:     %.6f seconds\n', time_inv);
fprintf(fid, '(d) Cholesky (SPD):   %.6f seconds\n', time_chol);
fclose(fid);

fprintf('Text results saved to %s\n', output_file);

% --- Save Plot ---
times = [time_gauss, time_lu, time_inv, time_chol];
labels = {'A\b (Gauss)', 'LU Decomp.', 'inv(A)', 'Cholesky'};

fig = figure('Visible', 'off');
bar(times);
set(gca, 'xticklabel', labels);
ylabel('Total Time (seconds)');
title(['Direct Solver Benchmark (n=' num2str(n) ', ' num2str(iterations) ' iterations)']);
grid on;
saveas(fig, output_plot);

fprintf('Plot saved to %s\n', output_plot);
fprintf('Problem 1, Part 1 complete.\n');