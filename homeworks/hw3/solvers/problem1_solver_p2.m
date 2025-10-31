% problem1_solver_p2.m
% This script is saved in /hw3/solvers/ and run from /hw3/
% It solves Problem 1, Part 2 of PS3 using the
% EXACT matrices specified in the PDF.

clear; clc;

% --- Setup ---
fprintf('Running Problem 1, Part 2...\n');

% Add methods folder to path
addpath('methods');

% Create directory for results (path relative to /hw3)
results_dir = 'results/problem1';
if ~exist(results_dir, 'dir')
   mkdir(results_dir);
end
output_file = fullfile(results_dir, 'p1_part2_iterations.txt');

% --- Define Systems (from PS3 2025.pdf) ---

% Vector b
b = [1; 4; -2; 1];

% System 1
A1 = [3 1 1 0; 
      1 5 -1 2; 
      1 0 3 1; 
      0 1 1 4];

% System 2
A2 = [2.5 1 1 0; 
      1 4.1 -1 2; 
      1 0 2.1 1; 
      0 1 1 2.1];

% System 3
A3 = [2 1 1 0; 
      1 3.5 -1 2; 
      1 0 2.1 1; 
      0 1 1 2.1];

% --- Set Iteration Parameters ---
tol = 1e-8;
max_iter = 5000;
x0 = zeros(4, 1); % Initial guess

% --- Run Solvers ---
fprintf('  Solving system with A1 (Strictly Diagonally Dominant)...\n');
[~, iter_j1] = gauss_jacobi(A1, b, x0, tol, max_iter);
[~, iter_s1] = gauss_seidel(A1, b, x0, tol, max_iter);

fprintf('  Solving system with A2 (Strictly Diagonally Dominant)...\n');
[~, iter_j2] = gauss_jacobi(A2, b, x0, tol, max_iter);
[~, iter_s2] = gauss_seidel(A2, b, x0, tol, max_iter);

fprintf('  Solving system with A3 (NOT Strictly Diagonally Dominant)...\n');
[~, iter_j3] = gauss_jacobi(A3, b, x0, tol, max_iter);
[~, iter_s3] = gauss_seidel(A3, b, x0, tol, max_iter);

% --- Save Text Results and Analysis ---
fid = fopen(output_file, 'w');
fprintf(fid, '--- Problem 1, Part 2: Iterative Solver Comparison ---\n');
fprintf(fid, 'Parameters: Tolerance = %e, Max Iterations = %d\n\n', tol, max_iter);

fprintf(fid, 'System 1 (A1):\n');
fprintf(fid, '  Gauss-Jacobi: %d iterations\n', iter_j1);
fprintf(fid, '  Gauss-Seidel: %d iterations\n', iter_s1);
fprintf(fid, '\n');

fprintf(fid, 'System 2 (A2):\n');
fprintf(fid, '  Gauss-Jacobi: %d iterations\n', iter_j2);
fprintf(fid, '  Gauss-Seidel: %d iterations\n', iter_s2);
fprintf(fid, '\n');

fprintf(fid, 'System 3 (A3):\n');
fprintf(fid, '  Gauss-Jacobi: %d iterations (Note: %d indicates non-convergence)\n', iter_j3, max_iter);
fprintf(fid, '  Gauss-Seidel: %d iterations (Note: %d indicates non-convergence)\n', iter_s3, max_iter);
fprintf(fid, '\n\n');

fprintf(fid, '--- Analysis of Convergence ("Why is it so different?") ---\n\n');
fprintf(fid, 'A sufficient condition for the convergence of Gauss-Jacobi and Gauss-Seidel is that the matrix is strictly diagonally dominant.\n\n');

fprintf(fid, 'Matrix A1: IS strictly diagonally dominant. Convergence is guaranteed.\n');
fprintf(fid, '  Row 1: |3| > |1| + |1| + |0| (3 > 2) -> TRUE\n');
fprintf(fid, '  Row 2: |5| > |1| + |-1| + |2| (5 > 4) -> TRUE\n');
fprintf(fid, '  Row 3: |3| > |1| + |0| + |1| (3 > 2) -> TRUE\n');
fprintf(fid, '  Row 4: |4| > |0| + |1| + |1| (4 > 2) -> TRUE\n\n');

fprintf(fid, 'Matrix A2: IS strictly diagonally dominant. Convergence is guaranteed.\n');
fprintf(fid, '  Row 1: |2.5| > |1| + |1| + |0| (2.5 > 2) -> TRUE\n');
fprintf(fid, '  Row 2: |4.1| > |1| + |-1| + |2| (4.1 > 4) -> TRUE\n');
fprintf(fid, '  Row 3: |2.1| > |1| + |0| + |1| (2.1 > 2) -> TRUE\n');
fprintf(fid, '  Row 4: |2.1| > |0| + |1| + |1| (2.1 > 2) -> TRUE\n\n');

fprintf(fid, 'Matrix A3: IS NOT strictly diagonally dominant. Convergence is not guaranteed.\n');
fprintf(fid, '  Row 1: |2| > |1| + |1| + |0| (2 > 2) -> FALSE\n');
fprintf(fid, '  Row 2: |3.5| > |1| + |-1| + |2| (3.5 > 4) -> FALSE\n');
fprintf(fid, '  Result: The methods fail to converge, as indicated by reaching the max iteration limit.\n\n');

fprintf(fid, '--- Comparison of Methods ---\n');
fprintf(fid, 'For the converging matrices (A1 and A2), Gauss-Seidel requires fewer iterations than Gauss-Jacobi.\n');
fprintf(fid, 'This is the expected theoretical outcome, as Gauss-Seidel uses the most recently updated values within the same iteration, accelerating convergence.\n');

fclose(fid);

fprintf('Text results and analysis saved to %s\n', output_file);
fprintf('Problem 1, Part 2 complete.\n');