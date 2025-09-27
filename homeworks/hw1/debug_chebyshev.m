%% ========================================================================
% Diagnostic Script for 2D Chebyshev Approximation
% Author: Ebrahim Ebrami
% Date: Fall 2025
%
% This script performs a step-by-step implementation of the 2D Chebyshev
% approximation method to diagnose a persistent bug. Each step's output
% is printed to the console for verification against theory.
% =========================================================================

clear; clc; close all;

fprintf('=== DIAGNOSTIC SCRIPT: 2D CHEBYSHEV APPROXIMATION ===\n\n');

%% SETUP: Define Test Function, Degree, and Grids
test_func = @(x,y) x.^2 + 2*y.^2 - 3*x.*y;
N = 4; % Degree of polynomial approximation.
m = N + 1;
TOLERANCE = 1e-12;

% Fine grid for final error evaluation
n_fine = 21;
grid_fine_1d = linspace(0, 1, n_fine);
[x_fine, y_fine] = ndgrid(grid_fine_1d, grid_fine_1d);
z_true = test_func(x_fine, y_fine);


%% STEP 1: Generate Chebyshev Nodes
% From theory (Eq. 13.30), nodes on [-1,1] are z_j = cos(j*pi/N).
% These are then mapped to the problem domain [0,1] using z = 2x - 1.
fprintf('--- STEP 1: Node Generation (N=%d) ---\n', N);
k = (0:N)';
z_nodes_1d = cos(k * pi / N); % Canonical nodes on [-1, 1]
x_nodes_1d = 0.5 * (1 + z_nodes_1d); % Mapped to [0, 1] (Note: 1-z or 1+z depends on z ordering)
[x_nodes, y_nodes] = ndgrid(x_nodes_1d, x_nodes_1d);

fprintf('Nodes on canonical domain [-1, 1]:\n');
disp(z_nodes_1d');
fprintf('Nodes on problem domain [0, 1]:\n');
disp(x_nodes_1d');


%% STEP 2: Evaluate Function at Nodes
% The true function is evaluated at the grid points generated in Step 1.
fprintf('\n--- STEP 2: Evaluate Function at Chebyshev Nodes ---\n');
z_at_nodes = test_func(x_nodes, y_nodes);

fprintf('Function values on the %dx%d grid of nodes:\n', m, m);
disp(z_at_nodes);


%% STEP 3: Calculate Chebyshev Coefficients
% This is the most complex step and the source of the bug.
% The formula is C = (4/N^2) * T * F_scaled * T', where T is the basis matrix.
% A final scaling of the boundary coefficients is also required.
fprintf('\n--- STEP 3: Calculate Coefficients ---\n');

% 3a: Construct the 1D Chebyshev basis matrix, T_ij = T_i(z_j)
T_basis = cos( (0:N)' * acos(z_nodes_1d') );
fprintf('Basis Matrix T (T_ij = T_i(z_j)) is %dx%d:\n', size(T_basis,1), size(T_basis,2));
% disp(T_basis); % This is large, uncomment to view

% 3b: Scale function values at the boundaries (for the sum)
F_scaled = z_at_nodes;
F_scaled([1, m], :) = F_scaled([1, m], :) * 0.5;
F_scaled(:, [1, m]) = F_scaled(:, [1, m]) * 0.5;
fprintf('Function values have been scaled at the boundaries.\n');

% 3c: Calculate the intermediate coefficients using the matrix formula
coeffs_intermediate = (4 / (N * N)) * (T_basis * F_scaled * T_basis');

% --- BUG FIX ---
% 3d: The coefficients themselves must ALSO be scaled at the boundaries.
% The reconstruction formula requires the first and last terms of the sum
% to be halved, which is equivalent to halving the coefficients a_0 and a_N.
coeffs_final = coeffs_intermediate;
coeffs_final([1, m], :) = coeffs_final([1, m], :) * 0.5;
coeffs_final(:, [1, m]) = coeffs_final(:, [1, m]) * 0.5;
fprintf('*** CRITICAL FIX APPLIED: Boundary coefficients have been scaled. ***\n');
fprintf('Final Coefficient Matrix C (%dx%d):\n', m, m);
disp(coeffs_final);


%% STEP 4: Evaluate the Approximation on the Fine Grid
% The polynomial is reconstructed using the formula: Z = T_eval * C * T_eval'
fprintf('\n--- STEP 4: Evaluate Approximation on Fine Grid ---\n');

% 4a: Get unique evaluation points and map to [-1, 1]
x_vec_fine = x_fine(:,1);
y_vec_fine = y_fine(1,:);
z_x_fine = 2 * x_vec_fine - 1;
z_y_fine = 2 * y_vec_fine' - 1;

% 4b: Build the evaluation basis matrices
Tx_eval = cos( acos(z_x_fine) * (0:N) );
Ty_eval = cos( acos(z_y_fine) * (0:N) );
fprintf('Evaluation basis matrices Tx (%dx%d) and Ty (%dx%d) created.\n',...
    size(Tx_eval,1), size(Tx_eval,2), size(Ty_eval,1), size(Ty_eval,2));

% 4c: Perform the evaluation using the final, scaled coefficients
% The reconstruction formula is f(x,y) = sum( C(i,j) * T_i(x) * T_j(y) ),
% where the sum is double-primed. Our scaling in step 3d accounts for this.
z_approx = Tx_eval * coeffs_final * Ty_eval';


%% STEP 5: Final Error Calculation
fprintf('\n--- STEP 5: Calculate Final Error ---\n');
error = abs(z_true - z_approx);
max_err = max(error(:));

fprintf('\nMaximum absolute error: %.4e\n', max_err);
if max_err < TOLERANCE
    fprintf('RESULT: PASS. The Chebyshev approximation is now correct.\n');
else
    fprintf('RESULT: FAIL. Error still exceeds tolerance of %.1e.\n', TOLERANCE);
end

fprintf('\n=== DIAGNOSTIC COMPLETE ===\n');