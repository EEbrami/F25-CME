%% ========================================================================
% Test Script for Problem 1 Helper Functions
% Author: Ebrahim Ebrami
% Date: Fall 2025
%
% This script tests the helper functions for 2D approximation, including
% Chebyshev and polynomial evaluations. It verifies the accuracy of the
% approximations against a known test function:
%   - cheb_eval_2d.m
%   - poly_eval_2d.m
% =========================================================================

clear; clc; close all;
addpath('functions');
fprintf('=== Testing Helper Functions for 2D Approximation ===\n\n');

%% 1. Setup
test_func = @(x,y) x.^2 + 2*y.^2 - 3*x.*y;
n_fine = 101;
grid_fine_1d = linspace(0, 1, n_fine);
[x_fine, y_fine] = ndgrid(grid_fine_1d, grid_fine_1d);
z_true = test_func(x_fine, y_fine);
TOLERANCE = 1e-12;

%% ========================================================================
% 2. Test: `cheb_eval_2d.m`
% =========================================================================
fprintf('--- Testing `cheb_eval_2d` ---\n');

% Parameters
N_cheb = 4;
m_cheb = N_cheb + 1;

% Grid Generation
k = (0:N_cheb)';
nodes_1d_cheb = 0.5 * (1 + cos(k * pi / N_cheb)); % Mapped to [0,1]
[x_nodes_cheb, y_nodes_cheb] = ndgrid(nodes_1d_cheb, nodes_1d_cheb);
z_nodes_cheb = test_func(x_nodes_cheb, y_nodes_cheb);

% --- BUG FIX ---
% Replace the entire coefficient calculation block with the proven logic
% from the debug script.
z_domain = cos(k * pi / N_cheb); % Canonical nodes on [-1, 1]

% Construct the basis matrix using the direct definition
T_basis = cos( (0:N_cheb)' * acos(z_domain') );

% Scale function values at the boundaries
F_scaled = z_nodes_cheb;
F_scaled([1, m_cheb], :) = F_scaled([1, m_cheb], :) * 0.5;
F_scaled(:, [1, m_cheb]) = F_scaled(:, [1, m_cheb]) * 0.5;

% Calculate intermediate coefficients
coeffs_cheb = (4 / (N_cheb * N_cheb)) * (T_basis * F_scaled * T_basis');

% Scale the boundary coefficients
coeffs_cheb([1, m_cheb], :) = coeffs_cheb([1, m_cheb], :) * 0.5;
coeffs_cheb(:, [1, m_cheb]) = coeffs_cheb(:, [1, m_cheb]) * 0.5;
% --- END OF FIX ---

% Evaluation (This will now receive the correct coefficients)
z_approx_cheb = cheb_eval_2d(coeffs_cheb, x_fine, y_fine);

% Verification
error_cheb = abs(z_true - z_approx_cheb);
max_err_cheb = max(error_cheb(:));
fprintf('Maximum absolute error for Chebyshev approximation: %.4e\n', max_err_cheb);

if max_err_cheb < TOLERANCE
    fprintf('RESULT: PASS. The `cheb_eval_2d` function is working correctly.\n\n');
else
    fprintf('RESULT: FAIL. The error exceeds the tolerance of %.1e.\n\n', TOLERANCE);
end

%% ========================================================================
% 3. Test: `poly_eval_2d.m`
% This part was already correct and remains unchanged.
% =========================================================================
fprintf('--- Testing `poly_eval_2d` ---\n');

N_poly = 2; m_poly = N_poly + 1;
nodes_1d_poly = linspace(0, 1, m_poly);
[x_nodes_poly, y_nodes_poly] = ndgrid(nodes_1d_poly, nodes_1d_poly);
z_nodes_poly = test_func(x_nodes_poly, y_nodes_poly);
X_mat = zeros(m_poly^2, m_poly^2);
y_vec = reshape(z_nodes_poly', [], 1);
idx = 1;
for i = 1:m_poly
    for j = 1:m_poly
        px = x_nodes_poly(i,j); py = y_nodes_poly(i,j);
        col_idx = 1;
        for p = 0:N_poly
            for q = 0:N_poly
                X_mat(idx, col_idx) = px^p * py^q;
                col_idx = col_idx + 1;
            end
        end
        idx = idx + 1;
    end
end
coeffs_poly = X_mat \ y_vec;
z_approx_poly = poly_eval_2d(coeffs_poly, N_poly, x_fine, y_fine);
error_poly = abs(z_true - z_approx_poly);
max_err_poly = max(error_poly(:));
fprintf('Maximum absolute error for Ordinary Polynomial approximation: %.4e\n', max_err_poly);
if max_err_poly < TOLERANCE
    fprintf('RESULT: PASS. The `poly_eval_2d` function is working correctly.\n\n');
else
    fprintf('RESULT: FAIL. The error exceeds the tolerance of %.1e.\n\n', TOLERANCE);
end

fprintf('=== Testing Complete ===\n');