%% ========================================================================
% Test Script for Problem 2 Helper Functions
% Author: Ebrahim
% Date: Fall 2025
%
% This script verifies the correctness of the key approximation functions
% used in solving the stochastic growth model:
%   - create_grid.m
%   - chebyshev_basis.m
%   - chebyshev_eval_2d_p2.m (the improved version)
% =========================================================================

clear; clc; close all;
addpath('p2_functions'); % Make sure this matches your folder name

fprintf('=== Testing Helper Functions for Problem 2 ===\n\n');

%% 1. Setup
% Define a simple known polynomial function to approximate
test_func = @(k,z) k.^2 - 2*k.*z + 3*z;

% Define an arbitrary domain for testing
domain.k_min = 0.8;
domain.k_max = 1.2;
domain.z_min = 0.9;
domain.z_max = 1.1;

% Parameters for the test
N_k = 10; % Number of nodes for k
N_z = 8;  % Number of nodes for z
cheb_deg = 4; % Degree of approximation

TOLERANCE = 1e-12;

%% ========================================================================
% 2. Test the Combination of Chebyshev Functions
% This test implicitly verifies create_grid, chebyshev_basis, and
% chebyshev_eval_2d_p2 by checking if they can perfectly replicate a known
% polynomial, which they should.
% =========================================================================
fprintf('--- Testing Chebyshev Approximation Logic ---\n');

% 2.1. Create Chebyshev nodes on the specified domain
k = (0:N_k-1)';
z = (0:N_z-1)';
k_nodes = domain.k_min + (domain.k_max-domain.k_min)/2 * (1 - cos(k*pi/(N_k-1)));
z_nodes = domain.z_min + (domain.z_max-domain.z_min)/2 * (1 - cos(z*pi/(N_z-1)));

% 2.2. Evaluate the true function on the grid of nodes
[k_grid_nodes, z_grid_nodes] = ndgrid(k_nodes, z_nodes);
z_at_nodes = test_func(k_grid_nodes, z_grid_nodes);

% 2.3. Construct the basis matrix using the function to be tested
basis_mat = chebyshev_basis(k_nodes, z_nodes, cheb_deg);

% 2.4. Solve for the coefficients via regression (OLS)
% These coefficients should allow us to perfectly reconstruct the function
coeffs = basis_mat \ z_at_nodes(:);
fprintf('Calculated %d Chebyshev coefficients.\n', length(coeffs));

% 2.5. Evaluate the approximation on a new, finer grid
n_fine = 51;
k_fine_vec = linspace(domain.k_min, domain.k_max, n_fine)';
z_fine_vec = linspace(domain.z_min, domain.z_max, n_fine);
[k_fine_grid, z_fine_grid] = ndgrid(k_fine_vec, z_fine_vec);

% Use the function we are testing
z_approx = chebyshev_eval_2d_p2(coeffs, cheb_deg, k_fine_grid, z_fine_grid, domain);

% 2.6. Calculate the true function values on the fine grid for comparison
z_true = test_func(k_fine_grid, z_fine_grid);

% 2.7. Verify the result
max_err = max(abs(z_true - z_approx), [], 'all');
fprintf('Maximum absolute error for Chebyshev approximation: %.4e\n', max_err);

if max_err < TOLERANCE
    fprintf('RESULT: PASS. The Chebyshev functions are working correctly together.\n\n');
else
    fprintf('RESULT: FAIL. The error exceeds the tolerance of %.1e.\n\n', TOLERANCE);
end

fprintf('=== Testing Complete ===\n');