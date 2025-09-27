%% ========================================================================
% Test Script for Function Approximation Helper Functions
% Author: Problem Set 1 Solution
% Date: Fall 2025
% This script tests the correctness of cheb_eval_2d and poly_eval_2d
% =========================================================================

clear; clc; close all;

% Add functions directory to path
addpath('functions');

fprintf('=== Testing Helper Functions ===\n\n');

%% Test 1: Simple function evaluation
fprintf('Test 1: Simple function evaluation...\n');

% Define a simple test function: f(x,y) = x^2 + 2*x*y + y^2
test_func = @(x,y) x.^2 + 2*x.*y + y.^2;

% Test points
x_test = [0, 0.5, 1];
y_test = [0, 0.5, 1];
[X_test, Y_test] = meshgrid(x_test, y_test);
true_values = test_func(X_test, Y_test);

fprintf('True function values at test points:\n');
disp(true_values);

%% Test 2: Polynomial evaluation with known coefficients
fprintf('\nTest 2: Polynomial evaluation with known coefficients...\n');

% For f(x,y) = x^2 + 2*x*y + y^2, the coefficients should be:
% [0, 0, 1, 0, 2, 0, 1, 0, 0] for degree 2
% Ordering: [c_00, c_01, c_02, c_10, c_11, c_12, c_20, c_21, c_22]
coeffs_known = [0; 0; 1; 0; 2; 0; 1; 0; 0];

% Test our polynomial evaluation
poly_values = poly_eval_2d(coeffs_known, 2, X_test, Y_test);
poly_error = abs(true_values - poly_values);

fprintf('Polynomial approximation values:\n');
disp(poly_values);
fprintf('Polynomial approximation error (should be ~0):\n');
disp(poly_error);
fprintf('Max polynomial error: %.2e\n', max(poly_error(:)));

%% Test 3: Chebyshev approximation consistency
fprintf('\nTest 3: Chebyshev approximation consistency...\n');

% Test with a smooth function: f(x,y) = sin(2*pi*x) * cos(2*pi*y)
smooth_func = @(x,y) sin(2*pi*x) .* cos(2*pi*y);

% Create Chebyshev grid
N = 7;
k = (0:N)';
nodes_1d = 0.5 * (1 - cos(k * pi / N));
[x_nodes, y_nodes] = ndgrid(nodes_1d, nodes_1d);
z_nodes = smooth_func(x_nodes, y_nodes);

% Calculate Chebyshev coefficients
m = N + 1;
z_domain = cos(k * pi / N);
T_Km = cos((0:N) .* acos(z_domain));
c = ones(m, 1); c([1, m]) = 2;
D = diag(2 ./ (c * m));
D(1,1) = D(1,1)/2;
F = z_nodes;
Gamma = D * T_Km' * F * T_Km * D;

% Test evaluation at the nodes (should be exact)
z_cheb_nodes = cheb_eval_2d(Gamma, x_nodes, y_nodes);
node_error = abs(z_nodes - z_cheb_nodes);

fprintf('Chebyshev interpolation error at nodes (should be ~0):\n');
fprintf('Max error at nodes: %.2e\n', max(node_error(:)));

% Test evaluation at intermediate points
x_test_fine = linspace(0, 1, 50);
y_test_fine = linspace(0, 1, 50);
[X_fine, Y_fine] = meshgrid(x_test_fine, y_test_fine);
z_true_fine = smooth_func(X_fine, Y_fine);
z_cheb_fine = cheb_eval_2d(Gamma, X_fine, Y_fine);
interp_error = abs(z_true_fine - z_cheb_fine);

fprintf('Chebyshev interpolation error at intermediate points:\n');
fprintf('Max error: %.2e, Mean error: %.2e\n', ...
        max(interp_error(:)), mean(interp_error(:)));

%% Test 4: Function dimension compatibility
fprintf('\nTest 4: Function dimension compatibility...\n');

% Test with different input dimensions
x_scalar = 0.5;
y_scalar = 0.3;
coeffs_test = ones(4,4);  % 4x4 coefficient matrix

try
    z_scalar = cheb_eval_2d(coeffs_test, x_scalar, y_scalar);
    fprintf('Scalar input test: PASSED (z = %.4f)\n', z_scalar);
catch e
    fprintf('Scalar input test: FAILED (%s)\n', e.message);
end

% Test with vector input
x_vec = [0.2, 0.7];
y_vec = [0.3, 0.8];
try
    z_vec = cheb_eval_2d(coeffs_test, x_vec, y_vec);
    fprintf('Vector input test: PASSED (size = %dx%d)\n', size(z_vec,1), size(z_vec,2));
catch e
    fprintf('Vector input test: FAILED (%s)\n', e.message);
end

%% Summary
fprintf('\n=== Test Summary ===\n');
fprintf('1. Simple function evaluation: Manual verification required\n');
if max(poly_error(:)) < 1e-10
    fprintf('2. Polynomial evaluation: PASSED (error < 1e-10)\n');
else
    fprintf('2. Polynomial evaluation: FAILED (error = %.2e)\n', max(poly_error(:)));
end
if max(node_error(:)) < 1e-10
    fprintf('3. Chebyshev node interpolation: PASSED (error < 1e-10)\n');
else
    fprintf('3. Chebyshev node interpolation: FAILED (error = %.2e)\n', max(node_error(:)));
end
if max(interp_error(:)) < 0.1  % Less strict for intermediate points
    fprintf('4. Chebyshev intermediate interpolation: PASSED (reasonable error)\n');
else
    fprintf('4. Chebyshev intermediate interpolation: WARNING (large error = %.2e)\n', max(interp_error(:)));
end

fprintf('\nFunction tests complete!\n');