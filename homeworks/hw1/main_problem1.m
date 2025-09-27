%% Problem Set 1: 2D Function Approximation Using the Franke Function
% Course: Computational Methods for Economists, Econ-81360
% Author: Problem Set 1 Solution
% Date: Fall 2025

%% Clear workspace and create results directory
clear; clc; close all;

% Create results directory structure
if ~exist('results/figures', 'dir')
    mkdir('results/figures');
end

%% Define the Franke function
franke = @(x, y) 0.75 * exp(-((9*x-2).^2 + (9*y-2).^2)/4) + ...
                 0.75 * exp(-((9*x+1).^2)/49 - (9*y+1)/10) + ...
                 0.5 * exp(-((9*x-7).^2 + (9*y-3).^2)/4) - ...
                 0.2 * exp(-(9*x-4).^2 - (9*y-7).^2);

%% Setup fine grid for evaluation and plotting
n_fine = 201;
x_fine = linspace(0, 1, n_fine);
y_fine = linspace(0, 1, n_fine);
[X_fine, Y_fine] = meshgrid(x_fine, y_fine);
Z_true = franke(X_fine, Y_fine);

%% Problem 1a: Chebyshev Approximation
fprintf('=== Problem 1a: Chebyshev Approximation ===\n');

% Parameters
N_cheb = 8;  % Polynomial degree

% Generate 2D Chebyshev nodes (extrema)
i_nodes = (0:N_cheb)';
j_nodes = (0:N_cheb)';
x_cheb_canon = cos(pi * i_nodes / N_cheb);  % [-1,1] domain
y_cheb_canon = cos(pi * j_nodes / N_cheb);  % [-1,1] domain

% Map to [0,1] domain
x_cheb = (x_cheb_canon + 1) / 2;
y_cheb = (y_cheb_canon + 1) / 2;

% Create tensor product grid
[X_cheb, Y_cheb] = meshgrid(x_cheb, y_cheb);
Z_cheb_nodes = franke(X_cheb, Y_cheb);

% Plot true function with Chebyshev nodes
fig1 = figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
surf(X_fine, Y_fine, Z_true, 'EdgeColor', 'none');
hold on;
scatter3(X_cheb(:), Y_cheb(:), Z_cheb_nodes(:), 50, 'ro', 'filled');
title('Franke Function with Chebyshev Nodes');
xlabel('x'); ylabel('y'); zlabel('z');
colorbar;
view(45, 30);

% Compute Chebyshev coefficients using the regression formula
% For 2D: Gamma_ij = (4/(N^2)) * sum_k sum_l f(x_k,y_l) * T_i(x_k) * T_j(y_l)
% with appropriate normalization factors

coeffs_cheb = zeros(N_cheb+1, N_cheb+1);
for i = 0:N_cheb
    for j = 0:N_cheb
        % Normalization factors
        c_i = (i == 0) + (i > 0);  % c_0 = 1, c_i = 2 for i > 0
        c_j = (j == 0) + (j > 0);  % c_0 = 1, c_j = 2 for j > 0
        
        % Compute coefficient
        sum_val = 0;
        for k = 0:N_cheb
            for l = 0:N_cheb
                T_i_k = cos(i * pi * k / N_cheb);
                T_j_l = cos(j * pi * l / N_cheb);
                sum_val = sum_val + Z_cheb_nodes(l+1, k+1) * T_i_k * T_j_l;
            end
        end
        coeffs_cheb(i+1, j+1) = (4 / (c_i * c_j * N_cheb^2)) * sum_val;
    end
end

% Evaluate Chebyshev approximation on fine grid
Z_cheb_approx = cheb_eval_2d(coeffs_cheb, X_fine, Y_fine);

% Plot Chebyshev approximation
subplot(1,3,2);
surf(X_fine, Y_fine, Z_cheb_approx, 'EdgeColor', 'none');
title('Chebyshev Approximation');
xlabel('x'); ylabel('y'); zlabel('z');
colorbar;
view(45, 30);

% Plot error
error_cheb = abs(Z_true - Z_cheb_approx);
subplot(1,3,3);
surf(X_fine, Y_fine, error_cheb, 'EdgeColor', 'none');
title('Chebyshev Approximation Error');
xlabel('x'); ylabel('y'); zlabel('|Error|');
colorbar;
view(45, 30);

saveas(fig1, 'results/figures/chebyshev_approximation.png');

% Report accuracy metrics
max_error_cheb = max(error_cheb(:));
mean_error_cheb = mean(error_cheb(:));
rmse_cheb = sqrt(mean(error_cheb(:).^2));

fprintf('Chebyshev Approximation Results:\n');
fprintf('  Max Error: %.6f\n', max_error_cheb);
fprintf('  Mean Error: %.6f\n', mean_error_cheb);
fprintf('  RMSE: %.6f\n\n', rmse_cheb);

%% Problem 1b: Ordinary Polynomial Approximation
fprintf('=== Problem 1b: Ordinary Polynomial Approximation ===\n');

% Generate uniform grid for polynomial approximation
N_poly = 8;
x_poly = linspace(0, 1, N_poly+1);
y_poly = linspace(0, 1, N_poly+1);
[X_poly, Y_poly] = meshgrid(x_poly, y_poly);
Z_poly_nodes = franke(X_poly, Y_poly);

% Plot true function with polynomial nodes
fig2 = figure('Position', [150, 150, 1200, 400]);

subplot(1,3,1);
surf(X_fine, Y_fine, Z_true, 'EdgeColor', 'none');
hold on;
scatter3(X_poly(:), Y_poly(:), Z_poly_nodes(:), 50, 'bo', 'filled');
title('Franke Function with Polynomial Nodes');
xlabel('x'); ylabel('y'); zlabel('z');
colorbar;
view(45, 30);

% Construct Vandermonde-like matrix for 2D polynomial
n_nodes = (N_poly+1)^2;
n_coeffs = (N_poly+1)^2;
X_matrix = zeros(n_nodes, n_coeffs);

% Fill the design matrix
node_idx = 1;
for i = 1:N_poly+1
    for j = 1:N_poly+1
        coeff_idx = 1;
        for p = 0:N_poly
            for q = 0:N_poly
                X_matrix(node_idx, coeff_idx) = X_poly(i,j)^p * Y_poly(i,j)^q;
                coeff_idx = coeff_idx + 1;
            end
        end
        node_idx = node_idx + 1;
    end
end

% Solve for coefficients: X * coeffs = z
coeffs_poly = X_matrix \ Z_poly_nodes(:);

% Evaluate polynomial on fine grid
Z_poly_approx = poly_eval_2d(coeffs_poly, N_poly, X_fine, Y_fine);

% Plot polynomial approximation
subplot(1,3,2);
surf(X_fine, Y_fine, Z_poly_approx, 'EdgeColor', 'none');
title('Polynomial Approximation');
xlabel('x'); ylabel('y'); zlabel('z');
colorbar;
view(45, 30);

% Plot error
error_poly = abs(Z_true - Z_poly_approx);
subplot(1,3,3);
surf(X_fine, Y_fine, error_poly, 'EdgeColor', 'none');
title('Polynomial Approximation Error');
xlabel('x'); ylabel('y'); zlabel('|Error|');
colorbar;
view(45, 30);

saveas(fig2, 'results/figures/polynomial_approximation.png');

% Report accuracy metrics
max_error_poly = max(error_poly(:));
mean_error_poly = mean(error_poly(:));
rmse_poly = sqrt(mean(error_poly(:).^2));

fprintf('Polynomial Approximation Results:\n');
fprintf('  Max Error: %.6f\n', max_error_poly);
fprintf('  Mean Error: %.6f\n', mean_error_poly);
fprintf('  RMSE: %.6f\n\n', rmse_poly);

%% Problem 1c: MATLAB Toolbox Methods
fprintf('=== Problem 1c: MATLAB Toolbox Methods ===\n');

% Method 1: griddedInterpolant with cubic splines
interp_cubic = griddedInterpolant(X_poly, Y_poly, Z_poly_nodes, 'cubic');
Z_cubic_approx = interp_cubic(X_fine, Y_fine);

% Method 2: scatteredInterpolant with random points
n_scatter = 100;
rng(42); % For reproducibility
x_scatter = rand(n_scatter, 1);
y_scatter = rand(n_scatter, 1);
z_scatter = franke(x_scatter, y_scatter);

interp_scattered = scatteredInterpolant(x_scatter, y_scatter, z_scatter, 'linear');
Z_scattered_approx = interp_scattered(X_fine, Y_fine);

% Plot toolbox methods
fig3 = figure('Position', [200, 200, 1200, 800]);

% Cubic spline
subplot(2,3,1);
surf(X_fine, Y_fine, Z_cubic_approx, 'EdgeColor', 'none');
title('Cubic Spline Interpolation');
xlabel('x'); ylabel('y'); zlabel('z');
colorbar;
view(45, 30);

% Cubic spline error
error_cubic = abs(Z_true - Z_cubic_approx);
subplot(2,3,2);
surf(X_fine, Y_fine, error_cubic, 'EdgeColor', 'none');
title('Cubic Spline Error');
xlabel('x'); ylabel('y'); zlabel('|Error|');
colorbar;
view(45, 30);

% Scattered data interpolation
subplot(2,3,4);
surf(X_fine, Y_fine, Z_scattered_approx, 'EdgeColor', 'none');
hold on;
scatter3(x_scatter, y_scatter, z_scatter, 30, 'ro', 'filled');
title('Scattered Data Interpolation');
xlabel('x'); ylabel('y'); zlabel('z');
colorbar;
view(45, 30);

% Scattered data error
error_scattered = abs(Z_true - Z_scattered_approx);
subplot(2,3,5);
surf(X_fine, Y_fine, error_scattered, 'EdgeColor', 'none');
title('Scattered Data Error');
xlabel('x'); ylabel('y'); zlabel('|Error|');
colorbar;
view(45, 30);

% Comparison plot
subplot(2,3,[3,6]);
methods = {'Chebyshev', 'Polynomial', 'Cubic Spline', 'Scattered'};
max_errors = [max_error_cheb, max_error_poly, max(error_cubic(:)), max(error_scattered(:))];
rmse_errors = [rmse_cheb, rmse_poly, sqrt(mean(error_cubic(:).^2)), sqrt(mean(error_scattered(:).^2))];

bar_x = 1:length(methods);
bar_width = 0.35;
bar(bar_x - bar_width/2, max_errors, bar_width, 'DisplayName', 'Max Error');
hold on;
bar(bar_x + bar_width/2, rmse_errors, bar_width, 'DisplayName', 'RMSE');
set(gca, 'XTickLabel', methods);
ylabel('Error');
title('Method Comparison');
legend('show');
yscale('log');

saveas(fig3, 'results/figures/toolbox_methods.png');

% Report toolbox method results
fprintf('Cubic Spline Results:\n');
fprintf('  Max Error: %.6f\n', max(error_cubic(:)));
fprintf('  RMSE: %.6f\n\n', sqrt(mean(error_cubic(:).^2)));

fprintf('Scattered Data Results:\n');
fprintf('  Max Error: %.6f\n', max(error_scattered(:)));
fprintf('  RMSE: %.6f\n\n', sqrt(mean(error_scattered(:).^2)));

%% Summary
fprintf('=== Summary ===\n');
fprintf('All approximation methods have been tested on the Franke function.\n');
fprintf('Results saved to results/figures/\n');
fprintf('- chebyshev_approximation.png\n');
fprintf('- polynomial_approximation.png\n');
fprintf('- toolbox_methods.png\n\n');

% Save workspace for potential further analysis
save('results/problem1_results.mat');
fprintf('Workspace saved to results/problem1_results.mat\n');