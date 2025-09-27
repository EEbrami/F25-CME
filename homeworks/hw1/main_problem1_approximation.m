%% ========================================================================
% ECON-81360 - Problem 1: Function Approximation
% Author: Ebrahim Ebrami
% Date: Fall 2025
% This script solves Problem 1 of the problem set 1: benchmarking 2D function
% approximation methods using the Franke function.
% =========================================================================

clear; clc; close all;

% Add functions directory to path
addpath('p1_functions');

% Create results directory
if ~exist('results/figures_problem1', 'dir')
    mkdir('results/figures_problem1');
end

fprintf('=== Problem 1: Function Approximation Benchmarking ===\n\n');

%% 1. Setup: Define Franke Function and Benchmark Surface

% The Franke function is a standard benchmark for testing approximation methods
franke = @(x,y) (0.75 * exp(-((9*x-2).^2 + (9*y-2).^2)/4) + ...
                 0.75 * exp(-((9*x+1).^2/49 + (9*y+1).^2/10)) + ...
                 0.5 * exp(-((9*x-7).^2 + (9*y-3).^2)/4) - ...
                 0.2 * exp(-((9*x-4).^2 + (9*y-7).^2)));

% Create a fine grid for plotting the "true" function and calculating error
n_fine = 201;
grid_fine_1d = linspace(0, 1, n_fine);
[x_fine, y_fine] = ndgrid(grid_fine_1d, grid_fine_1d);
z_true = franke(x_fine, y_fine);

% Visualize the true function
figure('Name', 'Benchmark: Franke Function', 'Position', [100, 100, 800, 600]);
surf(x_fine, y_fine, z_true);
title('True Franke Function');
xlabel('x'); ylabel('y'); zlabel('z');
colormap(parula);
view(145, 30);
set(gca, 'FontSize', 12);
grid on;
saveas(gcf, 'results/figures_problem1/01_franke_function.png');
fprintf('Step 1 Complete: True Franke function plotted and saved.\n');

%% ========================================================================
% Problem 1a: 2D Chebyshev Polynomial Approximation
% =========================================================================
fprintf('\n--- Running Problem 1a: Chebyshev Approximation ---\n');

% Parameters
N = 15; % Order of Chebyshev polynomial (number of nodes = N+1)
m = N + 1;

% 2.1. Grid Generation (Chebyshev Nodes)
% Using Chebyshev extrema (Gauss-Lobatto nodes)
k = (0:N)';
nodes_1d_cheb = 0.5 * (1 + cos(k * pi / N)); % Mapped to [0,1]
[x_nodes_cheb, y_nodes_cheb] = ndgrid(nodes_1d_cheb, nodes_1d_cheb);

% 2.2. Evaluate function on the grid
z_nodes_cheb = franke(x_nodes_cheb, y_nodes_cheb);

% Visualization: Show interpolation nodes on the true surface
figure('Name', 'Chebyshev Nodes', 'Position', [150, 150, 800, 600]);
surf(x_fine, y_fine, z_true, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
scatter3(x_nodes_cheb(:), y_nodes_cheb(:), z_nodes_cheb(:), 50, 'r', 'filled');
title(['Franke Function with ' num2str(m*m) ' Chebyshev Nodes (N=' num2str(N) ')']);
xlabel('x'); ylabel('y'); zlabel('z');
legend('True Surface', 'Interpolation Nodes', 'Location', 'northwest');
view(145, 30);
set(gca, 'FontSize', 12);
grid on;
saveas(gcf, 'results/figures_problem1/02_chebyshev_nodes.png');
fprintf('Generated %d Chebyshev nodes and plotted them.\n', m*m);

% 2.3. Coefficient Calculation
fprintf('Calculated Chebyshev coefficients.\n');

% Basis matrix T_ij = T_i(node_j)
% Note: The z_domain from your code is the same as nodes_1d in the diagnostic
z_domain = cos(k * pi / N); % Defines the nodes on the canonical [-1, 1] domain
T_basis = cos( (0:N)' * acos(z_domain') );

% Scale the function values at the nodes for the transform
F_scaled = z_nodes_cheb;
F_scaled([1, m], :) = F_scaled([1, m], :) * 0.5;
F_scaled(:, [1, m]) = F_scaled(:, [1, m]) * 0.5;

% Calculate coefficients via the correct formula
Gamma_cheb = (4 / (N * N)) * (T_basis * F_scaled * T_basis');
% CRITICAL FIX: Scale the boundary coefficients for the reconstruction formula
Gamma_cheb([1, m], :) = Gamma_cheb([1, m], :) * 0.5;
Gamma_cheb(:, [1, m]) = Gamma_cheb(:, [1, m]) * 0.5;

% 2.4. Evaluation/Interpolation
% Create an evaluation function to compute the approximation on the fine grid
% OLD
% z_approx_cheb = cheb_eval_2d(Gamma_cheb, x_coarse, y_coarse);

% NEW
z_approx_cheb = cheb_eval_2d(Gamma_cheb, x_fine, y_fine);

% 2.5. Accuracy and Visualization
error_cheb = abs(z_true - z_approx_cheb);
max_err_cheb = max(error_cheb(:));
mean_err_cheb = mean(error_cheb(:));
fprintf('Chebyshev (N=%d): Max absolute error = %.4e, Mean absolute error = %.4e\n', ...
        N, max_err_cheb, mean_err_cheb);

% Plot the approximated surface and error
figure('Name', 'Chebyshev Approximation Result', 'Position', [200, 200, 1200, 500]);
subplot(1, 2, 1);
surf(x_fine, y_fine, z_approx_cheb);
title(['Chebyshev Approx. (N=' num2str(N) ')']);
xlabel('x'); ylabel('y'); zlabel('z');
view(145, 30); set(gca, 'FontSize', 12); grid on;

subplot(1, 2, 2);
surf(x_fine, y_fine, error_cheb);
title(['Absolute Error (Max = ' sprintf('%.2e', max_err_cheb) ')']);
xlabel('x'); ylabel('y'); zlabel('Error');
view(145, 30); set(gca, 'FontSize', 12); grid on;
saveas(gcf, 'results/figures_problem1/03_chebyshev_result.png');
fprintf('Plotted Chebyshev approximation and error.\n');

%% ========================================================================
% Problem 1b: Ordinary Polynomial Approximation
% =========================================================================
fprintf('\n--- Running Problem 1b: Ordinary Polynomial Approximation ---\n');

% Parameters
N_poly = 8; % Degree of polynomial. Note: Higher degrees become unstable fast!
m_poly = N_poly + 1;

% 3.1. Grid Generation (Uniform Nodes)
nodes_1d_poly = linspace(0, 1, m_poly);
[x_nodes_poly, y_nodes_poly] = ndgrid(nodes_1d_poly, nodes_1d_poly);

% 3.2. Evaluate function on the grid
z_nodes_poly = franke(x_nodes_poly, y_nodes_poly);

% Visualization: Show interpolation nodes on the true surface
figure('Name', 'Uniform Nodes', 'Position', [250, 250, 800, 600]);
surf(x_fine, y_fine, z_true, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
scatter3(x_nodes_poly(:), y_nodes_poly(:), z_nodes_poly(:), 50, 'b', 'filled');
title(['Franke Function with ' num2str(m_poly^2) ' Uniform Nodes (N=' num2str(N_poly) ')']);
xlabel('x'); ylabel('y'); zlabel('z');
legend('True Surface', 'Interpolation Nodes', 'Location', 'northwest');
view(145, 30);
set(gca, 'FontSize', 12); grid on;
saveas(gcf, 'results/figures_problem1/04_uniform_nodes.png');
fprintf('Generated %d uniform nodes and plotted them.\n', m_poly^2);

% 3.3. Coefficient Calculation
% Construct the Vandermonde-like matrix X
X_mat = zeros(m_poly^2, m_poly^2);
y_vec = reshape(z_nodes_poly', [], 1); % Corrected for row-major order
idx = 1;
for i = 1:m_poly
    for j = 1:m_poly
        % Current point (x_i, y_j)
        px = x_nodes_poly(i,j);
        py = y_nodes_poly(i,j);
        
        % Create row for this point with basis functions x^p * y^q
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

% Solve the linear system X*c = y for coefficients c
coeffs_poly = X_mat \ y_vec;

% Check condition number to see potential for instability
cond_X = cond(X_mat);
fprintf('Condition number of Vandermonde matrix is %.2e. High values indicate instability.\n', cond_X);
fprintf('Calculated Ordinary Polynomial coefficients.\n');

% 3.4. Evaluation/Interpolation
z_approx_poly = poly_eval_2d(coeffs_poly, N_poly, x_fine, y_fine);

% 3.5. Accuracy and Visualization
error_poly = abs(z_true - z_approx_poly);
max_err_poly = max(error_poly(:));
mean_err_poly = mean(error_poly(:));
fprintf('Ordinary Poly (N=%d): Max absolute error = %.4e, Mean absolute error = %.4e\n', ...
        N_poly, max_err_poly, mean_err_poly);

% Plot the approximated surface and error
figure('Name', 'Ordinary Polynomial Result', 'Position', [300, 300, 1200, 500]);
subplot(1, 2, 1);
surf(x_fine, y_fine, z_approx_poly);
title(['Ordinary Poly Approx. (N=' num2str(N_poly) ')']);
xlabel('x'); ylabel('y'); zlabel('z');
view(145, 30); set(gca, 'FontSize', 12); grid on;

subplot(1, 2, 2);
surf(x_fine, y_fine, error_poly);
title(['Absolute Error (Max = ' sprintf('%.2e', max_err_poly) ')']);
xlabel('x'); ylabel('y'); zlabel('Error');
view(145, 30); set(gca, 'FontSize', 12); grid on;
saveas(gcf, 'results/figures_problem1/05_polynomial_result.png');
fprintf('Plotted Ordinary Polynomial approximation and error.\n');

%% ========================================================================
% Problem 1c: Toolbox Implementations
% =========================================================================
fprintf('\n--- Running Problem 1c: Toolbox Implementations ---\n');

% MATLAB griddedInterpolant for Splines
F_spline = griddedInterpolant(x_nodes_poly, y_nodes_poly, z_nodes_poly, 'cubic');
z_approx_spline = F_spline(x_fine, y_fine);
error_spline = abs(z_true - z_approx_spline);
max_err_spline = max(error_spline(:));
fprintf('Cubic Spline (grid %dx%d): Max absolute error = %.4e\n', ...
        m_poly, m_poly, max_err_spline);

figure('Name', 'Cubic Spline Approximation', 'Position', [350, 350, 1200, 500]);
subplot(1,2,1); 
surf(x_fine, y_fine, z_approx_spline); 
title('Cubic Spline Approx.'); 
xlabel('x'); ylabel('y'); zlabel('z');
view(145, 30); set(gca, 'FontSize', 12); grid on;
subplot(1,2,2); 
surf(x_fine, y_fine, error_spline); 
title(['Spline Error (Max = ' sprintf('%.2e', max_err_spline) ')']); 
xlabel('x'); ylabel('y'); zlabel('Error');
view(145, 30); set(gca, 'FontSize', 12); grid on;
saveas(gcf, 'results/figures_problem1/06_spline_result.png');
fprintf('Cubic spline approximation complete.\n');

% MATLAB scatteredInterpolant
num_scattered = m_poly^2;
x_scatter = rand(num_scattered, 1);
y_scatter = rand(num_scattered, 1);
z_scatter = franke(x_scatter, y_scatter);

F_scatter = scatteredInterpolant(x_scatter, y_scatter, z_scatter, 'natural');
z_approx_scatter = F_scatter(x_fine, y_fine);
error_scatter = abs(z_true - z_approx_scatter);
max_err_scatter = max(error_scatter(:));
fprintf('Scattered Interpolant (%d points): Max absolute error = %.4e\n', ...
        num_scattered, max_err_scatter);

figure('Name', 'Scattered Data Interpolation', 'Position', [400, 400, 1200, 500]);
subplot(1,2,1); 
surf(x_fine, y_fine, z_approx_scatter); 
title('Scattered Data Approx.'); 
xlabel('x'); ylabel('y'); zlabel('z');
view(145, 30); set(gca, 'FontSize', 12); grid on;
subplot(1,2,2); 
surf(x_fine, y_fine, error_scatter); 
title(['Scattered Error (Max = ' sprintf('%.2e', max_err_scatter) ')']); 
xlabel('x'); ylabel('y'); zlabel('Error');
view(145, 30); set(gca, 'FontSize', 12); grid on;
saveas(gcf, 'results/figures_problem1/07_scattered_result.png');
fprintf('Scattered data interpolation complete.\n');

%% ========================================================================
% Summary Results
% =========================================================================
fprintf('\n=== SUMMARY RESULTS ===\n');
fprintf('Method Comparison (Maximum Absolute Error):\n');
fprintf('%-25s %12s\n', 'Method', 'Max Error');
fprintf('%-25s %12s\n', '------', '---------');
fprintf('%-25s %12.4e\n', 'Chebyshev (N=15)', max_err_cheb);
fprintf('%-25s %12.4e\n', 'Ordinary Poly (N=8)', max_err_poly);
fprintf('%-25s %12.4e\n', 'Cubic Spline', max_err_spline);
fprintf('%-25s %12.4e\n', 'Scattered Interpolant', max_err_scatter);

% Create summary table figure
figure('Name', 'Method Comparison', 'Position', [450, 450, 800, 400]);
methods = {'Chebyshev', 'Ord. Poly', 'Cubic Spline', 'Scattered'};
errors = [max_err_cheb, max_err_poly, max_err_spline, max_err_scatter];
bar(log10(errors));
set(gca, 'XTickLabel', methods);
title('Approximation Method Comparison (Log10 Scale)');
ylabel('Log10(Maximum Absolute Error)');
grid on; set(gca, 'FontSize', 12);
for i = 1:length(errors)
    text(i, log10(errors(i))+0.1, sprintf('%.2e', errors(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end
saveas(gcf, 'results/figures_problem1/08_method_comparison.png');

fprintf('\nAll figures saved to results/figures_problem1/\n');
fprintf('Problem 1 complete!\n');