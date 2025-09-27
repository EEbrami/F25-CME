function plot_policy_function(results, params, method_name)
% PLOT_POLICY_FUNCTION visualizes the solved policy and consumption functions.
%
%   Inputs:
%       results       - The struct returned by the solver.
%       params        - The struct with model parameters.
%       method_name   - A string (e.g., 'Chebyshev') for titles and filenames.

fprintf('Generating plots for the %s method...\n', method_name);

%% 1. Setup for Plotting
% Create a fine grid for smooth plots
n_fine = 101;
k_plot_vec = linspace(params.k_ss * (1 - params.k_range), params.k_ss * (1 + params.k_range), n_fine)';
z_plot_vec = linspace(exp(-params.z_range), exp(params.z_range), n_fine);
[k_plot_grid, z_plot_grid] = ndgrid(k_plot_vec, z_plot_vec);

% Get the solved policy function from the results struct
if isfield(results, 'policy_coeffs') % Chebyshev
    coeffs = results.policy_coeffs;
    degree = results.method_params.degree;
    domain.k_min = min(k_plot_vec); domain.k_max = max(k_plot_vec);
    domain.z_min = min(z_plot_vec); domain.z_max = max(z_plot_vec);
    k_prime_grid = chebyshev_eval_2d_p2(coeffs, degree, k_plot_grid, z_plot_grid, domain);
else % Spline
    policy_fun = results.policy_fun;
    k_prime_grid = policy_fun(k_plot_grid, z_plot_grid);
end

% Derive the consumption policy function
c_grid = params.f(k_plot_grid, z_plot_grid) + (1-params.delta)*k_plot_grid - k_prime_grid;

%% 2. Plot 1: 3D Capital Policy Function
figure('Name', ['Capital Policy Function (' method_name ')']);
surf(k_plot_grid, z_plot_grid, k_prime_grid);
xlabel('Current Capital (k_t)');
ylabel('Productivity (z_t)');
zlabel("Next Period's Capital (k_{t+1})");
title(['Capital Policy Function: k'' = K(k,z) (' method_name ')']);
saveas(gcf, ['results/figures_problem2/' method_name '_Capital_Policy.png']);

%% 3. Plot 2: 3D Consumption Policy Function
figure('Name', ['Consumption Policy Function (' method_name ')']);
surf(k_plot_grid, z_plot_grid, c_grid);
xlabel('Current Capital (k_t)');
ylabel('Productivity (z_t)');
zlabel('Consumption (c_t)');
title(['Consumption Policy Function: c = C(k,z) (' method_name ')']);
saveas(gcf, ['results/figures_problem2/' method_name '_Consumption_Policy.png']);

%% 4. Plot 3: 2D Policy Function Slices
figure('Name', ['Policy Function Slices (' method_name ')']);
hold on;
% Define z levels for slices (low, steady-state, high)
z_slice_indices = [1, round(n_fine/2), n_fine];
colors = {'b', 'k', 'r'};
legend_entries = cell(length(z_slice_indices), 1);

for i = 1:length(z_slice_indices)
    idx = z_slice_indices(i);
    plot(k_plot_vec, k_prime_grid(:, idx), 'Color', colors{i}, 'LineWidth', 2);
    legend_entries{i} = sprintf('z = %.2f', z_plot_vec(idx));
end

% Add 45-degree line for reference (where k' = k)
plot(k_plot_vec, k_plot_vec, 'k--', 'LineWidth', 1);
legend_entries{end+1} = '45-degree line (k''=k)';

xlabel('Current Capital (k_t)');
ylabel("Next Period's Capital (k_{t+1})");
title(['Policy Function Slices (' method_name ')']);
legend(legend_entries, 'Location', 'northwest');
grid on;
hold off;
saveas(gcf, ['results/figures_problem2/' method_name '_Policy_Slices.png']);
end