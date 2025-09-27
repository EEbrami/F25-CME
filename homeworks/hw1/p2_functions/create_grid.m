function [k_grid, z_grid, k_nodes, z_nodes] = create_grid(params, grid_params, method_name)
% CREATE_GRID generates the state space grid for k and z.

% 1. Create 1D capital grid
% --- BUG FIX: Read k_range from the main params struct ---
k_min = params.k_ss * (1 - params.k_range);
k_max = params.k_ss * (1 + params.k_range);

if strcmp(method_name, 'chebyshev')
    % Chebyshev nodes on [-1, 1]
    cheb_nodes = -cos((0:grid_params.k_points-1)' * pi / (grid_params.k_points-1));
    % Map to [k_min, k_max]
    k_nodes = k_min + (k_max - k_min) * (cheb_nodes + 1) / 2;
else % Spline uses a uniform grid
    k_nodes = linspace(k_min, k_max, grid_params.k_points)';
end

% 2. Create 1D productivity grid (log-uniform)
% --- BUG FIX: Read z_range from the main params struct ---
z_min = exp(-params.z_range);
z_max = exp(params.z_range);
z_nodes = linspace(z_min, z_max, grid_params.z_points);

% 3. Create 2D grid from 1D node vectors
[k_grid, z_grid] = ndgrid(k_nodes, z_nodes);
end