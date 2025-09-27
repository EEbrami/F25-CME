function basis_matrix = chebyshev_basis(k_nodes, z_nodes, degree)
% CHEBYSHEV_BASIS creates the 2D Chebyshev basis matrix for regression.

[k_grid, z_grid] = ndgrid(k_nodes, z_nodes);
k_flat = k_grid(:);
z_flat = z_grid(:);

% Map grid points from their domain to the canonical domain [-1, 1]
k_min = min(k_nodes); k_max = max(k_nodes);
z_min = min(z_nodes); z_max = max(z_nodes);
k_mapped = 2 * (k_flat - k_min) ./ (k_max - k_min) - 1;
z_mapped = 2 * (z_flat - z_min) ./ (z_max - z_min) - 1;

% Construct the basis matrix
num_points = length(k_flat);
num_coeffs = (degree + 1)^2;
basis_matrix = zeros(num_points, num_coeffs);

col_idx = 1;
for i = 0:degree
    for j = 0:degree
        % T_i(k) * T_j(z)
        basis_matrix(:, col_idx) = cos(i * acos(k_mapped)) .* cos(j * acos(z_mapped));
        col_idx = col_idx + 1;
    end
end
end