function z_approx = chebyshev_eval_2d_p2(coeffs, degree, k_eval_grid, z_eval_grid, domain)
% CHEBYSHEV_EVAL_2D_P2 evaluates the 2D Chebyshev polynomial on a given domain.
% (CORRECTED to accept domain bounds as input)

k_flat = k_eval_grid(:);
z_flat = z_eval_grid(:);

% Map evaluation points from their domain to the canonical domain [-1, 1]
k_mapped = 2 * (k_flat - domain.k_min) ./ (domain.k_max - domain.k_min) - 1;
z_mapped = 2 * (z_flat - domain.z_min) ./ (domain.z_max - domain.z_min) - 1;

% Clamp mapped values to [-1, 1] to avoid complex numbers from acos
k_mapped(k_mapped > 1) = 1; k_mapped(k_mapped < -1) = -1;
z_mapped(z_mapped > 1) = 1; z_mapped(z_mapped < -1) = -1;

num_points = length(k_flat);
z_approx_flat = zeros(num_points, 1);
coeff_idx = 1;

% Reconstruct the polynomial value from coefficients
for i = 0:degree
    for j = 0:degree
        basis_term = cos(i * acos(k_mapped)) .* cos(j * acos(z_mapped));
        z_approx_flat = z_approx_flat + coeffs(coeff_idx) * basis_term;
        coeff_idx = coeff_idx + 1;
    end
end

z_approx = reshape(z_approx_flat, size(k_eval_grid));
end