function z_approx = poly_eval_2d(coeffs, N, x_eval, y_eval)
% POLY_EVAL_2D Evaluates a 2D ordinary polynomial on a grid. (CORRECTED)
%
%   Inputs:
%       coeffs  - Vector of polynomial coefficients.
%       N       - The maximum degree of the polynomial in each dimension.
%       x_eval  - Grid of x-coordinates for evaluation.
%       y_eval  - Grid of y-coordinates for evaluation.
%
%   Output:
%       z_approx - The approximated values on the (x_eval, y_eval) grid.

z_approx = zeros(size(x_eval));
coeff_idx = 1;

% Loop through all combinations of powers for x and y up to degree N.
% This loop structure must match the coefficient generation script.
for p = 0:N % Power of x
    for q = 0:N % Power of y
        basis_term = (x_eval .^ p) .* (y_eval .^ q);
        z_approx = z_approx + coeffs(coeff_idx) * basis_term;
        coeff_idx = coeff_idx + 1;
    end
end

end