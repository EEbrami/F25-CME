function z_out = poly_eval_2d(coeffs, N, x_eval, y_eval)
%POLY_EVAL_2D Evaluates a 2D ordinary polynomial.
%   INPUTS:
%     coeffs: Vector of polynomial coefficients.
%     N:      Degree of the polynomial.
%     x_eval: Matrix of x-coordinates for evaluation.
%     y_eval: Matrix of y-coordinates for evaluation.
%   OUTPUT:
%     z_out:  Matrix of approximated z-values.

z_out = zeros(size(x_eval));
m = N + 1;

for i = 1:numel(x_eval)
    px = x_eval(i);
    py = y_eval(i);
    basis_vec = zeros(1, m^2);
    idx = 1;
    for p = 0:N
        for q = 0:N
            basis_vec(idx) = px^p * py^q;
            idx = idx + 1;
        end
    end
    z_out(i) = basis_vec * coeffs;
end
end