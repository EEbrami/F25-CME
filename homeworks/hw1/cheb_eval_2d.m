function z_out = cheb_eval_2d(coeffs, x_eval, y_eval)
%CHEB_EVAL_2D Evaluates a 2D Chebyshev polynomial.
%   INPUTS:
%     coeffs: (N+1)x(N+1) matrix of Chebyshev coefficients.
%     x_eval: Matrix of x-coordinates for evaluation (domain [0,1]).
%     y_eval: Matrix of y-coordinates for evaluation (domain [0,1]).
%   OUTPUT:
%     z_out:  Matrix of approximated z-values.

% Map evaluation points from [0,1] to canonical domain [-1,1]
x_canon = 2*x_eval - 1;
y_canon = 2*y_eval - 1;

N = size(coeffs, 1) - 1;

% Evaluate basis matrices for x and y dimensions
Tx = cos((0:N)' .* acos(x_canon(:)'));
Ty = cos((0:N)' .* acos(y_canon(:)'));

% The approximation is the matrix product: T_y' * Gamma * T_x
z_out = Ty' * coeffs * Tx;
end