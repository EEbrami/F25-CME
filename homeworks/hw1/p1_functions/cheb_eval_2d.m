function z_approx = cheb_eval_2d(coeffs, x_eval_grid, y_eval_grid)
% CHEB_EVAL_2D Evaluates a 2D Chebyshev polynomial on a grid. (FINAL CORRECTED)
%
%   Inputs:
%       coeffs        - (N+1)x(N+1) matrix of final, scaled Chebyshev coeffs.
%       x_eval_grid   - Grid of x-coordinates for evaluation (from ndgrid).
%       y_eval_grid   - Grid of y-coordinates for evaluation (from ndgrid).
%
%   Output:
%       z_approx - The approximated values on the evaluation grid.

N = size(coeffs, 1) - 1;

% Extract unique 1D vectors for the grid points and map to [-1, 1]
z_x = 2 * x_eval_grid(:,1) - 1;
z_y = 2 * y_eval_grid(1,:)' - 1; % Needs to be a column vector

% Construct the 1D Chebyshev basis matrices for the evaluation points
% using the cosine definition, which is robust and matches the debug script.
Tx_eval = cos( acos(z_x) * (0:N) );
Ty_eval = cos( acos(z_y) * (0:N) );

% Evaluate the 2D polynomial using matrix multiplication.
% This is the direct and correct formula confirmed by the debug script.
z_approx = Tx_eval * coeffs * Ty_eval';

end