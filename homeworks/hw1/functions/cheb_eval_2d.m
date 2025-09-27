function z_out = cheb_eval_2d(coeffs, x_eval, y_eval)
% CHEB_EVAL_2D Evaluates a 2D Chebyshev polynomial given coefficients and evaluation points
%
% SYNTAX:
%   z_out = cheb_eval_2d(coeffs, x_eval, y_eval)
%
% INPUTS:
%   coeffs  - (N+1)x(N+1) matrix of Chebyshev coefficients
%   x_eval  - Matrix of x evaluation points (domain [0,1])
%   y_eval  - Matrix of y evaluation points (domain [0,1])
%
% OUTPUTS:
%   z_out   - Matrix of evaluated function values
%
% DESCRIPTION:
%   Evaluates a 2D Chebyshev polynomial using the tensor product approach.
%   The domain of x_eval and y_eval is assumed to be [0,1], which is mapped
%   internally to the canonical domain [-1,1].
%
% EXAMPLE:
%   % Create coefficients matrix and evaluation points
%   coeffs = rand(6,6);
%   [x,y] = meshgrid(linspace(0,1,10), linspace(0,1,10));
%   z = cheb_eval_2d(coeffs, x, y);
%
% See also: poly_eval_2d

% Author: Problem Set 1, Econ-81360
% Date: Fall 2025

    % Map evaluation points from [0,1] to canonical domain [-1,1]
    x_canon = 2*x_eval - 1;
    y_canon = 2*y_eval - 1;
    
    N = size(coeffs, 1) - 1;

    % Evaluate basis matrices for x and y dimensions
    Tx = cos((0:N)' .* acos(x_canon(:)'));
    Ty = cos((0:N)' .* acos(y_canon(:)'));

    % The approximation is the matrix product: T_y' * Gamma * T_x
    % This efficiently computes the double summation
    z_out = Ty' * coeffs * Tx;
end