function z_out = poly_eval_2d(coeffs, N, x_eval, y_eval)
% POLY_EVAL_2D Evaluates a 2D ordinary polynomial given coefficients and evaluation points
%
% SYNTAX:
%   z_out = poly_eval_2d(coeffs, N, x_eval, y_eval)
%
% INPUTS:
%   coeffs  - Vector of polynomial coefficients (length (N+1)^2)
%   N       - Degree of the polynomial
%   x_eval  - Matrix of x evaluation points
%   y_eval  - Matrix of y evaluation points
%
% OUTPUTS:
%   z_out   - Matrix of evaluated function values
%
% DESCRIPTION:
%   Evaluates a 2D ordinary (monomial) polynomial of the form:
%   sum_{i=0}^N sum_{j=0}^N c_{ij} * x^i * y^j
%   
%   The coefficients are assumed to be ordered as:
%   [c_00, c_01, ..., c_0N, c_10, c_11, ..., c_NN]
%
% EXAMPLE:
%   % Create coefficients for a quadratic polynomial
%   coeffs = [1; 0.5; -0.3; 0.2; 0.1; 0; -0.1; 0; 0.05];
%   [x,y] = meshgrid(linspace(0,1,10), linspace(0,1,10));
%   z = poly_eval_2d(coeffs, 2, x, y);
%
% See also: cheb_eval_2d

% Author: Problem Set 1, Econ-81360
% Date: Fall 2025

    z_out = zeros(size(x_eval));
    m = N + 1;
    
    % Loop through each point in the evaluation grid
    for i = 1:numel(x_eval)
        px = x_eval(i);
        py = y_eval(i);
        
        % Construct the basis vector for this point
        basis_vec = zeros(1, m^2);
        idx = 1;
        for p = 0:N
            for q = 0:N
                basis_vec(idx) = px^p * py^q;
                idx = idx + 1;
            end
        end
        
        % Compute the value as the dot product of the basis and coefficients
        z_out(i) = basis_vec * coeffs;
    end
end