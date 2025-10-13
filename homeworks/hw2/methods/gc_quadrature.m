function I = gc_quadrature(f, a, b, n)
% gc_quadrature Computes the definite integral of a function f over [a, b] 
% using the Gauss-Chebyshev quadrature formula, adapted for a standard 
% interval with no weight function in the integrand.
%
% Inputs:
%   f: Function handle for the integrand, f(x).
%   a: Lower limit of integration.
%   b: Upper limit of integration.
%   n: Order of the polynomial (number of nodes).
%
% Output:
%   I: The approximate value of the definite integral.

% --- 1. Compute Chebyshev Nodes (x_i) over [-1, 1] ---
% Formula: x_i = cos((2*i - 1) * pi / (2*n)) for i = 1,...,n[cite: 6497].
i = (1:n)';
x_i = cos((2*i - 1) * pi / (2*n));

% --- 2. Transform Nodes to the Integration Domain [a, b] (y_i) ---
% Formula: y_i = a + (b - a) * (x_i + 1) / 2[cite: 6507].
% Note: The formula for the transformed nodes is embedded in the overall quadrature sum.
% The expression inside the function call is y_i.
y_i = ((x_i + 1) * (b - a) / 2) + a;

% --- 3. Compute Weights and Apply Quadrature Formula ---
% The formula for arbitrary domain [a, b] with $f(y)$ (not weighted $f(x)w(x)$):
% I = (pi * (b-a) / (2*n)) * sum( f(y_i) * sqrt(1 - x_i^2) )[cite: 6519].
% Note: The term sqrt(1 - x_i^2) serves as the correction factor derived from
% applying the change of variables and integrating unweighted f(y)[cite: 6517, 6519].

% Compute the intermediate factor sqrt(1 - x_i^2)
sqrt_term = sqrt(1 - x_i.^2);

% Evaluate the function at the transformed nodes
f_y_i = f(y_i);

% Calculate the summation term
sum_term = sum(f_y_i .* sqrt_term);

% Final integral approximation
I = (pi * (b - a) / (2 * n)) * sum_term;

end