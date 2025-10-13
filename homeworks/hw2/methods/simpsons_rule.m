function I = simpsons_rule(f, a, b, M)
% SIMPSONS_RULE Computes the definite integral of a function f over [a, b] 
% using the Compound Simpson's Rule with M subintervals. M MUST be even.
%
% Inputs:
%   f: Function handle for the integrand, f(x).
%   a: Lower limit of integration.
%   b: Upper limit of integration.
%   M: Number of subintervals (subperiods). M MUST BE EVEN.
%
% Output:
%   I: The approximate value of the definite integral.

    % 1. Check constraint: M must be even
    if mod(M, 2) ~= 0
        error('Simpsons_Rule: M must be an even number of subintervals.');
    end

    % 2. Calculate step size
    h = (b - a) / M;
    
    % 3. Define the nodes (x_0 to x_M)
    x_nodes = linspace(a, b, M + 1);
    
    % 4. Evaluate the function at all nodes
    f_nodes = f(x_nodes);
    
    % 5. Apply the Compound Simpson's Rule formula weights:
    % 1 * f(x0) + 4 * sum(odd nodes) + 2 * sum(even nodes) + 1 * f(xM)
    
    % Indices for the internal nodes:
    % Odd indices (x1, x3, x5, ...): f_nodes(2), f_nodes(4), ...
    % Even indices (x2, x4, x6, ...): f_nodes(3), f_nodes(5), ...
    
    % Summation of odd-indexed nodes (weighted by 4)
    sum_odd = sum(f_nodes(2:2:M));
    
    % Summation of even-indexed nodes (weighted by 2)
    % Only calculate if M >= 4 (i.e., there is at least one even internal node x2)
    if M >= 4
        sum_even = sum(f_nodes(3:2:M-1));
    else
        sum_even = 0;
    end
    
    % Final integration value
    I = (h / 3) * (f_nodes(1) + 4 * sum_odd + 2 * sum_even + f_nodes(M + 1));
end