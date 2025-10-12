function I = trapezoid_rule(f, a, b, M)
% TRAPEZOID_RULE Computes the definite integral of a function f over [a, b] 
% using the Compound Trapezoid Rule with M subintervals.
%
% Inputs:
%   f: Function handle for the integrand, f(x).
%   a: Lower limit of integration.
%   b: Upper limit of integration.
%   M: Number of subintervals (subperiods).
%
% Output:
%   I: The approximate value of the definite integral.

    % 1. Calculate step size (width of each subinterval)
    h = (b - a) / M;
    
    % 2. Define the nodes (x_0 to x_M)
    x_nodes = linspace(a, b, M + 1);
    
    % 3. Evaluate the function at all nodes
    f_nodes = f(x_nodes);
    
    % 4. Apply the Compound Trapezoid Rule formula:
    % (f(x0) + f(xM)) + 2 * sum(f(x1) to f(xM-1))
    
    % Sum of the internal nodes (from x_2 to x_M, which is f_nodes(2) to f_nodes(M))
    if M > 1
        sum_internal = sum(f_nodes(2:M));
    else
        sum_internal = 0; % No internal nodes if M=1
    end
    
    % Final integration value
    I = (h / 2) * (f_nodes(1) + 2 * sum_internal + f_nodes(M + 1));
end