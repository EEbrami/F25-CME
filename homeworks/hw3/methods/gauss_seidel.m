function [x, k] = gauss_seidel(A, b, x0, tol, max_iter)
% Solves the linear system Ax=b using the Gauss-Seidel iterative method.
%
% Inputs:
%   A         - n x n coefficient matrix
%   b         - n x 1 vector
%   x0        - n x 1 initial guess
%   tol       - convergence tolerance
%   max_iter  - maximum number of iterations
%
% Outputs:
%   x         - n x 1 solution vector
%   k         - number of iterations taken

    n = length(b);
    x = x0;
    
    % Splitting A = N + P, where N is the lower triangle of A
    N = tril(A);
    P = A - N;
    
    for k = 1:max_iter
        % x_new = N \ (-P * x + b);
        % We solve N*x_new = -P*x + b
        % This can be solved more efficiently element-by-element
        
        x_new = x; % Start with old x to update in-place
        
        for i = 1:n
            % Sum of terms from P (upper triangle)
            sum_p = P(i, i+1:n) * x(i+1:n);
            
            % Sum of terms from N (lower triangle, *already updated*)
            sum_n = N(i, 1:i-1) * x_new(1:i-1);
            
            % Solve for x_new(i)
            x_new(i) = (b(i) - sum_p - sum_n) / A(i, i);
        end
        
        % Check for convergence
        if norm(x_new - x, inf) < tol
            x = x_new;
            return; % Solution converged
        end
        
        x = x_new;
    end
    
    % If loop finishes, it did not converge
    warning('Gauss-Seidel did not converge within %d iterations.', max_iter);
    k = max_iter;
end