function [x, k] = gauss_jacobi(A, b, x0, tol, max_iter)
% Solves the linear system Ax=b using the Gauss-Jacobi iterative method.
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
    
    % Splitting A = N + P, where N is the diagonal of A
    N = diag(diag(A));
    P = A - N;
    
    % Pre-compute N_inv * b and -N_inv * P for efficiency
    % N_inv is just the reciprocal of the diagonal elements
    N_inv = diag(1 ./ diag(A));
    N_inv_b = N_inv * b;
    T = -N_inv * P; % Iteration matrix
    
    for k = 1:max_iter
        x_new = T * x + N_inv_b;
        
        % Check for convergence
        if norm(x_new - x, inf) < tol
            x = x_new;
            return; % Solution converged
        end
        
        x = x_new;
    end
    
    % If loop finishes, it did not converge
    warning('Gauss-Jacobi did not converge within %d iterations.', max_iter);
    k = max_iter;
end