function I_2d = gc_quadrature_2d(f_2d, a_x, b_x, a_y, b_y, N_nodes)
% GC_QUADRATURE_2D Computes the definite integral of a function f(x, y) 
% over a rectangular domain [a_x, b_x] x [a_y, b_y] using the 
% 2D Gauss-Chebyshev Type 1 Product Rule.
%
% Inputs:
%   f_2d: Function handle for the integrand, f(x, y). Must be vectorized: f(X, Y).
%   a_x, b_x: Limits for the x-dimension.
%   a_y, b_y: Limits for the y-dimension.
%   N_nodes: Number of nodes to use in *each* dimension (N_x = N_y = N).
%
% Output:
%   I_2d: The approximate value of the definite integral.

    N = N_nodes;

    % --- 1. Compute 1D Nodes and Weights (Chebyshev Type 1 Adapted) ---

    i = (1:N)';
    
    % Chebyshev Nodes over [-1, 1] for x and y
    x_i = cos((2*i - 1) * pi / (2*N));
    
    % Compensation factor (sqrt(1 - x_i^2))
    sqrt_term = sqrt(1 - x_i.^2);
    
    % The constant factor of the adapted 1D rule: pi / (2*N)
    constant_factor = pi / (2 * N);
    
    % --- 2. X-Dimension Transformation and Weights ---
    % Transform nodes to [a_x, b_x]
    y_x = ((x_i + 1) * (b_x - a_x) / 2) + a_x;
    % Weights for X (W_x) = Constant_Factor * (b_x - a_x) * Compensation_Factor
    W_x = constant_factor * (b_x - a_x) * sqrt_term;

    % --- 3. Y-Dimension Transformation and Weights ---
    % Transform nodes to [a_y, b_y]
    y_y = ((x_i + 1) * (b_y - a_y) / 2) + a_y;
    % Weights for Y (W_y) = Constant_Factor * (b_y - a_y) * Compensation_Factor
    W_y = constant_factor * (b_y - a_y) * sqrt_term;

    % --- 4. Apply 2D Product Rule ---
    
    % Generate the 2D grid of nodes (X_ij, Y_ij) using meshgrid
    [X_grid, Y_grid] = meshgrid(y_x, y_y);
    
    % Generate the 2D grid of weights (W_ij = W_x_i * W_y_j)
    % Outer product of the two 1D weight vectors
    W_grid = W_y * W_x'; 
    
    % Evaluate the 2D function at ALL N*N grid points
    F_grid = f_2d(X_grid, Y_grid);
    
    % The 2D integral is the sum of (Weights * Function_Evaluations)
    % I_2d = sum_i(sum_j( W_ij * F_ij ))
    I_2d = sum(sum(W_grid .* F_grid));

end
