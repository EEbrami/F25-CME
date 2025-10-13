function I_2d = trapz_rule_2d(f_2d, a_x, b_x, a_y, b_y, M)
% TRAPZ_RULE_2D Computes the integral of f(x, y) over [a_x, b_x] x [a_y, b_y]
% using the 2D Compound Trapezoid Rule (Product Rule).
%
% Inputs:
%   f_2d: Function handle for the integrand, f(x, y). Must be vectorized: f(X, Y).
%   a_x, b_x: Limits for the x-dimension.
%   a_y, b_y: Limits for the y-dimension.
%   M: Number of subintervals in *each* dimension (M_x = M_y = M).
%
% Output:
%   I_2d: The approximate value of the definite integral.

    M_x = M;
    M_y = M;
    
    % --- 1. Generate 1D Nodes and Weights (for Trapezoid Rule) ---
    h_x = (b_x - a_x) / M_x;
    h_y = (b_y - a_y) / M_y;
    
    % Nodes in X and Y (M+1 points)
    x_nodes = linspace(a_x, b_x, M_x + 1);
    y_nodes = linspace(a_y, b_y, M_y + 1);
    
    % Trapezoid 1D Weights Pattern: [h/2, h, h, ..., h, h/2]
    % W_1D = [0.5, 1, 1, ..., 1, 0.5] * h
    W_x_pattern = ones(1, M_x + 1);
    W_y_pattern = ones(1, M_y + 1);
    
    % Set endpoint weights to 0.5
    W_x_pattern(1) = 0.5; W_x_pattern(end) = 0.5;
    W_y_pattern(1) = 0.5; W_y_pattern(end) = 0.5;
    
    % Apply step size (h) to get final 1D weight vectors
    W_x = h_x * W_x_pattern;
    W_y = h_y * W_y_pattern;
    
    % --- 2. Apply 2D Product Rule ---
    
    % Create 2D Grid of Nodes
    [X_grid, Y_grid] = meshgrid(x_nodes, y_nodes);
    
    % Evaluate the 2D function at all M^2 grid points
    F_grid = f_2d(X_grid, Y_grid);
    
    % The 2D weight matrix W_ij = W_y_i * W_x_j
    % Since W_x and W_y are row vectors, use outer product: W_y' * W_x
    W_grid = W_y' * W_x; 
    
    % The 2D integral is the sum of (Weights * Function_Evaluations)
    I_2d = sum(sum(W_grid .* F_grid));

end
