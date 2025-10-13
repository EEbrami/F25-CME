function [I_mc, Samples] = monte_carlo_2d(f_2d, a_x, b_x, a_y, b_y, T)
% MONTE_CARLO_2D Computes the 2D definite integral of f(x, y) over a rectangular 
% domain [a_x, b_x] x [a_y, b_y] using the Monte Carlo method.
%
% Inputs:
%   f_2d: Function handle for the integrand, f(x, y).
%   a_x, b_x: Limits for the x-dimension.
%   a_y, b_y: Limits for the y-dimension.
%   T: Total number of sampled points (T = 10,000 for P2c mandate).
%
% Outputs:
%   I_mc: The Monte Carlo approximation of the integral.
%   Samples: Matrix containing the T sampled (x, y) points.

    % 1. Calculate the Area of the integration domain
    Area = (b_x - a_x) * (b_y - a_y);
    
    % 2. Generate T uniformly distributed random points
    % rand(T, 1) generates T random numbers in [0, 1]. Scale to domain [a, b].
    X_samples = a_x + (b_x - a_x) * rand(T, 1);
    Y_samples = a_y + (b_y - a_y) * rand(T, 1);
    
    % Combine samples for debugging/storage (Optional output)
    Samples = [X_samples, Y_samples]; 
    
    % 3. Evaluate the integrand at all sampled points
    F_samples = f_2d(X_samples, Y_samples);
    
    % 4. Apply the Monte Carlo formula: I = Area * (1/T) * sum(f(x_i))
    % The core MC estimator for expected value is (1/T) * sum(f).
    % The integral approximation is I = Area * E[f(x,y)]
    I_mc = Area * (sum(F_samples) / T);

end
