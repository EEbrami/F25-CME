function ExpectedRHS = integrate_euler_quad(k_next, z_current, K_coef, D, params)
% INTEGRATE_EULER_QUAD: Calculates E_t[u'(c')R'] using MATLAB's 'quad' (Adaptive Simpson's Rule).

    beta = params.beta;
    
    % Define the integrand function handle, binding all required parameters
    integrand = @(epsilon) integrand_function(epsilon, k_next, z_current, K_coef, D, params);
    
    % Use 'quad' over a practically infinite domain (e.g., +/- 10 standard deviations)
    % Note: 'quad' is older and less suited for infinite bounds than 'integral'.
    % Setting high tolerance for accuracy.
    Integral_Value = quad(integrand, -10*params.sigma, 10*params.sigma, 1e-10);
    
    ExpectedRHS = beta * Integral_Value;
end