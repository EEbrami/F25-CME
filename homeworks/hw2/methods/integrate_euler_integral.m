function ExpectedRHS = integrate_euler_integral(k_next, z_current, K_coef, D, params)
% INTEGRATE_EULER_INTEGRAL: Calculates E_t[u'(c')R'] using MATLAB's modern 'integral'.

    beta = params.beta;
    
    % Define the integrand function handle, binding all required parameters
    integrand = @(epsilon) integrand_function(epsilon, k_next, z_current, K_coef, D, params);
    
    % Use 'integral' over the infinite domain (-Inf to Inf) for robustness.
    Integral_Value = integral(integrand, -Inf, Inf, 'AbsTol', 1e-10, 'RelTol', 1e-8);
    
    ExpectedRHS = beta * Integral_Value;
end