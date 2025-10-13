function h_epsilon = integrand_function(epsilon, k_next, z_current, K_coef, D, params)
% INTEGRAND_FUNCTION: Computes the full integrand h(epsilon) = g(epsilon) * phi(epsilon) 
% for the Euler Equation, where phi(epsilon) is the N(0, sigma^2) PDF.
    
    % --- Unpack parameters ---
    delta = params.delta;
    alpha = params.alpha;
    A     = params.A;
    rho   = params.rho;
    sigma = params.sigma;

    % 1. Future state variables (z')
    z_next = z_current.^rho .* exp(epsilon);
    
    % 2. Future policy function k'' = K(k', z')
    % Assumes Polynomial_2d(X, D) is accessible for basis evaluation
    X1 = Polynomial_2d([k_next, z_next], D); 
    k_next_next = X1 * K_coef; % k''
    
    % 3. Future consumption c'
    output_next = A .* z_next .* (k_next.^alpha); 
    c_next = (1 - delta) .* k_next + output_next - k_next_next;
    c_next(c_next <= 1e-12) = 1e-12; % numerical stability
    
    % 4. Marginal Utility u'(c') and Marginal Return R' (The g(epsilon) part)
    mu_next = 1 ./ c_next; % u'(c') for gamma=1
    R_prime = 1 - delta + A .* z_next .* (alpha .* (k_next.^(alpha-1))); 
    RHS_term_g_eps = mu_next .* R_prime;
    
    % 5. Normal PDF phi(epsilon)
    pdf_epsilon = normpdf(epsilon, 0, sigma);
    
    % The full integrand: h(epsilon) = g(epsilon) * phi(epsilon)
    h_epsilon = RHS_term_g_eps .* pdf_epsilon;
end