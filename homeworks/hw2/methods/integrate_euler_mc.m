function ExpectedRHS = integrate_euler_mc(k_next, z_current, K_coef, D, params)
% INTEGRATE_EULER_MC: Calculates E_t[u'(c')R'] using Monte Carlo integration (T=10 nodes).
% Note: The total number of Monte Carlo nodes T is hardcoded here (T=10).
    
    T_nodes = 10;
    beta    = params.beta;
    sigma   = params.sigma;
    
    % 1. Draw 10 random shocks from the shock distribution epsilon ~ N(0, sigma^2)
    % randn() draws from N(0, 1), so scale by sigma.
    epsilon_j = params.sigma * randn(T_nodes, 1);
    weight_j = (1/T_nodes) * ones(T_nodes, 1); % Uniform weights 
    
    % --- Core Integrand Evaluation (Vectorized over 10 nodes) ---
    
    % Future state variables (z') for all T nodes
    z_next = z_current.^params.rho .* exp(epsilon_j);
    
    % Policy evaluation k'' = K(k', z') 
    X1 = Polynomial_2d([k_next*ones(T_nodes, 1), z_next], D); % Basis matrix (T x N_basis)
    k_next_next = X1 * K_coef; % k'' (T x 1 vector)
    
    % Future consumption c'
    output_next = params.A .* z_next .* (k_next.^params.alpha); 
    c_next = (1 - params.delta) .* k_next + output_next - k_next_next;
    c_next(c_next <= 1e-12) = 1e-12; % Numerical stability (clamping)
    
    % Marginal Utility u'(c') and Marginal Return R'
    mu_next = 1 ./ c_next; 
    R_prime = 1 - params.delta + params.A .* z_next .* (params.alpha .* (k_next.^(params.alpha-1))); 
    RHS_term_g_eps = mu_next .* R_prime; % g(epsilon_j)
    
    % --- Discrete Summation (Monte Carlo Approximation) ---
    % E[g(epsilon)] = sum( weight_j * g(epsilon_j) )
    Integral_Value = sum(weight_j .* RHS_term_g_eps);
    
    ExpectedRHS = beta * Integral_Value;
end