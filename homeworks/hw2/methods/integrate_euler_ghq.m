function ExpectedRHS = integrate_euler_ghq(k_next, z_current, K_coef, D, params, epsilon_nodes, weight_nodes)
% INTEGRATE_EULER_GHQ: Replicates the original discrete integration using 
% Gauss-Hermite Quadrature (GHQ) nodes/weights.

    % --- Unpack parameters ---
    beta  = params.beta;
    delta = params.delta;
    alpha = params.alpha;
    A     = params.A;
    rho   = params.rho;
    
    % --- Compute Policy and Consumption in all GHQ Nodes (Vectorized) ---
    n_nodes = length(epsilon_nodes);
    
    % Scale nodes for N(0, sigma^2) 
    % Note: GH nodes are typically for N(0, 1) and must be scaled by sigma.
    % Assuming the provided GH_Quadrature.m outputs pre-scaled nodes (epsilon).
    epsilon_j = epsilon_nodes; 

    % Future state z'
    z_next = z_current.^rho .* exp(epsilon_j);

    % Policy evaluation k'' = K(k', z') at nodes
    X1 = Polynomial_2d([k_next*ones(n_nodes, 1), z_next], D); % Create basis matrix
    k_next_next = X1 * K_coef; % k''

    % Consumption c' at nodes
    output_next = A .* z_next .* (k_next.^alpha); 
    k_next_dupl = k_next * ones(n_nodes, 1);
    c_next = (1 - delta) .* k_next_dupl + output_next - k_next_next;
    c_next(c_next <= 1e-12) = 1e-12;

    % Marginal Utility u'(c') and Marginal Return R' (Vectorized g(epsilon))
    mu_next = 1 ./ c_next;
    R_prime = 1 - delta + A .* z_next .* (alpha .* (k_next_dupl.^(alpha-1)));
    RHS_term_g_eps = mu_next .* R_prime;
    
    % --- Discrete Summation (Integration) ---
    % E[g(epsilon)] = sum( w_j * g(epsilon_j) )
    Integral_Value = RHS_term_g_eps' * weight_nodes;
    
    ExpectedRHS = beta * Integral_Value;
end