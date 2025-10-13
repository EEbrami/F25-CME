function [log10_max_error, log10_mean_error] = evaluate_accuracy(K_coef, D, params, k_path, z_path, eps_res_nodes, w_res_nodes)
% EVALUATE_ACCURACY: Computes Euler residuals (R_EE) along a simulated path 
% using the 10-node GHQ quadrature for the expectation term.

    T_simul_path = length(k_path);
    Residuals = zeros(T_simul_path, 1);
    
    % Unpack parameters
    delta = params.delta;
    alpha = params.alpha;
    A = params.A;
    rho = params.rho;

    % Policy function evaluator
    policy_eval = @(k, z) Polynomial_2d([k, z], D) * K_coef;
    
    % Policy derivative
    f_prime = @(k, z) A * alpha * z .* k.^(alpha-1);

    % Loop over the entire path 
    for t = 1:T_simul_path
        
        k_t = k_path(t);
        z_t = z_path(t);
        
        % 1. Calculate c_t and u'(c_t)
        k_next_guess = policy_eval(k_t, z_t); 
        c_t = A * z_t * (k_t^alpha) + (1 - delta) * k_t - k_next_guess;
        c_t = max(c_t, 1e-12); % Clamp consumption
        mu_t = 1 / c_t; % u'(c_t) for gamma=1
        
        % 2. Calculate the Expected RHS (E[u'(c')R']) using 10-node GHQ
        Expected_RHS = 0;
        k_next_guess = k_next_guess(1); % Ensure k_next is scalar for feval policy_eval

        % Loop over the 10 GHQ nodes (j=1...10)
        for j = 1:length(eps_res_nodes) 
            eps_j = eps_res_nodes(j);
            w_j = w_res_nodes(j);
            
            % Future state z' (for E[.])
            z_prime = z_t^rho * exp(eps_j); 
            
            % Future policy k'' = K(k', z')
            k_next_next = policy_eval(k_next_guess, z_prime); 
            
            % Future consumption c'
            c_prime = A * z_prime * (k_next_guess^alpha) + (1 - delta) * k_next_guess - k_next_next;
            c_prime = max(c_prime, 1e-12); % Clamp consumption
            
            % Marginal Utility and Return
            mu_prime = 1 / c_prime; 
            R_prime = (1 - delta + f_prime(k_next_guess, z_prime)); 
            
            Expected_RHS = Expected_RHS + w_j * mu_prime * R_prime;
        end
        
        % 3. Compute Residual: R(k, z) = (beta * E_RHS) / u'(c_t) - 1
        % R = (beta * E_RHS) * c_t - 1
        residual = (params.beta * Expected_RHS) * c_t - 1;
        
        Residuals(t) = residual;
    end
    
    % 4. Final log10 reporting (Max and Mean Absolute Residuals)
    log10_max_error = log10(max(abs(Residuals)));
    log10_mean_error = log10(mean(abs(Residuals)));
end