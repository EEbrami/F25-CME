function [k_path, z_path] = run_stochastic_simulation(K_coef, D, params, T_simul)
% RUN_STOCHASTIC_SIMULATION: Simulates the k and z paths for T_simul periods 
% starting at the steady state. Discards T_burn periods.
    
    T_burn = 200; % Mandated burn-in
    
    k_path = zeros(T_simul, 1);
    z_path = zeros(T_simul, 1);
    
    k_ss = params.k_ss;
    z_ss = 1.0; % z_ss=1.0 is implied by log(z_ss)=0 at steady state

    k_path(1) = k_ss;
    z_path(1) = z_ss;
    
    % Generate T_simul-1 random shocks from N(0, sigma^2)
    rng('default'); % Ensures reproducibility of the stochastic path
    shocks = params.sigma * randn(T_simul, 1);
    
    % Policy function evaluator (reusing Polynomial_2d)
    policy_eval = @(k, z) Polynomial_2d([k, z], D) * K_coef;

    for t = 1:(T_simul - 1)
        k_t = k_path(t);
        z_t = z_path(t);
        
        % z_{t+1} = z_t^rho * exp(epsilon_t+1)
        z_path(t+1) = z_t^params.rho * exp(shocks(t));
        
        % k_{t+1} = K(k_t, z_t)
        k_path(t+1) = policy_eval(k_t, z_t);
        
        % State Variable Clamps (Practical Stabilization for simulation)
        k_path(t+1) = max(min(k_path(t+1), 1.2 * k_ss), 0.8 * k_ss); 
        z_path(t+1) = max(min(z_path(t+1), 1.5 * z_ss), 0.5 * z_ss); 
    end
    
    % Discard burn-in periods (T_burn=200 is used here)
    k_path = k_path((T_burn+1):end);
    z_path = z_path((T_burn+1):end);
end