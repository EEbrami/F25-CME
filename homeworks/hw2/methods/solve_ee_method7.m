function [K_coef, iter_count, k_mean_final] = solve_ee_method7(X0, k0, z0, K_coef_initial, D, params, integration_method)
% Generic fixed-point solver for Method 7.
% integration_method is a function handle (@integrate_euler_quad or @integrate_euler_integral)

    K_coef = K_coef_initial;
    M_grid = size(k0, 1);
    max_iter = params.max_iter;
    tol = params.tol;
    
    for iter_count = 1:max_iter
        K_coef_prev = K_coef;
        k1_new = zeros(M_grid, 1); % Target k' (LHS of regression)

        % --- Step 1: Compute Target k' for each grid point ---
        for m = 1:M_grid 
            k_m = k0(m, 1);
            z_m = z0(m, 1);
            X0_m = X0(m, :); 
            
            % Step 1.a: Current guess for k' = K_hat(k, z; V(i))
            k_next_guess = X0_m * K_coef; 
            
            % Step 1.c: Compute Expected Value E_t[u'(c')R']
            Expected_RHS = feval(integration_method, k_next_guess, z_m, K_coef, D, params);
            
            % Compute current consumption c_m (gamma=1: c_m = 1 / (beta * E_RHS))
            c_m = 1 / (params.beta * Expected_RHS);
                         
            % Step 1.d: Compute new capital target k'_m
            k1_new(m) = (1 - params.delta) * k_m + params.A * z_m * (k_m^params.alpha) - c_m;    
        end
        
        % --- Step 2: Regression and Update v ---
        % --- In solve_ee_method7, Step 2 ---
        
        % 1. Define the Tikhonov Regularization parameter (lambda)
        % This must be defined globally, perhaps in params.lambda_reg
        lambda_reg = 1e-9; 
        
        % 2. Check for Numerical Collapse (The Tikhonov intervention)
        if any(isnan(k1_new)) || any(isinf(k1_new)) % k1_new is the target vector k_prime_updated
            % If a crash occurred, DO NOT UPDATE THE COEFFICIENTS. 
            % Simply retain the previous coefficients (K_coef_prev) and increase damping.
            % This is a conservative measure to climb out of the sink region.
            K_coef = K_coef_prev;
            params.kdamp = max(params.kdamp * 0.5, 0.01); % Damp the policy update aggressively
            warning('Numerical instability detected. Retaining previous coefficients and increasing damping.');
        else
            % --- Stabilized Tikhonov Regression (The PS1 equivalent fix) ---
            
            X_matrix = X0; % basis_matrix (X0) 
            XTX = X_matrix' * X_matrix;
            XTY = X_matrix' * k1_new;
            
            % Add the penalty term to the diagonal of X'X
            XTX_reg = XTX + lambda_reg * eye(size(XTX));
            
            % Solve the Regularized Linear System (Tikhonov Regularization)
            K_coef_new = XTX_reg \ XTY;
            
            % Update the coefficients using damping
            K_coef = params.kdamp * K_coef_new + (1 - params.kdamp) * K_coef; 
        end

        % --- Check Convergence ---
        norm_diff = max(abs(K_coef - K_coef_prev));
        if norm_diff < tol
            break;
        end
    end
    
    % Final mean capital value on the grid (for outputting a descriptive result)
    k_mean_final = mean(k1_new); 
    
    if iter_count == max_iter
        warning('Solver did not converge within %d iterations.', max_iter);
    end
end