function results = solve_growth_model(params, grid_params, solver_params, method)
% SOLVE_GROWTH_MODEL solves the stochastic growth model via Algorithm 7.

% 1. Initialization (Algorithm 7, Step a-d)
[k_grid, z_grid, k_nodes, z_nodes] = create_grid(params, grid_params, method.name);
[z_shocks, weights] = qnwnorm(5, 0, params.sigma^2); % Gauss-Hermite nodes

% --- START OF FIX ---
% Define the domain struct needed by the Chebyshev evaluation function
domain.k_min = min(k_nodes);
domain.k_max = max(k_nodes);
domain.z_min = min(z_nodes);
domain.z_max = max(z_nodes);
% --- END OF FIX ---

% Initial guess for the policy function
k_prime_guess = (1 - params.delta) * k_grid;

if strcmp(method.name, 'chebyshev')
    basis_matrix = chebyshev_basis(k_nodes, z_nodes, method.params.degree);
    coeffs = basis_matrix \ k_prime_guess(:);
else % Spline
    policy_fun = griddedInterpolant({k_nodes, z_nodes}, k_prime_guess, 'cubic', 'linear');
end

% Iteration loop
dist = 1;
iter = 0;
while dist > solver_params.tol && iter < solver_params.max_iter
    iter = iter + 1;

    % 2. Computation of k' on the grid (Algorithm 7, Step 1)
    if strcmp(method.name, 'chebyshev')
        k_prime_on_grid = basis_matrix * coeffs;
        k_prime_on_grid = reshape(k_prime_on_grid, grid_params.k_points, grid_params.z_points);
    else % Spline
        k_prime_on_grid = policy_fun(k_grid, z_grid);
    end

    % RHS of Euler Equation (Expectation Term)
    RHS_euler = zeros(size(k_grid));
    for i_shock = 1:length(z_shocks)
        z_prime = z_grid.^params.rho * exp(z_shocks(i_shock));
        
        if strcmp(method.name, 'chebyshev')
            % --- BUG FIX: Pass the domain struct as the 5th argument ---
            k_prime_prime = chebyshev_eval_2d_p2(coeffs, method.params.degree, k_prime_on_grid, z_prime, domain);
        else % Spline
            k_prime_prime = policy_fun(k_prime_on_grid, z_prime);
        end
        
        c_prime = params.f(k_prime_on_grid, z_grid) + (1-params.delta)*k_prime_on_grid - k_prime_prime;
        c_prime(c_prime < 0) = 1e-8;
        
        marginal_utility = params.u_prime(c_prime);
        marginal_return = (1 - params.delta + params.f_prime(k_prime_on_grid, z_prime));
        
        RHS_euler = RHS_euler + weights(i_shock) * marginal_utility .* marginal_return;
    end
    RHS_euler = params.beta * RHS_euler;

    % Invert Euler to find updated c_t and k_t+1
    c_updated = params.u_prime_inv(RHS_euler);
    k_prime_updated = params.f(k_grid, z_grid) + (1-params.delta)*k_grid - c_updated;

    % 3. Update Policy Function (Algorithm 7, Step 2)
    if strcmp(method.name, 'chebyshev')
        coeffs_new = basis_matrix \ k_prime_updated(:);
        dist = max(abs(coeffs_new - coeffs));
        coeffs = (1 - solver_params.damping) * coeffs + solver_params.damping * coeffs_new; % Damping
    else % Spline
        dist = max(abs(policy_fun(k_grid, z_grid) - k_prime_updated), [], 'all');
        policy_fun.Values = (1-solver_params.damping)*policy_fun.Values + solver_params.damping*k_prime_updated;
    end
end

% Store results
results.iterations = iter;
if strcmp(method.name, 'chebyshev')
    results.policy_coeffs = coeffs;
    results.method_params = method.params;
    results.domain = domain; % Store domain for later use
else
    results.policy_fun = policy_fun;
end

% Evaluate accuracy of the solution
[results.log10_max_error, results.log10_mean_error] = evaluate_accuracy(results, params);
end