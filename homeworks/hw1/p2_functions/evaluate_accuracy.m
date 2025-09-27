function [log10_max_err, log10_mean_err] = evaluate_accuracy(results, params)
% EVALUATE_ACCURACY simulates the model to find Euler equation errors.

T = 10200; % Total simulation length
T_burn = 200; % Burn-in period

% Draw shocks
rng(2025); % for reproducibility
shocks = randn(T, 1) * params.sigma;

% Initialize simulation vectors
k_sim = zeros(T, 1);
z_sim = zeros(T, 1);
k_sim(1) = params.k_ss;
z_sim(1) = params.z_ss;

% Get policy function
if isfield(results, 'policy_coeffs') % Chebyshev
    coeffs = results.policy_coeffs;
    degree = results.method_params.degree;
    % --- START OF FIX ---
    % The domain struct is now required as the 5th argument
    domain = results.domain;
    policy_fun = @(k,z) chebyshev_eval_2d_p2(coeffs, degree, k, z, domain);
    % --- END OF FIX ---
else % Spline
    policy_fun = results.policy_fun;
end

% Simulate the model
for t = 1:T-1
    z_sim(t+1) = z_sim(t)^params.rho * exp(shocks(t));
    k_sim(t+1) = policy_fun(k_sim(t), z_sim(t));
end

% Calculate Euler Errors on the simulated path (excluding burn-in)
k = k_sim(T_burn:end-1);
z = z_sim(T_burn:end-1);
k_prime = k_sim(T_burn+1:end);

c = params.f(k, z) + (1-params.delta)*k - k_prime;

% Integrate to find RHS of Euler equation
[z_shocks_acc, weights_acc] = qnwnorm(10, 0, params.sigma^2); % More accurate quadrature
RHS_euler = zeros(size(k));
for i_shock = 1:length(z_shocks_acc)
    z_prime = z.^params.rho * exp(z_shocks_acc(i_shock));
    k_prime_prime = policy_fun(k_prime, z_prime);
    
    c_prime = params.f(k_prime, z_prime) + (1-params.delta)*k_prime - k_prime_prime;
    c_prime(c_prime < 0) = 1e-8;
    
    marginal_utility = params.u_prime(c_prime);
    marginal_return = (1 - params.delta + params.f_prime(k_prime, z_prime));
    
    RHS_euler = RHS_euler + weights_acc(i_shock) * marginal_utility .* marginal_return;
end
RHS_euler = params.beta * RHS_euler;

% Euler error is the proportional difference
euler_error = abs(RHS_euler ./ params.u_prime(c) - 1);

log10_max_err = log10(max(euler_error));
log10_mean_err = log10(mean(euler_error));
end