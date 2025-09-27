%% ========================================================================
% PS2, ECON-81360: Solving the Stochastic Growth Model
% Author: Ebrahim Ebrami
% Date: Fall 2025
%
% This script solves the standard neoclassical stochastic growth model using
% an Euler equation method. It compares two approximation methods and
% generates plots of the resulting policy functions.
% =========================================================================

clear; clc; close all;
addpath('p2_functions');

% --- ADD THIS BLOCK ---
% Add the CompEcon Toolbox to the MATLAB Path
addpath(genpath('CompEcon'));
fprintf('CompEcon Toolbox added to path.\n\n');
% --- END OF BLOCK ---

% Create results directory for figures
if ~exist('results/figures_problem2', 'dir')
    mkdir('results/figures_problem2');
end

fprintf('=== Problem Set 2: Solving the Stochastic Growth Model ===\n\n');

%% 1. Model Setup and Steady State
params.alpha   = 0.36;      % Capital share
params.beta    = 0.99;      % Discount factor
params.delta   = 0.025;     % Depreciation rate
params.gamma   = 1;         % CRRA parameter (gamma=1 -> log utility)
params.rho     = 0.95;      % TFP persistence
params.sigma   = 0.01;      % TFP shock std. dev.

% Anonymous functions for utility and production
if params.gamma == 1
    params.u = @(c) log(c);
    params.u_prime = @(c) 1./c;
    params.u_prime_inv = @(x) 1./x;
else
    params.u = @(c) c.^(1-params.gamma) / (1-params.gamma);
    params.u_prime = @(c) c.^(-params.gamma);
    params.u_prime_inv = @(x) x.^(-1/params.gamma);
end

% Find normalizing constant A for k_ss = 1
params.k_ss = 1;
params.z_ss = 1; % Steady state TFP is 1
params.A = (1/params.beta - (1-params.delta)) / params.alpha;

% Production function and its derivative
params.f = @(k,z) params.A * z .* k.^params.alpha;
params.f_prime = @(k,z) params.A * params.alpha * z .* k.^(params.alpha-1);
fprintf('Step 1: Model parameters defined and steady state calibrated.\n');
fprintf('         Normalizing constant A = %.4f for k_ss = 1.\n\n', params.A);

%% 2. Define Grids and Approximation Parameters
% Grid setup
grid_params.k_points = 20; % Number of grid points for capital
grid_params.z_points = 10; % Number of grid points for TFP
params.k_range = 0.2; % +/- 20% deviation from k_ss for grid (now in params)
params.z_range = 3 * params.sigma; % +/- 3 std devs for log(z) (now in params)

% Chebyshev approximation setup
cheb_params.degree = 5;

% Solver parameters
solver_params.max_iter = 1000;
solver_params.tol = 1e-7;
solver_params.damping = 0.1;

%% 3. Solve the Model for Each Approximation Method
fprintf('--- Solving the model using different approximation methods ---\n');

% --- Method 1: Chebyshev Approximation ---
fprintf('\nRunning solver for Chebyshev polynomials (Degree %d)...\n', cheb_params.degree);
method_cheb.name = 'chebyshev';
method_cheb.params = cheb_params;
results_cheb = solve_growth_model(params, grid_params, solver_params, method_cheb);
fprintf('Solver converged in %d iterations.\n', results_cheb.iterations);
fprintf('Log10 Max Euler Error: %.4f\n', results_cheb.log10_max_error);
fprintf('Log10 Mean Euler Error: %.4f\n', results_cheb.log10_mean_error);
plot_policy_function(results_cheb, params, 'Chebyshev'); % Call plotting function

% --- Method 2: Spline Approximation ---
fprintf('\nRunning solver for Cubic Splines...\n');
method_spline.name = 'spline';
results_spline = solve_growth_model(params, grid_params, solver_params, method_spline);
fprintf('Solver converged in %d iterations.\n', results_spline.iterations);
fprintf('Log10 Max Euler Error: %.4f\n', results_spline.log10_max_error);
fprintf('Log10 Mean Euler Error: %.4f\n', results_spline.log10_mean_error);
plot_policy_function(results_spline, params, 'Spline'); % Call plotting function

%% 4. Compare Results
fprintf('\n\n--- SUMMARY OF RESULTS ---\n');
fprintf('%-20s | %-15s | %-15s\n', 'Approximation Method', 'Log10 Max Error', 'Log10 Mean Error');
fprintf('------------------------------------------------------------\n');
fprintf('%-20s | %-15.4f | %-15.4f\n', 'Chebyshev', results_cheb.log10_max_error, results_cheb.log10_mean_error);
fprintf('%-20s | %-15.4f | %-15.4f\n', 'Cubic Spline', results_spline.log10_max_error, results_spline.log10_mean_error);

if results_cheb.log10_max_error < results_spline.log10_max_error
    fprintf('\nConclusion: Chebyshev approximation performed best.\n');
else
    fprintf('\nConclusion: Spline approximation performed best.\n');
end

fprintf('\nProblem 2 complete!\n');