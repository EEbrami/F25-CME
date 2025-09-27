% Example Optimization Script for Computational Economics
% This script demonstrates basic numerical optimization techniques
% Course: Econ-81360, Fall 2025

%% Clear workspace
clear; clc; close all;

%% Example 1: Unconstrained Optimization
% Minimize f(x) = (x-2)^2 + 3
fprintf('Example 1: Unconstrained Optimization\n');
fprintf('Minimizing f(x) = (x-2)^2 + 3\n');

% Define objective function
f = @(x) (x-2)^2 + 3;

% Find minimum using fminunc
x0 = 0;  % Initial guess
options = optimoptions('fminunc', 'Display', 'iter');
[x_opt, f_opt] = fminunc(f, x0, options);

fprintf('Optimal x: %.4f\n', x_opt);
fprintf('Optimal value: %.4f\n', f_opt);

%% Example 2: Plot the function
x = linspace(-2, 6, 100);
y = arrayfun(f, x);

figure(1);
plot(x, y, 'b-', 'LineWidth', 2);
hold on;
plot(x_opt, f_opt, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
xlabel('x');
ylabel('f(x)');
title('Optimization Example: f(x) = (x-2)^2 + 3');
legend('f(x)', 'Minimum', 'Location', 'best');
grid on;

%% Example 3: Economic Application - Utility Maximization
fprintf('\nExample 3: Utility Maximization\n');
fprintf('Maximize U(x,y) = x^0.5 * y^0.5 subject to px*x + py*y = m\n');

% Parameters
px = 2;    % Price of good x
py = 1;    % Price of good y
m = 100;   % Income

% Utility function (we minimize negative utility)
utility = @(vars) -(vars(1)^0.5 * vars(2)^0.5);

% Constraint function: px*x + py*y - m = 0
constraint = @(vars) px*vars(1) + py*vars(2) - m;

% Initial guess
x0 = [20, 20];

% Solve using fmincon
options = optimoptions('fmincon', 'Display', 'iter');
[x_opt, u_opt] = fmincon(utility, x0, [], [], [], [], [0, 0], [], ...
                         @(vars) deal(constraint(vars), []), options);

fprintf('Optimal consumption:\n');
fprintf('x* = %.2f, y* = %.2f\n', x_opt(1), x_opt(2));
fprintf('Maximum utility: %.4f\n', -u_opt);
fprintf('Budget check: %.2f (should equal %.2f)\n', ...
        px*x_opt(1) + py*x_opt(2), m);

%% Save results
save('optimization_results.mat', 'x_opt', 'f_opt', 'u_opt');
fprintf('\nResults saved to optimization_results.mat\n');