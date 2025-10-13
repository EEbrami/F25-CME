% Clear the workspace and command window
clear; clc;

%% problem1_solver.m

%%
% -----------------------------------------------------------------
% --- Part 1(a) Code: Integration by Parts ---
% -----------------------------------------------------------------
% Integral: I = integral_0^2pi (exp(-x) * sin(10x)) dx

% --- Setup: Ensure method functions are accessible ---
% The working directory is /hw2, and the function is in /hw2/methods/.
addpath('./methods'); 

% --- Part 1(a) Code (Exact Solution) ---
syms x
f_sym = exp(-x) * sin(10*x);
I_sym_expression = int(f_sym, x, 0, 2*pi);
I_exact = double(I_sym_expression);

% 4. Display Results (Console Output)
fprintf('---------------------------------------------------\n');
fprintf('Problem 1(a) Results: Integration by Parts\n');
fprintf('---------------------------------------------------\n');
disp(['Symbolic Expression (I_sym_expression): ', char(I_sym_expression)]);
fprintf('Exact Numerical Value (I_exact): %.15f\n\n', I_exact);


% 5. Output to Text File (PATH ADJUSTMENT)
output_dir = './results/problem1/'; % Path relative to the execution directory (/hw2)
p1a_base = 'p1a_integration_by_parts.txt';
filepath_p1a = fullfile(output_dir, p1a_base); % Constructs the full path robustly

fid = fopen(filepath_p1a, 'wt'); % Open file for writing text at the specified path

if fid == -1
    % Use mkdir to create the directory if it doesn't exist (robust coding)
    % NOTE: This assumes the 'results' folder exists, otherwise use mkdir(output_dir, 'p')
    [success, msg, msgID] = mkdir(output_dir);
    if success
        fid = fopen(filepath_p1a, 'wt'); % Try opening again after creation
        if fid == -1
            error('Could not create directory and open file %s. Error: %s', filepath_p1a, msg);
        end
    else
        error('Could not open file %s for writing. Check permissions or path existence. Error: %s', filepath_p1a, msg);
    end
end

% Write content
fprintf(fid, '--- HW2 Problem 1(a) Results: Integration by Parts ---\n');
fprintf(fid, 'Integral: I = integral_0^2pi (exp(-x) * sin(10x)) dx\n\n');
fprintf(fid, '1. Exact Symbolic Expression (I_sym_expression):\n');
fprintf(fid, '%s\n\n', char(I_sym_expression));
fprintf(fid, '2. Exact Numerical Value (I_exact, double precision):\n');
fprintf(fid, '%.15f\n', I_exact);

fclose(fid); % Close the file handle

fprintf('Results saved successfully to: %s\n', filepath_p1a);

%%
% -----------------------------------------------------------------
% --- Part 1(b) Code: Gauss-Chebyshev Quadrature ---
% -----------------------------------------------------------------

% 6. Method Parameters (REVISED to set D and derive N)
method_name = 'Gauss-Chebyshev Quadrature';
gc_type = 'Type 1 (Adapted for Unweighted Function)'; % Type 1 from slide 33
weight_function_source = 'w(x) = (1-x^2)^(-1/2) (Source Weight)'; % Source weight from slide 33
integral_type = 'Unweighted Integral on Arbitrary Domain [a, b]'; % Problem type from slide 34
% --- CONTROLLING PARAMETER ---
D_poly_order = 40; % The constraint from the HW: exact up to 5th order polynomial
% --- DERIVED PARAMETER ---
N_nodes = (D_poly_order + 1) / 2; % N = (D+1)/2 = 3 nodes
a = 0;        % Lower integration limit
b = 2*pi;     % Upper integration limit
degree_exactness = 2 * N_nodes - 1; % D = 2N-1, must match D_poly_order
f_handle = @(x) exp(-x) .* sin(10*x); 

% 7. Function Call and Result Storage
% NOTE: I_gc stores the real value for comparison later.
try
    % Placeholder: replace with the actual function call when gc_quadrature is available
    % I_gc = gc_quadrature(f_handle, a, b, N_nodes);
    I_gc = gc_quadrature(f_handle, a, b, N_nodes); % N_nodes is correctly set to 3.
catch ME
    warning('gc_quadrature.m call failed. Storing NaN for I_gc.');
    I_gc = NaN; 
end

% 8. Console Output for P1b
fprintf('\n\n---------------------------------------------------\n');
fprintf('Problem 1(b) Results: %s\n', method_name);
fprintf('---------------------------------------------------\n');
fprintf('Approximation (I_gc): %.15f\n', I_gc);


% Write detailed content to the file

output_dir = './results/problem1/'; % Path relative to the execution directory (/hw2)
p1b_base = 'p1b_Gauss-Chebyshev.txt';
filepath_p1b = fullfile(output_dir, p1b_base); % Constructs the full path robustly


fid = fopen(filepath_p1b, 'wt'); 

if fid == -1
    error('Could not open file %s for writing.', filepath_p1b);
end


fprintf(fid, '--- HW2 Problem 1(b) Results: Gauss-Chebyshev Quadrature ---\n');
fprintf(fid, 'Method Used: %s\n', method_name);
fprintf(fid, 'Method Type: %s\n', gc_type);
fprintf(fid, 'Integral: I = integral_0^2pi (exp(-x) * sin(10x)) dx\n\n');

fprintf(fid, '--- Parameters Used ---\n');
fprintf(fid, 'Integration Interval: [%.1f, %.3f]\n', a, b);
fprintf(fid, 'Source Weight Function: %s\n', weight_function_source); 
fprintf(fid, 'Integral Form Solved: %s\n', integral_type); 
fprintf(fid, 'Required Polynomial Order (D): %d\n', D_poly_order); % INPUT PARAMETER
fprintf(fid, 'Number of Quadrature Nodes (N): %d\n', N_nodes); % DERIVED PARAMETER
fprintf(fid, 'Theoretical Degree of Exactness: %d (2N-1)\n\n', degree_exactness); % DERIVED PARAMETER

fprintf(fid, '--- Results ---\n');
fprintf(fid, '1. Numerical Approximation (I_gc):\n');
fprintf(fid, '%.15f\n', I_gc);

fclose(fid); 
fprintf('Results saved successfully to: %s\n', filepath_p1b);

% -----------------------------------------------------------------
% --- Part 1(b) EXTENSION: Dynamic Convergence Study and Plotting ---
% -----------------------------------------------------------------

% --- Setup Parameters for Convergence Loop ---
TOLERANCE = 1e-3; % Define "good enough" as 0.1% relative error
D_start = 3;      % Start with the minimum required degree of D=3
D_poly_orders = [];
Relative_Errors = [];
D_last_satisfied = D_start;

% Initial parameters are already defined: I_exact, f_handle, a, b

% --- 1. Iteratively find N and D until tolerance is met ---
D = D_start;
MAX_ITER = 100; 

for iter = 1:MAX_ITER
    % Derive N (Number of Nodes) from the current odd Degree (D)
    N_current = (D + 1) / 2; 

    try
        % Execute the Gauss-Chebyshev calculation
        I_gc_current = gc_quadrature(f_handle, a, b, N_current); 
    catch ME
        fprintf('Error: gc_quadrature failed for D=%d (N=%d). Halting loop.\n', D, N_current);
        break; 
    end
    
    % Calculate Relative Error: Abs(Exact - Approx) / Abs(Exact)
    Rel_Error = abs(I_exact - I_gc_current) / abs(I_exact);

    % Store results for plotting
    D_poly_orders = [D_poly_orders, D];
    Relative_Errors = [Relative_Errors, Rel_Error];

    % Check for tolerance condition
    if Rel_Error <= TOLERANCE
        D_last_satisfied = D;
        fprintf('\nConvergence met at D=%d (N=%d). Relative Error: %.2e < %.2e\n', D, N_current, Rel_Error, TOLERANCE);
        break;
    end

    % Increment D by 2 to maintain the odd polynomial degree sequence (3, 5, 7, ...)
    D = D + 2; 

    if iter == MAX_ITER
        fprintf('Warning: Max iterations reached without meeting tolerance of %.2e\n', TOLERANCE);
    end
end

% --- 2. Plotting the Convergence (D vs. Relative Error) ---

% Use existing output directory structure
plot_dir = './results/problem1/'; 
plot_filename = 'p1b_convergence_D_vs_Error.png';
filepath_plot = fullfile(plot_dir, plot_filename);

figure('Renderer', 'painters', 'Position', [100 100 800 500]); % Setup figure
semilogy(D_poly_orders, Relative_Errors, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');

% Final caption for the plot title
plot_caption = sprintf('Tolerance (%.2e) Satisfied at D = %d', TOLERANCE, D_last_satisfied);

title({'Gauss-Chebyshev Convergence: Required D vs. Relative Error'; plot_caption});
xlabel('Polynomial Degree (D)');
ylabel('Relative Error (log_{10} scale)');
grid on;

% Add a horizontal line for the tolerance threshold
hold on;
plot([D_poly_orders(1), D_poly_orders(end)], [TOLERANCE, TOLERANCE], 'r--', 'DisplayName', 'Tolerance Threshold');
legend('Convergence Path', 'Tolerance Threshold', 'Location', 'southwest');
hold off;

% 3. Saving the Plot
saveas(gcf, filepath_plot);
close(gcf);

fprintf('Convergence study complete. Plot saved to: %s\n', filepath_plot);

%%
% -----------------------------------------------------------------
% --- Part 1(c) Block 1: Trapezoid Rule Convergence Study ---
% -----------------------------------------------------------------
TOLERANCE = 1e-3; % Tolerance defined as 0.1% relative error
M_start = 2;      % Start with 2 subperiods
M_periods_trapz = [];
Error_trapz_rel = [];
M_last_satisfied_trapz = M_start;

fprintf('\n---------------------------------------------------\n');
fprintf('Problem 1(c) | Trapezoid Rule Convergence\n');
fprintf('Target Tolerance: %.2e\n', TOLERANCE);
fprintf('---------------------------------------------------\n');

M = M_start;
MAX_M = 2^15; % Max M set to 32768 for safety

while M <= MAX_M
    % Execute the Trapezoid Rule calculation
    try
        I_trapz_current = trapezoid_rule(f_handle, a, b, M);
    catch ME
        fprintf('Error: Trapezoid Rule failed for M=%d. Halting loop.\n', M);
        break; 
    end
    
    % Calculate Relative Error: Abs(Exact - Approx) / Abs(Exact)
    Rel_Error = abs(I_exact - I_trapz_current) / abs(I_exact);

    % Store results for plotting and documentation
    M_periods_trapz = [M_periods_trapz, M];
    Error_trapz_rel = [Error_trapz_rel, Rel_Error];

    % Check for tolerance condition
    if Rel_Error <= TOLERANCE
        M_last_satisfied_trapz = M;
        fprintf('M=%-5d | Relative Error: %10.3e (CONVERGED)\n', M, Rel_Error);
        break;
    end
    
    fprintf('M=%-5d | Relative Error: %10.3e\n', M, Rel_Error);

    % Double M for exponential increase
    M = M * 2; 

    if M > MAX_M
        fprintf('Warning: Max iterations reached without meeting tolerance.\n');
    end
end

% --- Plotting the Convergence (Log-Log Scale to show Order) ---
plot_dir = './results/problem1/'; 
plot_filename_trapz = 'p1c_trapezoid_convergence.png';
filepath_plot_trapz = fullfile(plot_dir, plot_filename_trapz);

figure('Renderer', 'painters', 'Position', [100 100 800 500]); 
loglog(M_periods_trapz, Error_trapz_rel, 'rs-', 'LineWidth', 2, 'MarkerFaceColor', 'r');

plot_caption_trapz = sprintf('Tolerance (%.2e) Satisfied at M = %d', TOLERANCE, M_last_satisfied_trapz);

title({'Trapezoid Rule Convergence (Algebraic: O(h^2))'; plot_caption_trapz});
xlabel('Number of Subperiods (M) [log scale]');
ylabel('Relative Error (log_{10} scale)');

hold on;
% Add tolerance line
plot([M_periods_trapz(1), M_periods_trapz(end)], [TOLERANCE, TOLERANCE], 'k--', 'DisplayName', 'Tolerance Threshold');
legend('Convergence Path', 'Tolerance Threshold', 'Location', 'southwest');
grid on;
hold off;

saveas(gcf, filepath_plot_trapz);
close(gcf); 
fprintf('Trapezoid plot saved to: %s\n', filepath_plot_trapz);

% -----------------------------------------------------------------
% --- Part 1(c) Block 2: Simpson's Rule Convergence Study ---
% -----------------------------------------------------------------
TOLERANCE = 1e-3; 
M_start = 2;      
M_periods_simpson = [];
Error_simpson_rel = [];
M_last_satisfied_simpson = M_start;

fprintf('\n---------------------------------------------------\n');
fprintf('Problem 1(c) | Simpson''s Rule Convergence\n');
fprintf('Target Tolerance: %.2e\n', TOLERANCE);
fprintf('---------------------------------------------------\n');

M = M_start;
MAX_M = 2^15; 

while M <= MAX_M
    % M is guaranteed to be even (starts at 2, doubles)

    % Execute the Simpson's Rule calculation
    try
        I_simpson_current = simpsons_rule(f_handle, a, b, M);
    catch ME
        fprintf('Error: Simpson''s Rule failed for M=%d. Halting loop.\n', M);
        break; 
    end
    
    % Calculate Relative Error
    Rel_Error = abs(I_exact - I_simpson_current) / abs(I_exact);

    % Store results for plotting and documentation
    M_periods_simpson = [M_periods_simpson, M];
    Error_simpson_rel = [Error_simpson_rel, Rel_Error];

    % Check for tolerance condition
    if Rel_Error <= TOLERANCE
        M_last_satisfied_simpson = M;
        fprintf('M=%-5d | Relative Error: %10.3e (CONVERGED)\n', M, Rel_Error);
        break;
    end
    
    fprintf('M=%-5d | Relative Error: %10.3e\n', M, Rel_Error);

    % Double M for exponential increase
    M = M * 2; 

    if M > MAX_M
        fprintf('Warning: Max iterations reached without meeting tolerance.\n');
    end
end

% --- Plotting the Convergence (Log-Log Scale to show Order) ---
plot_dir = './results/problem1/'; 
plot_filename_simpson = 'p1c_simpsons_convergence.png';
filepath_plot_simpson = fullfile(plot_dir, plot_filename_simpson);

figure('Renderer', 'painters', 'Position', [100 100 800 500]); 
loglog(M_periods_simpson, Error_simpson_rel, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');

plot_caption_simpson = sprintf('Tolerance (%.2e) Satisfied at M = %d', TOLERANCE, M_last_satisfied_simpson);

title({'Simpson''s Rule Convergence (Algebraic: O(h^4))'; plot_caption_simpson});
xlabel('Number of Subperiods (M) [log scale]');
ylabel('Relative Error (log_{10} scale)');

hold on;
% Add tolerance line
plot([M_periods_simpson(1), M_periods_simpson(end)], [TOLERANCE, TOLERANCE], 'k--', 'DisplayName', 'Tolerance Threshold');
legend('Convergence Path', 'Tolerance Threshold', 'Location', 'southwest');
grid on;
hold off;

saveas(gcf, filepath_plot_simpson);
close(gcf); 
fprintf('Simpson''s plot saved to: %s\n', filepath_plot_simpson);

% -----------------------------------------------------------------
% --- Part 1(c) Block 3: Text File Output for Convergence Profile
% -----------------------------------------------------------------

output_dir = './results/problem1/'; % Path relative to the execution directory (/hw2)
p1c_base = 'p1c_compound_rules.txt';
filepath_p1c = fullfile(output_dir, p1c_base); 

fid = fopen(filepath_p1c, 'wt'); 
if fid == -1
    % Create directory if needed (robust coding check)
    [success, msg, msgID] = mkdir(output_dir);
    if success
        fid = fopen(filepath_p1c, 'wt'); 
        if fid == -1
            error('Could not create directory and open file %s for analysis report.', filepath_p1c);
        end
    else
        error('Could not open file %s for analysis report. Error: %s', filepath_p1c, msg);
    end
end

% --- 1. Write General Study Parameters ---
fprintf(fid, '--- HW2 Problem 1(c) Results: Trapezoid and Simpson Convergence Study ---\n');
fprintf(fid, 'Integral: I = integral_0^2pi (exp(-x) * sin(10x)) dx\n\n');
fprintf(fid, '--- Study Parameters ---\n');
fprintf(fid, 'Integration Interval: [%.1f, %.3f]\n', a, b);
fprintf(fid, 'Exact Value (I_exact): %.15f\n', I_exact);
fprintf(fid, 'Required Relative Tolerance: %.2e\n\n', TOLERANCE);

% --- 2. Write Final Convergence Results ---
fprintf(fid, '--- Final Convergence Points ---\n');
fprintf(fid, 'Trapezoid Rule converged at M (Subperiods): %d\n', M_last_satisfied_trapz);
fprintf(fid, 'Simpson''s Rule converged at M (Subperiods): %d\n\n', M_last_satisfied_simpson);

% --- 3. Write Tabulated Convergence Profile for Analysis ---
fprintf(fid, '--- Convergence Profile (M vs. Relative Error) ---\n');
fprintf(fid, 'This table documents the M-h relationship for order of convergence analysis.\n');
fprintf(fid, 'M (Subperiods) | h (Step Size) | Trapz Rel. Error | Simpson Rel. Error\n');
fprintf(fid, '----------------------------------------------------------------------------\n');

% NOTE: The loop uses the full length of the Trapezoid array for completeness.
n_rows = length(M_periods_trapz); 
h_step = (b - a) ./ M_periods_trapz;

for k = 1:n_rows
    % Use M_periods_trapz as the master M list
    
    % Simpson's may have converged earlier, so its array might be shorter.
    if k <= length(Error_simpson_rel)
        simpson_error_str = sprintf('%16.3e', Error_simpson_rel(k));
    else
        % If Simpson's converged earlier, its subsequent error is effectively zero (or below machine epsilon)
        simpson_error_str = '   < TOLERANCE   '; 
    end
    
    fprintf(fid, '%-14d | %13.5e | %16.3e | %s\n', ...
        M_periods_trapz(k), h_step(k), Error_trapz_rel(k), simpson_error_str);
end

fprintf('----------------------------------------------------------------------------\n');
fclose(fid); 
fprintf('Convergence profile saved successfully to: %s\n', filepath_p1c);

%%
% -----------------------------------------------------------------
% --- Part 1(d) Code: Built-in Comparison and Final Report ---
% -----------------------------------------------------------------

% 1. Define Standard Inputs (re-used from previous blocks)
a = 0;
b = 2*pi;
f_handle = @(x) exp(-x) .* sin(10*x); 

% 2. Recalculate Custom Method Results (Robustness Check)

% A. Gauss-Chebyshev (GC): Compare HW Mandate vs. Optimal Accuracy (from P1b extension)
D_gc_optimal = D_last_satisfied; % The degree found by the P1b extension loop (e.g., 41)
N_gc_optimal = (D_gc_optimal + 1) / 2; % The required N for that D (e.g., 21)
D_gc_hw_mandate = 5;
N_gc_hw_mandate = 3;

I_gc_optimal_final = gc_quadrature(f_handle, a, b, N_gc_optimal);
I_gc_hw_mandate_final = gc_quadrature(f_handle, a, b, N_gc_hw_mandate);

% B. Compound Rules (P1c Tolerance)
I_trapz_manual_final = trapezoid_rule(f_handle, a, b, M_last_satisfied_trapz);
I_simpson_manual_final = simpsons_rule(f_handle, a, b, M_last_satisfied_simpson);

% C. Built-in Adaptive Quadrature (High Accuracy Standard)
tol_high = 1e-9;
I_integral_high = integral(f_handle, a, b, 'AbsTol', tol_high);
I_quad = quad(f_handle, a, b); 

% --- 3. Consolidate ALL Results for Final Table ---

% D. Built-in Trapezoid (Comparison with your manual method)
M_trapz_final = M_last_satisfied_trapz; 
x_nodes_final = linspace(a, b, M_trapz_final + 1);
I_trapz_builtin = trapz(x_nodes_final, f_handle(x_nodes_final)); 

% Retrieve existing custom results (variables must be accessible)
I_gc_optimal_final = gc_quadrature(f_handle, a, b, N_gc_optimal);
I_gc_hw_mandate_final = gc_quadrature(f_handle, a, b, N_gc_hw_mandate);
I_trapz_manual_final = trapezoid_rule(f_handle, a, b, M_last_satisfied_trapz);
I_simpson_manual_final = simpsons_rule(f_handle, a, b, M_last_satisfied_simpson);
I_integral_high = integral(f_handle, a, b, 'AbsTol', 1e-9);
I_quad = quad(f_handle, a, b);

% Define the full array of results
Final_Results = struct('Method', {}, 'I_approx', {}, 'Relative_Error', {}, 'Params', {});

idx = 1;

% Row 1: Exact Value (Benchmark)
Final_Results(idx) = struct('Method', 'Exact Solution (P1a)', 'I_approx', I_exact, 'Relative_Error', 0.0, 'Params', 'N/A');
idx = idx + 1;

% Row 2: Your GC (P1b HW Mandate - LOW ACCURACY)
Final_Results(idx) = struct('Method', 'GC Quadrature (D=5 poly.)', ...
                          'I_approx', I_gc_hw_mandate_final, ...
                          'Relative_Error', abs(I_exact - I_gc_hw_mandate_final) / abs(I_exact), ...
                          'Params', 'D=5, N=3. No direct MATLAB counterpart.');
idx = idx + 1;

% Row 3: Your GC (P1b Extension - OPTIMAL ACCURACY)
Final_Results(idx) = struct('Method', 'GC Quadrature (Optimized)', ...
                          'I_approx', I_gc_optimal_final, ...
                          'Relative_Error', abs(I_exact - I_gc_optimal_final) / abs(I_exact), ...
                          'Params', sprintf('D=41, N=21 (Tol=1e-3). No direct MATLAB counterpart.'));
idx = idx + 1;

% Row 4: Your Trapezoid (P1c Tolerance)
Final_Results(idx) = struct('Method', 'Trapz Rule (Custom)', 'I_approx', I_trapz_manual_final, 'Relative_Error', abs(I_exact - I_trapz_manual_final) / abs(I_exact), 'Params', sprintf('M=%d (Tol=1e-3). Direct MATLAB counterpart exists.', M_last_satisfied_trapz));
idx = idx + 1;

% Row 5: Built-in Trapz (Direct Comparison)
Final_Results(idx) = struct('Method', 'Trapz (Built-in)', 'I_approx', I_trapz_builtin, 'Relative_Error', abs(I_exact - I_trapz_builtin) / abs(I_exact), 'Params', sprintf('M=%d (Same Nodes as Custom). Verifies Custom Code.', M_last_satisfied_trapz));
idx = idx + 1;

% Row 6: Your Simpson's (P1c Tolerance)
Final_Results(idx) = struct('Method', 'Simpson''s Rule (Custom)', 'I_approx', I_simpson_manual_final, 'Relative_Error', abs(I_exact - I_simpson_manual_final) / abs(I_exact), 'Params', sprintf('M=%d (Tol=1e-3). Compared to adaptive solvers.', M_last_satisfied_simpson));
idx = idx + 1;

% Row 7: MATLAB's 'integral' (Numerical Gold Standard)
Final_Results(idx) = struct('Method', 'Integral (AbsTol=1e-9)', 'I_approx', I_integral_high, 'Relative_Error', abs(I_exact - I_integral_high) / abs(I_exact), 'Params', 'Adaptive Gauss-Kronrod Method. Modern MATLAB Standard.');
idx = idx + 1;

% Row 8: MATLAB's 'quad' (Old Adaptive)
Final_Results(idx) = struct('Method', 'Quad (Default)', 'I_approx', I_quad, 'Relative_Error', abs(I_exact - I_quad) / abs(I_exact), 'Params', 'Adaptive Simpson''s Rule. Older Method (Tol=1e-6 Default).');
% End of replacement


% --- 4. Write Final Text Report (p1d_final_comparison.txt) ---

output_dir = './results/problem1/';
p1d_base = 'p1d_final_comparison.txt';
filepath_p1d = fullfile(output_dir, p1d_base);

fid = fopen(filepath_p1d, 'wt'); 
if fid == -1
    error('Could not open file %s.', filepath_p1d);
end

fprintf(fid, '--- HW2 Problem 1(d) Results: Cross-Method Comparison ---\n');
fprintf(fid, 'Integral: I = integral_0^2pi (exp(-x) * sin(10x)) dx\n\n');

fprintf(fid, 'Baseline Exact Value (P1a): %18.15f\n\n', I_exact);
fprintf(fid, 'Comparison Metric: Relative Error = |I_exact - I_approx| / |I_exact|\n\n');

fprintf(fid, 'Method                      | Approximation Value    | Relative Error (%%) | Parameters/Constraints\n');
fprintf(fid, '----------------------------|------------------------|--------------------|----------------------------------------\n');

for i = 1:length(Final_Results)
    Rel_Err_Percent = Final_Results(i).Relative_Error * 100;

    fprintf(fid, '%-27s | %22.15f | %18.5e | %s\n', ...
        Final_Results(i).Method, ...
        Final_Results(i).I_approx, ...
        Rel_Err_Percent, ...
        Final_Results(i).Params);
end

fprintf(fid, '----------------------------|------------------------|--------------------|----------------------------------------\n');
fclose(fid); 
fprintf('Final comparison report saved to: %s\n', filepath_p1d);
