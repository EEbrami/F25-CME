% Clear the workspace and command window
clear; clc;

%% problem2_solver.m
% Solves Problem 2, parts (a) and (b) of HW2: 2D Quadrature (Product Rules).
% Integral: I = integral_0^1 integral_0^1 exp(-xy) * (sin(6*pi*x) + cos(8*pi*y)) dx dy

% --- Setup: Ensure method functions are accessible and constants defined ---
addpath('./methods'); 
% 1. Define the 2D Integrand Function Handle
f_2d = @(x, y) exp(-x .* y) .* (sin(6 * pi * x) + cos(8 * pi * y));
% 2. Define Integration Domain
a_x = 0; b_x = 1; % Domain for x
a_y = 0; b_y = 1; % Domain for y
output_dir = './results/problem2/'; 

% -----------------------------------------------------------------
% --- Part 2(a) Code: 2D Gauss-Chebyshev Quadrature (Calculation) ---
% -----------------------------------------------------------------
% Constraint: up to degree 5 (D=5 -> N=3). We calculate multiple D values.
D_test_values = [3, 5, 7, 9, 11]; % Fixed D values (e.g., up to D=11)
Results_P2a = struct('D', {}, 'N', {}, 'I_approx', {});

fprintf('---------------------------------------------------\n');
fprintf('Problem 2(a) | 2D Gauss-Chebyshev Results\n');
fprintf('---------------------------------------------------\n');
for k = 1:length(D_test_values)
    D = D_test_values(k);
    N = (D + 1) / 2; % Derived N nodes required for degree D
    
    try
        I_approx = gc_quadrature_2d(f_2d, a_x, b_x, a_y, b_y, N);
    catch ME
        warning('Quadrature failed for D=%d (N=%d). Error: %s', D, N, ME.message);
        I_approx = NaN; 
    end

    Results_P2a(k).D = D;
    Results_P2a(k).N = N;
    Results_P2a(k).I_approx = I_approx;
    fprintf('D=%-3d (N=%-2d) | Approximation: %.15f\n', D, N, I_approx);
end

% --- Output Results to File (p2a_gc_2d_results.txt) ---
mkdir(output_dir); 
filepath_p2a = fullfile(output_dir, 'p2a_gc_2d_results.txt');
fid = fopen(filepath_p2a, 'wt'); 
if fid == -1, error('Could not open file %s.', filepath_p2a); end

fprintf(fid, '--- HW2 Problem 2(a) Results: 2D Gauss-Chebyshev Approximations ---\n');
fprintf(fid, 'Method: 2D Gauss-Chebyshev Type 1 Product Rule\n');
fprintf(fid, 'Domain: [%.1f, %.1f] x [%.1f, %.1f]\n\n', a_x, b_x, a_y, b_y);
fprintf(fid, 'Polynomial Degree (D) | Nodes per Dimension (N) | Total Nodes (N^2) | Approximation (I_approx)\n');
fprintf(fid, '----------------------|-------------------------|-------------------|---------------------------\n');
for k = 1:length(Results_P2a)
    fprintf(fid, '%-21d | %-23d | %-17d | %25.15f\n', ...
        Results_P2a(k).D, Results_P2a(k).N, Results_P2a(k).N^2, Results_P2a(k).I_approx);
end
fprintf(fid, '----------------------|-------------------------|-------------------|---------------------------\n');
fclose(fid); 
fprintf('\nResults saved successfully to: %s\n', filepath_p2a);


% -----------------------------------------------------------------
% --- Part 2(b) Code: 2D Trapezoid Rule (Calculation) ---
% -----------------------------------------------------------------
% Requirement: Use "different number of subperiods." We choose fixed values.
M_test_values = [4, 16, 64, 256, 512]; % Fixed M values (must be non-zero and reasonable)
Results_P2b = struct('M', {}, 'I_approx', {});

fprintf('\n---------------------------------------------------\n');
fprintf('Problem 2(b) | 2D Trapezoid Rule Results\n');
fprintf('---------------------------------------------------\n');

for k = 1:length(M_test_values)
    M = M_test_values(k); 
    
    try
        I_approx = trapz_rule_2d(f_2d, a_x, b_x, a_y, b_y, M);
    catch ME
        warning('Quadrature failed for M=%d. Error: %s', M, ME.message);
        I_approx = NaN; 
    end
    
    Results_P2b(k).M = M;
    Results_P2b(k).I_approx = I_approx;
    fprintf('M=%-5d | Total Nodes: %-7d | Approximation: %19.15f\n', M, M^2, I_approx);
end

% --- Output Results to File (p2b_trapz_results.txt) ---
filepath_p2b_txt = fullfile(output_dir, 'p2b_trapz_results.txt'); 
fid = fopen(filepath_p2b_txt, 'wt'); 
if fid == -1, error('Could not open file %s for writing.', filepath_p2b_txt); end

fprintf(fid, '--- HW2 Problem 2(b) Results: 2D Compound Trapezoid Rule ---\n');
fprintf(fid, 'Method: 2D Compound Trapezoid Rule (Product Rule)\n');
fprintf(fid, 'Domain: [%.1f, %.1f] x [%.1f, %.1f]\n\n', a_x, b_x, a_y, b_y);
fprintf(fid, 'M (Subperiods) | Total Nodes (M^2) | h (Step Size) | Approximation (I_approx)\n');
fprintf(fid, '---------------|-------------------|---------------|--------------------------\n');

for k = 1:length(Results_P2b)
    M = Results_P2b(k).M;
    h_step = (b_x - a_x) / M;
    
    fprintf(fid, '%-14d | %-17d | %13.5e | %26.15f\n', ...
        M, M^2, h_step, Results_P2b(k).I_approx);
end

fprintf(fid, '--------------------------------------------------------------------------------\n');
fclose(fid); 
fprintf('Trapezoid results saved to: %s\n', filepath_p2b_txt);

%%
% -----------------------------------------------------------------
% --- Part 2(c) Code: Monte Carlo Integration ---
% -----------------------------------------------------------------

% --- Monte Carlo Parameters ---
T_SAMPLES = 10000; % Mandated number of samples (T)
N_EXPERIMENTS = 5; % Experiment runs to observe variance (answering the question)
I_mc_estimates = zeros(N_EXPERIMENTS, 1);

fprintf('\n---------------------------------------------------\n');
fprintf('Problem 2(c) | Monte Carlo Integration Results\n');
fprintf('Mandated Samples (T): %d\n', T_SAMPLES);
fprintf('---------------------------------------------------\n');

% Execute the Monte Carlo simulation N_EXPERIMENTS times
for k = 1:N_EXPERIMENTS
    % The Area=1 is handled inside monte_carlo_2d.m
    try
        [I_mc_estimates(k), ~] = monte_carlo_2d(f_2d, a_x, b_x, a_y, b_y, T_SAMPLES);
        fprintf('Experiment %d | Approximation: %.15f\n', k, I_mc_estimates(k));
    catch ME
        warning('Monte Carlo failed on experiment %d. Error: %s', k, ME.message);
        I_mc_estimates(k) = NaN;
    end
end

% The final output is defined as the mean of the experiments
I_mc_final = mean(I_mc_estimates(~isnan(I_mc_estimates))); 

% --- Output Results to File (p2c_monte_carlo.txt) ---
filepath_p2c_txt = fullfile(output_dir, 'p2c_monte_carlo.txt'); 
fid = fopen(filepath_p2c_txt, 'wt'); 
if fid == -1, error('Could not open file %s for writing.', filepath_p2c_txt); end

fprintf(fid, '--- HW2 Problem 2(c) Results: Monte Carlo Integration ---\n');
fprintf(fid, 'Method: Crude Monte Carlo (2D)\n');
fprintf(fid, 'Domain Area: 1.0\n');
fprintf(fid, 'Samples per Run (T): %d\n\n', T_SAMPLES);

fprintf(fid, '--- Observational Experiment ---\n');
fprintf(fid, 'The integral was run %d times with T=10000 random points to facilitate observation.\n', N_EXPERIMENTS);
fprintf(fid, 'Run Index | Monte Carlo Approximation (I_MC)\n');
fprintf(fid, '----------|-----------------------------------\n');
for k = 1:N_EXPERIMENTS
    fprintf(fid, '%-9d | %28.15f\n', k, I_mc_estimates(k));
end

fprintf(fid, '----------|-----------------------------------\n');
fprintf(fid, 'Mean Estimate (I_MC_final): %28.15f\n\n', I_mc_final);

fprintf(fid, '--- Observation Note for P2d Analysis ---\n');
fprintf(fid, 'The primary observation data point is the variation in results across runs (variance), which will be analyzed in P2d.\n');
fprintf(fid, '--------------------------------------------------------------------------------------------------\n');

fclose(fid); 
fprintf('Monte Carlo results saved to: %s\n', filepath_p2c_txt);

%%
% -----------------------------------------------------------------
% --- Part 2(d) Code: Built-in Comparison and Specific Analysis ---
% -----------------------------------------------------------------

% --- 1. Baseline Calculation and Comparison Targets ---
a_x = 0; b_x = 1; 
a_y = 0; b_y = 1;
f_2d_handle = @(x, y) exp(-x .* y) .* (sin(6 * pi * x) + cos(8 * pi * y));

% Calculate BASELINE: Use dblquad as mandated by the prompt.
try
    I_dblquad_baseline = dblquad(@(x, y) f_2d_handle(x, y), a_x, b_x, a_y, b_y);
catch
    I_dblquad_baseline = NaN;
end

% Calculate Numerical Gold Standard (for reference and final row)
I_integral2_gold = integral2(f_2d_handle, a_x, b_x, a_y, b_y, 'AbsTol', 1e-12);

% --- 2. Recalculate Custom Results at MANDATED PARAMETERS ---
% A. GC at D=5 (N=3) - Mandated constraint from P2a
I_gc_D5 = gc_quadrature_2d(f_2d_handle, a_x, b_x, a_y, b_y, 3);
% B. Trapezoid at M=4 - Mandated subperiod minimum
I_trapz_M4 = trapz_rule_2d(f_2d_handle, a_x, b_x, a_y, b_y, 4);
% C. Monte Carlo (Final Mean Estimate) - From P2c output variable
I_mc_final_p2c = mean(I_mc_estimates(~isnan(I_mc_estimates)));

% --- 3. Generate Specific Comparison Table ---
Comparison_Params = struct('Method', ...
    {'dblquad (Base)', 'GCQ (D=5)', 'Trapezoid (M=4)', 'Monte(T=10k)', 'integral2'}, ...
    'I_approx', ...
    {I_dblquad_baseline, I_gc_D5, I_trapz_M4, I_mc_final_p2c, I_integral2_gold}, ...
    'Nodes_Used', ...
    {NaN, 9, 16, 10000, NaN});

output_dir = './results/problem2/';
filepath_p2d_txt = fullfile(output_dir, 'p2d_specific_comparison.txt'); 
fid = fopen(filepath_p2d_txt, 'wt'); 
if fid == -1, error('Could not open file %s.', filepath_p2d_txt); end

fprintf(fid, '--- HW2 Problem 2(d) Results: Accuracy Comparison at Mandated Parameters ---\n');
fprintf(fid, 'Comparison Baseline (dblquad Value): %22.15f\n\n', I_dblquad_baseline);
fprintf(fid, 'NOTE: Relative Error is calculated against the dblquad result, as mandated.\n\n');

fprintf(fid, 'Method            | I_approximate Value    | Nodes/Samples | Relative Error (%%)\n');
fprintf(fid, '------------------|------------------------|---------------|--------------------\n');

for i = 1:length(Comparison_Params)
    I_current = Comparison_Params(i).I_approx;
    
    if i == 1
        Rel_Error_Percent = NaN; % Baseline has 0% error against itself
        Node_str = '    N/A   ';
    else
        % Calculate Relative Error against dblquad (the baseline)
        Rel_Error = abs(I_dblquad_baseline - I_current) / abs(I_dblquad_baseline) * 100;
        Rel_Error_Percent = Rel_Error;
        Node_str = sprintf('%13d', Comparison_Params(i).Nodes_Used);
    end
    
    fprintf(fid, '%-17s | %22.15f | %-13s | %13e\n', ...
        Comparison_Params(i).Method, I_current, Node_str, Rel_Error_Percent);
end
fprintf(fid, '------------------|------------------------|---------------|--------------------\n');
fclose(fid); 
fprintf('\nSpecific comparison report saved to: %s\n', filepath_p2d_txt);

% -----------------------------------------------------------------
% --- Part 2(d) EXTENSION: Convergence Study and Separate Plots ---
% -----------------------------------------------------------------

% --- Baseline: Use integral2 as the TRUE Gold Standard for all plots ---
I_exact_2d = integral2(f_2d_handle, a_x, b_x, a_y, b_y, 'AbsTol', 1e-12); 
TOLERANCE = 1e-3; % Target tolerance for satisfaction marker

% Define Convergence Parameters (same as previous plan)
M_values_plot = [4, 8, 16, 32, 64, 128, 256, 512, 1024]; % Trapezoid Subperiods (M must be same for Simpson's later)
D_values_plot = 3:2:39; % GC Degrees (up to D=39, N=20)

% --- A. GC Quadrature Data Collection (Exponential Convergence) ---
Nodes_GC = [];
Rel_Error_GC = [];
D_satisfied_gc = NaN;

for D = D_values_plot
    N = (D + 1) / 2; 
    try
        I_gc_current = gc_quadrature_2d(f_2d_handle, a_x, b_x, a_y, b_y, N);
        Rel_Error = abs(I_exact_2d - I_gc_current) / abs(I_exact_2d);
        Nodes_GC = [Nodes_GC, N^2];
        Rel_Error_GC = [Rel_Error_GC, Rel_Error];
        
        if Rel_Error <= TOLERANCE && isnan(D_satisfied_gc)
            D_satisfied_gc = D;
        end
    catch
        % Skip failed GC points
    end
end

% --- B. Trapezoid Rule Data Collection (Algebraic Convergence) ---
Nodes_Trapz = [];
Rel_Error_Trapz = [];
M_satisfied_trapz = NaN;

for M = M_values_plot
    try
        I_trapz_current = trapz_rule_2d(f_2d_handle, a_x, b_x, a_y, b_y, M);
        Rel_Error = abs(I_exact_2d - I_trapz_current) / abs(I_exact_2d);
        Nodes_Trapz = [Nodes_Trapz, M^2];
        Rel_Error_Trapz = [Rel_Error_Trapz, Rel_Error];
        
        if Rel_Error <= TOLERANCE && isnan(M_satisfied_trapz)
            M_satisfied_trapz = M;
        end
    catch
        % Skip failed Trapz points
    end
end

% -----------------------------------------------------------------
% --- C. Plot 1: GC Convergence (Exponential/Low Nodes) ---
% -----------------------------------------------------------------
filepath_plot_gc = fullfile(output_dir, 'p2d_gc_convergence.png');

figure('Renderer', 'painters', 'Position', [100 100 850 600]); 
semilogy(Nodes_GC, Rel_Error_GC, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;

% Add tolerance line
plot([Nodes_GC(1), Nodes_GC(end)], [TOLERANCE, TOLERANCE], 'k--', 'DisplayName', 'Tolerance Threshold');

if ~isnan(D_satisfied_gc)
    % Find the exact point of convergence for the annotation
    N_conv_gc = (D_satisfied_gc + 1) / 2;
    Error_conv_gc = abs(I_exact_2d - gc_quadrature_2d(f_2d_handle, a_x, b_x, a_y, b_y, N_conv_gc)) / abs(I_exact_2d);
    
    % Plot and annotate the satisfactory step
    loglog(N_conv_gc^2, Error_conv_gc, 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', sprintf('Satisfied D=%d', D_satisfied_gc));
    plot_caption_gc = sprintf('Tolerance (%.2e) Satisfied at D=%d (N=%d)', TOLERANCE, D_satisfied_gc, N_conv_gc);
else
    plot_caption_gc = sprintf('Tolerance (%.2e) NOT satisfied up to D=%d', TOLERANCE, D_values_plot(end));
end

title({'2D GC Quadrature Convergence (Exponential)'; plot_caption_gc});
xlabel('Total Nodes (N^2)');
ylabel('Relative Error (log_{10} scale)');
legend('Location', 'southwest');
grid on;
hold off;
saveas(gcf, filepath_plot_gc);
close(gcf); 
fprintf('GC convergence plot saved to: %s\n', filepath_plot_gc);

% -----------------------------------------------------------------
% --- D. Plot 2: Trapezoid Convergence (Algebraic/High Nodes) ---
% -----------------------------------------------------------------
filepath_plot_trapz = fullfile(output_dir, 'p2d_trapz_convergence.png');

figure('Renderer', 'painters', 'Position', [100 100 850 600]); 
loglog(Nodes_Trapz, Rel_Error_Trapz, 'rs-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on;

% Add tolerance line
plot([Nodes_Trapz(1), Nodes_Trapz(end)], [TOLERANCE, TOLERANCE], 'k--', 'DisplayName', 'Tolerance Threshold');

if ~isnan(M_satisfied_trapz)
    % Find the exact point of convergence for the annotation
    I_conv_trapz = trapz_rule_2d(f_2d_handle, a_x, b_x, a_y, b_y, M_satisfied_trapz);
    Error_conv_trapz = abs(I_exact_2d - I_conv_trapz) / abs(I_exact_2d);
    
    % Plot and annotate the satisfactory step
    loglog(M_satisfied_trapz^2, Error_conv_trapz, 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', sprintf('Satisfied M=%d', M_satisfied_trapz));
    plot_caption_trapz = sprintf('Tolerance (%.2e) Satisfied at M=%d (Total Nodes=%d)', TOLERANCE, M_satisfied_trapz, M_satisfied_trapz^2);
else
    plot_caption_trapz = sprintf('Tolerance (%.2e) NOT satisfied up to M=%d', TOLERANCE, M_values_plot(end));
end

title({'2D Trapezoid Rule Convergence (Algebraic)'; plot_caption_trapz});
xlabel('Total Nodes (M^2) [log scale]');
ylabel('Relative Error (log_{10} scale)');
legend('Location', 'southwest');
grid on;
hold off;
saveas(gcf, filepath_plot_trapz);
close(gcf); 
fprintf('Trapezoid convergence plot saved to: %s\n', filepath_plot_trapz);

% --- Final data saving block (D. Save Convergence Data Table) remains unchanged ---