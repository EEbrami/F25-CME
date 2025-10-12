
## Modular Code Structure for HW2 (Integration)

This blueprint eliminates the monolithic file issue experienced in HW1 by logically grouping files into layers based on their responsibility.

### I. Root Directory: `~/homeworks/hw2/`

This is the central execution directory. It contains only the master execution script and high-level documentation.

| File Name | Responsibility (Tier 1: Orchestration) |
| :--- | :--- |
| **`run_hw2.m`** | **Master Script.** Clears workspace, adds paths, loads parameters, and calls the two main solver functions (`solve_prob1_deterministic`, `solve_prob3_stochastic`) in sequence. |
| `25F - PS2.md`, `hw2_overview.md`, `surface_plot_description.png` | Assignment Documentation. |

***

### II. Problem Logic & Solvers: `~/homeworks/hw2/solvers/`

This folder holds the high-level logic and flow for solving each problem set objective.

| File Name | Functional Role (Tier 2: Problem Flow) |
| :--- | :--- |
| **`solve_prob1_deterministic.m`** | Driver script for Problems 1 and 2 (1D and 2D deterministic integration benchmarks). Loops through different $N$ values and methods (Trapezoid, Simpson, Chebyshev) for comparison. |
| **`solve_prob3_stochastic.m`** | **Driver for Problem 3 (NSGM).** Implements the main GSSA/EEM convergence loop. Calls primitive functions for approximation (from HW1's legacy functions) and the integration routines from `/lib_primitives/`. |
| `model_equations_nsgm.m` | Defines the actual economic functions (e.g., utility, production, Euler equation residual $R(k, z, \alpha)$). |

***

### III. Numerical Primitives: `~/homeworks/hw2/lib_primitives/`

This is your toolbox of reusable mathematical functions.

| File Name | Functional Role (Tier 3: Reusable Tools) |
| :--- | :--- |
| **`setup_params.m`** | Defines all model and numerical constants ($\beta, \delta, N, \sigma$). **(Moved from root)** |
| **`calculate_expected_value.m`** | **Abstraction Layer:** Contains the logic to perform the conditional expectation integral, calling one of the three helper files based on a method input flag (e.g., `method='M1'`). |
| `GH_Quadrature.m` | **(Provided Helper, Moved)** Gauss-Hermite quadrature routine. |
| `Monomials_1.m` | **(Provided Helper, Moved)** Monomial Rule 1 for stochastic integration. |
| `Monomials_2.m` | **(Provided Helper, Moved)** Monomial Rule 2 for stochastic integration. |
| `nodes_chebyshev.m` | Calculates the Gauss-Lobatto nodes and weights (needed for both P1 approximation checks and P3 basis function setup). |
| `quad_trapezoid_1d.m` | Custom implementation of the 1D composite trapezoid rule (for P1 accuracy checks). |

***

### IV. Results and Visualization: `~/homeworks/hw2/results/`

This folder is created explicitly for output and visualization code, fulfilling your requirement to see every step.

| File Name | Functional Role (Tier 4: Visualization Hooks) |
| :--- | :--- |
| **`plot_node_distribution.m`** | **Visual Step 1:** Plots the distribution of nodes (e.g., showing non-uniform Chebyshev clustering or the sparse Monomial rule points). |
| **`plot_error_convergence.m`** | Plots the $\log_{10}(\text{Error})$ vs. $\log_{10}(N)$ curve for Problem 1's quadrature comparison. |
| **`plot_error_surface_p3.m`** | Plots the final Euler error surface to compare the spatial accuracy of the GH, M1, and M2 integration methods across the state space. |
| `output_results_table.m` | Generates the final LaTeX or text output tables required for submission. |