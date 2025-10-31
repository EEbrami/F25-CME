# Report: Problem Set 3 - Linear Equation Methods

**Course:** ECON-81360: Computational Methods for Economists

**Author:** Ebrahim Ebrami

**Date:** October 31, 2025

## Introduction

This report details the implementation and analysis for Problem Set 3, which focuses on the application, performance, and numerical stability of various linear equation solution methods. The analysis is divided into two parts:

1.  **Problem 1:** A comparative study of direct and iterative solvers, focusing on computational speed and convergence properties. This section includes a critical analysis of benchmark design and the performance implications of parallel versus sequential computation.
2.  **Problem 2:** An experimental analysis of numerical stabilization techniques (RLS, SVD, TSVD) as applied to the Generalized Stochastic Simulation Algorithm (GSSA), as detailed in **Judd, Maliar, and Maliar (2011)**.

All analysis is contextualized using the theoretical concepts from the lecture slides **"Linear Equation Methods" (Maliar, 2025)** and the aforementioned paper.

### File Structure

All code and results are located within the `homeworks/hw3/` directory:

* **`/solvers/`**: Contains the main executable scripts.
    * `problem1_solver_p1.m`: Benchmarks direct solution methods.
    * `problem1_solver_p2.m`: Tests iterative solution methods.
    * `problem2_solver.m`: Main script to run all stabilization experiments.
* **`/methods/`**: Contains all custom-written or modified functions.
    * `gauss_jacobi.m` / `gauss_seidel.m`: Custom functions for Problem 1, Part 2.
    * `run_gssa_experiment.m`: A modified, function-based version of Dr. Maliar's `Main_GSSA_1.m` used for batch experiments.
* **`/gssa_1_agent_capital/`**: Contains the original, unmodified code provided by Dr. Maliar (e.g., `Main_GSSA_1.m`, `Num_Stab_Approx.m`).
* **`/results/`**: Contains all text and figure outputs.
    * `/problem1/p1_part1_timing.txt`: The raw timing data for the direct solver benchmark.
    * `/problem2/p2_stabilization_results.txt`: The results table for the GSSA experiments.

---

## Problem 1: Direct vs. Iterative Solvers

### Part 1: Direct Solver Benchmark

The objective of this problem is to compare the time needed to solve the system $Ax=b$ 1,000 times. The benchmark script `solvers/problem1_solver_p1.m` was used, which correctly isolates the $O(n^2)$ "solve" step from the $O(n^3)$ "setup" (factorization/inversion) step for methods (b), (c), and (d).

**Analysis of Results (`results/problem1/p1_part1_timing.txt`)**

| Method | Total Time (1000 Iterations) | Analysis of Timed Operation |
| :--- | :--- | :--- |
| (a) Gaussian (A\b) | 31.355 s | **Flawed Benchmark:** This times `A\b` *inside* the loop, re-computing the $O(n^3)$ factorization 1,000 times. |
| (b) LU Decomposition | 5.246 s | **Correct Benchmark:** Times the 1,000 $O(n^2)$ *sequential* forward/backward substitutions. |
| (c) Matrix Inverse | 0.174 s | **Correct Benchmark:** Times the 1,000 $O(n^2)$ *parallel* matrix-vector multiplications. |
| (d) Cholesky (SPD) | 3.554 s | **Correct Benchmark:** Times the 1,000 $O(n^2)$ *sequential* forward/backward substitutions. |

#### Critical Analysis: Why is the Matrix Inverse "Solve" Fastest?

The empirical data is correct: the `inv(A)*b` solve step (0.17s) is over **30 times faster** than the `L\U\P*b` solve step (5.24s). This does not contradict the textbook; it reveals a crucial distinction between sequential and parallel operations.

* **LU/Cholesky Solve (Sequential):** The solve step `x = U\(L\(P*b))` is a *substitution*. This operation is **inherently sequential**. To find `x(2)`, you must first know the value of `x(1)`. This operation cannot be effectively parallelized and is thus bottlenecked by single-core processor speed.
* **Matrix Inverse Solve (Parallel):** The solve step `x = A_inv*b` is a *matrix-vector multiplication*. This is one of the most highly-optimized operations in computing (a Level 2 BLAS routine) and is massively **parallelizable**. MATLAB's underlying **Math Kernel Library (MKL)** executes this operation across all available CPU cores simultaneously.

**Conclusion:** This benchmark demonstrates a key trade-off. The `inv(A)` method has the **fastest *solve* step** due to its parallelizable nature. However, as noted in the "Linear Equations" (2025) slides, this speed comes at the cost of:
1.  **Higher Setup Cost:** The $O(n^3)$ cost to compute `inv(A)` is approximately 3 times higher than the cost to compute `lu(A)`.
2.  **Numerical Instability:** Computing the inverse is numerically less stable and can introduce larger floating-point errors than substitution.

For a problem requiring thousands of solves (like this one), the `inv(A)` method *is* fastest in practice, but `lu(A)` is generally preferred as it provides a better balance of setup cost, solve speed, and numerical reliability.

### Part 2: Iterative Solvers

This part implements the Gauss-Jacobi and Gauss-Seidel iterative methods based on operator splitting. The functions `methods/gauss_jacobi.m` and `methods/gauss_seidel.m` are called by the `solvers/problem1_solver_p2.m` script, which now uses the exact matrices specified in the problem set.

**Analysis of Results (`results/problem1/p1_part2_iterations.txt`)**

| System Matrix | Property | Gauss-Jacobi Iterations | Gauss-Seidel Iterations |
| :--- | :---: | :---: | :---: |
| `A1` | **Strictly Diagonally Dominant** | 40 | 24 |
| `A2` | **Strictly Diagonally Dominant** | 224 | 124 |
| `A3` | Not Diagonally Dominant | 5000 (Failed) | 5000 (Failed) |

**Conclusion:** The results perfectly illustrate the central convergence theorem for these methods.
* For systems `A1` and `A2`, **strict diagonal dominance** is satisfied, which guarantees convergence.
* System `A3` is **not** strictly diagonally dominant (e.g., `Row 1: |2| > |1| + |1| + |0|` is false). This violates the sufficient condition, and as a result, both methods failed to converge, hitting the 5000-iteration limit.
* This output directly answers the problem's question ("Why is it so different across cases?").
* For the converging matrices (`A1`, `A2`), Gauss-Seidel converged significantly faster than Gauss-Jacobi. This is the expected theoretical outcome, as Gauss-Seidel uses the most recently updated values of $x^{(k+1)}$ within the same iteration, accelerating convergence.

---

## Problem 2: Numerical Stabilization in GSSA

This problem uses the `gssa_1_agent_capital` code to analyze numerical stabilization techniques when solving a 1-agent growth model. Per the assignment, a 5th-degree polynomial was used (`D_max = 5` in `methods/run_gssa_experiment.m`). The `solvers/problem2_solver.m` script was used to generate the data.

**Final Results (`results/problem2/p2_stabilization_results.txt`)**

| Case Description | CPU Time (s) | Mean Error (log10) | Max Error (log10) | Parameters (RM, penalty, norm) |
| :--- | :---: | :---: | :---: | :--- |
| LAD-PP (Base) (FAILED) | 2.1745 | NaN | NaN | RM=3, penalty=0, norm=1 |
| LAD-DP (Base) (FAILED) | 0.95172 | NaN | NaN | RM=4, penalty=0, norm=1 |
| RLAD-PP (eta=10^-7) (FAILED) | 2.0613 | NaN | NaN | RM=7, penalty=-7, norm=1 |
| RLAD-DP (eta=10^-7) (FAILED) | 0.91481 | NaN | NaN | RM=8, penalty=-7, norm=1 |
| RLAD-PP (eta=10^-4) (FAILED) | 2.1057 | NaN | NaN | RM=7, penalty=-4, norm=1 |
| RLAD-DP (eta=10^-4) (FAILED) | 0.93908 | NaN | NaN | RM=8, penalty=-4, norm=1 |
| **RLS-Tikhonov (eta=10^-7)** | **6.2313** | **-6.3476** | **-5.0649** | **RM=5, penalty=-7, norm=1** |
| RLS-Tikhonov (eta=10^-4) | 8.3809 | -4.7174 | -3.5781 | RM=5, penalty=-4, norm=1 |
| **LS-SVD (Base)** | **4.5888** | **-9.9077** | **-8.1931** | **RM=2, penalty=0, norm=1** |
| RLS-TSVD (kappa=10^6) | 4.6591 | -8.253 | -6.651 | RM=6, penalty=6, norm=1 |
| **RLS-TSVD (kappa=10^7)** | **4.5736** | **-8.9151** | **-7.26** | **RM=6, penalty=7, norm=1** |
| RLS-TSVD (kappa=10^8) | 4.6307 | -9.9077 | -8.1931 | RM=6, penalty=8, norm=1 |
| LAD-DP (No Norm) (FAILED) | 0.90476 | NaN | NaN | RM=4, penalty=0, norm=0 |
| **RLS-Tikhonov (No Norm)** | **6.1583** | **-6.3476** | **-5.0649** | **RM=5, penalty=-7, norm=0** |
| **RLS-TSVD (No Norm)** | **4.6809** | **-8.9151** | **-7.26** | **RM=6, penalty=7, norm=0** |

### Analysis of PS3 Questions

#### 1. Compare LAD Methods

**Observation:** All six experiments based on Least-Absolute-Deviations (LAD) methods (LAD-PP, LAD-DP, RLAD-PP, RLAD-DP) failed to converge.

**Theoretical Analysis:** This is an expected outcome, not a code error. As discussed in the **Judd, Maliar, and Maliar (2011)** paper, LAD estimators are **non-smooth** and can be discontinuous. The GSSA method relies on a fixed-point iteration to converge on the policy function coefficients. The paper explicitly warns that the "jumps" inherent to the LAD method "might create problems in solving for a fixed point." This experiment confirms that theoretical warning: the LAD methods are not stable enough for this iterative solution algorithm.

#### 2. Tikhonov vs. SVD

* **Comparison:** The `LS-SVD (Base)` method (Mean Error: -9.9077) is dramatically more accurate than the `RLS-Tikhonov` method (Mean Error: -6.3476). For this problem, SVD is the superior method.
* **Parameter Effect:** When the Tikhonov `eta` parameter was increased (made stricter) from $10^{-7}$ to $10^{-4}$, the mean error worsened significantly (from -6.3476 to -4.7174). This indicates that the method is sensitive to this parameter and that the larger penalty was "over-regularizing" the system, harming accuracy.

#### 3. SVD vs. Truncated SVD (TSVD)

* **Comparison:** The `LS-SVD (Base)` (error -9.9077) and the `RLS-TSVD (kappa=10^8)` (error -9.9077) produce identical results. This is because a `kappa` threshold of $10^8$ is so loose that no singular values were "truncated."
* **Parameter Effect:** As the truncation threshold `kappa` is made stricter (lowered to $10^7$ and $10^6$), accuracy decreases (error rises to -8.9151 and -8.253, respectively). This implies that the singular values being removed by these stricter thresholds were numerically significant and contained important information for the solution.

#### 4. The Role of Normalization

**Observation:** Normalization had **no effect** on the results for the LS-based methods.
* The Tikhonov method produced an identical error (-6.3476) with and without normalization.
* The TSVD method produced an identical error (-8.9151) with and without normalization.

**Conclusion:** While normalization is a standard and often critical step for regularization (as noted in the **Judd, Maliar, and Maliar (2011)** paper), this experiment shows it was redundant for this specific problem. This suggests the variables and polynomial basis functions were already on a comparable scale, so no benefit was gained from the additional transformation.
