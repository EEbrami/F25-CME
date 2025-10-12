## Analysis of Homework 1: Function Approximation in Dynamic Economics

The provided materials document a computational economics problem set (ECON-81360) focused on **function approximation** techniques. The context is established by the lecture slides on *Approximation and Interpolation of Functions*, which motivates these methods as essential tools for solving dynamic economic models, specifically for approximating the unobservable **policy functions** (e.g., the capital accumulation function $k' = K(k, \theta)$) within a stochastic setting. The primary goal is to compare the performance of different basis functions and node choices.

### Problem 1: Function Approximation Benchmarking (The Franke Function)

The first problem benchmarks four approximation methods against the known, smooth **Franke function** in two dimensions. This setup isolates the approximation error from the economic model's complexity.

| Method | Nodes/Grid | Basis Type | Condition Number | Max Absolute Error |
| :--- | :--- | :--- | :--- | :--- |
| **Chebyshev** ($N=15$) | **Gauss-Lobatto** (Chebyshev Extrema) | Orthogonal Polynomial (Tensor Product) | N/A (Orthogonality Used) | **$3.8670 \times 10^{-3}$** |
| **Ordinary Polynomial** ($N=8$) | Uniform Grid | Monomials (Tensor Product) | $4.04 \times 10^{12}$ | $1.7380 \times 10^{-1}$ |
| **Cubic Spline** ($9 \times 9$) | Uniform Grid | Piecewise Polynomial | N/A (Toolbox Impl.) | $5.4613 \times 10^{-2}$ |
| **Scattered Interpolant** | Random Points | Natural Neighbors | N/A (Toolbox Impl.) | $1.6727 \times 10^{-1}$ |

#### Critical Interpretation of Problem 1 Results

The results unequivocally validate a core theoretical tenet of numerical methods: the **superiority of orthogonal polynomial bases on specialized grids** for smooth functions.

* The **Chebyshev** method, leveraging the **discrete orthogonality** property on the carefully chosen **Gauss-Lobatto nodes**, yielded an error that was at least one order of magnitude smaller than all other methods. This performance aligns with the Chebyshev theory's promise to achieve near **minmax ($L^{\infty}$) approximation**.
* The **Ordinary Polynomial** approach suffered catastrophic failure, exhibiting a **condition number of $4.04 \times 10^{12}$**, leading to the highest observed error. This demonstrates the numerical instability of the **Vandermonde matrix** with a uniform node choice and high polynomial degree, confirming the warning articulated in the lecture notes.
* The **Cubic Spline**, while numerically stable, performed significantly worse than Chebyshev for this smooth function, achieving only $5.46 \times 10^{-2}$ versus $3.87 \times 10^{-3}$. Splines are fundamentally designed to enforce local smoothness at the knots; for functions without kinks or non-smooth derivatives, their piecewise nature limits their global accuracy relative to a high-order global polynomial method like Chebyshev.

***

### Problem 2: Solving the Stochastic Growth Model (NSGM)

The second problem applies the numerical methods to a practical problem in computational economics: solving the two-dimensional **Neoclassical Stochastic Growth Model** using a **Projection (Euler Equation) Method**. The objective is to approximate the capital policy function $K(k, z)$ and find the coefficients that minimize the error in the Euler equation.

| Approximation Method | $\log_{10}$ Max Euler Error | $\log_{10}$ Mean Euler Error |
| :--- | :--- | :--- |
| **Chebyshev** (Degree 5) | $-1.8653$ | $-2.7060$ |
| **Cubic Spline** (20x10 Grid) | **$-2.0246$** | **$-2.7879$** |

#### Critical Interpretation of Problem 2 Results

The observed outcome in Problem 2 presents a notable counterpoint to the findings of Problem 1:

* **Result:** The **Cubic Spline** approximation achieved the smaller overall maximum Euler error ($\log_{10}$ Max Error of **$-2.0246$** vs. $-1.8653$), leading the final log to conclude it "performed best".
* **Forward-Thinking Analysis:** This reversal, where the less accurate global method from Problem 1 (Spline) outperforms the globally optimal method (Chebyshev) in the economic model, suggests that the policy function being approximated exhibits **local complexities** that are better captured by a piecewise approach.
    * In the context of the EEM, the goal is not general function approximation, but satisfying a first-order condition (Euler equation) that involves future-period variables and expectations. **Splines** are known to be advantageous in economic models that inherently introduce **non-smoothness** (e.g., occasionally binding constraints, non-concave value functions), as highlighted in the notes.
    * A standard NSGM policy function is smooth, yet numerical imprecision, integration quadrature errors, or local flattening/steepening around the steady state can introduce numerical "kinks" that the piecewise continuity of a spline handles more gracefully than a low-degree global Chebyshev polynomial that might overshoot or oscillate locally. The performance in this dynamic model suggests the **practical robustness of splines** when local accuracy is paramount, even if Chebyshev remains the theoretical champion for non-constrained, perfectly smooth global approximation.

In substance, the homework successfully demonstrates that the choice of approximation scheme must be guided not merely by theoretical convergence properties, but by the **numerical stability and local regularity** demanded by the specific economic problem. The failure of the poorly parameterized Ordinary Polynomial and the superior practical performance of the Spline in the NSGM provide critical instructional insights beyond the simple validation of Chebyshev theory.
