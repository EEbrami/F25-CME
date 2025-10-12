## Computational Methods for Economists (Econ-81360) - Problem Set 2 Analysis

### General Submission Mandate

Problem Set 2 shifts the focus from **function approximation** (HW1) to **numerical integration** methods. The submission systematically requires: (1) **Implementing** deterministic and stochastic integration rules, (2) **Benchmarking** their accuracy and convergence properties across different dimensions, and (3) **Applying** these various methods to the core task in computational economics: numerically solving the **conditional expectation** in a dynamic stochastic model.

| Problem | Domain | Methods Implemented/Compared | Key Analytical Goal |
| :--- | :--- | :--- | :--- |
| **1: 1D Integral** | $[0, 2\pi]$ | Gauss-Chebyshev, Trapezoid Rule, Simpson's Rule, Analytical Integration, MATLAB built-ins (`quad`, `integral`, `trapz`). | Compare accuracy and convergence rate (dependency on number of nodes $n$). |
| **2: 2D Integral** | $[0, 1] \times [0, 1]$ | 2D Gauss-Chebyshev (Product Rule), 2D Trapezoid Rule (Product Rule), Monte Carlo Integration, MATLAB built-in (`dblquad`). | Examine the computational cost and accuracy when moving to higher dimensions. |
| **3: Economic Model** | Stochastic (Normal) | Gauss-Hermite Quadrature (original), MATLAB's `quad`, Monte Carlo (using 10 nodes), Monomial Rules. | Evaluate the impact of expectation-approximating routines on the **accuracy of the economic policy function solution**. |

***

### Direct Relationship to Integration Lecture Slides

The structure of HW2 is a direct, practical application and critique of the methods presented in the "Deterministic Integration" slides.

1.  **Foundational Methods (Problem 1):**
    * The problem requires coding and comparing the **Trapezoid Rule** and **Simpson's Rule**, both of which are classified explicitly as **Newton-Cotes Formulas**. Students must analyze how accuracy depends on the number of subperiods, which directly relates to the concept of **truncation error** and the fact that Newton-Cotes methods exhibit polynomial convergence.
    * The implementation of the **Gauss-Chebyshev integration** rule directly applies the theory of **Gaussian formulas**, which are superior because they achieve exact integration for polynomials of degree $2n-1$ (for $n$ points) by optimizing the node and weight locations. This comparison highlights the practical trade-off between arbitrary vs. specialized node placement.

2.  **The Curse of Dimensionality (Problem 2 & 3):**
    * Problems 2 and 3 pivot to the central challenge in computational economics: **Multidimensional Integration**, specifically noting that the problem is ubiquitous due to multiple assets or error terms.
    * The request for 2D Gauss-Chebyshev and Trapezoid rules (Problem 2) represents the simplest extension using **Product Rules**, an approach critiqued in the lecture for suffering from the **curse of dimensionality** ($m^d$ functional evaluations).
    * **Problem 3** connects the mathematical tools directly to the NSGM, where the integral is the **Expected utility** term $E_{t}[\dots]$ in the Euler equation.
        * The problem's use of a shock $\epsilon_t \sim N(0, \sigma^2)$ mandates the use of **Gauss-Hermite quadrature**, a specific Gaussian formula tailored for integrals with the **Normal density weight function** $e^{-x^2}$ over the domain $[-\infty, \infty]$.
        * It then tests two major **forward-thinking** alternatives to the curse of dimensionality: **Monte Carlo Integration** and **Monomial Formulas**. The slides note that Monomial formulas are a **cheap, non-product alternative** that grows linearly in $d$, positioning them as a critical tool for large-scale economic models.

In essence, HW2 is a rigorous exercise forcing the confrontation between classic quadrature techniques and the practical demands of solving high-dimensional, stochastic economic models.
