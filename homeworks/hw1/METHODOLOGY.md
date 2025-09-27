# Methodology: Problem Set 1 Solution Steps

This document outlines the computational steps taken to solve Problem Set 1.

### 1. Environment Setup
- Define the Franke function as a MATLAB anonymous function.
- Create a fine grid (201x201) for the domain `[0,1]x[0,1]` to serve as the benchmark for plotting and error calculation.

### 2. Problem 1a: Chebyshev Approximation
- **Grid Generation:** Generate 2D Chebyshev nodes (extrema) using a tensor product.
- **Visualization:** Plot the true Franke function surface and overlay the Chebyshev interpolation nodes as a scatter plot.
- **Coefficient Calculation:** Implement the 2D Chebyshev regression formula to compute the coefficient matrix.
- **Interpolation & Analysis:** Evaluate the interpolant on the fine grid, plot the approximated surface and the error surface, and report accuracy metrics.

### 3. Problem 1b: Ordinary Polynomial Approximation
- **Grid Generation:** Create a 2D uniform grid. Visualize these nodes on the Franke function surface.
- **Coefficient Calculation:** Construct the Vandermonde-like matrix `X` and solve the linear system `X*c = y` for the coefficients `c`.
- **Interpolation & Analysis:** Evaluate the polynomial on the fine grid, visualize the result and error, and report accuracy.

### 4. Problem 1c: Toolbox Implementation
- **MATLAB Splines:** Use the `griddedInterpolant` function with the `'cubic'` method on a uniform grid.
- **MATLAB Scattered Data:** Generate random points in the domain and use `scatteredInterpolant`.
- **Fackler Toolbox (Optional):** Use `funfitxy` to compute Chebyshev coefficients and `funeval` to evaluate the interpolant, demonstrating a streamlined workflow.

### 5. Problem 2: Economic Model Solution
- **Baseline Setup:** Run the provided baseline code (Method 7) to replicate the original results.
- **Method Replacement:** Sequentially replace the ordinary polynomial interpolation step within the model's solution loop with a Chebyshev regression and then a spline-based approach.
- **Accuracy Comparison:** For each method, run the full algorithm to convergence and compute the Euler equation residuals on a long stochastic simulation to compare accuracy.