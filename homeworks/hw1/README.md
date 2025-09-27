# Homework 1: Function Approximation

**Course**: Econ-81360 - Computational Methods for Economists  
**Term**: Fall 2025  
**Points**: 100

## Overview

This problem set builds and tests proficiency in function approximation, a foundational technique in computational economics. It consists of two main problems:

1. **Function Approximation Benchmarking**: Compare various 2D interpolation methods using the Franke function
2. **Economic Model Application**: Apply these methods to solve a neoclassical stochastic growth model

## Files Structure

```
hw1/
├── README.md                           # This file
├── main_problem1_approximation.m       # Main script for Problem 1
├── main_problem2_growth_model.m        # Main script for Problem 2
├── functions/                          # Helper functions
│   ├── cheb_eval_2d.m                 # 2D Chebyshev evaluation
│   ├── poly_eval_2d.m                 # 2D polynomial evaluation
│   └── solve_growth_model_baseline.m  # Growth model solver
├── documentation/                      # Project documentation
│   ├── PROJECT_OVERVIEW.md           # Problem overview
│   └── METHODOLOGY.md                # Detailed methodology
└── results/                           # Output directory (auto-created)
    ├── figures_problem1/              # Problem 1 plots
    └── figures_problem2/              # Problem 2 plots
```

## How to Run

### Problem 1: Function Approximation Benchmarking

1. Navigate to the `hw1/` directory in MATLAB
2. Run the main script:
   ```matlab
   main_problem1_approximation
   ```
3. The script will:
   - Generate and display the benchmark Franke function
   - Implement and compare 4 approximation methods:
     - 2D Chebyshev polynomials
     - Ordinary (monomial) polynomials  
     - Cubic splines via `griddedInterpolant`
     - Scattered data interpolation
   - Save all figures to `results/figures_problem1/`
   - Display accuracy comparison table

### Problem 2: Economic Model Application

1. Run the framework script:
   ```matlab
   main_problem2_growth_model
   ```
2. This will show the implementation plan and placeholder results
3. To complete Problem 2, you need to:
   - Obtain the baseline EGM policy iteration code
   - Implement the enhanced Chebyshev and Spline versions
   - Compare accuracy via Euler equation residuals

## Methods Implemented

### Problem 1 Methods

| Method | Description | Strengths | Weaknesses |
|--------|-------------|-----------|------------|
| **Chebyshev** | 2D tensor product using extrema nodes | Excellent convergence, avoids Runge phenomenon | Complex implementation |
| **Ordinary Polynomial** | Monomial basis on uniform grid | Simple concept | Unstable for high degrees |
| **Cubic Splines** | Piecewise cubic via MATLAB toolbox | Good accuracy, stable | Limited to grid data |
| **Scattered Data** | Natural neighbor interpolation | Handles irregular data | Less accurate than grid methods |

### Problem 2 Framework

- **Algorithm 7**: Euler equation method parameterizing capital decision function
- **Baseline**: Ordinary polynomial regression on uniform grid
- **Enhanced**: Replace with Chebyshev or spline approximations
- **Accuracy**: Measured via Euler equation residuals

## Expected Results

### Problem 1
The Chebyshev method should achieve the highest accuracy (lowest maximum error), followed by splines, with ordinary polynomials showing instability at higher degrees.

### Problem 2  
Enhanced methods (Chebyshev/Splines) should produce smaller Euler equation residuals compared to the baseline ordinary polynomial method.

## Key Learning Objectives

1. **Understand function approximation trade-offs**: Accuracy vs. computational cost vs. implementation complexity
2. **Experience the curse of dimensionality**: See how approximation difficulty increases with dimensionality
3. **Apply methods to economic problems**: Bridge abstract numerical methods to real economic applications
4. **Compare methodological approaches**: Develop intuition for when to use which method

## Technical Notes

- All methods use the domain [0,1] × [0,1] for consistency
- Chebyshev nodes are mapped from canonical [-1,1] domain
- Error analysis uses a fine 201×201 evaluation grid
- Figures are automatically saved as PNG files
- Code includes comprehensive documentation and help blocks

## Submission Requirements

When complete, your submission should include:
- All MATLAB files with working implementations
- Generated figures showing method comparisons
- Written report discussing results and method trade-offs
- Completed Euler equation analysis for Problem 2

---

*For detailed methodology and mathematical background, see the files in the `documentation/` folder.*