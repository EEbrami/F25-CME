# Project Overview: Problem Set 1

This problem set is designed to build and test proficiency in **function approximation**, a foundational technique in computational economics. It is divided into two distinct but related problems.

## Problem 1: Interpolation Methods Benchmark

The first problem focuses on the practical application and comparison of various two-dimensional function approximation methods. We interpolate a known benchmark, the **Franke function**, using three different approaches:
1.  **Chebyshev Polynomials:** This method is noted for its stability and efficiency, avoiding the oscillatory behavior common with ordinary polynomials.
2.  **Ordinary (Monomial) Polynomials:** This serves as a comparison, highlighting potential instability and Runge's phenomenon.
3.  **Toolbox Routines (Splines and Scattered Data):** This part emphasizes practical efficiency by using pre-built, optimized functions for spline and scattered data interpolation.

## Problem 2: Application to a Dynamic Economic Model

The second problem applies these approximation techniques to solve the standard **neoclassical stochastic growth model**. The core task is to replace the existing, less robust interpolation method (ordinary polynomials) within the provided Euler equation solution algorithm with the superior methods analyzed in Problem 1. The objective is to assess how the choice of approximation method impacts the accuracy of the economic model's solution, measured by the magnitude of the Euler equation residuals.