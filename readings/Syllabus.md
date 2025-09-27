# Computational Methods for Economists  
**Econ-81360**  
**Fall 2025**  
The City University of New York  
The Graduate Center  

---

**Instructor:** Lilia Maliar, Office 5313.04, lmaliar@gc.cuny.edu  
**Lectures:** Tuesday, 2:00 pm - 4:00 pm  
**Office hours:** Tuesday, 4:00 pm – 5:00 pm  

---

## Course Description

### Overview
This course studies computational approaches for solving dynamic economic models. The objectives of the course are threefold:

1. It provides background in numerical analysis (approximation, integration, optimization, error analysis), and describes local and global numerical methods (perturbation, Smolyak, endogenous grid, stochastic simulation, cluster grid methods).
2. It shows applications from recent economic literature representing challenges to computational methods (new Keynesian models with a zero lower bound, default risk models, Krusell-Smith models, international trade models, overlapping-generations models, nonstationary growth models, dynamic games).
3. It surveys recent developments in software and hardware (Python, Julia, GPUs, parallel computing, supercomputers), as well as machine learning techniques.

### Learning Goals and Outcomes

1. Demonstrate a strong understanding of numerical techniques for solving dynamic economic problems.
2. Develop strong programming skills.
3. Build solution algorithms suitable for given economic problems.
4. Identify new research questions in light of existing theories.
5. Design theoretical models to answer relevant research questions in the field.
6. Develop writing skills consistent with professional publications.

### Assessment

- **Individual problem sets (40%)**: Apply lecture material in practice. Relates to learning goals 1, 2, and 3.
- **Group project (20%)**: Each group presents an existing computational technique. Relates to learning goals 1 and 2.
- **Final individual project (40%)**: A short research paper using computational methods to address economic questions, presented in class. Relates to learning goals 4, 5, and 6.

| ACTIVITIES               | PERCENTAGES |
|--------------------------|-------------|
| Problem Sets             | 40%         |
| Group Project            | 20%         |
| Final Individual Project | 40%         |

### Prerequisites

Familiarity with basic dynamic optimization theory and undergraduate mathematics is recommended. Students who have passed the first-year core econometrics, microeconomics, and macroeconomics courses will be able to follow the course.

### Computing Languages

Students should know some computational language. MATLAB and Python are suggested. MATLAB is easy to learn and popular in economic literature; Python is open-source and widely used for machine learning. MATLAB will be used for most topics, Python for deep learning. Other languages (C, C++, Fortran, Julia) may also be used.

---

## Course Schedule

**Note:** Content may change depending on student progress and feedback.

### INTRODUCTION
1. Introduction to the course: richer and more complex economic models, analytical solutions, challenges of economic dynamics, model estimation, methodology, objectives, and outline.
2. Key ingredients of global solution methods: neoclassical growth model, Bellman equation, Euler equation, and an example of a global solution method.

### PART I. BACKGROUND IN NUMERICAL ANALYSIS WITH AN APPLICATION TO THE ONE-AGENT GROWTH MODEL

1. Approximation and interpolation: Lagrange and Hermite interpolation, orthogonal polynomials, Chebyshev regression, splines, complete polynomials.
2. Deterministic and Monte Carlo integration: monomial and quasi-Monte Carlo integration, Newton-Cotes, Monte Carlo and Gaussian integration, nonproduct deterministic methods, random and quasi-random sequences.
3. Linear equations and numerically stable methods: Cholesky decomposition, condition numbers, Gauss-Jacobi and Gauss-Seidel, SVD, LAD, regularization, principal component analysis.
4. Nonlinear equations and solvers: search methods, bisection, Newton method, BFGS and DFP updates, continuation and homotopy methods, fitting methods.
5. Optimization: iterative schemes, ECM and EGM methods.
6. Projection methods: for differential, integral, and functional equations.
7. Convergence rates and accuracy: residuals, errors, and economically meaningful accuracy measures.

### PART II. NUMERICAL ANALYSIS FOR HIGH DIMENSIONAL APPLICATIONS (WITH AN EXAMPLE OF A MULTI-AGENT GROWTH MODEL)

1. Smolyak technique: Smolyak grid, Smolyak polynomials, Lagrange interpolation.
2. Epsilon-distinguishable set and cluster grid techniques; ergodic-set methods.
3. Precomputation of intertemporal-choice functions and integrals; multivariate expectations.

### PART III. CHALLENGING ECONOMIC APPLICATIONS IN THE LITERATURE

1. Nonstationary and unbalanced growth models: time-varying parameters, seasonal changes, volatility, unbalanced growth.
2. New Keynesian models: stylized models, Taylor rule and zero lower bound, central banking models.
3. Models with many state variables: Aiyagari model, Krusell-Smith model, approximate aggregation.
4. Deep learning for solving dynamic models: consumption-saving, Krusell-Smith, discrete-choice models, HANK models.
5. Models with default risk: application of envelope condition method to sovereign default.
6. Dynamic games: Markov perfect equilibria, stochastic games, time inconsistency.

---

## Main Texts

- Kenneth L. Judd (1998). *Numerical Methods in Economics*. The MIT Press.
- Lilia Maliar and Serguei Maliar (2014). “Numerical Methods for Large Scale Dynamic Economic Models” in: Schmedders, K. and K. Judd (Eds.), *Handbook of Computational Economics*, Volume 3, Chapter 7, Elsevier Science.

## Supplementary Texts

- Jerome Adda & Russell W. Cooper (2003). *Dynamic Economics: Quantitative Methods and Applications*. MIT Press.
- Paolo Brandimarte (2006). *Numerical Methods in Finance and Economics*. Wiley.
- Burkhard Heer & Alfred Maussner (2004). *Dynamic General Equilibrium Modeling*. Springer.
- Kenneth L. Judd, Lilia Maliar & Serguei Maliar (2017). *Ergodic set methods for solving dynamic economic models*. MIT Press, in progress.
- Jianjun Miao (2014). *Economic Dynamics in Discrete Time*. MIT Press.
- Mario J. Miranda & Paul L. Fackler (2002). *Applied Computational Economics and Finance*. MIT Press.
- G. C. Lim & Paul D. McNelis (2008). *Computational Macroeconomics for the Open Economy*. MIT Press.
- David A. Kendrick, P. Ruben Mercado, Hans M. Amman (2006). *Computational Economics*. Princeton University Press.
- John Stachurski (2009). *Economic Dynamics. Theory and Computation*. MIT Press.

---

## Selected Readings

### Part I Readings

- Arellano, C., L. Maliar, S. Maliar & V. Tsyrennikov (2016). Envelope condition method with an application to default models. *Journal of Economic Dynamics and Control* 69, 436-459.
- Barillas, F. & J. Fernández-Villaverde (2007). A generalization of the endogenous grid method. *Journal of Economic Dynamics and Control* 31, 2698-2712.
- Carroll, K. (2005). The method of endogenous grid points for solving dynamic stochastic optimal problems. *Economic Letters* 91, 312-320.
- Judd, K., Maliar, L. & S. Maliar (2011). Numerically stable and accurate stochastic simulation approaches for solving dynamic economic models. *Quantitative Economics* 22, 173-210.
- Judd, K. (1998). *Numerical Methods in Economics*. MIT Press.
- Judd, K., L. Maliar & S. Maliar (2017). Lower bounds on approximation errors to numerical solutions of dynamic economic models. *Econometrica* 85(3), 991-1020.
- Maliar L. & S. Maliar (2013). Envelope condition method versus endogenous grid method for solving dynamic programming problems. *Economic Letters* 120, 262-266.
- Maliar, L. & S. Maliar (2014). “Numerical Methods for Large Scale Dynamic Economic Models” in: Schmedders, K. and K. Judd (Eds.), *Handbook of Computational Economics*, Volume 3, Chapter 7, Elsevier Science.
- Marcet, A. & G. Lorenzoni (1999). The parameterized expectation approach: some practical issues. In: R. Marimon and A. Scott (Eds.) *Computational Methods for Study of Dynamic Economies*. Oxford University Press.
- Rust, J. (1997). Using randomization to break the curse of dimensionality. *Econometrica* 65, 487-516.
- Stroud A. (1971). *Approximate Integration of Multiple Integrals*. Prentice Hall.

### Part II, Section 1 Readings

- Brumm, J., Scheidegger, S. (2017). Using adaptive sparse grids to solve high-dimensional dynamic models. *Econometrica*, forthcoming.
- Judd, K, Maliar, L., Maliar, S. & R. Valero (2014). Smolyak method for solving dynamic economic models: Lagrange interpolation, anisotropic grid and adaptive domain. *Journal of Economic Dynamics and Control* 44, 92-103.
- Krueger, D. & F. Kubler (2004). Computing equilibrium in OLG models with production. *Journal of Economic Dynamics and Control* 28, 1411-1436.
- Malin, B., Krueger, D., & F. Kubler (2011). Solving the multi-country real business cycle model using a Smolyak-collocation method. *Journal of Economic Dynamics and Control* 35(2), 229-239.

### Part II, Section 2 Readings

- Eldén, L. (2007). *Matrix Methods in Data Mining and Pattern Recognition*. SIAM.
- Hastie, T., R. Tibshirani & J. Friedman (2009). *The Elements of Statistical Learning*. Springer.
- Judd, K., Maliar, L. & S. Maliar (2010). A cluster-grid projection algorithm: solving problems with high dimensionality, NBER 15965.
- Maliar, S., Maliar, L. & K. Judd (2011). Solving the multi-country real business cycle model using ergodic set methods. *Journal of Economic Dynamics and Control* 35(2), 207-228.
- Maliar, L. & S. Maliar (2015). Merging simulation and projection approaches to solve high-dimensional problems with an application to a new Keynesian model. *Quantitative Economics* 6, 1-47.
- Niederreiter, H. (1992). *Random Number Generation and Quasi-Monte Carlo Methods*. SIAM.
- Temlyakov, V. (2011). *Greedy approximation*. Cambridge University Press.

### Part II, Section 3 Readings

- Judd, K., L. Malia, S. Malia & I. Tsener (2017). How to solve dynamic stochastic models computing expectations just once. *Quantitative Economics* 8(3), 851-893.
- Maliar, L. & S. Maliar (2005). Solving nonlinear stochastic growth models: an algorithm computing value function by simulations. *Economics Letters* 87, 135-140.

### Part III, Section 1 Readings

- Adjemian, S. & M. Juillard (2013). Stochastic extended path approach. Manuscript.
- Fair, R. & J. Taylor (1983). Solution and maximum likelihood estimation of dynamic nonlinear rational expectations models. *Econometrica* 51, 1169-1185.
- Schmitt-Grohé, S. & M. Uribe (2012). What's news in business cycles? *Econometrica* 80, 2733-2764.
- Lepetyuk, V., L. Maliar, S. Maliar & J. B. Taylor (2019). Extended function path perturbation for nonstationary and unbalanced growth models. Manuscript.
- Maliar, L., S. Maliar, J. Taylor & I. Tsener (2015). A tractable framework for analyzing a class of nonstationary Markov Models, NBER 21155.
- Maliar, L., Maliar S. & I. Tsener (2018). Capital-skill complementarity: twenty years after. Manuscript.

### Part III, Section 2 Readings

- Aruoba, S., Cuba-Borda, P. & F. Schorfheide (2017). Macroeconomic dynamics near the ZLB: a tale of two countries. *Review of Economic Studies* 85(1), 87-118.
- Christiano, L., M. Eichenbaum & S. Rebelo (2011). When is the government spending multiplier large? *Journal of Political Economy* 119(1), 78-121.
- Fernández-Villaverde, J., Gordon, G., Guerrón-Quintana, P. & J. Rubio-Ramírez (2012, 2015). Nonlinear adventures at the zero lower bound. NBER 18058 and *Journal of Economic Dynamics and Control*, 182-204.
- Lepetuyk, V., L. Maliar & S. Maliar (2019). When the U.S. catches a cold, Canada sneezes: a lower-bound tale told by deep learning. CEPR working paper DP 14025.
- Maliar, L. & S. Maliar (2014). Merging simulation and projection approaches to solve high-dimensional problems with an application to a new Keynesian model. *Quantitative Economics* 6, 1-47.
- Mertens, K. & M. Ravn (2011). Credit channels in a liquidity trap. CEPR discussion paper 8322.
- Smets, F. & R. Wouters (2007). Shocks and frictions in US business cycles: a Bayesian DSGE approach. *American Economic Review* 97(3), 586-606.

### Part III, Section 3 Readings

- Ahn, S., G. Kaplan, B. Moll, T. Winberry & C. Wolf (2017). When Inequality matters for macro and macro matters for inequality. Manuscript.
- Auclert, A., B. Bardóczy, M. Rognlie & L. Straub (2019). Using the sequence-space Jacobian to solve and estimate heterogeneous-agent models. Manuscript.
- Bayer, C. & R. Lutticke (2018). Solving heterogeneous agent models with aggregate uncertainty and many idiosyncratic states in discrete time by perturbation methods. Manuscript.
- Boppart, T., Krusell, P. & K. Mitman (2018). Exploiting MIT shocks in heterogeneous-agent: The impulse-response as a numerical derivative. *Journal of Economic Dynamics and Control* 89, 68-92.
- Den Haan (2010a, b). Assessing the accuracy of the aggregate law of motion and comparison of solutions to the incomplete markets model with aggregate uncertainty. *Journal of Economic Dynamics and Control* 34, 4-99.
- Kim, S., R. Kollmann & J. Kim (2010). Solving the incomplete markets model with aggregate uncertainty using a perturbation method. *Journal of Economic Dynamics and Control* 34, 50-58.
- Krusell, P. & A. Smith (1998). Income and wealth heterogeneity in the macroeconomy. *Journal of Political Economy* 106:5, 868-96.
- Ljungqvist, L. & T. Sargent (2000). *Recursive Macroeconomic Theory*. MIT Press.
- Maliar, L., Maliar, S. & F. Valli (2009). Solving the incomplete markets model with aggregate uncertainty using the Krusell-Smith algorithm. *Journal of Economic Dynamics and Control* 34, 42-49.
- Maliar, L., S. Maliar & P. Winant (2019). Will AI replace computational economists any time soon? CEPR working paper DP 14024.
- Reiter, M. (2010, 2018). Solving the incomplete markets economy with aggregate uncertainty by backward induction and HetSol toolkit. *Journal of Economic Dynamics and Control* 34, 28-35; Slides.
- Winberry, T. (2016). A toolbox for solving and estimating heterogeneous agent macro models. Manuscript.

### Part III, Section 4 Readings

- Azinovic, M., G. Luca & S. Scheidegger (2019). Deep equilibrium nets. *International Economic Review* 63(4), 1471-1525.
- Duarte, V. (2018). Machine learning for continuous-time economics. [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3012602)
- Fernández-Villaverde, J., S. Hurtado & G. Nuño (2018). Financial frictions and the wealth distribution. Manuscript.
- Goodfellow, I., Y. Bengio & A. Courville (2016). *Deep Learning*. MIT Press.
- Hastie, T., R. Tibshirani & J. Friedman (2009). *The Elements of Statistical Learning*. Springer.
- Maliar, L., S. Maliar & P. Winant (2019). Deep learning for solving dynamic economic models. *Journal of Monetary Economics* 122.
- Maliar, L. & S. Maliar (2021). Deep learning classification: modeling discrete labor choice. *Journal of Economic Dynamics and Control* 135, 104295.
- Villa, A. & V. Valaitis (2024). Machine learning projection methods for macro-finance models. *Quantitative Macroeconomics* 15(1), 145-173.

### Part III, Section 5 Readings

- Aguiar, M. & G. Gopinath (2006). Defaultable debt, interest rates and the current account. *Journal of International Economics* 69(1), 64-83.
- Arellano, C. (2008). Default risk and income fluctuations in emerging economies. *American Economic Review* 98(3), 690-712.
- Arellano, C., Maliar, L., Maliar S. & V. Tsyrennikov (2014). Envelope condition method with an application to default risk models. Manuscript.
- Gordon, G. & S. Qui (2015). A divide and conquer algorithm for exploiting policy function monotonicity. Manuscript, Indiana University.
- Villemot, S. (2012). Accelerating the resolution of sovereign debt models using an endogenous grid method. Dynare working paper 17, CEPREMAP.

### Part III, Section 6 Readings

- Cao, D. & I. Werning (2018). Saving and dissaving with hyperbolic discounting. *Econometrica* 86(3), 805-857.
- Harris, C. & D. Laibson (2001). Dynamic choices of hyperbolic consumers. *Econometrica* 69(4), 935-959.
- Ferris, M., K.L. Judd & K. Schmedders (2011). Solving dynamic games with Newton’s method. Manuscript.
- Judd, K. (2004). Existence, uniqueness, and computational theory for time consistent equilibria: a hyperbolic discounting example. Manuscript.
- Krusell, P. & A. Smith (2003). Consumption-savings decisions with quasi-geometric discounting. *Econometrica* 71, 365-375.
- Krusell, P., B. Kuruscu & A. Smith (2002). Equilibrium welfare and government policy with quasi-geometric discounting. *Journal of Economic Theory* 105, 42-72.
- Laibson, D., A. Repetto & J. Tobacman (1998). Self-control and saving for retirement. *Brookings Papers on Economic Activity* 1, 91-172.
- Maliar, L. & S. Maliar (2005). Solving the neoclassical growth model with quasi-geometric discounting: A grid-based Euler-equation method. *Computational Economics* 26, 163-172.
- Maliar, L. & S. Maliar (2006). The neoclassical growth model with heterogeneous quasi-geometric consumers. *Journal of Money, Credit, and Banking* 38(3), 635-654.
- Maliar, L. & S. Maliar (2016). Ruling out multiplicity of smooth equilibria in dynamic games: a hyperbolic discounting example. *Dynamic Games and Applications* 6(2), 243-261.

---
