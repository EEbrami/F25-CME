---
Model Progression in Validation: A Proof of Concept
---

The paper's three-step validation is a "proof of concept" designed to build trust in the methodology. It starts with the simplest possible case and incrementally adds the two key complexities the paper aims to solve: **aggregate nonlinearity** and **agent heterogeneity**. 

---

### Model 1: Linearized Small-Scale DSGE (Linear RANK)

This is the baseline model, which is a standard, "off-the-shelf" three-equation New Keynesian ($\text{NK}$) model. It is fully linear and has a known analytical solution, making it a perfect benchmark.

#### Key Formulas

The model is defined by a system of log-linearized equations:
* **IS Curve (Euler Eq.)**:
    $$\hat{Y}_{t}=E_{t}\hat{Y}_{t+1}-\sigma^{-1}(\hat{R}_{t}-E_{t}\hat{\Pi}_{t+1})$$
* **New Keynesian Phillips Curve**:
    $$\hat{\Pi}_{t}=\kappa(\hat{Y}_{t}-\hat{Y}_{t}^{*})+\beta E_{t}\hat{\Pi}_{t+1}$$
* **Monetary Policy Rule (Linear)**:
    $$\hat{R}_{t}=\phi_{\Pi}\hat{\Pi}_{t}+\phi_{Y}\hat{X}_{t}$$

There is **no ZLB constraint** and **no agent heterogeneity** (it is a **Representative Agent New Keynesian**, or $\text{RANK}$, model).

#### Key Results

* **NN-Solution Accuracy**: The paper tests the **NN-Solution** by comparing its approximated policy functions ($\text{PFs}$) against the model's true analytical solution. The NN approximation and the analytical solution "almost perfectly coincide." The mean squared residual error converges to a very low value of $10^{-10}$.
* **NN-Evaluation Accuracy**: The paper tests the **NN-Evaluation** ($\text{NN}$ particle filter) by checking if it can accurately approximate the likelihood function. The NN (blue line) successfully "cuts through the cloud" of noisy likelihood points from a standard particle filter (orange dots) and "peaks closely" at the true parameter value.

---

### Model 2: Nonlinear RANK with ZLB

This model takes the $\text{RANK}$ framework from Model 1 and introduces the first layer of complexity.

#### What's Added: Aggregate Nonlinearity ($\text{ZLB}$)

The model is no longer linearized. The linear Taylor rule is replaced with a **recurrently binding Zero Lower Bound ($\text{ZLB}$) constraint**. The model still features a single "representative household."

#### Key Formulas

The household and firm equations are now written in their full nonlinear form. The critical change is the monetary policy rule, which introduces the $\max[\cdot]$ operator:
* **Monetary Policy Rule (Nonlinear with ZLB)**:
    $$R_{t}=\max[1,R(\frac{\Pi_{t}}{\Pi})^{\theta_{t1}}(\frac{Y_{t}}{Y})^{\theta_{Y}}]$$

#### Key Results

The test's goal is to compare the paper's full $\text{NN}$-based estimation ($\text{NN}$-Solution + $\text{NN}$-Evaluation) against a different "state-of-the-art" (non-NN) global solution method.

* **Posterior Accuracy**: The posterior distributions estimated by the $\text{NN}$ method (blue line) and the alternative, state-of-the-art method (red dashed line) are "**remarkably similar**."
* **Parameter Recovery**: The paper concludes that its $\text{NN}$ method performs **equally well as the benchmark**, with the posterior medians falling "close to the true value" and within the $90\%$ credible intervals.

---

### Model 3: Nonlinear HANK with ZLB

This model adds the second and final layer of complexity.

#### What's Added: Agent Heterogeneity & Idiosyncratic Risk

The model replaces the "representative agent" from Model 2 with a "**continuum of households**" (i.e., **Heterogeneous Agent New Keynesian**, or $\text{HANK}$). These households face:

* **Idiosyncratic Income Risk** ($s_t^i$).
* An **individual borrowing limit** ($B_t \ge \underline{B}$), which is a new **micro-level nonlinearity**.

The model retains the **aggregate $\text{ZLB}$ nonlinearity** from Model 2.

#### Key Formulas

The model is now fully defined, integrating both nonlinearities:

* **Aggregate Nonlinearity ($\text{ZLB}$)**: The $\text{ZLB}$ from Model 2 is still present:
    $$R_{t}=\max[1,R_{t}^{n}]$$
* **Heterogeneity (Individual Budget Constraint)**: The household budget constraint is now indexed by agent $i$ and includes the idiosyncratic productivity shock $s_t^i$:
    $$C_{t}^{i}+B_{t}^{i}=\tau_{t}(\frac{W_{t}}{A_{t}}\exp(s_{t}^{i})H_{t}^{i})^{1-\gamma_{\tau}}+\frac{R_{t-1}}{\Pi_{t}}B_{t-1}^{i}+\text{Div}_{t}^{i}$$
* **Micro Nonlinearity (Borrowing Constraint)**: Agents also face an individual borrowing limit:
    $$B_{t}\ge\underline{B}$$

#### Key Results

The test's goal is to see if the full $\text{NN}$ methodology can accurately recover the true parameters of this highly complex model using **simulated data**.

* **Parameter Recovery**: The paper reports that the method successfully estimates all 10 structural parameters.
* **Accuracy**: The posterior median for each parameter is "**very close to the true value**," and the true value is "**mostly contained in the 90% credible interval**." This is a crucial result, as it proves the method can successfully estimate parameters that affect the model's **Deterministic Steady State** (e.g., $\sigma_s$ and $\underline{B}$), which is a key bottleneck for older methods.