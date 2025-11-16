The **"Unmanageable" Nested Loop (Traditional Method)** is fundamentally characterized by the fact that the **Solution Step** (solving the model's policy functions) and the **Evaluation Step** (calculating likelihood) are repeated *inside* the high-level **Estimation Loop** for every single parameter draw.

In contrast, the **Decoupled NN Method** breaks this computationally expensive dependency by executing the **Solution** and **Evaluation** steps *once* during an **Upfront Training Phase** to create fast, pre-trained Neural Network surrogates. The final **Estimation Loop** then runs extremely fast, relying only on these pre-computed NN approximations of the likelihood function.

---

### 1. The "Unmanageable" Nested Loop (Traditional Method)

This is the traditional approach to estimating a structural model. It consists of a high-level **Estimation Loop** (the optimizer) that, at **every single step**, must call a **Nested Computational Loop** to solve and evaluate the model for a new parameter guess.

This nested structure is the source of the computational bottleneck.

**Algorithm**:
* Start **High-Level Estimation Loop** (e.g., **RWMH** algorithm, for $N$ draws).
    * Iteration 1: Generate a parameter guess, $\Theta_1$.
    * Start **Nested Computational Loop (The Bottleneck)**:
        * A. **Solution Step**: Run a traditional, **slow global solver** to find the model's policy functions $\psi(\cdot | \Theta_1)$. This must be done **from scratch**.
        * B. **Evaluation Step**: Use the policy functions $\psi(\cdot | \Theta_1)$ to run a traditional, **slow Monte Carlo filter** (e.g., particle filter) to get the likelihood value, $\mathcal{L}(\mathbb{Y}_{1:T} | \Theta_1)$.
    * End Nested Computational Loop.
    * Use the value $\mathcal{L}(\mathbb{Y}_{1:T} | \Theta_1)$ to decide whether to accept or reject the draw $\Theta_1$.
    * Iteration 2: Generate a new parameter guess, $\Theta_2$.
    * Start **Nested Computational Loop (The Bottleneck)**:
        * A. **Solution Step**: Run the slow solver **all over again** for $\Theta_2$.
        * B. **Evaluation Step**: Run the slow particle filter **all over again** for $\Theta_2$.
    * End Nested Computational Loop.
    * Use the value $\mathcal{L}(\mathbb{Y}_{1:T} | \Theta_2)$ to accept/reject $\Theta_2$.
    * ...Repeat $N$ times.

**Visualizing the Loop**:
Estimation Loop (Optimizer):
$\dots \rightarrow$ [Generate $\Theta_n$] $\rightarrow$ (START Nested Loop) $\rightarrow$ Solve Model for $\Theta_n$ $\rightarrow$ Evaluate Likelihood for $\Theta_n$ $\rightarrow$ (END Nested Loop) $\rightarrow$ [Accept/Reject $\Theta_n$] $\rightarrow$ [Generate $\Theta_{n+1}$] $\rightarrow \dots$

---

### 2. Breaking the Loop: The Decoupled NN Method

The paper's innovation is to **break this nested loop**. The **entire nested computational loop** (Steps A and B) is removed from the estimation process and replaced by NNs that are trained **once** upfront.

This **"decouples"** the process into two distinct, sequential phases.

#### Phase 1: Upfront Training (The New "Computational" Loop)

This is a new, one-time process where all the heavy computation is performed.

**Algorithm**:
1.  Train a suite of NNs (**"NN-Solution"**) to approximate the model's key components over the **entire** parameter space.
    * **NN-DSS** approximates the **Deterministic Steady State**.
    * **NN-PolicyFunctions** approximates the **"extended policy functions"** by treating parameters as inputs (**"pseudo-state variables"**).
    * *This replaces: The need to run the **Solution Step (A)** at every iteration of the estimation.*
2.  Use the **"NN-Solution"** to quickly generate a large training sample (e.g., 10,000s) of parameter-likelihood pairs: $\{(\Theta_i, \mathcal{L}_i), \dots\}$.
3.  Train a final NN (**"NN-Evaluation"**, or the **"NN particle filter"**) on this dataset to learn the mapping from any $\Theta$ to its likelihood $\mathcal{L}$.
    * *This replaces: The need to run the **Evaluation Step (B)** at every iteration of the estimation.*

#### Phase 2: Fast Estimation (The New "Estimation" Loop)

With the NNs now fully trained, the estimation loop becomes incredibly simple and fast.

**Algorithm**:
* Start **High-Level Estimation Loop** (e.g., **RWMH** algorithm, for $N$ draws).
    * Iteration 1: Generate a parameter guess, $\Theta_1$.
    * Start **FAST Evaluation**:
        * Feed $\Theta_1$ into the pre-trained **"NN-Evaluation"**.
        * The NN **instantaneously** returns the approximated likelihood $\mathcal{L}(\Theta_1)$.
    * End FAST Evaluation.
    * Use the value $\mathcal{L}(\Theta_1)$ to decide whether to accept or reject the draw $\Theta_1$.
    * Iteration 2: Generate a new parameter guess, $\Theta_2$.
    * Start **FAST Evaluation**:
        * Feed $\Theta_2$ into the pre-trained **"NN-Evaluation"**.
        * The NN **instantaneously** returns the approximated likelihood $\mathcal{L}(\Theta_2)$.
    * End FAST Evaluation.
    * ...Repeat $N$ times.

**Visualizing the Decoupled Loops**:
Phase 1 (One Time):
[Train NN-Solution] $\rightarrow$ [Train NN-Evaluation] $\rightarrow$ (Trained NNs are now "frozen")

Phase 2 (Estimation Loop):
$\dots \rightarrow$ [Generate $\Theta_n$] $\rightarrow$ (Call NN-Evaluation on $\Theta_n$) $\rightarrow$ [Get Likelihood $\mathcal{L}(\Theta_n)$] $\rightarrow$ [Accept/Reject $\Theta_n$] $\rightarrow$ [Generate $\Theta_{n+1}$] $\rightarrow \dots$