# Presentation Flow and Transitions

This guide provides a narrative script for a presentation based on the sorted markdown files. It outlines the goal of each slide and suggests language for transitioning smoothly to the next topic.

---

### Slide 1: `01_general_understanding.md`

*   **Goal:** Start with the 30,000-foot view of the problem.
*   **Script:** "Today, we're tackling a major challenge in modern macroeconomics: how to actually *estimate* complex structural models, specifically Heterogeneous Agent New Keynesian or 'HANK' models. The core issue boils down to two massive computational bottlenecks: first, solving the model for what agents will do, and second, evaluating how well the model fits the data. Our goal is to find a way to do this efficiently."
*   **Transition:** "But before we dive into the technical weeds, it's crucial to understand *why* this is such a timely and important problem in the field."

---

### Slide 2: `02_literature_direction.md`

*   **Goal:** Place the work in its academic context.
*   **Script:** "For decades, macroeconomics had two separate conversations. One group focused on agent heterogeneity but had to simplify their models by linearizing them. The other group studied important nonlinearities, like the Zero Lower Bound, but had to assume everyone was the same—the 'representative agent.' This paper sits at the very frontier, attempting to fuse these two branches. We are among the first to not just *solve*, but fully *estimate* a HANK model that is also nonlinear."
*   **Transition:** "Now that we see the academic frontier we're pushing, let's get specific about the computational hurdles that have made this so difficult until now. This brings us back to the two bottlenecks."

---

### Slide 3: `03_inner_loops.md`

*   **Goal:** Detail the technical nature of the problem.
*   **Script:** "The bottlenecks exist within what we call 'nested computational loops.' For every single parameter set we want to test, we first have an 'inner loop'—the Solution Step—where we iterate endlessly to find the agents' policy functions. Then, we have a *second* 'inner loop'—the Evaluation Step—where we run thousands of Monte Carlo simulations to calculate the model's likelihood. Doing both of these for every guess is what makes the traditional approach computationally unmanageable."
*   **Transition:** "So, how do we escape this nested loop nightmare? The core innovation is to completely restructure the problem."

---

### Slide 4: `04_loop_transformation.md`

*   **Goal:** Present the conceptual solution.
*   **Script:** "The solution is to 'decouple' the process. Instead of solving and evaluating *inside* the main estimation loop, we do all the heavy lifting in a one-time, upfront training phase. We train neural networks to become surrogates for the model's solution and its likelihood function. The main estimation then becomes incredibly fast, as it just needs to query these pre-trained networks."
*   **Transition:** "This naturally leads to the question: what do these neural networks actually look like?"

---

### Slide 5: `05_NN_design.md`

*   **Goal:** Explain the specific NN architecture.
*   **Script:** "We use a suite of specialized networks. The 'NN-Solution' is a set of networks trained to approximate the policy functions. The key trick here is that we treat the model's structural parameters as *inputs* to the network. Then, the 'NN-Evaluation' network learns the mapping from any set of parameters directly to the final likelihood value, effectively replacing the entire particle filter."
*   **Transition:** "This all sounds great in theory, but does it actually work? To build confidence, we validate our method step-by-step."

---

### Slide 6: `06_model_progression_validation.md`

*   **Goal:** Show the proof of concept and build trust in the method.
*   **Script:** "We start with a simple, linear textbook model where we know the right answer. Our method nails it. Then, we add the first layer of complexity: the Zero Lower Bound. We show our method's results are 'remarkably similar' to other state-of-the-art solvers. Finally, we go to the full-blown Nonlinear HANK model and show that we can accurately recover the true parameters from simulated data. This progression proves our method is reliable."
*   **Transition:** "Now that we've established the method is sound, we can ask the big question: what do we learn about the economy from applying it?"

---

### Slide 7: `07_ZLB_importance.md`

*   **Goal:** Present the key economic finding and the payoff for the complex method.
*   **Script:** "The payoff is a crucial economic insight. The model shows that the interaction between household income risk and the ZLB is a massive amplifier of economic volatility, accounting for 22 percentage points of GDP volatility. Simpler models completely miss this channel. This justifies the complexity: only the full nonlinear HANK model can connect the micro data on inequality to the macro data on business cycles, and it gets the wealth Gini coefficient nearly right without even being trained on it."
*   **Transition:** "The engine of this whole story is risk. Let's briefly define the specific shocks that drive our model."

---

### Slide 8: `08_shock_sd.md`

*   **Goal:** Define the key parameters that drive the model's results.
*   **Script:** "The model is driven by four shocks, but the most important one is σ_s, the standard deviation of idiosyncratic income risk. This is the parameter that drives precautionary savings, pushes the economy closer to the ZLB, and is a key source of the amplification we just saw. Our method successfully identifies this crucial parameter using only macro data."
*   **Concluding Thought:** "This demonstrates the power of the NN approach: it not only solves a major computational problem but also allows us to uncover economically significant mechanisms that were previously beyond our reach."
