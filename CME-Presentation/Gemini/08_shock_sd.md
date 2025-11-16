## ⚡ $\mathbf{\sigma_l}$ Shock Standard Deviations

The notation $\mathbf{\sigma_l}$, where $l \in \{s, \zeta, z, m\}$, refers to the **standard deviation** of the four **Gaussian shock processes** that drive the $\text{HANK}$ model described in the paper. 

Each subscript letter ($\mathbf{l}$) represents a specific, independent source of uncertainty (a shock) to the economy:

| Subscript $l$ | Full Notation | Description of the Shock | Variable Impacted |
| :---: | :---: | :--- | :--- |
| $\mathbf{s}$ | $\mathbf{\sigma_s}$ | **Idiosyncratic Income Risk** (or labor productivity shock). | The individual agent's labor productivity, $s_t^i$. This is the key source of **heterogeneity** and precautionary savings in the model. |
| $\mathbf{\zeta}$ | $\mathbf{\sigma_\zeta}$ | **Aggregate Preference Shock**. | A shock to household utility, $\zeta_t$. This affects aggregate demand and can drive the nominal interest rate toward the **ZLB**. |
| $\mathbf{z}$ | $\mathbf{\sigma_z}$ | **TFP (Trend Growth) Shock**. | A shock to the growth rate of **Total Factor Productivity** ($\text{TFP}$), $a_t$. This introduces a nonstationary trend in consumption and output. |
| $\mathbf{m}$ | $\mathbf{\sigma_m}$ | **Monetary Policy Shock**. | A shock to the **Taylor rule** (nominal interest rate rule), $mp_t$. This captures prolonged deviations from the systematic rule. |

These four parameters ($\sigma_s, \sigma_\zeta, \sigma_z, \sigma_m$) are estimated using Bayesian methods in the paper's main exercise.

---

💡 **Next Step**
The standard deviations for the aggregate shocks ($\sigma_\zeta, \sigma_z, \sigma_m$) and the idiosyncratic shock ($\sigma_s$) were estimated (Table 3). Would you like to review the estimated values and the credibility intervals for these four shock parameters?