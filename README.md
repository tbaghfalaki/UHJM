## Unified Zero-Inflated Hurdle Joint Modeling
Perform a Gibbs sampler for hurdle joint models focusing on estimating joint models for zero-inflated longitudinal measurements and time-to-event data. This versatile package accommodates various distributional assumptions, including Gaussian, Gamma, inverse Gaussian, Weibull, exponential, beta, Poisson, negative binomial, logarithmic, Bell, generalized Poisson, and binomial. These models, as described in Ganjali et al. (2024), are implemented through two key functions: "ZIJM" for joint modeling with a proportional hazard sub-model and a piecewise constant baseline hazard, and "ZISRE" for joint modeling with a Weibull sub-model incorporating a shared random effects model.

### Installation
To acquire the latest development version of UHJM, you may utilize the following code snippet to install it directly from GitHub:

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/UHJM")
```
This will seamlessly fetch and install the most up-to-date version of UHJM for your use.

### Reference 
Ganjali, M., Baghfalaki, T. & Balakrishnan, N. (2024). A Unified Joint Modeling of Zero-Inflated Longitudinal Measurements and Time-to-Event Outcome. *Submitted*.



