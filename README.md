## Unified Zero-Inflated Hurdle Joint Modeling
Performing a Gibbs sampler for hurdle joint models involves estimating joint models for zero-inflated longitudinal measurements and time-to-event data. This versatile package accommodates various distributional assumptions, including Gaussian, Gamma, inverse Gaussian, Weibull, exponential, beta, Poisson, negative binomial, logarithmic, Bell, generalized Poisson, and binomial. As described in Ganjali et al. (2025), these models are implemented through two key functions:

- "ZIJMCV" facilitates joint modeling with a proportional hazard sub-model and a piecewise constant baseline hazard, considering associations based on the current values.
- "ZISRE" enables joint modeling with a Weibull sub-model by incorporating a shared random effects model.


### Installation
To acquire the latest development version of UHJM, you may utilize the following code snippet to install it directly from GitHub:

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/UHJM")
```
This will seamlessly fetch and install the most up-to-date version of UHJM for your use.


### Example Usage

This analysis is presented [here](/Exam1.md)

### Reference 
Ganjali, M., Baghfalaki, T., Balakrishnan, N. and Jacqmin-Gadda, H. (2025). A Unified Joint Modeling of Zero-Inflated Longitudinal Measurements and Time-to-Event Outcomes with Applications to HIV and Colorectal Cancer Data. *Journal of Statistical Computation and Simulation*, 10.1080/00949655.2025.2588591.



