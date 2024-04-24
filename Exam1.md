Getting Started
---------------

```
library(UHJM)
```
Loading the data from the package includes both longitudinal data in long format and survival data. It's essential to ensure that the same subject (ID) is present in both datasets.

```
data(long_data_p)
data(surv_data_p)
```

Dividing data to 70% training data and 30% validation set:

```
set.seed(2)
  INDTRAIN <- sample(surv_data_nb$id, 0.7 * (dim(surv_data_nb)[1]))
  INDVALID <- surv_data_nb$id[-INDTRAIN]
  dataLong_t <- subset(
    long_data_nb,
    long_data_nb$id %in% INDTRAIN
  )
  dataSurv_t <- subset(
    surv_data_nb,
    surv_data_nb$id %in% INDTRAIN
  )
  names(dataSurv_t)

  dataLong_v <- subset(
    long_data_nb,
    long_data_nb$id %in% INDVALID
  )
  dataSurv_v <- subset(
    surv_data_nb,
    surv_data_nb$id %in% INDVALID
  )
```
We are considering one marker; therefore, we require one fixed effects model for rate and probability, one random effects model for , and another for the survival model.

```
FixedY = Y1 ~ obstime + x1 + x2
RandomY = ~obstime
GroupY = ~ id
FixedZ = ~ obstime + x1
RandomZ = ~obstime
GroupZ = ~id
formSurv = Surv(survtime, death) ~ w1 + w2,
```

This versatile package supports a wide range of distributional assumptions, encompassing Gaussian, Gamma, inverse Gaussian, Weibull, exponential, beta, Poisson, negative binomial, logarithmic, Bell, generalized Poisson, and binomial distributions. Further elaboration on the model's specifics can be found in Ganjali et al. (2024).

These joint models are operationalized through two pivotal functions: (1) "ZIJMCV", facilitating joint modeling with a proportional hazard sub-model and a piecewise constant baseline hazard, considering associations based on the current values, and (2) "ZISRE," enabling joint modeling with a Weibull sub-model by incorporating a shared random effects model. At first, we consider the first one:

Finally, we have to use the VS function with the following arguments:

-  FixedY formula for fixed part of longitudinal count model
-  RandomY formula for random part of longitudinal count model
-  GroupY formula specifying the cluster variable for Y (e.g. = ~ subject)
-  FixedZ formula for fixed part of longitudinal probability model
-  RandomZ formula for random part of longitudinal probability model
-  GroupZ formula specifying the cluster variable for Z (e.g. = ~ subject)
-  formSurv formula for survival model
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  obstime the observed time in longitudinal data
-  id the id variable in longitudinal data
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.
-  K Number of nodes and weights for calculating Gaussian quadrature
-  family Family objects provide a convenient way to specify the details of the models. They cover various distributions like "Gaussian", "Exponential", "Weibull", "Gamma", "Beta", "inverse.gaussian", "Poisson", "NB", "Logarithmic", "Bell", "GP", and "Binomial". Specifically, "NB" and "GP" are tailored for hurdle negative binomial and hurdle generalized Poisson joint models, respectively, while the others are utilized for the corresponding models based on their names.

As an example, consider the following command, where this implementation has been performed on training data:

```
Z2 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1 + w2,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    obstime = "obstime", id = "id", n.chains = 2,
    n.iter = 200, n.burnin = 100, n.thin = 1, K = 15, family = "NB"
  )
```

A part of its output of the function is as follows:

```
$Estimation
$Estimation$Count_model
                   Est         SD         L_CI       U_CI     Rhat
(Intercept) -0.9393380 0.15128469 -1.224993616 -0.6484760 1.621542
obstime      0.2625674 0.16164251  0.005090672  0.6192264 1.290469
x1          -0.6079522 0.19236562 -1.060006716 -0.2798539 1.072495
x2          -0.6832724 0.08898479 -0.852994882 -0.5144040 1.133840
Dispersion   5.2285879 0.98801470  3.546213139  7.4822129 1.035696

$Estimation$Zero_inflated_model
                   Est        SD       L_CI       U_CI     Rhat
(Intercept) -1.0435121 0.1543119 -1.3578442 -0.7637013 1.025776
obstime     -0.8803696 0.1796026 -1.2730814 -0.5391646 1.056521
x1           1.0631190 0.1720780  0.7092521  1.3984815 1.017474

$Estimation$Survival_model
                    Est         SD       L_CI       U_CI     Rhat
w1            1.0193130 0.09347805  0.8465882  1.2062719 1.019889
w2           -1.1007597 0.16804981 -1.4314663 -0.7755796 1.005280
gamma_lambda -0.5437266 0.08571440 -0.7250236 -0.3866456 1.113134
gamma_pi     -0.4982187 0.09159170 -0.6773262 -0.3245060 1.002459
h1            0.9093232 0.21471229  0.5623703  1.4385323 1.171186
h2            0.8910990 0.17250013  0.5998604  1.2731592 1.126302
h3            0.8472350 0.14652661  0.5958458  1.1625943 1.099862
h4            0.9301828 0.17185870  0.6405757  1.3017021 1.037230
h5            0.8675282 0.39187958  0.3110176  1.7871397 1.001031

$Estimation$D
$Estimation$D$D11
          Intercept     Slope
Intercept 0.9908800 0.5480709
Slope     0.5480709 1.1476204

$Estimation$D$D22
          Intercept     Slope
Intercept 0.8144997 0.2615043
Slope     0.2615043 0.7440443



$DIC
[1] 5073.027

$LPML
[1] -2296.774
```


For the share random effects, the function ZISRE has the following arguments:

-  FixedY formula for fixed part of longitudinal count model
-  RandomY formula for random part of longitudinal count model
-  GroupY formula specifying the cluster variable for Y (e.g. = ~ subject)
-  FixedZ formula for fixed part of longitudinal probability model
-  RandomZ formula for random part of longitudinal probability model
-  GroupZ formula specifying the cluster variable for Z (e.g. = ~ subject)
-  offset the offset or library size for discrete response. If offset=NULL, it is considered without an offset.
-  obstime the observed time in longitudinal data
-  formSurv formula for survival model
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.
-  family Family objects provide a convenient way to specify the details of the models. They cover various distributions like "Gaussian", "Exponential", "Weibull", "Gamma", "Beta", "inverse.gaussian", "Poisson", "NB", "Logarithmic", "Bell", "GP", and "Binomial". Specifically, "NB" and "GP" are tailored for hurdle negative binomial and hurdle generalized Poisson joint models, respectively, while the others are utilized for the corresponding models based on their names.

As an example, consider the following command, where this implementation has been performed on training data:
```
set.seed(2)
  INDTRAIN <- sample(surv_data_e$id, 0.7 * (dim(surv_data_e)[1]))
  INDVALID <- surv_data_e$id[-INDTRAIN]
  dataLong_t <- subset(
    long_data_e,
    long_data_e$id %in% INDTRAIN
  )
  dataSurv_t <- subset(
    surv_data_e,
    surv_data_e$id %in% INDTRAIN
  )
  names(dataSurv_t)

  dataLong_v <- subset(
    long_data_e,
    long_data_e$id %in% INDVALID
  )
  dataSurv_v <- subset(
    surv_data_e,
    surv_data_e$id %in% INDVALID
  )


  Z1 <- ZISRE(
    FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    obstime = "obstime", offset = NULL,
    n.chains = 2,
    n.iter = 2000, n.burnin = 1000, n.thin = 1, family = "Exponential"
  )

```
A part of the output of this function is as follows: 
```
$Estimation
$Estimation$Y_model
                   Est         SD       L_CI       U_CI     Rhat
(Intercept)  0.9610801 0.10669602  0.7610596  1.1562755 2.890635
x1           0.5646248 0.07778790  0.4059200  0.7160297 1.016472
x2           0.5314088 0.04833851  0.4320189  0.6207531 1.500519
obstime     -0.7546451 0.08694646 -0.9253188 -0.5879938 1.351733

$Estimation$Zero_inflated_model
                    Est         SD       L_CI         U_CI     Rhat
(Intercept) -0.17257617 0.09292774 -0.3634907  0.001927923 1.022518
x1          -0.84544404 0.10879313 -1.0567273 -0.642239328 1.178647
x2          -0.01148273 0.06148162 -0.1317004  0.108729183 1.233322
obstime      1.03868094 0.08752284  0.8678618  1.208407196 1.126586

$Estimation$Survival_model
                     Est         SD       L_CI       U_CI     Rhat
(Intercept) -0.003074646 0.09910825 -0.1885954  0.1863271 2.512048
w1           0.909595438 0.14609170  0.6651618  1.1905946 3.119033
(Intercept)  4.081171087 1.68499032  1.3666209  6.7323928 5.548363
obstime      0.076677448 1.67419538 -2.3656178  3.3874757 4.281811
(Intercept)  2.939603247 2.37930402 -1.5429257  6.3966420 4.364750
obstime     -4.787436631 1.80208901 -7.9888266 -1.3260451 2.283047
Scale        0.731763580 0.10188681  0.5833568  0.9318901 4.087750

$Estimation$D
            (Intercept)     obstime  (Intercept)      obstime
(Intercept)  0.18291289  0.16848338 -0.011761110  0.076022046
obstime      0.16848338  0.17668661 -0.010728000  0.076949456
(Intercept) -0.01176111 -0.01072800  0.017657126 -0.001174521
obstime      0.07602205  0.07694946 -0.001174521  0.046796535


$DIC
[1] 6326.464

$LPML
[1] -2397.208
```




Dynamic prediction
---------------
To reduce estimation biases resulting from variable selection, we propose incorporating an additional stage to calculate dynamic predictions. After variable selection using CS or DS prior, we recommend re-estimating the proportional hazard model by substituting CS or Ds with non-informative normal priors for the association parameters of the selected markers and the regression coefficients of the selected covariates. This has been done by considering *VS2* function in the package. The main arguments in this function are:

- object an object inheriting from class VS function.
- Method the method for variable selection including "LBFDR" for LBFDR and "BF" for Bayes factor.

The following command is considered for this aim:

```
Step2 <- VS2(VS,
  Method = "LBFDR", n.chains = 2, n.iter = 2000, n.burnin = 1000,
  n.thin = 1, dataLong = dataLong_t, dataSurv = dataSurv_t
)
```
Finally, for dynamic prediction, we should utilize the *DP* function, specifying the following arguments:


- object an object inheriting from class VS
- object2 an object inheriting from class VS2
- Method the method for variable selection including "LBFDR" for LBFDR and "BF" for Bayes factor.
- s the landmark time for prediction
- t the window of prediction for prediction
- cause_main the main cause for prediction
- dataLong data set of observed longitudinal variables (validation set).
- dataSurv data set of observed survival variables (validation set).



```
DP <- DP(VS, Step2,
  Method = "LBFDR", s = 0.1, t = 0.5, n.chains = 1, n.iter = 3000, n.burnin = 2000,
  n.thin = 1, cause_main = 1,
  DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)
```

The outputs of this function are as follows:

```
$DP
     id         est
1     2 0.710044463
2     5 0.121999798
3     7 0.155945232
4    10 0.059012615
5    11 0.031432560
6    15 0.059668201
7    16 0.029397494
8    18 0.104735078
9    19 0.134540067
10   21 0.223344760
11   22 0.405120681
12   24 0.009959658
.
.
.
243 488 0.054126771
244 492 0.051628731
245 493 0.089307181
246 494 0.007899571
247 495 0.058865561
248 497 0.125211111
249 498 0.045489998
250 500 0.054972112

$s
[1] 0.1

$t
[1] 0.5
```
Computing AUC and BS for the predictions
---------------
For this purpose, we use DPCri package <https://github.com/tbaghfalaki/DPCri>.

Computing the criteria using this package is straightforward, as demonstrated by the following commands:

- s the landmark time for prediction
- t the window of prediction for prediction
- Survt the survival time
- CR the indicator for competing risks or censoring
- P the risk predictions
- cause the main cause for prediction


Consider the following command: 

```
library(survival)
library(DPCri)

Criteria(
  s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$CR, P = DP$DP$est, cause = 1
)
```
with the following outputs:

```
$Cri
           est         sd
AUC 0.74090038 0.05008035
BS  0.09049356 0.01382442
```

Other beneficial functions
---------

Although all the requirements for computing risk prediction are completed at this stage, we propose some other beneficial functions.


- The first one is *DP0*, which computes dynamic prediction without reestimation and the use of *VS2* functions.

```
DP0 <- DP0(VS, s = 0.1, t = 0.5, n.chains = 1, n.iter = 3000, n.burnin = 2000,
         n.thin = 1, cause_main = 1,
         DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)

$Cri
          est         sd
AUC 0.7344633 0.04313503
BS  0.1251009 0.01484688
```

- The second one is dynamic prediction for one-marker by *DPOM* as follows:

```
D1 <- DPOM(VS,
  N_marker = 1, s = 0.1, t = 0.5, cause_main = 1, n.chains = 1,
  n.iter = 2000, n.burnin = 1000,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE, dataLong_v, dataSurv_v
)

Criteria(
  s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$CR, P = D1$DP$est, cause = 1
)$Cri
```
with the following outputs:

```
$DP
     id        est
1     2 0.35079126
2     5 0.13189275
3     7 0.08357039
4    10 0.12484472
5    11 0.14788052
.
.
.
247 495 0.12207898
248 497 0.08566377
249 498 0.09357884
250 500 0.08983822

$s
[1] 0.1

$t
[1] 0.5


    est         sd
AUC 0.68654215 0.05618483
BS  0.09340079 0.01390426
```

- The third one is dynamic prediction for some markers JM by *DPSM* as follows which markers 1 and 2 are considered:
 
```
D1 <- DPSM(VS, Step2,
  N_markers = c(1, 2), s = 0.1, t = 0.5, cause_main = 1, n.chains = 1,
  n.iter = 2000, n.burnin = 1000,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE, dataLong_v, dataSurv_v
)
```
with the following outputs:

```
$DP
     id         est
1     2 0.636124878
2     5 0.108507783
3     7 0.133040704
4    10 0.062460364
5    11 0.030595321
6    15 0.054849917
7    16 0.029366528
8    18 0.102721955
9    19 0.126131630
10   21 0.216189878
11   22 0.370319961
12   24 0.010574606
.
.
.
247 495 0.054092313
248 497 0.121026926
249 498 0.047753484
250 500 0.053321847

     est         sd
AUC 0.74880268 0.04967864
BS  0.08925278 0.01368505
````


- The last one is the Monte Carlo approximation of dynamic prediction, which is the Monte Carlo version of the *DP* function with a new argument as follows:

- mi the number of multiple imputation for Monte-Carlo approximation; default is 10.


Using this function facilitates the computation of credible intervals for each prediction.

```
MCDP <- MCDP(VS, Step2,
  Method = "LBFDR", s = 0.1, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
  n.thin = 1, cause_main = 1, mi = 10,
  DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)

Criteria(
  s = 0.1, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$CR, P = MCDP$DP$est, cause = 1
)
```
with the following outputs:
```
$DP
     id         est       lower       upper
1     2 0.716386507 0.675167218 0.755475954
2     5 0.114635792 0.107778261 0.124362366
3     7 0.149360907 0.144127688 0.160760439
4    10 0.061768634 0.058670334 0.064408627
5    11 0.031963961 0.030191980 0.034225293
.
.
.
247 495 0.056800414 0.051971578 0.061358109
248 497 0.119515472 0.114530735 0.125109349
249 498 0.045570271 0.042719707 0.048772580
250 500 0.058485596 0.054425111 0.063311557

$Cri
           est         sd
AUC 0.74209770 0.05029602
BS  0.08990641 0.01377804
```
