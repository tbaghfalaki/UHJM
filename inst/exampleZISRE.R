library(survival)


\donttest{
  data(long_data_nb)
 data(surv_data_nb)
 Z1 <- ZISRE(
   FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
   FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
   formSurv = Surv(survtime, death) ~ w1,
   obstime = "obstime", offset = NULL,
   dataLong = long_data_nb, dataSurv = surv_data_nb,
   n.chains = 2,
   n.iter = 200, n.burnin = 100, n.thin = 1, family = "Poisson"
 )


 Z2 <- ZISRE(
   FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
   FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
   formSurv = Surv(survtime, death) ~ w1,
   obstime = "obstime", offset = NULL,
   dataLong = long_data_nb, dataSurv = surv_data_nb,
   n.chains = 2,
   n.iter = 200, n.burnin = 100, n.thin = 1, family = "NB"
 )

 Z3 <- ZISRE(
   FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
   FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
   formSurv = Surv(survtime, death) ~ w1,
   obstime = "obstime", offset = NULL,
   dataLong = long_data_nb, dataSurv = surv_data_nb,
   n.chains = 2,
   n.iter = 200, n.burnin = 100, n.thin = 1, family = "GP"
 )


 Z4 <- ZISRE(
   FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
   FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
   formSurv = Surv(survtime, death) ~ w1,
   obstime = "obstime", offset = NULL,
   dataLong = long_data_nb, dataSurv = surv_data_nb,
   n.chains = 2,
   n.iter = 200, n.burnin = 100, n.thin = 1, family = "Bell"
 )

 Z5 <- ZISRE(
   FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
   FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
   formSurv = Surv(survtime, death) ~ w1,
   obstime = "obstime", offset = NULL,
   dataLong = long_data_nb, dataSurv = surv_data_nb,
   n.chains = 2,
   n.iter = 200, n.burnin = 100, n.thin = 1, family = "Logarithmic"
 )


 Z6 <- ZISRE(
   FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
   FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
   formSurv = Surv(survtime, death) ~ w1,
   obstime = "obstime", offset = NULL,
   dataLong = long_data_nb, dataSurv = surv_data_nb,
   n.chains = 2,
   n.iter = 200, n.burnin = 100, n.thin = 1, family = "binomial"
 )
}





\donttest{

  data(long_data_n)
data(surv_data_n)
Z1 <- ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  obstime = "obstime", offset = NULL,
  dataLong = long_data_n, dataSurv = surv_data_n,
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, family = "Gaussian"
)

data(long_data_b)
data(surv_data_b)
Z2 <- ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  obstime = "obstime", offset = NULL,
  dataLong = long_data_b, dataSurv = surv_data_b,
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, family = "Beta"
)

data(long_data_e)
data(surv_data_e)
Z3 <- ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  obstime = "obstime", offset = NULL,
  dataLong = long_data_e, dataSurv = surv_data_e,
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, family = "Exponential"
)


Z4 <- ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  obstime = "obstime", offset = NULL,
  dataLong = long_data_e, dataSurv = surv_data_e,
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, family = "Weibull"
)

Z5 <- ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  obstime = "obstime", offset = NULL,
  dataLong = long_data_e, dataSurv = surv_data_e,
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, family = "Gamma"
)


Z6 <- ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  obstime = "obstime", offset = NULL,
  dataLong = long_data_e, dataSurv = surv_data_e,
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, family = "inverse.gaussian"
)
}
