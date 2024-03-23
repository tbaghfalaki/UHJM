library(survival)


\donttest{
  data(long_data_nb)
  data(surv_data_nb)

Z1 <- ZIJMCV(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = survival::Surv(survtime, death) ~ w1,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 1,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Poisson"
)


Z2 <- ZIJMCV(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "NB"
)



Z3 <- ZIJMCV(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 1,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "GP"
)



Z4 <- ZIJMCV(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Bell"
)


Z5 <- ZIJMCV(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Logarithmic"
)

Z6 <- ZIJMCV(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "binomial"
)

}
