library(survival)


\donttest{
  data(long_data_p)
  data(surv_data_p)

Z1 <- ZIJMCV(
  FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = survival::Surv(survtime, death) ~ w1+w2,
  dataLong = long_data_p, dataSurv = surv_data_p,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 30, n.burnin = 20, n.thin = 1, K = 15, family = "Poisson"
)



Z2 <- ZIJMCV(
  FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = survival::Surv(survtime, death) ~ w1+w2,
  dataLong = long_data_p, dataSurv = surv_data_p,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 1000, n.burnin = 500, n.thin = 1, K = 15, family = "NB"
)

data(long_data_nb)
data(surv_data_nb)

Z1 <- ZIJMCV(
  FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = survival::Surv(survtime, death) ~ w1+w2,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Poisson"
)



  Z2 <- ZIJMCV(
  FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1+w2,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id", n.chains = 2,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "NB"
)




Z3 <- ZIJMCV(
  FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1+w2,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 1,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "GP"
)



Z4 <- ZIJMCV(
  FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1+w2,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Bell"
)


Z5 <- ZIJMCV(
  FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1+w2,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Logarithmic"
)

Z6 <- ZIJMCV(
  FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1+w2,
  dataLong = long_data_nb, dataSurv = surv_data_nb,
  obstime = "obstime", id = "id",
  n.chains = 2,
  n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "binomial"
)

}



\donttest{
  data(long_data_n)
  data(surv_data_n)

  Z1 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = survival::Surv(survtime, death) ~ w1 ,
    dataLong = long_data_n, dataSurv = surv_data_n,
    obstime = "obstime", id = "id",
    n.chains = 2,
    n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Gaussian"
  )

  data(long_data_e)
  data(surv_data_e)

  Z2 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = survival::Surv(survtime, death) ~ w1+w2,
    dataLong = long_data_e, dataSurv = surv_data_e,
    obstime = "obstime", id = "id",
    n.chains = 2,
    n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "inverse.gaussian"
  )



    Z2 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = survival::Surv(survtime, death) ~ w1+w2,
    dataLong = long_data_e, dataSurv = surv_data_e,
    obstime = "obstime", id = "id",
    n.chains = 2,
    n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Exponential"
  )


    Z3 <- ZIJMCV(
      FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
      FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
      formSurv = survival::Surv(survtime, death) ~ w1+w2,
      dataLong = long_data_e, dataSurv = surv_data_e,
      obstime = "obstime", id = "id",
      n.chains = 2,
      n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Weibull"
    )


    Z3 <- ZIJMCV(
      FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
      FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
      formSurv = survival::Surv(survtime, death) ~ w1+w2,
      dataLong = long_data_e, dataSurv = surv_data_e,
      obstime = "obstime", id = "id",
      n.chains = 2,
      n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Gamma"
    )


    Z4 <- ZIJMCV(
      FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
      FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
      formSurv = survival::Surv(survtime, death) ~ w1+w2,
      dataLong = long_data_b, dataSurv = surv_data_b,
      obstime = "obstime", id = "id",
      n.chains = 2,
      n.iter = 1000, n.burnin = 500, n.thin = 1, K = 15, family = "Beta"
    )




  }
