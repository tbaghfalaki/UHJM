library(survival)
rm(list=ls())

\donttest{

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


  Z2 <- ZIJMCV(
  FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1+w2,
  dataLong = dataLong_t, dataSurv = dataSurv_t,
  obstime = "obstime", id = "id", n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, K = 15, family = "NB"
)


DD=DP_CV(object=Z2, s = 0.5, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
        n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v)


  Criteria(
    s = 0.5, t = 0.5, Survt = dataSurv_v$survtime,
    CR = dataSurv_v$death, P = DD$DP$est, cause = 1
  )$Cri



}



\donttest{
  rm(list=ls())
  INDTRAIN <- sample(surv_data_n$id, 0.7 * (dim(surv_data_n)[1]))
  INDVALID <- surv_data_n$id[-INDTRAIN]
  dataLong_t <- subset(
    long_data_n,
    long_data_n$id %in% INDTRAIN
  )
  dataSurv_t <- subset(
    surv_data_n,
    surv_data_n$id %in% INDTRAIN
  )
  names(dataSurv_t)

  dataLong_v <- subset(
    long_data_n,
    long_data_n$id %in% INDVALID
  )
  dataSurv_v <- subset(
    surv_data_n,
    surv_data_n$id %in% INDVALID
  )

  Z1 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2 , RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = survival::Surv(survtime, death) ~ w1 ,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    obstime = "obstime", id = "id",
    n.chains = 2,
    n.iter = 2000, n.burnin = 1000, n.thin = 1, K = 15, family = "Gaussian"
  )


  DD=DP_CV(object=Z1, s = 0.5, t = 0.5, n.chains = 1, n.iter = 1000, n.burnin = 500,
        n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v)


  Criteria(
    s = 0.5, t = 0.5, Survt = dataSurv_v$survtime,
    CR = dataSurv_v$death, P = DD$DP$est, cause = 1
  )$Cri




  }