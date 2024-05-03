library(survival)
rm(list=ls())

\donttest{
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
    formSurv = Surv(survtime, death) ~ w1, IStructure=TRUE,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    obstime = "obstime", offset = NULL,
    n.chains = 2,
    n.iter = 20, n.burnin = 10, n.thin = 1, family = "Exponential"
  )


  DD <- DP_SRE(Z1,
    s = 0.1, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )


  DPplot2(Z1,
  s = 0.4, id_new = 498, by = 0.1, mi = 2,
  Marker_lab="Biomarker", Time_lab="Time (week)",
  n.chains = 1, n.iter = 20, n.burnin = 10,
  dataLong = dataLong_v, dataSurv = dataSurv_v
)

  ########

  set.seed(2)
  INDTRAIN <- sample(surv_data_b$id, 0.7 * (dim(surv_data_b)[1]))
  INDVALID <- surv_data_b$id[-INDTRAIN]
  dataLong_t <- subset(
    long_data_b,
    long_data_b$id %in% INDTRAIN
  )
  dataSurv_t <- subset(
    surv_data_b,
    surv_data_b$id %in% INDTRAIN
  )
  names(dataSurv_t)

  dataLong_v <- subset(
    long_data_b,
    long_data_b$id %in% INDVALID
  )
  dataSurv_v <- subset(
    surv_data_b,
    surv_data_b$id %in% INDVALID
  )


  Z1 <- ZISRE(
    FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ x1 + x2 + obstime, RandomZ = ~1, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1, IStructure=TRUE,
    obstime = "obstime", offset = NULL,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    n.chains = 2,
    n.iter = 200, n.burnin = 100, n.thin = 1, family = "Beta"
  )



  DD <- DP_SRE(Z1,
    s = 0.1, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 500,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )



  ##################
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


  Z2 <- ZISRE(
    FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1, IStructure=TRUE,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    obstime = "obstime", offset = NULL,
    n.chains = 2,
    n.iter = 200, n.burnin = 100, n.thin = 1, family = "Bell"
  )

  DD <- DP_SRE(Z2,
    s = 0.1, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 500,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )



##################
set.seed(2)
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


Z2 <- ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~1, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1, IStructure=FALSE,
  dataLong = dataLong_t, dataSurv = dataSurv_t,
  obstime = "obstime", offset = NULL,
  n.chains = 2,
  n.iter = 20, n.burnin = 10, n.thin = 1, family = "Gaussian"
)

DD <- DP_SRE(Z2,
             s = 0.1, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 500,
             n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)
}

