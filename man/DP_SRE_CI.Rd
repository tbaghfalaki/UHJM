% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DP_SRE_CI.R
\name{DP_SRE_CI}
\alias{DP_SRE_CI}
\title{Dynamic prediction with credible interval}
\usage{
DP_SRE_CI(
  object,
  s = s,
  t = t,
  mi = 10,
  offset = NULL,
  n.chains = n.chains,
  n.iter = n.iter,
  n.burnin = floor(n.iter/2),
  n.thin = max(1, floor((n.iter - n.burnin)/1000)),
  dataLong,
  dataSurv
)
}
\arguments{
\item{object}{an object inheriting from class VS}

\item{s}{the landmark time for prediction}

\item{t}{the window of prediction for prediction}

\item{mi}{the number of multiple imputation for Monte-Carlo approximation; default is 10.}

\item{offset}{the offset or library size for discrete response. If offset=NULL, it is considered without an offset.}

\item{n.chains}{the number of parallel chains for the model; default is 1.}

\item{n.iter}{integer specifying the total number of iterations; default is 1000.}

\item{n.burnin}{integer specifying how many of n.iter to discard as burn-in ; default is 5000.}

\item{n.thin}{integer specifying the thinning of the chains; default is 1.}

\item{dataLong}{data set of observed longitudinal variables.}

\item{dataSurv}{data set of observed survival variables.}
}
\value{
\itemize{
\item mu.vect list of posterior mean for each parameter
\item sd.vect list of standard error for each parameter
\item 2.5\% list of posterior mode for each parameter
\item 97.5\% list of posterior median for each parameter
\item Rhat Gelman and Rubin diagnostic for all parameter
}
}
\description{
Dynamic prediction for ZISRE with CI
}
\details{
Estimate DP for joint modeling based on VS
}
\examples{
library(survival)
rm(list=ls())

\donttest{
  set.seed(2)
  INDTRAIN <- sample(surv_data_e$id, 0.7 * (dim(surv_data_e)[1]))
  INDVALID <- surv_data_e$id[-INDTRAIN]
  dataLong_t <- subset(
    long_data_e,
    long_data_e$id \%in\% INDTRAIN
  )
  dataSurv_t <- subset(
    surv_data_e,
    surv_data_e$id \%in\% INDTRAIN
  )
  names(dataSurv_t)

  dataLong_v <- subset(
    long_data_e,
    long_data_e$id \%in\% INDVALID
  )
  dataSurv_v <- subset(
    surv_data_e,
    surv_data_e$id \%in\% INDVALID
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
    s = 0.1, t = 0.5, n.chains = 1, offset = NULL, n.iter = 2000, n.burnin = 1000,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )


  DD <- DP_SRE_CI(Z1,
               s = 0.1, t = 0.5, offset = NULL, mi=5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
               n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )

  DPplot2(Z1,
  s = 0.4, id_new = 498, by = 0.1, mi = 2,digits=1,
  offset = NULL,
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
    long_data_b$id \%in\% INDTRAIN
  )
  dataSurv_t <- subset(
    surv_data_b,
    surv_data_b$id \%in\% INDTRAIN
  )
  names(dataSurv_t)

  dataLong_v <- subset(
    long_data_b,
    long_data_b$id \%in\% INDVALID
  )
  dataSurv_v <- subset(
    surv_data_b,
    surv_data_b$id \%in\% INDVALID
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
    s = 0.1, t = 0.5, offset = NULL, n.chains = 1, n.iter = 2000, n.burnin = 500,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )



  ##################
  set.seed(2)
  INDTRAIN <- sample(surv_data_nb$id, 0.7 * (dim(surv_data_nb)[1]))
  INDVALID <- surv_data_nb$id[-INDTRAIN]
  dataLong_t <- subset(
    long_data_nb,
    long_data_nb$id \%in\% INDTRAIN
  )
  dataSurv_t <- subset(
    surv_data_nb,
    surv_data_nb$id \%in\% INDTRAIN
  )
  names(dataSurv_t)

  dataLong_v <- subset(
    long_data_nb,
    long_data_nb$id \%in\% INDVALID
  )
  dataSurv_v <- subset(
    surv_data_nb,
    surv_data_nb$id \%in\% INDVALID
  )


  Z2 <- ZISRE(
    FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1, IStructure=TRUE,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    obstime = "obstime", offset = NULL,
    n.chains = 2,
    n.iter = 200, n.burnin = 100, n.thin = 1, family = "NB"
  )

  DD <- DP_SRE(Z2,
    s = 0.1, t = 0.5, offset = NULL, n.chains = 1, n.iter = 2000, n.burnin = 500,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )


  DD <- DP_SRE_CI(Z2,
                  s = 0.1, t = 0.5, offset = NULL, mi=5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
                  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )

  DPplot2(Z2,
          s = 0.4, id_new = 168, by = 0.1, mi = 2,digits=1,
          offset = NULL,
          Marker_lab="Biomarker", Time_lab="Time (week)",
          n.chains = 1, n.iter = 20, n.burnin = 10,
          dataLong = dataLong_v, dataSurv = dataSurv_v
  )



##################
set.seed(2)
INDTRAIN <- sample(surv_data_n$id, 0.7 * (dim(surv_data_n)[1]))
INDVALID <- surv_data_n$id[-INDTRAIN]
dataLong_t <- subset(
  long_data_n,
  long_data_n$id \%in\% INDTRAIN
)
dataSurv_t <- subset(
  surv_data_n,
  surv_data_n$id \%in\% INDTRAIN
)
names(dataSurv_t)

dataLong_v <- subset(
  long_data_n,
  long_data_n$id \%in\% INDVALID
)
dataSurv_v <- subset(
  surv_data_n,
  surv_data_n$id \%in\% INDVALID
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
             s = 0.1, t = 0.5, offset = NULL, n.chains = 1, n.iter = 2000, n.burnin = 500,
             n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)
}

}
\author{
Taban Baghfalaki \email{t.baghfalaki@gmail.com}
}
