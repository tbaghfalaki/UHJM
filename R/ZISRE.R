#' Zero-inflation joint modeling
#'
#' @description
#' Fits zero-inflated hurdle shared random effects models under Gaussian, Gamma, inverse Gaussian, Weibull, exponential, beta, Poisson, negative binomial, logarithmic, Bell, generalized Poisson, and binomial.
#'
#' @details
#' Function using 'JAGS' software to estimate the hurdle joint modeling
#'
#' @param FixedY formula for fixed part of longitudinal count model
#' @param RandomY formula for random part of longitudinal count model
#' @param GroupY formula specifying the cluster variable for Y (e.g. = ~ subject)
#' @param FixedZ formula for fixed part of longitudinal probability model
#' @param RandomZ formula for random part of longitudinal probability model
#' @param GroupZ formula specifying the cluster variable for Z (e.g. = ~ subject)
#' @param offset the offset or library size for discrete response. If offset=NULL, it is considered without an offset.
#' @param obstime the observed time in longitudinal data
#' @param formSurv formula for survival model
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#' @param n.thin integer specifying the thinning of the chains; default is 1.
#' @param family Family objects provide a convenient way to specify the details of the models. They cover various distributions like "Gaussian", "Exponential", "Weibull", "Gamma", "Beta", "inverse.gaussian", "Poisson", "NB", "Logarithmic", "Bell", "GP", and "Binomial". Specifically, "NB" and "GP" are tailored for hurdle negative binomial and hurdle generalized Poisson joint models, respectively, while the others are utilized for the corresponding models based on their names.
#'
#'
#' @importFrom
#'
#' @return
#' - mu.vect list of posterior mean for each parameter
#' - sd.vect list of standard error for each parameter
#' - 2.5% list of posterior mode for each parameter
#' - 97.5% list of posterior median for each parameter
#' - Rhat Gelman and Rubin diagnostic for all parameter
#' - betaL1 the regression coefficients of rate of the hurdle model
#' - betaL2 the regression coefficients of probability of the hurdle model
#' - betaS the regression coefficients of the survival model
#' - Sigma the variance of random effects
#' - gamma the association parameters
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}, Mojtaba Ganjali \email{m-ganjali@sbu.ac.ir}
#'
#'
#' @example inst/exampleZISRE.R
#'
#' @md
#' @export

ZISRE <- function(FixedY, RandomY, GroupY, FixedZ, RandomZ, GroupZ, formSurv, dataLong, dataSurv,
                  n.chains = n.chains, obstime = "obstime", offset = NULL,
                  n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, family = "Poisson") {
  data_long <- dataLong[unique(c(
    all.vars(GroupY), all.vars(FixedY), all.vars(RandomY),
    all.vars(GroupZ), all.vars(FixedZ), all.vars(RandomZ), offset
  ))]
  y <- data_long[all.vars(FixedY)][, 1]
  mfX <- stats::model.frame(FixedY, data = data_long)
  X1 <- stats::model.matrix(FixedY, mfX)
  mfU <- stats::model.frame(RandomY, data = data_long)
  Z1 <- stats::model.matrix(RandomY, mfU)
  # id_prim <- as.integer(data_long[all.vars(GroupY)][, 1])
  id <- as.integer(data_long[all.vars(GroupY)][, 1])


  M <- table(id)
  id_prim <- rep(1:length(M), M)

  mfX2 <- stats::model.frame(FixedZ, data = data_long)
  X2 <- stats::model.matrix(FixedZ, mfX2)
  mfU2 <- stats::model.frame(RandomZ, data = data_long)
  Z2 <- stats::model.matrix(RandomZ, mfU2)

  n2 <- length(unique(id_prim))

  n <- length(dataLong[, 1])
  z <- rep(0, n)
  z[y == 0] <- 1
  ## Number of patients and number of longitudinal observations per patient
  n1 <- length(y)
  Nbeta1 <- dim(X1)[2]
  Nbeta2 <- dim(X2)[2]



  if (length(offset) == 0) {
    offset <- rep(0, n1)
  } else {
    offset <- data_long[offset][, 1]
  }


  #####################
  # formSurv=survival::formSurv
  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  death <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
  nTime <- length(Time) # number of subject having Time
  # design matrice
  mfZ <- stats::model.frame(formSurv, data = tmp)
  XS <- stats::model.matrix(formSurv, mfZ)



  surt.cen <- rep(0, n2)
  for (i in 1:n2) {
    if (death[i] == 1) (surt.cen[i] <- 10000)
    if (death[i] == 0) ((surt.cen[i] <- Time[i]) & (Time[i] <- NA))
  }
  ###########
  Bell_mean <- "model{


    for(i in 1:n){
          cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

       ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

     log(theta1[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


for(jj in 1:10){

tab[i,jj]=pow(-jj,(jj-1))/exp(logfact(jj))*pow(theta1[i],jj)
}

theta[i]<-sum(tab[i,])
    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"
  ###########
  Bellwc_mean <- "model{


    for(i in 1:n){
    cpoinvy[i]<- 1/(0.0001+exp(ll[i]))
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

       ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

    log(theta1[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


for(jj in 1:10){

tab[i,jj]=pow(-jj,(jj-1))/exp(logfact(jj))*pow(theta1[i],jj)
}

theta[i]<-sum(tab[i,])
    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])

    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"




  ##########
  NB <- "model{

  for(i in 1:n){
  cpoinvy[i]<- 1/(0.0001+exp(ll[i]))
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+loggam(r+y[i])-loggam(r)-loggam(y[i]+1)+r*log(r/(r+lambda[i]))+y[i]*log(lambda[i]/(lambda[i]+r))-log(1-pow(r/(r+lambda[i]),r)))

      log(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])
  }
    #####
   for(k in 1:n2){
   cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])


    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
  r~ dgamma(0.1,0.1)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }

 Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))

   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

}"

  ###########
  Poisson <- "model{


    for(i in 1:n){
    cpoinvy[i]<- 1/(0.0001+exp(ll[i]))
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

      ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))


     log(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"

  ###########
  Bell <- "model{


    for(i in 1:n){
    cpoinvy[i]<- 1/(0.0001+exp(ll[i]))
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

       ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

     log(theta[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])

#theta1[i]<-theta[i]*exp(theta[i])
    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])

    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"
  ###########
  Bellwc <- "model{


    for(i in 1:n){
    cpoinvy[i]<- 1/(0.0001+exp(ll[i]))
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

       ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

     log(theta[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])

    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"

  ###########
  logar <- "model{


    for(i in 1:n){
    cpoinvy[i]<- 1/(0.0001+exp(ll[i]))
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

ll[i]<-(1-z[i])*(log(-1/log(1-pi[i]))+y[i]*log(pi[i])-log(y[i]))+z[i]*log(lambda[i])+
  (1-z[i])*log(1-lambda[i])

     logit(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])




    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])

    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"

  ###########
  binomial <- "model{

  m<-max(y)

    for(i in 1:n){
    cpoinvy[i]<- 1/(0.0001+exp(ll[i]))
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

ll[i]<-(1-z[i])*(loggam(m+1)-loggam(y[i]+1)-loggam(m-y[i]+1)+y[i]*log(pi[i])+(m-y[i])*log(1-pi[i])-
                   log(1-pow(1-pi[i],m)))+z[i]*log(lambda[i])+(1-z[i])*log(1-lambda[i])


     logit(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])



    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"
  ###########
  GP <- "model{

    for(i in 1:n){
          cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

      ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-y[i]*log(1+phiz*lambda[i])+(y[i]-1)*log(1+phiz*y[i])- loggam(y[i]+1)-
        lambda[i]*(1+phiz*y[i])/(1+phiz*lambda[i])-log(1-exp(-lambda[i]/(1+phiz*lambda[i]))))

      log(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])

    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)

    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }
      phiz~dnorm(0,.001)

  Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))

   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"

  Exp <- "model{

    for(i in 1:n){
          cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

ll[i] <- (1-z[i])*(logdensity.exp(y[i],lambda[i]))+z[i]*log(pi[i])+(1-z[i])*log(1-pi[i])


     log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"


  ######
  Gamma <- "model{

    for(i in 1:n){
          cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.gamma(y[i],sigma1, mu1[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

 mu1[i]<-sigma1/mu[i]
     log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}


sigma1~dgamma(.1,.1)
  sigma<-1/sigma1
  }"
  #####

  Beta <- "model{

    for(i in 1:n){
          cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.beta(y[i],a1[i], b1[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

a1[i]<-mu[i]*phi11
    b1[i]<-(1-mu[i])*phi11

     logit(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}
  phi11~dgamma(.1,.1)
  phi1<-1/phi11

  }"
  #####
  Weibull <- "model{

    for(i in 1:n){
          cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

ll[i] <- (1-z[i])*(logdensity.weib(y[i],kappa, mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


     log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])



    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}


  kappa~dgamma(.1,.1)

  }"

  Gaussian <- "model{

    for(i in 1:n){
          cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.norm(y[i], mu[i],tau))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

    mu[i] <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])



    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}


 tau~dgamma(.01,.01)
  sigma<-1/tau
  }"




  IGauss <- "model{

    for(i in 1:n){
          cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

 ll[i]<- (1-z[i])*(0.5*(log(lambda)-log(2*3.14)-3*log(0.000000001+y[i]))-
 0.5*lambda*pow((y[i]-mu[i]),2)/(pow(mu[i],2)*(0.000000001+y[i]))+log(1-muz[i]))+z[i]*log(muz[i])

     log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
cpoinvt[k]<- (1/(0.0001+dweib(Time[k],p,mut[k])))*death[k]+
(1/(0.0001+1-pweib(surt.cen[k],p,mut[k])))*(1-death[k])
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

  }"











  ########################################################################################


  if (family == "Exponential") {
    model.file <- textConnection(Exp)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
      sigma = sim1$sims.list$sigma
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"




      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- rbind(cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2))
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"




      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))





      ZPM <- rbind(cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2))

      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
  }
  ##
  #####
  if (family == "Beta") {
    model.file <- textConnection(Beta)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p", "phi1")

    y[y == 1] <- 0.99999

    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
      phi = sim1$sims.list$phi
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      names(sim1$mean$phi) <-
        names(sim1$sd$phi) <-
        names(sim1$q2.5$phi) <-
        names(sim1$q97.5$phi) <-
        names(sim1$Rhat$phi) <- "Dispersion"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$phi, sim1$sd$phi, sim1$q2.5$phi, sim1$q97.5$phi, sim1$Rhat$phi)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      names(sim1$mean$phi) <-
        names(sim1$sd$phi) <-
        names(sim1$q2.5$phi) <-
        names(sim1$q97.5$phi) <- "Dispersion"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))





      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)

      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$phi, sim1$sd$phi, sim1$q2.5$phi, sim1$q97.5$phi)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
  }
  ###
  if (family == "Weibull") {
    model.file <- textConnection(Weibull)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p", "kappa")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
      kappa = sim1$sims.list$kappa
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      names(sim1$mean$kappa) <-
        names(sim1$sd$kappa) <-
        names(sim1$q2.5$kappa) <-
        names(sim1$q97.5$kappa) <-
        names(sim1$Rhat$kappa) <- "Scale"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$kappa, sim1$sd$kappa, sim1$q2.5$kappa, sim1$q97.5$kappa, sim1$Rhat$kappa)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      names(sim1$mean$kappa) <-
        names(sim1$sd$kappa) <-
        names(sim1$q2.5$kappa) <-
        names(sim1$q97.5$kappa) <- "Scale"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))





      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)

      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$kappa, sim1$sd$kappa, sim1$q2.5$kappa, sim1$q97.5$kappa)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
  }
  ##
  if (family == "inverse.gaussian") {
    model.file <- textConnection(IGauss)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p", "sigma")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
      sigma = sim1$sims.list$sigma
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Shape"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Shape"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))





      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)

      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
  }
  ###
  if (family == "Gamma") {
    model.file <- textConnection(Gamma)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p", "sigma")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
      sigma = sim1$sims.list$sigma
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Scale"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Scale"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))





      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)

      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
    ###############################################
  }
  ####
  if (family == "Gaussian") {
    model.file <- textConnection(Gaussian)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p", "sigma")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
      sigma = sim1$sims.list$sigma
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Variance"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Variance"

      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))





      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)

      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Y_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
    ###############################################
  }

  if (family == "Poisson") {
    model.file <- textConnection(Poisson)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, offset = offset, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
  }
  ####################### $$$$$$$$$$$$$
  if (family == "Bell") {
    model.file <- textConnection(Bell)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p")
    C <- c()
    for (i in 1:n) {
      C[i] <- log(numbers::bell(y[i])) - lfactorial(y[i])
    }


    if (is.infinite(numbers::bell(max(y))) == TRUE) {
      model.file <- textConnection(Bellwc)
    }


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, C = C, offset = offset, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
  }


  ####################### $$$$$$$$$$$$$
  if (family == "Logarithmic") {
    model.file <- textConnection(logar)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, offset = offset, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
    ###############################################
  }

  if (family == "binomial") {
    model.file <- textConnection(binomial)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "p")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, offset = offset, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
    ###############################################
  }
  if (family == "NB") {
    model.file <- textConnection(NB)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS), r = 1,
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "Sigma", "gamma", "r", "p", "cpoinvy", "cpoinvt")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, offset = offset, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )


    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
      r = sim1$sims.list$r
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <-
        names(sim1$Rhat$r) <- "Dispersion"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r, sim1$Rhat$r)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")


      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <-
        names(sim1$Rhat$r) <- "Dispersion"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }

    #########
  }

  if (family == "GP") {
    model.file <- textConnection(GP)
    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS), phiz = 1,
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigma", "gamma", "phiz", "p")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 100000, is.censored=1-death,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS, offset = offset, death = death,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
      XS = XS
    )

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
      r = sim1$sims.list$phiz
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <-
        names(sim1$Rhat$phiz) <- "Dispersion"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz, sim1$Rhat$phiz)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")


      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <-
        names(sim1$Rhat$phiz) <- "Dispersion"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
  }

  ############

  K <- 100000
  DIC <- sim1$DIC - 2 * n * K
  LPML <- -sum(log(sim1$mean$cpoinvy)) - sum(log(sim1$mean$cpoinvt))




  if (family == "Bell") {
    if (is.infinite(numbers::bell(max(y))) == TRUE) {
      DIC <- "It is not possible to compute the DIC because the values of the Bell number are not finite."
      LPML <- "It is not possible to compute the LPML because the values of the Bell number are not finite."
    }
  }

  list(
    FixedY = FixedY, FixedZ = FixedZ, formSurv = formSurv,
    RandomY = RandomY, RandomZ = RandomZ, GroupY = GroupY, GroupZ = GroupZ, id = id,
    obstime = obstime,
    family = family, MCMC = MCMC, Estimation = results,
    DIC = DIC, LPML = LPML
  )
}
