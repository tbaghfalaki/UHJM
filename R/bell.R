Bell1 <- "model{


  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


Bell <- "model{


  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

    ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])


      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])



    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"



if (is.matrix(XS) == FALSE) {
  if (family == "Bell") {
    model.file <- textConnection(Poisson1)

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }


    parameters <- c("betaL1", "betaL2", "betaS", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")

    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 10000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
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
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
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
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      rownames(sim1$mean$Sigmab) <-
        rownames(sim1$sd$Sigmab) <-
        rownames(sim1$q2.5$Sigmab) <-
        rownames(sim1$q97.5$Sigmab) <-
        rownames(sim1$Rhat$Sigmab) <-
        colnames(sim1$mean$Sigmab) <-
        colnames(sim1$sd$Sigmab) <-
        colnames(sim1$q2.5$Sigmab) <-
        colnames(sim1$q97.5$Sigmab) <-
        colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- sim1$mean$Sigmaa
      D22 <- sim1$mean$Sigmab


      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
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
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      rownames(sim1$mean$Sigmab) <-
        rownames(sim1$sd$Sigmab) <-
        rownames(sim1$q2.5$Sigmab) <-
        rownames(sim1$q97.5$Sigmab) <-
        colnames(sim1$mean$Sigmab) <-
        colnames(sim1$sd$Sigmab) <-
        colnames(sim1$q2.5$Sigmab) <-
        colnames(sim1$q97.5$Sigmab) <- c("Intercept", "Slope")



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- sim1$mean$Sigmaa
      D22 <- sim1$mean$Sigmab


      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }else {
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
          Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
          gamma_lambda = stats::rnorm(1)
        )
      }

      parameters <- c("betaL1", "betaL2", "betaS", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")


      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 10000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
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
        Sigmaa = sim1$sims.list$Sigmaa,
        Sigmab = sim1$sims.list$Sigmab,
        gamma_lambda = sim1$sims.list$gamma_lambda,
        gamma_pi = sim1$sims.list$gamma_pi
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


        names(sim1$mean$h) <-
          names(sim1$sd$h) <-
          names(sim1$q2.5$h) <-
          names(sim1$q97.5$h) <-
          names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


        names(sim1$mean$gamma_lambda) <-
          names(sim1$sd$gamma_lambda) <-
          names(sim1$q2.5$gamma_lambda) <-
          names(sim1$q97.5$gamma_lambda) <-
          names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


        names(sim1$mean$gamma_pi) <-
          names(sim1$sd$gamma_pi) <-
          names(sim1$q2.5$gamma_pi) <-
          names(sim1$q97.5$gamma_pi) <-
          names(sim1$Rhat$gamma_pi) <- "gamma_pi"


        rownames(sim1$mean$Sigmaa) <-
          rownames(sim1$sd$Sigmaa) <-
          rownames(sim1$q2.5$Sigmaa) <-
          rownames(sim1$q97.5$Sigmaa) <-
          rownames(sim1$Rhat$Sigmaa) <-
          colnames(sim1$mean$Sigmaa) <-
          colnames(sim1$sd$Sigmaa) <-
          colnames(sim1$q2.5$Sigmaa) <-
          colnames(sim1$q97.5$Sigmaa) <-
          colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")



        ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
        colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

        D11 <- sim1$mean$Sigmaa
        D22 <- sim1$mean$Sigmab


        MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
        colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



        SM <- rbind(
          cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
          cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
          cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
          cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
        )

        colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



        results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
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


        names(sim1$mean$h) <-
          names(sim1$sd$h) <-
          names(sim1$q2.5$h) <-
          names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


        names(sim1$mean$gamma_lambda) <-
          names(sim1$sd$gamma_lambda) <-
          names(sim1$q2.5$gamma_lambda) <-
          names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


        names(sim1$mean$gamma_pi) <-
          names(sim1$sd$gamma_pi) <-
          names(sim1$q2.5$gamma_pi) <-
          names(sim1$q97.5$gamma_pi) <- "gamma_pi"


        rownames(sim1$mean$Sigmaa) <-
          rownames(sim1$sd$Sigmaa) <-
          rownames(sim1$q2.5$Sigmaa) <-
          rownames(sim1$q97.5$Sigmaa) <-
          colnames(sim1$mean$Sigmaa) <-
          colnames(sim1$sd$Sigmaa) <-
          colnames(sim1$q2.5$Sigmaa) <-
          colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <- c("Intercept", "Slope")



        ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
        colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

        D11 <- sim1$mean$Sigmaa
        D22 <- sim1$mean$Sigmab


        MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
        colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



        SM <- rbind(
          cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
          cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
          cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
          cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
        )


        colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



        results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
      }

  #########################
