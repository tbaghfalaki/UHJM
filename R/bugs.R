Exp1 <- "model{


  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.exp(y[i],mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
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


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"


Exp <- "model{


  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.exp(y[i],mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])



      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
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


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"
