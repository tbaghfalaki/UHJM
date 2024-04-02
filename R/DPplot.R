DPplot=function(object, s = s, id=id, n.chains = n.chains, n.iter = n.iter, n.burnin = floor(n.iter / 2),
                n.thin = max(1, floor((n.iter - n.burnin) / 1000)), dataLong, dataSurv){
  Z1=object




  DD <- DP_SRE(Z1,
               s = s, t = 0.5, n.chains = 1, n.iter = n.iter, n.burnin = n.burnin,
               dataLong = dataLong, dataSurv = dataSurv
  )


}
