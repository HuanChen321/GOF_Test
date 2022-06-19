#estimating phi which is a monotone function in covariate Z while Z is time-dependent
#allow ties

library(Iso)
library(survival)

isoph_td <- function(Start, Stop, Z, Status,K = median(Z), maxiter = 10^4, eps = 10^-3){
  
  n = length(Status)
  Data0 = data.frame(Start, Stop, Z, Status)
  Data0 = Data0[order(Data0$Z),]
  
  #start from the first observed failure
  location_start <- min(seq(n)[Data0$Status > 0])
  Data <- Data0[location_start:n,]
  
  z = Data$Z
  status = Data$Status
  start = Data$Start
  stop = Data$Stop
  n = length(status)
    
  z.obs = unique(z[status == 1])
  m = length(z.obs)
  
  t.obs = sort(unique(stop))
  nt = length(t.obs)
  
  #anchor
  k = sum(z.obs < K)
  if(k == 0) k=1
  zk = z.obs[k]
  
  #counting process
  Y=dN=matrix(0,n,nt)
  for (i in 1:n) {
    rank.t = which(stop[i]==t.obs)
    Y[i,][1:rank.t] = 1
    if(status[i] == 1) dN[i,][rank.t] = 1
  }
  
  #initial
  z.bar = z - zk
  
  coxfit <- coxph(Surv(start,stop,status) ~ z.bar, control = coxph.control(timefix = FALSE))
  
  beta.hat = max(0.01, coxfit$coefficients)
  psi = beta.hat*(z.obs - zk)
  betaZ = psi
  
  #picm
  Y2 = dN2 = matrix(0, m, nt)#Ystar
  
  z.obs = c(z.obs,Inf)
  
  for(h in 1:m){
    #risk set (z)
    idx = which(z.obs[h] <= z & z < z.obs[h+1])#idx in {1,...,n}
    
    for(ht in idx){
      #risk set (t)
      Y2[h, ] = Y2[h, ] + Y[ht, ]*(start[ht]< t.obs)*(stop[ht]>= t.obs)  
      dN2[h, ] = dN2[h, ] + dN[ht, ]
    }
    
  }
  
  dNsum = colSums(dN2)
  Delta = rowSums(dN2)
  
  iter = 0
  d.e = 1
  while (d.e > eps) {
    iter = iter + 1
    if(iter > maxiter) break
    
    den = colSums(Y2*exp(psi))
    index.zero = which(den > 0)
    
    weight = c()
    for (s in 1:m) 
      weight[s] = sum((Y2[s,]*dNsum/den)[index.zero])
      
    exp.psi.new = pava(Delta/weight, weight)
    psi.new = log(exp.psi.new)
    
    d.e = sum(abs(exp(psi.new) - exp(psi)))
    psi = psi.new
    
  }
  
  psi.new = psi.new - psi.new[k]
  
  conv = 0
  if(d.e < eps) conv = 1
  
  #compute psi.full
  psi.full = rep(NA, length(Data0$Z))
  for (j in 1:m)
    psi.full[which(Data0$Z %in% z.obs[j])] = psi.new[j]
  
  if(is.na(psi.full[1]))
    psi.full[1] = -Inf
  
  which.na = which(is.na(psi.full))
  if(length(which.na) > 0){
    for (j in which.na) {
      psi.full[j] = psi.full[j-1]
    }
  }
  
  
  
  ##negative log partial likelihood
  
  #monotone
  ll_mono1 <- -sum(psi.new*Delta)
  ll_mono2 <- sum(log(1e-3+colSums(Y2*exp(psi.new)))*dNsum)
  ll_mono <- ll_mono1 + ll_mono2
  
  if(is.na(ll_mono)){
    cat('bib')
  }
  
  #cox
  ll_linear1 <- -sum(betaZ*Delta)
  ll_linear2 <- sum(log(colSums(Y2*exp(betaZ)))*dNsum)
  ll_linear <- ll_linear1 + ll_linear2
  if(is.na(ll_linear)){
    cat('bib')
  }
  
  return(list(Z= Data0$Z,psi.hat = psi.full, zk = zk, betahat = beta.hat, conv = conv, ll_linear = ll_linear, ll_mono = ll_mono))
  
}



