#Negative Log Partial Likelihood without Censoring for Time-independent Covariate

nlpl_noCensor_ind<- function(mono_f = NA, est_beta = NA, zk, uni_cov, surv_time){
  # flag = get("flag", envir = .GlobalEnv)
  # if(flag){
  #   return (1)
  # }

  #mono_f : phi, the monotone estimator (vector). Already sorted if it comes from isoph.
  #est_beta : parameter for linearity (scalar)
  #uni_cov : covariate (vector)
  #surv_time : survival time (vector)
  
  #Sort covariate and survival time by covariate
  X <- surv_time[order(uni_cov)] # sorted survival time
  z <- sort(uni_cov) #sorted covariate
  
  #then sort survival time, covariate and psi by time
  phi <- mono_f[order(X)]
  z <- z[order(X)]
  
  #partial likelihood
  if(is.na(est_beta)){
    #monotone
    ll1 <- -sum(phi)
    ll2 <- sum(log(rev(cumsum(rev(exp(phi))))))
    ll <- ll1 + ll2
  }
  else if(is.na(mono_f)){
    #linear
    linear <- est_beta*(z - zk)
    ll1 <- -sum(linear)
    ll2 <- sum(log(rev(cumsum(rev(exp(linear))))))
    ll <- ll1 + ll2
  }
  return(ll)
}