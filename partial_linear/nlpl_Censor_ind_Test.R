#Only for time-independent covariate currently.

#Calculate negative log partial likelihood for Cox model or monotone hazard model.
############Not work for failures at the same time
nlpl_Censor_ind <- function(events, surv_time, mono_f){
  # flag = get("flag", envir = .GlobalEnv)
  # if(flag){
  #   return (1)
  # }
  #events : 0,1 indicating whether an observed failure happened (vector)
  #surv_time : observed survival time (vector)
  #mono_f : estimated phi value using isoph in isoSurv package (vector). It is sorted by covariate when output by isoph!
  
  star <- events[order(surv_time)]
  mono_f <- mono_f[order(surv_time)]
  X <- sort(surv_time)
  n_star <- sum(star)
  
  #partial likelihood
  denominator <- c()
  risk_set <- cumsum(star)

  ll1 <- -sum(mono_f*star)
  for(i in 1:n_star){
    denominator[i] <- (risk_set == i)%*%exp(mono_f)
  }
  ll2 <- sum(log(rev(cumsum(rev(denominator)))))

  return(ll1+ll2)
}

