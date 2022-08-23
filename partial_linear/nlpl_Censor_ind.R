############Not work for failures at the same time
nlpl_Censor_ind <- function(data, zk, mono_f = NA, est_beta = NA, W){
                            # events = MySimu$event, 
                            # surv_time = MySimu$stop,
                            # uni_cov = MySimu$covs,
  
  # flag = get("flag", envir = .GlobalEnv)
  # if(flag){
  #   return (1)
  # }
  
  
  #events : 0,1 indicating whether an observed failure happened (vector)
  #surv_time : observed survival time (vector)
  #z : monotone covariate (vector)
  #W: linear covariate (matrix)
  #mono_f : estimated phi value using isoph in isoSurv package (vector). It is sorted by covariate when output by isoph!
  #est_beta : estimated parameter for linear hazard if applicable (scalar)
  surv_time <- data$stop
  W <- as.data.frame(W)[order(surv_time), ]
  z <- data$covs[order(surv_time)]
  star <- data$event[order(surv_time)]
  mono_f <- mono_f[order(surv_time)]
  X <- sort(surv_time)
  n_star <- sum(star)

  #partial likelihood
  denominator <- c()
  risk_set <- cumsum(star)
  
  #linear
  if(is.na(mono_f)[1]){
    #denominator
    ll1 <- -sum((est_beta[1]*(z - zk) + as.matrix(W)%*%as.matrix(est_beta[2:length(est_beta)]))*star)
    for(i in 1:n_star){
      denominator[i] = (risk_set >= i)%*%exp(est_beta[1]*(z - zk) + as.matrix(W)%*%as.matrix(est_beta[2:length(est_beta)]))#sum for all T_j > = T_i(exp(beta*z_j))
    }
    ll2 <- sum(log(denominator))
  }
  
  #monotone
  else{
    ll1 <- -sum((mono_f + as.matrix(W)%*%as.matrix(est_beta))*star)
    for(i in 1:n_star){
      denominator[i] <- (risk_set == i)%*%exp(mono_f + as.matrix(W)%*%as.matrix(est_beta))
    }
    ll2 <- sum(log(rev(cumsum(rev(denominator)))))
  }
  
  
  return(ll1+ll2)
}
