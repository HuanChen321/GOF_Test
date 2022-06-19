#Generating censoring time using Kaplan Meier estimator with uniform noise.

CenTime <- function(X, Status, M){
  #X: observed survival time
  #Status: uncensored=1, censored=0
  
  n <- length(X)
  res <- X
  obs <- seq(from = 1, to  = n)[Status == 1]
  Status0 <- Status
  X0 <- X
  
  #sort the data by observed survival time
  Status <- Status[order(X)]
  X <- sort(X)
  
  ind <- 1:n
  LHS <- ((n-ind)/(n-ind+1))^Status
  LHS <- cumprod(LHS)
  KMEst <- LHS[Status == 0]#LHS
  
  Nc <- n - sum(Status)#number of censored cases
  
  Tcen <- X[Status == 0]#vector of censoring time
  num_largeu <- 0#record the time when the inequality doesn't hold
  
  
  for (i in obs) {
    Ci <- min(seq(Nc)[KMEst <= runif(1)])
    if(Ci == Inf){
        res[i] <- Tcen[Nc]#in case min(KMEst) > u[i]
        num_largeu <- num_largeu + 1
      }
    else{
        res[i] <- Tcen[Ci]
      }
      
    #conditiong on generated censoring time is larger than the survival time.
    if (res[i] <= X0[i]){
        res[i] <- runif(n = 1, min = X0[i], max = M)
        #the observed survival time must be smaller than the upperbound of censoring distribution
      }
    
  }
  
  return(list("res" = jitter(res, factor = 2),
              "urate" = num_largeu/n))
}