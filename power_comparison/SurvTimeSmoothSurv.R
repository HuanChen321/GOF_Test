#Generating survival time by Breslow-type estimator using monotone hazard model or cox model
#using beta_hat, covariate, observed survival time and Status.
#the output is ordered by corresponding covariate.

SurvTime <- function(psi, data){
  #Z = data$covs, X= data$stop, Status = data$event
  ##Z: monotone covariate
  ##W: linear covariates
  ##X: observed survival time
  ##Status: uncensored=1, censored=0
  
  
  #sort by observed survival time to evaluate the integral
  psi <- psi[order(data$stop)]
  data <- data[order(data$stop), ]
  Z <- data$covs
  Status <- data$event
  X <- data$stop
  
  num_largeu <- 0
  
  samp_size <- length(Z)
  nstar <- sum(Status)
  Tstar <- X[Status == 1]
 
  Baseline <- -log(runif(samp_size))/exp(psi)
  
  risk <- cumsum(Status)
  denom <- c()
  for (i in 1:nstar) {
    denom[i] <- sum(exp(psi)[risk >= i])
  }
  BSEst <- cumsum(1/denom)
  
  res <- c()
  for (i in 1:samp_size) {
    Ti <- min(seq(nstar)[BSEst >=Baseline[i]])
    if(Ti == Inf){
      res[i] <- Tstar[nstar]##in case max(BSEst)<Baseline[i]
      num_largeu <- num_largeu + 1
    }
    else{
      res[i] <- Tstar[Ti]
    }
  }
  
  data$X <- res
  data <- data[order(data$covs),]
  
  #smoothing survival time
  observed <- seq(samp_size)[data$event > 0]
  for (i in 1:(nstar-1)) {
    nums <- observed[i] - observed[i+1]
    if(nums >= 2){
      data$X[(observed[i]+1) : (observed[i+1]-1)] <- sort(runif(nums - 1, 
                                                       min = data$X[observed[i]], 
                                                       max = data$X[observed[i+1]]), decreasing = TRUE)
    }
  }
  data$X[observed[i+1]:samp_size] <- abs(jitter(data$X[observed[i+1]:samp_size], factor = 2))
  return(list("res" = data$X,
              "urates" = num_largeu/samp_size))
}
