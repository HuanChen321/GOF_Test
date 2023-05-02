#Generating survival time by Breslow-type estimator 
#using beta_hat, covariate, observed survival time and Status.
#the output is ordered by corresponding covariate.

SurvTime <- function(beta, MySimu){
  #Z: covariate
  #X: observed survival time.
  #Status: uncensored=1, censored=0.
  #Z = MySimu_ori$covs, X= MySimu_ori$stop, Status = MySimu_ori$event
  
  #sort by observed survival time to evaluate the integral
  MySimu <- MySimu[order(MySimu$stop),]
  Z <- MySimu$covs
  Status <- MySimu$event
  X <- MySimu$stop
  num_largeu <- 0
  
  samp_size <- length(Z)
  nstar <- sum(Status)
  Tstar <- X[Status == 1]
 
  Baseline <- -log(runif(samp_size))/exp(Z*beta)
  
  risk <- cumsum(Status)
  denom <- c()
  for (i in 1:nstar) {
    denom[i] <- sum(exp(beta*Z)[risk >= i])
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
  
  MySimu$X <- res
  MySimu <- MySimu[order(MySimu$covs),]
  
  # #smoothing survival time
  # observed <- seq(samp_size)[MySimu$event > 0]
  # for (i in 1:(nstar-1)) {
  #   nums <- observed[i] - observed[i+1]
  #   if(nums >= 2){
  #     MySimu$X[(observed[i]+1) : (observed[i+1]-1)] <- sort(runif(nums - 1, 
  #                                                      min = MySimu$X[observed[i]], 
  #                                                      max = MySimu$X[observed[i+1]]), decreasing = TRUE)
  #   }
  # }
  # MySimu$X[observed[i+1]:samp_size] <- abs(jitter(MySimu$X[observed[i+1]:samp_size], factor = 2))
  
  return(list("res" = MySimu$X,
              "urates" = num_largeu/samp_size))
}
