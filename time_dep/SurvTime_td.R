#Generating survival time by Breslow-type estimator 
#using beta_hat, covariate, observed survival time and Status.
#the output has the same order as Z

SurvTime_td <- function(beta, OST, Status, Z, X){
  #beta: beta_hat
  #OST: Observed survival time, pmin(T, C)
  #Status: uncensored=1, censored=0.
  #Z: time-dependent covariate, a list
  #X: jump points for Z
  X <- X[-1]#drop the first 0 since we only need the right endpoint
  
  #define Z^star_(i)
  n <- length(Status)
  loc <- seq(n)[Status > 0]
  
  zn_star <- sum(Status)
  Z_star <- c()
  
  for(i in 1:zn_star){
      Z_star[i] <- Z[[loc[i]]][length(Z[[loc[i]]])]
  }
  
  Z_star <- sort(unique(Z_star))#unique covariate for failed subjects at the time of failure
  zn_star <- length(Z_star)#number of risk set
  I <- c(Z_star, Inf)
  
  #Brewslow-type estimator for cumulative baseline hazard 
  
  #dN_star and T_star
  T_fail <- OST[Status > 0]
  T_star <- sort(unique(T_fail))#unique time for failure
  tn_star <- length(T_star)
  
  dN_star <- c()
  for (i in 1:tn_star) {
    dN_star[i] <- sum(T_fail == T_star[i])
  }
  
  #i. Lambda0hat
  dLambda0hat <- c()
  
  Y_star <- matrix(nrow = tn_star, ncol = zn_star)
  #Y_star[i,j] = Ystar_j(T_i)
  #column j constains R_j at T_star[1] to T_star[tn_star] 
  
  for (i in 1:tn_star) {
    
    #R = [R1(T_star_i),...,Rzn_star(T_star_i)]
    
    for (j in 1:zn_star) {
      
      #R_j(T_star_i)
      R <- c()
      
      for (h in 1:n) {
        #Is Z[[h]](T_star_i) in R_j(T_star_i)
        
        m <- length(Z[[h]])
        
        if(OST[h] >= T_star[i])#(1)
          {
          #(2)
          Z_hT <- Z[[h]][min(seq(m)[(X >= T_star[i])[1:m]])]
          #(3)
          if((Z_hT >= I[j]) & (Z_hT < I[j+1]))
            R <- c(R, h)
        }
      }
      Y_star[i,j] <- length(R)
    }
  }
  
  ##dLambda0hat
  #use Y_star, dN_star and Z_star
  num <- dN_star
  denom <- colSums(t(Y_star)*exp(beta*Z_star))
  
  dLambda0hat <- num/denom
  Lambda0hat <- cumsum(dLambda0hat)
  
  
  #TX
  TX = sort(unique(c(X, T_star)))
  nTX <- length(TX)
  Lambda0TX <- rep(Lambda0hat[1], nTX)
  
  for (i in 2:tn_star) {
    loc <- min(seq(nTX)[TX >= T_star[i]])
    Lambda0TX[loc:nTX] <- Lambda0hat[i]
  }
  
  
  
  
  #ii. w
  
  SurvTime <- c()
  urate <- 0
  
  for (j in 1:n) {
    Zj <- c()
    
    if(OST[j] < min(X, T_star))######################corrected
      TXj <- OST[j]
    else
      TXj <- sort(unique(c(X[X<=OST[j]],T_star[T_star <= OST[j]])))
    
    nj <- length(TXj)
    loc <- 1
    i <- 1
  
    while(i <= nj) {
      if(TXj[i] <= X[loc]){
        Zj[i] <- Z[[j]][loc]
        i <- i + 1
      }
      else{
        loc <- loc + 1
      }
    }
    
    #############################Survival Time##################################
    RHS <- -log(runif(1))
    LHS <- cumsum(exp(beta*Zj)*diff(Lambda0TX)[1:nj])
    #cat(RHS, LHS, nj)
    
    error_found <- 0 
    tryCatch({ some_result <- if(RHS > LHS[nj]){cat('')} }, error = function(e) { error_found <<- 1 } )
    if (error_found){
      cat("debug stop here")
    }
    
    
    if( RHS > LHS[nj] ){
      SurvTime[j] <- jitter(OST[j])
      urate <- urate + 1
    }
    else
      SurvTime[j] <- TXj[min(seq(nj)[LHS>=RHS])]
  }
  
  
  urate <- urate/n
  
  
  

  
  # MySimu$X <- res
  # MySimu <- MySimu[order(MySimu$covs),]
  # 
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
  # MySimu$X[observed[i+1]:samp_size] <- jitter(MySimu$X[observed[i+1]:samp_size], factor = 2)
  return(list("res" = SurvTime,
              "urates" = urate))
}
