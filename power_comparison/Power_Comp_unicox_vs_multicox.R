#Bootstrap Hypothesis Testing with censoring
#(smoothed survival time, KM estimator for censoring, est beta)
#H0: linear
#H1: monotone
#Test statistic T = l(psi_hat) - l(beta_hat)
#Survival time is generated from exp(1)

# setwd("C:/Users/Huan/Dropbox/0My research/codes/partial_linear")
source("nlpl_Censor_ind_Test.R")
source("SurvTimeSmoothSurv.R")
source("CenTimeFix.R")


Power_unicox_vs_multicox <- function(Z, phiZ, W, beta_trt, B = 100, samp_size = 500, k,qt){
  #Z (vector) : monotone covariate. Z is sorted.
  #phiZ (vector) : the value of phi(z).
  #hazard is MONOTONE INCREASING in Z.
  
  #W (dataframe) : linear covariate. Column names cannot be 'covs'.
  #beta_trt (vector) : the coefficient for W.
  
  
  
  
  #Generate the original sample
  sig_level <- 0.05
  tf <- FALSE #timefix
  
  set.seed(k)  
  T_survival <- -log(runif(samp_size))/exp(phiZ + as.matrix(W) %*% as.matrix(beta_trt))
  maxC <- quantile(T_survival, qt)
  T_censoring <- runif(samp_size, min = 0, max = maxC)
  Status <- as.numeric(T_survival <= T_censoring)
  # censorrate <- 1 - mean(Status)
  X <- pmin(T_survival, T_censoring)
  MySimu_ori <- data.frame("stop" = X, "covs" = Z, "event" = Status)
  MySimu_ori <- cbind(MySimu_ori, W)
  W_names = colnames(W)
  
  # save(T_survival, file = paste0("survivaltime_",k,".rda"))
  
  
  
  
  ##Start from the first observed failure (after sorting by Z)
  N <- length(Z)
  W <- as.data.frame(W[order(MySimu_ori$covs), ])
  colnames(W) <- W_names
  
  MySimu_ori <- MySimu_ori[order(MySimu_ori$covs), ]
  location_start <- min(seq(N)[MySimu_ori$event > 0])
  MySimu <- MySimu_ori[location_start:N,]
  W2 <- as.data.frame(as.matrix(W)[location_start:N, ])
  colnames(W2) = W_names
  
  
  
  
  #Calculate test statistic
  
    
  #Fit Cox model with only covariate Z(under H0)
  
  cox_H0 <- coxph(Surv(stop, event) ~ covs, data = MySimu, control = coxph.control(timefix = tf))
  
  mybeta_H0 <- max(cox_H0$coefficients, 0)#monotone increasing
  #mybeta[1] <- max(mybeta[1],0)
  nlpl_H0 <- nlpl_Censor_ind(events = MySimu$event, 
                             surv_time = MySimu$stop, 
                             uni_cov = MySimu$covs,
                             zk = median(MySimu$covs),
                             est_beta = mybeta_H0)
  
  
  #Fit Cox model with all covariates(under H1)
  
  cmd <- "cox_H1 <- coxph(Surv(stop, event) ~ covs"
  
  for (v in W_names) {
    cmd <- paste(cmd, v, sep = "+")
  }
  cmd <- paste(cmd, ", data = MySimu, control = coxph.control(timefix = tf))")
  
  # cox <- coxph(Surv(stop, event) ~ covs, data = MySimu, control = coxph.control(timefix = tf), ties = "breslow")
  eval(parse(text = cmd))
  
  mybeta <- cox_H1$coefficients
  mybeta[1] <- max(mybeta[1],0)#monotone increasing
  nlpl_H1 <- nlpl_Censor_ind_PL(data = MySimu,
                               # events = MySimu$event, 
                               # surv_time = MySimu$stop,
                               # uni_cov = MySimu$covs,
                               est_beta = mybeta, 
                               zk = median(MySimu$covs),
                               W = W2)
    
    
  
  
  Tobs <- nlpl_H1 - nlpl_H0
  
  
  #Generate bootstrap sample
  Tn <- c()
  # cenrate <- 0

  
  
  for (i in 1:B) {
    
    
    
    #Generating survival time under the null
    #generating survival time
    survgen <- SurvTime(psi = mybeta_H0*MySimu$covs, data = MySimu)
    survtime_B <- survgen$res
    
    #generating censoring time
    cengen <- CenTime(X= MySimu$stop, Status = MySimu$event, M=max(X))
    centime_B <-cengen$res
    
    Status_B <- as.numeric(survtime_B <= centime_B)
    X_B <- pmin(survtime_B, centime_B)
    # cenrate <- cenrate + 1 - mean(Status_B)
    
    
    MySimu_B <- data.frame("stop" = X_B, "covs" = MySimu$covs, "event" = Status_B)
    MySimu_B <- cbind(MySimu_B, W2)
    N <- length(MySimu_B$covs)
    location_start <- min(seq(N)[MySimu_B$event > 0])
    MySimu_B <- MySimu_B[location_start:N, ]
    W_B <- as.matrix(W2)[location_start:N, ]
    
    
    
    
    #Fit Cox model with only covariate Z(under H0)
    cox_H0 <- coxph(Surv(stop, event) ~ covs, data = MySimu_B, control = coxph.control(timefix = tf))
    
    mybeta <- max(cox_H0$coefficients, 0)#monotone increasing
    
    nlpl_H0 <- nlpl_Censor_ind(events = MySimu_B$event, 
                               surv_time = MySimu_B$stop, 
                               uni_cov = MySimu_B$covs,
                               zk = median(MySimu_B$covs),
                               est_beta = mybeta)
    
    
    
    
    #Fit Cox model with all covariates(under H1)
    cmd <- "cox <- coxph(Surv(stop, event) ~ covs"
    
    for (v in W_names) {
      cmd <- paste(cmd, v, sep = "+")
    }
    cmd <- paste(cmd, ", data = MySimu_B, control = coxph.control(timefix = tf))")
    
    # cox <- coxph(Surv(stop, event) ~ covs, data = MySimu, control = coxph.control(timefix = tf), ties = "breslow")
    eval(parse(text = cmd))
    
    mybeta <- cox$coefficients
    mybeta[1] <- max(mybeta[1],0)#monotone increasing
    nlpl_H1 <- nlpl_Censor_ind_PL(data = MySimu_B,
                                  # events = MySimu$event, 
                                  # surv_time = MySimu$stop,
                                  # uni_cov = MySimu$covs,
                                  est_beta = mybeta, 
                                  zk = median(MySimu_B$covs),
                                  W = W_B)
    
    
    Tn[i] <- nlpl_H1 - nlpl_H0
    
  }
  
  T_critical <- quantile(Tn, sig_level)
  
  
  results <- list(
    "T_obs" = Tobs,
    "T_critical" = T_critical,
    "p_value" = mean(as.numeric(Tn < Tobs)))
  
  return(results)
}


