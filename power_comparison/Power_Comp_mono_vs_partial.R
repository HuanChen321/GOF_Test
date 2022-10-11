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


Power_mono_vs_partial <- function(Z, phiZ, W, beta_trt, B = 100, samp_size = 500, k,qt){
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
  
    
  ##Fit monotoe hazard model(under H0)
  #BE CAREFUL the hazard should be increasing in Z
  iso_H0 <- isoph(Surv(stop, event) ~ iso(covs), data = MySimu)
  
  psi_hat_H0 <- iso_H0$iso.cov$psi.hat
  # conv <- iso_H0$conv
  # steps <- iso_H0$iter
  
  nlpl_H0 <- nlpl_Censor_ind(events = MySimu$event, 
                             surv_time = MySimu$stop, 
                             uni_cov = MySimu$covs, 
                             zk = iso$Zk,
                             mono_f = psi_hat_H0)
  
  
  ##Fit partial linear model(under H1)
  #BE CAREFUL the hazard should be increasing in Z
  cmd = "iso <- isoph(Surv(stop, event) ~ iso(covs) "
  
  for (v in W_names) {
    cmd <- paste(cmd, v, sep = "+")
  }
  cmd <- paste(cmd, ", data = MySimu)")
  eval(parse(text = cmd))
  
  psi_hat <- iso$iso.cov$psi.hat
  # save(psi_hat, file = paste0("psi_hat_",k,".rda"))
  nlpl_H1 <- nlpl_Censor_ind_PL(data = MySimu,
                               # events = MySimu$event, 
                               # surv_time = MySimu$stop, 
                               # uni_cov = MySimu$covs, 
                               zk = iso$Zk,
                               mono_f = psi_hat,
                               est_beta = iso$beta[,1],
                               W = W2)

  
  
  Tobs <- nlpl_H1 - nlpl_H0
  
  
  #Generate bootstrap sample
  Tn <- c()
  # cenrate <- 0

  
  
  for (i in 1:B) {
    
    
    
    #Generating survival time under the null
    #generating survival time
    survgen <- SurvTime(psi = psi_hat_H0, data = MySimu)
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
    
    
    
      
    ##Fit monotoe hazard model(under H0)
    #BE CAREFUL the hazard should be increasing in Z
    iso_H0 <- isoph(Surv(stop, event) ~ iso(covs), data = MySimu_B)
    
    psi_hat <- iso_H0$iso.cov$psi.hat
    
    nlpl_H0 <- nlpl_Censor_ind(events = MySimu_B$event, 
                               surv_time = MySimu_B$stop, 
                               uni_cov = MySimu_B$covs, 
                               zk = iso$Zk,
                               mono_f = psi_hat)
    
    
    
    
    ##Fit partial linear model(under H1)
    #BE CAREFUL the hazard should be increasing in Z
    cmd = "iso <- isoph(Surv(stop, event) ~ iso(covs) "
    
    for (v in W_names) {
      cmd <- paste(cmd, v, sep = "+")
    }
    cmd <- paste(cmd, ", data = MySimu_B)")
    eval(parse(text = cmd))
    
    psi_hat <- iso$iso.cov$psi.hat
    # save(psi_hat, file = paste0("psi_hat_",k,".rda"))
    nlpl_H1 <- nlpl_Censor_ind_PL(data = MySimu_B,
                                  zk = iso$Zk,
                                  mono_f = psi_hat,
                                  est_beta = iso$beta[,1],
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


