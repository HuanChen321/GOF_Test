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
  
  # Testing the coefficient for the first linear covariate, i.e. the first column in W
  
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
  
  # Anchor constraints
  W2 <- W2 - apply(W2, 2, median)[col(W2)]
  
  
  
  
  #Calculate test statistic
  
  
    
  ##Fit monotone or partial linear hazard model(under H0)
  #BE CAREFUL the hazard should be increasing in Z
  cmd = "iso_H0 <- isoph(Surv(stop, event) ~ iso(covs) "
  for (v in  W_names[!W_names %in% "V1"]) {
    cmd <- paste(cmd, v, sep = "+")
  }
  cmd <- paste(cmd, ", data = MySimu)")
  eval(parse(text = cmd))
  
  if(is.na(iso_H0$beta)){
    mono_ZW0 <- iso_H0$iso.cov$psi.hat
  }else{
    mono_ZW0 <- iso_H0$iso.cov$psi.hat + as.matrix(subset(W2, select = -c(V1)))%*%iso_H0$beta$est
  }
  
  nlpl_H0 <- nlpl_Censor_ind(events = MySimu$event, 
                             surv_time = MySimu$stop, 
                             mono_f = mono_ZW0)
  
  
  ##Fit partial linear model(under H1)
  #BE CAREFUL the hazard should be increasing in Z
  cmd = "iso <- isoph(Surv(stop, event) ~ iso(covs) "
  for (v in W_names) {
    cmd <- paste(cmd, v, sep = "+")
  }
  cmd <- paste(cmd, ", data = MySimu)")
  eval(parse(text = cmd))
  
  mono_ZW <- iso$iso.cov$psi.hat + as.matrix(W2)%*%iso$beta$est
 
  nlpl_H1 <- nlpl_Censor_ind(events = MySimu$event, 
                             surv_time = MySimu$stop, 
                             mono_f = mono_ZW)

  
  
  Tobs <- nlpl_H1 - nlpl_H0
  
  
  #Generate bootstrap sample
  Tn <- c()
  # cenrate <- 0

  
  
  for (i in 1:B) {
    
    
    #Generating survival time under the null
    #generating survival time
    survgen <- SurvTime(psi = mono_ZW0, data = MySimu)
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
    
    #do we need to apply anchor constraints again to W_B?
    # Anchor constraints
    W_B <- W_B - apply(W_B, 2, median)[col(W_B)]
    
    
    
      
    ##Fit monotoe hazard model(under H0)
    #BE CAREFUL the hazard should be increasing in Z
    
    cmd = "iso_H0 <- isoph(Surv(stop, event) ~ iso(covs) "
    for (v in  W_names[!W_names %in% "V1"]) {
      cmd <- paste(cmd, v, sep = "+")
    }
    cmd <- paste(cmd, ", data = MySimu_B)")
    eval(parse(text = cmd))
    
    if(is.na(iso_H0$beta)){
      mono_ZW <- iso_H0$iso.cov$psi.hat
    }else{
      mono_ZW <- iso_H0$iso.cov$psi.hat + as.matrix(subset(W_B, select = -c(V1)))%*%iso_H0$beta$est
    }
    
    nlpl_H0 <- nlpl_Censor_ind(events = MySimu_B$event, 
                               surv_time = MySimu_B$stop, 
                               mono_f = mono_ZW)
    
    
    ##Fit partial linear model(under H1)
    #BE CAREFUL the hazard should be increasing in Z
    cmd = "iso <- isoph(Surv(stop, event) ~ iso(covs) "
    for (v in W_names) {
      cmd <- paste(cmd, v, sep = "+")
    }
    cmd <- paste(cmd, ", data = MySimu_B)")
    eval(parse(text = cmd))
    
    mono_ZW <- iso$iso.cov$psi.hat + as.matrix(W_B)%*%iso$beta$est
    
    nlpl_H1 <- nlpl_Censor_ind(events = MySimu_B$event, 
                               surv_time = MySimu_B$stop, 
                               mono_f = mono_ZW)
    
    
  
    
    Tn[i] <- nlpl_H1 - nlpl_H0
    
  }
  
  T_critical <- quantile(Tn, sig_level)
  
  
  results <- list(
    "T_obs" = Tobs,
    "T_critical" = T_critical,
    "p_value" = mean(as.numeric(Tn < Tobs)))
  
  return(results)
}


