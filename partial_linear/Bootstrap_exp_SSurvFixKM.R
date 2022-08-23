#Bootstrap Hypothesis Testing with censoring
#(smoothed survival time, KM estimator for censoring, est beta)
#H0: linear
#H1: monotone
#Test statistic T = l(psi_hat) - l(beta_hat)
#Survival time is generated from exp(1)

source("nlpl_Censor_ind.R")
source("SurvTimeSmoothSurv.R")
source("CenTimeFix.R")


Bootstrap_KMBS <- function(Z, phiZ, W, beta_trt, B = 100, samp_size = 500, k,qt){
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
censorrate <- 1 - mean(Status)
X <- pmin(T_survival, T_censoring)
MySimu_ori <- data.frame("stop" = X, "covs" = Z, "event" = Status)
MySimu_ori <- cbind(MySimu_ori, W)
W_names = colnames(W)


##Start from the first observed failure (after sorting by Z)
N <- length(Z)
W <- as.data.frame(W[order(MySimu_ori$covs), ])
colnames(W) <- W_names
 
MySimu_ori <- MySimu_ori[order(MySimu_ori$covs), ]
location_start <- min(seq(N)[MySimu_ori$event > 0])
MySimu <- MySimu_ori[location_start:N,]
W2 <- as.matrix(W)[location_start:N, ]


##Fit partial linear model
#BE CAREFUL the hazard should be increasing in Z
cmd = "iso <- isoph(Surv(stop, event) ~ iso(covs) "

for (v in W_names) {
  cmd <- paste(cmd, v, sep = "+")
}
cmd <- paste(cmd, ", data = MySimu)")

eval(parse(text = cmd))

psi_hat <- iso$iso.cov$psi.hat
nlpl_iso0 <- nlpl_Censor_ind(data = MySimu,
                             # events = MySimu$event, 
                             # surv_time = MySimu$stop, 
                             # uni_cov = MySimu$covs, 
                             zk = iso$Zk,
                             mono_f = psi_hat,
                             est_beta = iso$beta[,1],
                             W = W2)


##Fit Cox model
cmd <- "cox <- coxph(Surv(stop, event) ~ covs"

for (v in W_names) {
  cmd <- paste(cmd, v, sep = "+")
}
cmd <- paste(cmd, ", data = MySimu, control = coxph.control(timefix = tf))")

# cox <- coxph(Surv(stop, event) ~ covs, data = MySimu, control = coxph.control(timefix = tf), ties = "breslow")
eval(parse(text = cmd))

mybeta <- cox$coefficients
mybeta[1] <- max(mybeta[1],0)#monotone increasing
nlpl_cox0 <- nlpl_Censor_ind(data = MySimu,
                             # events = MySimu$event, 
                             # surv_time = MySimu$stop,
                             # uni_cov = MySimu$covs,
                             est_beta = mybeta, 
                             zk = iso$Zk,
                             W = W2)

Tobs <- nlpl_iso0 - nlpl_cox0


#Generate bootstrap sample
Tn <- c()
cenrate <- 0
urateS <- 0
urateC <- 0

#sort Z since MySimu_ori is sorted by Z

Z_B <- sort(Z)

for (i in 1:B) {
 
  hzrate_B <- mybeta[1]*Z_B + as.matrix(W)%*%as.matrix(mybeta[2:length(mybeta)])
  
  #generating survival time
  survgen <- SurvTime(beta = mybeta, data = MySimu_ori, linear_covs = W)
  survtime_B <- survgen$res
  urateS <- urateS + survgen$urate
  
  #generating censoring time
  cengen <- CenTime(X= MySimu_ori$stop, Status = MySimu_ori$event, M=max(X))
  centime_B <-cengen$res
  urateC <- urateC + cengen$urate
  
  Status_B <- as.numeric(survtime_B <= centime_B)
  X_B <- pmin(survtime_B, centime_B)
  cenrate <- cenrate + 1 - mean(Status_B)
  
  
  MySimu_B <- data.frame("stop" = X_B, "covs" = Z_B, "event" = Status_B)
  MySimu_B <- cbind(MySimu_B, W)
  # MySimu_B <- MySimu_B[order(MySimu_B$covs),]#sort by covariate to drop the records with covariate smaller than the covariate of the first observed failure
  location_start <- min(seq(N)[MySimu_B$event > 0])
  MySimu_B <- MySimu_B[location_start:N, ]
  W_B <- as.matrix(W)[location_start:N, ]
  
  
  #fit partial linear model
  cmd = "isofit <- isoph(Surv(stop, event) ~ iso(covs) "
  
  for (v in W_names) {
    cmd <- paste(cmd, v, sep = "+")
  }
  cmd <- paste(cmd, ", data = MySimu_B)")
  
  eval(parse(text = cmd))
  # isofit <- isoph(Surv(stop, event) ~ iso(covs, shape = "inc"), data = MySimu_B)
  
  nlpl_iso <- nlpl_Censor_ind(data = MySimu_B,
                              # events = MySimu$event, 
                              # surv_time = MySimu$stop, 
                              # uni_cov = MySimu$covs, 
                              zk = isofit$Zk,
                              mono_f = isofit$iso.cov$psi.hat,
                              est_beta = iso$beta[,1],
                              W = W_B)
  
  
  #fit cox model
  cmd <- "coxfit <- coxph(Surv(stop, event) ~ covs"
  
  for (v in W_names) {
    cmd <- paste(cmd, v, sep = "+")
  }
  cmd <- paste(cmd, ", data = MySimu_B, control = coxph.control(timefix = tf))")
  
  # coxfit <- coxph(Surv(stop, event) ~ covs, data = MySimu_B, control = coxph.control(timefix = tf), ties = "breslow")
  eval(parse(text = cmd))
  
  beta_hat <- coxfit$coefficients
  beta_hat[1] <- max(beta_hat[1],0)#monotone increasing
  nlpl_cox <- nlpl_Censor_ind(data = MySimu_B,
                               # events = MySimu$event, 
                               # surv_time = MySimu$stop,
                               # uni_cov = MySimu$covs,
                               est_beta = beta_hat, 
                               zk = iso$Zk,
                               W = W_B)
  
  Tn[i] <- nlpl_iso - nlpl_cox
}

T_critical <- quantile(Tn, sig_level)


results <- list(
  "T_obs" = Tobs,
  "T_n" = Tn,
  "beta_hat" = mybeta,
  "T_critical" = T_critical,
  #"censorrate_B" = cenrate/B,
  #"censorrate_data" = censorrate,
  #"urateS" = urateS/B,
  #"urateC" = urateC/B,
  "p_value" = mean(as.numeric(Tn < Tobs)))

return(results)
}


