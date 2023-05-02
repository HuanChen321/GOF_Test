#Bootstrap Hypothesis Testing with censoring
#(smoothed survival time, KM estimator for censoring, est beta)
#H0: linear
#H1: monotone
#Test statistic T = l(psi_hat) - l(beta_hat)
#Survival time is generated from exp(1)

#library(isoSurv)
#library(survival)
source("nlpl_Censor_ind.R")
#source("SurvTime.R")
#source("CenTime.R")
source("SurvTimeSmoothSurv.R")
source("CenTimeFix.R")


Bootstrap_KMBS <- function(Z, hzrate, lambda, B = 100, samp_size = 500, k,qt){
  
#Generate the original sample
sig_level <- 0.05
tf <- FALSE #timefix

set.seed(k)  
T_survival <- -log(runif(samp_size))/exp(hzrate)
maxC <- quantile(T_survival, qt)
T_censoring <- runif(samp_size, min = 0, max = maxC)
Status <- as.numeric(T_survival <= T_censoring)
censorrate <- 1 - mean(Status)
X <- pmin(T_survival, T_censoring)
MySimu_ori <- data.frame("stop" = X, "covs" = Z, "event" = Status)

##Start from the first(in covariate) observed failure
N <- length(Z)
MySimu_ori <- MySimu_ori[order(MySimu_ori$covs),]
location_start <- min(seq(N)[MySimu_ori$event > 0])
MySimu <- MySimu_ori[location_start:N,]

#Estimate beta and phi
iso <- isoph(Surv(stop, event) ~ iso(covs, shape = "inc"), data = MySimu)
psi_hat <- iso$iso.cov$psi.hat
nlpl_iso0 <- nlpl_Censor_ind(events = MySimu$event, 
                             surv_time = MySimu$stop, 
                             mono_f = psi_hat, 
                             uni_cov = MySimu$covs, 
                             zk = iso$Zk)

cox <- coxph(Surv(stop, event) ~ covs, data = MySimu, control = coxph.control(timefix = tf), ties = "breslow")
mybeta <- max(cox$coefficients,0)
nlpl_cox0 <- nlpl_Censor_ind(events = MySimu$event, 
                             surv_time = MySimu$stop, 
                             est_beta = mybeta, 
                             uni_cov = MySimu$covs,
                             zk = iso$Zk)

Tobs <- nlpl_iso0 - nlpl_cox0


#Generate bootstrap sample
Tn <- c()
cenrate <- 0
urateS <- 0
urateC <- 0

for (i in 1:B) {
  covs_B <- Z
  hzrate_B <- mybeta*covs_B
  
  survgen <- SurvTime(beta = mybeta, MySimu_ori)
  survtime_B <- survgen$res
  urateS <- urateS + survgen$urate
  
  cengen <- CenTime(X= MySimu_ori$stop, Status = MySimu_ori$event, M=max(X))
  centime_B <-cengen$res
  urateC <- urateC + cengen$urate
  
  Status_B <- as.numeric(survtime_B <= centime_B)
  X_B <- pmin(survtime_B, centime_B)
  cenrate <- cenrate + 1 - mean(Status_B)
  
  MySimu_B <- data.frame("stop" = X_B, "covs" = covs_B, "event" = Status_B)
  MySimu_B <- MySimu_B[order(MySimu_B$covs),]#sort by covariate to drop the records with covariate smaller than the covariate of the first observed failure
  location_start <- min(seq(N)[MySimu_B$event > 0])
  MySimu_B <- MySimu_B[location_start:N,]
  
  isofit <- isoph(Surv(stop, event) ~ iso(covs, shape = "inc"), data = MySimu_B)
  nlpl_iso <- nlpl_Censor_ind(events = MySimu_B$event,
                              mono_f = isofit$iso.cov$psi.hat, 
                              uni_cov = MySimu_B$covs, 
                              surv_time = MySimu_B$stop, 
                              zk = isofit$Zk)
  
  coxfit <- coxph(Surv(stop, event) ~ covs, data = MySimu_B, control = coxph.control(timefix = tf), ties = "breslow")
  beta_hat <- max(0, coxfit$coefficients)
  nlpl_cox <- nlpl_Censor_ind(events = MySimu_B$event,
                              est_beta = beta_hat, 
                              uni_cov = MySimu_B$covs, 
                              surv_time = MySimu_B$stop, 
                              zk = isofit$Zk)
  
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


