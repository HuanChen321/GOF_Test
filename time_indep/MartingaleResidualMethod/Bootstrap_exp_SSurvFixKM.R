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

fit <- cox.aalen(Surv(stop, event)~prop(covs), data = MySimu, residuals = 1, rate.sim=0, n.sim=100)
mt_resids = cum.residuals(fit, data=MySimu, cum.resid = 1, n.sim = 100)




results <- list("p_value" = mt_resids$pval.test)

return(results)
}


