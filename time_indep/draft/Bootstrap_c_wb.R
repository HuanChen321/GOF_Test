#Bootstrap Hypothesis Testing with no censoring
#H0: linear
#H1: monotone
#Test statistic T = l(psi_hat) - l(beta_hat)

library(isoSurv)
library(survival)
source("D:/0My research/031821/nlpl_Censor_ind.R")
source("D:/0My research/031821/nlpl_noCensor_ind.R")


Bootstrap_c_wb <- function(Z, hzrate, nu, lambda, B = 100, samp_size = 500){
  #scale = lambda
  #shift= nu
  #bootstrap sample size = samp_size
  
#Generate the original sample
sig_level <- 0.05
tf <- FALSE #timefix

T_survival <- -log(runif(samp_size))/exp(hzrate)
T_censoring <- runif(samp_size, min = 0, max = quantile(T_survival, 0.96))
Status_c <- as.numeric(T_survival <= T_censoring)
MySimu <- data.frame("stop" = pmin(T_survival, T_censoring), "covs" = Z, "event" = Status_c)

##Start from the first observed failure
N <- length(Z)
location_start <- min(seq(N)[MySimu$event > 0])
MySimu_c <- MySimu[location_start:N,]

#Estimate beta_hat
iso_c <- isoph(Surv(stop, event) ~ iso(covs, shape = "inc"), data = MySimu_c)
psi_hat <- iso_c$iso.cov$psi.hat
nlpl_iso0 <- nlpl_Censor_ind(events = MySimu_c$event, 
                             surv_time = MySimu_c$stop, 
                             mono_f = psi_hat, 
                             uni_cov = MySimu_c$covs, 
                             zk = iso_c$Zk)

cox_c <- coxph(Surv(stop, event) ~ covs, data = MySimu_c, control = coxph.control(timefix = tf), ties = "breslow")
mybeta_c <- max(cox_c$coefficients,0)
nlpl_cox0 <- nlpl_Censor_ind(events = MySimu_c$event, 
                             surv_time = MySimu_c$stop, 
                             est_beta = mybeta_c, 
                             uni_cov = MySimu_c$covs,
                             zk = iso_c$Zk)

Tobs <- nlpl_iso0 - nlpl_cox0


#Generate bootstrap sample
Tn <- c()
for (i in 1:B) {
  covs_B <- sort(runif(samp_size, min = 0.1, max = 2))
  hzrate_B <- mybeta_c*covs_B
  survtime_B <- -nthroot((log(runif(samp_size))/(lambda*exp(hzrate_B))),nu)
  Status_B <- rep(1, samp_size)
  
  isofit <- isoph(Surv(survtime_B, Status_B) ~ iso(covs_B, shape = "inc"))
  nlpl_iso <- nlpl_noCensor_ind(mono_f = isofit$iso.cov$psi.hat, 
                                uni_cov = covs_B, 
                                surv_time = survtime_B, 
                                zk = isofit$Zk)
  
  coxfit <- coxph(Surv(survtime_B, Status_B) ~ covs_B, control = coxph.control(timefix = tf), ties = "breslow")
  beta_hat <- max(0, coxfit$coefficients)
  nlpl_cox <- nlpl_noCensor_ind(est_beta = beta_hat, 
                                uni_cov = covs_B, 
                                surv_time = survtime_B, 
                                zk = isofit$Zk)
  
  Tn[i] <- nlpl_iso - nlpl_cox

}

T_critical <- quantile(Tn, sig_level)

results <- list(
  "T_obs" = Tobs,
  "beta_hat" = mybeta_c,
  "T_critical" = T_critical,
  "p_value" = mean(as.numeric(Tn < Tobs)),
  "reject" = as.numeric(Tobs < T_critical))

return(results)
}