#Testing significance of the linear covariate
#Want to compare the power when using partial linear model and cox model
#i.e when the linear covariate is significant, i.e. beta for W is not zero, which test will have higher power

source("Power_Comp_mono_vs_partial.R")
source("Power_Comp_unicox_vs_multicox.R")

B_rep = 100
sampsize_ori = 100
B_size = 100

library(foreach)
library(parallel)
library(doSNOW)

cl <- makeCluster(24)
registerDoSNOW(cl)
#
pb <- txtProgressBar(max = B_rep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time = Sys.time()

TSLC <- foreach(k = 1:B_rep,
                  .multicombine=TRUE, .packages=c('isoSurv', 'survival', 'pracma'),
                  .options.snow=opts) %dopar% {
                    
                    
                    source('myisoph_initial.R')
                    assignInNamespace("isoph.initial", isoph.initial, pos = "package:isoSurv")
                    
                    #generating sample 
                    set.seed(k)
                    Z <- sort(runif(sampsize_ori, min = 0.1, max = 2))
                    W <- data.frame(
                      W1 = runif(sampsize_ori, min = 0, max = 1)
                    )
                    
                    hzrate1 = Z
                    hzrate_m1 = Z^2
                    
                    # library(isoSurv)
                    # library(survival)
                    # library(pracma)
                    
                    
                    #when true model for hazard is Cox model
                    #Univariate Cox model vs multivariate Cox model
                    test_cox_Tcox <- Power_unicox_vs_multicox(Z, phiZ = hzrate1, W = W, beta_trt = c(0), B = B_size, samp_size = sampsize_ori, k = k, qt = 0.93)#0.93
                    #Monotone hazard model vs partial linear model
                    test_pl_Tcox <- Power_mono_vs_partial(Z, phiZ = hzrate1, W = W, beta_trt = c(0), B = B_size, samp_size = sampsize_ori, k = k, qt = 0.93)#0.93
                    
                    #when true model for hazard is partial linear model
                    #Univariate Cox model vs multivariate Cox model
                    test_cox_Tpl <- Power_unicox_vs_multicox(Z, phiZ = hzrate_m1, W = W, beta_trt = c(0), B = B_size, samp_size = sampsize_ori, k = k, qt = 0.9)#0.9
                    #Monotone hazard model vs partial linear model
                    test_pl_Tpl <- Power_mono_vs_partial(Z, phiZ = hzrate_m1, W = W, beta_trt = c(0), B = B_size, samp_size = sampsize_ori, k = k, qt = 0.9)#0.9
                    
                    
                    return(list("test_cox_Tcox" = test_cox_Tcox,
                                "test_cox_Tpl" = test_cox_Tpl,
                                
                                "test_pl_Tcox" = test_pl_Tcox,
                                "test_pl_Tpl" = test_pl_Tpl))
                  }


close(pb)
stopCluster(cl)

end_time = Sys.time()
cat('Time used:\n')
print(end_time - start_time)

save(TSLC, file = "SizeComp100s100iter100bt.rda")
