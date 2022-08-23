# source("Bootstrap_c_exp.R")
# source("Bootstrap_c_gp.R")
# source("Bootstrap_c_wb.R")
#source("Bootstrap_exp_SSurvFixKM_moutofn.R")
source("Bootstrap_exp_SSurvFixKM.R")


B_rep = 1
sampsize_ori = 100
B_size = 100

library(foreach)
library(parallel)
library(doSNOW)
# library(doParallel)

# cl <- makeCluster(25)
# registerDoParallel(cl)
# 
# start_time = Sys.time()

# cl <- makeCluster(1)
# registerDoSNOW(cl)
# 
# pb <- txtProgressBar(max = B_rep, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)

start_time = Sys.time()

# BHT_CB <- foreach(k = 1:B_rep, 
#                          .multicombine=TRUE, .packages=c('isoSurv', 'survival', 'pracma'),
#                          .options.snow=opts) %dopar% {

test1_c <- list()
testm1_c <- list()

for (k in 1:B_rep) {
  
  
   source('myisoph_initial.R')
   assignInNamespace("isoph.initial", isoph.initial, pos = "package:isoSurv")
   
  set.seed(k)
  Z <- sort(runif(sampsize_ori, min = 0.1, max = 2))
  W <- data.frame(
    W1 = rnorm(sampsize_ori, mean = 0, sd = 1)
  )
  
  hzrate1 = Z
  hzrate_m1 = Z^2

  
  library(isoSurv)
  library(survival)
  library(pracma)
  
  test1_c[[k]] <- Bootstrap_KMBS(Z, phiZ = hzrate1, W = W, beta_trt = c(1), B = B_size, samp_size = sampsize_ori, k = k, qt = 0.93)

  testm1_c[[k]] <- Bootstrap_KMBS(Z, phiZ = hzrate_m1, W = W, beta_trt = c(1), B = B_size, samp_size = sampsize_ori, k = k, qt = 0.9)

  
   
   # return(list("test1_c" = test1_c,
   #            
   #             "testm1_c" = testm1_c
   #            ))
}



# close(pb)
# stopCluster(cl)

end_time = Sys.time()
cat('Time used:\n')
print(end_time - start_time)


# stopImplicitCluster()
# stopCluster(cl)
# 
# end_time = Sys.time()
# cat('Time used:\n')
# print(end_time - start_time)

# save(BHT_CB, file = "BHT02SSurvFixKMmaxX100s2b.rda")
