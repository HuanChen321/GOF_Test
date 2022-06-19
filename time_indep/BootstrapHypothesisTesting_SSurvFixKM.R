# source("Bootstrap_c_exp.R")
# source("Bootstrap_c_gp.R")
# source("Bootstrap_c_wb.R")
#source("Bootstrap_exp_SSurvFixKM_moutofn.R")
source("Bootstrap_exp_SSurvFixKM.R")


B_rep = 500
sampsize_ori = 1000
B_size = 500

library(foreach)
library(parallel)
library(doSNOW)
# library(doParallel)

# cl <- makeCluster(25)
# registerDoParallel(cl)
# 
# start_time = Sys.time()

cl <- makeCluster(24)
registerDoSNOW(cl)

pb <- txtProgressBar(max = B_rep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time = Sys.time()

BHT_CB <- foreach(k = 1:B_rep, 
                         .multicombine=TRUE, .packages=c('isoSurv', 'survival', 'pracma'),
                         .options.snow=opts) %dopar% {
  
   source('myisoph_initial.R')
   assignInNamespace("isoph.initial", isoph.initial, pos = "package:isoSurv")
   
  set.seed(k)
  Z <- sort(runif(sampsize_ori, min = 0.1, max = 2))
  hzrate0 = 0*Z
  hzrate1 = Z
  hzrate5 = 5*Z
  hzrate_m1 = sqrt(Z)
  hzrate_m2 = Z^2
  hzrate_m3 = exp(Z)
  hzrate_m4 = log(Z)
  
  # library(isoSurv)
  # library(survival)
  # library(pracma)
  
  test0_c_simple <- Bootstrap_KMBS(Z, hzrate0, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.96)
  test1_c_simple <- Bootstrap_KMBS(Z, hzrate1, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.93)
  test5_c_simple <- Bootstrap_KMBS(Z, hzrate5, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.8)
  testm1_c_simple <- Bootstrap_KMBS(Z, hzrate_m1, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.95)
  testm2_c_simple <- Bootstrap_KMBS(Z, hzrate_m2, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.9)
  testm3_c_simple <- Bootstrap_KMBS(Z, hzrate_m3, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.87)
  testm4_c_simple <- Bootstrap_KMBS(Z, hzrate_m4, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.91)
  
  #test0_nc_simple <- Bootstrap_BS(Z, hzrate0, B = B_size, samp_size = sampsize_ori)
  #test1_nc_simple <- Bootstrap_BS(Z, hzrate1, B = B_size, samp_size = sampsize_ori)
  #test5_nc_simple <- Bootstrap_BS(Z, hzrate5, B = B_size, samp_size = sampsize_ori)
  #testm1_nc_simple <- Bootstrap_BS(Z, hzrate_m1, B = B_size, samp_size = sampsize_ori)
  #testm2_nc_simple <- Bootstrap_BS(Z, hzrate_m2, B = B_size, samp_size = sampsize_ori)
  #testm3_nc_simple <- Bootstrap_BS(Z, hzrate_m3, B = B_size, samp_size = sampsize_ori)
  #testm4_nc_simple <- Bootstrap_BS(Z, hzrate_m4, B = B_size, samp_size = sampsize_ori)
  
   
   return(list(#"test0_nc_simple" = test0_nc_simple,
   #"test1_nc_simple" = test1_nc_simple,
   #           "test5_nc_simple" = test5_nc_simple,
   #           "testm1_nc_simple" = testm1_nc_simple,
   #           "testm2_nc_simple" = testm2_nc_simple,
   #           "testm3_nc_simple" = testm3_nc_simple,
   #           "testm4_nc_simple" = testm4_nc_simple,
              
              "test0_c_simple" = test0_c_simple,
              "test1_c_simple" = test1_c_simple,
              "test5_c_simple" = test5_c_simple,
               "testm1_c_simple" = testm1_c_simple,
               "testm2_c_simple" = testm2_c_simple,
               "testm3_c_simple" = testm3_c_simple,
               "testm4_c_simple" = testm4_c_simple
              ))
}



close(pb)
stopCluster(cl)

end_time = Sys.time()
cat('Time used:\n')
print(end_time - start_time)


# stopImplicitCluster()
# stopCluster(cl)
# 
# end_time = Sys.time()
# cat('Time used:\n')
# print(end_time - start_time)

save(BHT_CB, file = "BHT02SSurvFixKMmaxX1000s500b.rda")
