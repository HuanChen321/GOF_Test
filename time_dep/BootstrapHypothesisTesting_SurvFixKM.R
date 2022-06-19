
source("Bootstrap_exp_SurvFixKM.R")
B_rep = 100#500


library(foreach)
library(parallel)
library(doSNOW)


cl <- makeCluster(2)
registerDoSNOW(cl)

pb <- txtProgressBar(max = B_rep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time = Sys.time()

BHT_CB <- foreach(k = 1:B_rep, 
                         .multicombine=TRUE, .packages=c('Iso', 'survival', 'pracma'),
                         .options.snow=opts) %dopar% {
  
   
  

                           
  sampsize_ori = 100 #500
  num_intervals = 20
  B_size = 100 #100
  
  set.seed(k)                                     
  Z <- matrix(runif(sampsize_ori*num_intervals, min = 0, max = 2), nrow = sampsize_ori)
  #hzrate0 = 0*Z
  hzrate1 = Z
  hzrate5 = 5*Z
  hzrate_m1 = sqrt(Z)
  hzrate_m2 = Z^2
  hzrate_m3 = exp(Z)
  hzrate_m4 = log(Z)
  
  
  # test0_c <- Bootstrap_KMBS(Z, hzrate0, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.96)
  test1_c <- Bootstrap_KMBS(Z, hzrate1, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.93, num_intervals  = num_intervals)
  test5_c <- Bootstrap_KMBS(Z, hzrate5, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.8, num_intervals  = num_intervals)
  testm1_c <- Bootstrap_KMBS(Z, hzrate_m1, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.95, num_intervals  = num_intervals)
  testm2_c <- Bootstrap_KMBS(Z, hzrate_m2, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.9, num_intervals  = num_intervals)
  testm3_c <- Bootstrap_KMBS(Z, hzrate_m3, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.87, num_intervals  = num_intervals)
  testm4_c <- Bootstrap_KMBS(Z, hzrate_m4, B = B_size, samp_size = sampsize_ori, k = k, qt = 0.91, num_intervals  = num_intervals)

  
   
   return(list(#"test0_c" = test0_c
               "test1_c" = test1_c,
               "test5_c" = test5_c,
               "testm1_c" = testm1_c,
               "testm2_c" = testm2_c,
               "testm3_c" = testm3_c,
               "testm4_c" = testm4_c
              ))
}



close(pb)
stopCluster(cl)

end_time = Sys.time()
cat('Time used:\n')
print(end_time - start_time)



save(BHT_CB, file = "BHT02TD_SurvFixKM100s100bt100b.rda")
