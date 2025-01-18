source("Bootstrap_c_exp.R")
source("Bootstrap_c_gp.R")
source("Bootstrap_c_wb.R")



B_rep = 500
sampsize_ori = 500
B_size = 100

library(foreach)
library(parallel)
library(doParallel)

cl <- makeCluster(30)
registerDoParallel(cl)

start_time = Sys.time()

BHT_cnc <- foreach(k = 1:B_rep) %dopar% {
  
  set.seed(k)
  Z <- sort(runif(sampsize_ori, min = 0.1, max = 2))
  hzrate0 = 0*Z
  hzrate1 = Z
  hzrate5 = 5*Z
  hzrate_m1 = sqrt(Z)
  hzrate_m2 = Z^2
  
  library(isoSurv)
  library(survival)
  library(pracma)
  
   test0_nc_exp1 <- Bootstrap_c_exp(Z, lambda = 1, hzrate0, B = B_size, samp_size = sampsize_ori)
   test1_nc_exp1 <- Bootstrap_c_exp(Z, lambda = 1, hzrate1, B = B_size, samp_size = sampsize_ori)
   test5_nc_exp1 <- Bootstrap_c_exp(Z, lambda = 1, hzrate5, B = B_size, samp_size = sampsize_ori)
   testm1_nc_exp1 <- Bootstrap_c_exp(Z, lambda = 1, hzrate_m1, B = B_size, samp_size = sampsize_ori)
   testm2_nc_exp1 <- Bootstrap_c_exp(Z, lambda = 1, hzrate_m2, B = B_size, samp_size = sampsize_ori)
   
   test0_nc_exp05 <- Bootstrap_c_exp(Z, lambda = 0.5, hzrate0, B = B_size, samp_size = sampsize_ori)
   test1_nc_exp05 <- Bootstrap_c_exp(Z, lambda = 0.5, hzrate1, B = B_size, samp_size = sampsize_ori)
   test5_nc_exp05 <- Bootstrap_c_exp(Z, lambda = 0.5, hzrate5, B = B_size, samp_size = sampsize_ori)
   testm1_nc_exp05 <- Bootstrap_c_exp(Z, lambda = 0.5, hzrate_m1, B = B_size, samp_size = sampsize_ori)
   testm2_nc_exp05 <- Bootstrap_c_exp(Z, lambda = 0.5, hzrate_m2, B = B_size, samp_size = sampsize_ori)

   test0_nc_gp011 <- Bootstrap_c_gp(Z, hzrate0, alpha = 0.1, lambda = 1, B = B_size, samp_size = sampsize_ori)
   test1_nc_gp011 <- Bootstrap_c_gp(Z, hzrate1, alpha = 0.1, lambda = 1, B = B_size, samp_size = sampsize_ori)
   test5_nc_gp011 <- Bootstrap_c_gp(Z, hzrate5, alpha = 0.1, lambda = 1, B = B_size, samp_size = sampsize_ori)
   testm1_nc_gp011 <- Bootstrap_c_gp(Z, hzrate_m1, alpha = 0.1, lambda = 1, B = B_size, samp_size = sampsize_ori)
   testm2_nc_gp011 <- Bootstrap_c_gp(Z, hzrate_m2, alpha = 0.1, lambda = 1, B = B_size, samp_size = sampsize_ori)
   
   test0_nc_gp12 <- Bootstrap_c_gp(Z, hzrate0, alpha = 1, lambda = 2, B = B_size, samp_size = sampsize_ori)
   test1_nc_gp12 <- Bootstrap_c_gp(Z, hzrate1, alpha = 1, lambda = 2, B = B_size, samp_size = sampsize_ori)
   test5_nc_gp12 <- Bootstrap_c_gp(Z, hzrate5, alpha = 1, lambda = 2, B = B_size, samp_size = sampsize_ori)
   testm1_nc_gp12 <- Bootstrap_c_gp(Z, hzrate_m1, alpha = 1, lambda = 2, B = B_size, samp_size = sampsize_ori)
   testm2_nc_gp12 <- Bootstrap_c_gp(Z, hzrate_m2, alpha = 1, lambda = 2, B = B_size, samp_size = sampsize_ori)
   
   test0_nc_wb51 <- Bootstrap_c_wb(Z, hzrate0, nu = 5, lambda = 1, B = B_size, samp_size = sampsize_ori)
   test1_nc_wb51 <- Bootstrap_c_wb(Z, hzrate1, nu = 5, lambda = 1, B = B_size, samp_size = sampsize_ori)
   test5_nc_wb51 <- Bootstrap_c_wb(Z, hzrate5, nu = 5, lambda = 1, B = B_size, samp_size = sampsize_ori)
   testm1_nc_wb51 <- Bootstrap_c_wb(Z, hzrate_m1, nu = 5, lambda = 1, B = B_size, samp_size = sampsize_ori)
   testm2_nc_wb51 <- Bootstrap_c_wb(Z, hzrate_m2, nu = 5, lambda = 1, B = B_size, samp_size = sampsize_ori)
   
   test0_nc_wb11 <- Bootstrap_c_wb(Z, hzrate0, nu = 1, lambda = 1, B = B_size, samp_size = sampsize_ori)
   test1_nc_wb11 <- Bootstrap_c_wb(Z, hzrate1, nu = 1, lambda = 1, B = B_size, samp_size = sampsize_ori)
   test5_nc_wb11 <- Bootstrap_c_wb(Z, hzrate5, nu = 1, lambda = 1, B = B_size, samp_size = sampsize_ori)
   testm1_nc_wb11 <- Bootstrap_c_wb(Z, hzrate_m1, nu = 1, lambda = 1, B = B_size, samp_size = sampsize_ori)
   testm2_nc_wb11 <- Bootstrap_c_wb(Z, hzrate_m2, nu = 1, lambda = 1, B = B_size, samp_size = sampsize_ori)
    
   
  return(list("test0_nc_exp1" = test0_nc_exp1,
              "test1_nc_exp1" = test1_nc_exp1,
              "test5_nc_exp1" = test5_nc_exp1,
              "testm1_nc_exp1" = testm1_nc_exp1,
              "testm2_nc_exp1" = testm2_nc_exp1,
              
              "test0_nc_exp05" = test0_nc_exp05,
              "test1_nc_exp05" = test1_nc_exp05,
              "test5_nc_exp05" = test5_nc_exp05,
              "testm1_nc_exp05" = testm1_nc_exp05,
              "testm2_nc_exp05" = testm2_nc_exp05,
              
              "test0_nc_gp011" = test0_nc_gp011,
              "test1_nc_gp011" = test1_nc_gp011,
              "test5_nc_gp011" = test5_nc_gp011,
              "testm1_nc_gp011" = testm1_nc_gp011,
              "testm2_nc_gp011" = testm2_nc_gp011,
              
              "test0_nc_gp12" = test0_nc_gp12,
              "test1_nc_gp12" = test1_nc_gp12,
              "test5_nc_gp12" = test5_nc_gp12,
              "testm1_nc_gp12" = testm1_nc_gp12,
              "testm2_nc_gp12" = testm2_nc_gp12,
              
              "test0_nc_wb51" = test0_nc_wb51,
              "test1_nc_wb51" = test1_nc_wb51,
              "test5_nc_wb51" = test5_nc_wb51,
              "testm1_nc_wb51" = testm1_nc_wb51,
              "testm2_nc_wb51" = testm2_nc_wb51,
              
              "test0_nc_wb11" = test0_nc_wb11,
              "test1_nc_wb11" = test1_nc_wb11,
              "test5_nc_wb11" = test5_nc_wb11,
              "testm1_nc_wb11" = testm1_nc_wb11,
              "testm2_nc_wb11" = testm2_nc_wb11))
  
}

stopImplicitCluster()
stopCluster(cl)

end_time = Sys.time()
cat('Time used:\n')
print(end_time - start_time)

save(BHT_cnc, file = "BHT_cnc.rda")
