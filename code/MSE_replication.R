#loading packages
rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(optimx)
library(fdapace)
library(sme)
library(dplyr)
library(smoothmest)



source('code/functions.R')



#----
errorEstMean <- function(sampleSize = 10, muf = muf1, mint = 2, maxt = 8) {
  
  # generate data
  
  source('code/RBF.R')
  fdata1 <- fdagen_fs(n = sampleSize, mint= mint, maxt = maxt, muf = muf1, theta = c(1, 5, 2))
  
  
  # GP fit for muf1
  
  testTime <- seq(0, 1, length.out = 100)
  
  #Gaussian kernel
  source('code/RBF.R')
  fet.muf1.rbf <- feature(fdata1$id, fdata1$t, fdata1$y)
  
  source('code/laplace.R')
  
  fet.muf1.lap <- feature(fdata1$id, fdata1$t, fdata1$y)
  
  source('code/matern52.R')
  
  fet.muf1.m52 <- feature(fdata1$id, fdata1$t, fdata1$y)
  
  source('code/matern32.R')
  
  fet.muf1.m32 <- feature(fdata1$id, fdata1$t, fdata1$y)
  
  
  # PACE
  
  fpca.dt1 <- MakeFPCAInputs(IDs = fdata1$id, fdata1$t, fdata1$y, sort = FALSE) 
  fpca1 <- FPCA(fpca.dt1$Ly, fpca.dt1$Lt, optns = list(dataType = 'Sparse')) #fit model
  
  # SME
  
  sme1 <- sme(data.frame(tme = fdata1$t, ind = as.factor(fdata1$id), y = fdata1$y), 
              criteria = 'AIC', maxIter = 10000, initial.lambda.mu = 0, initial.lambda.v = 0)

  
  
  # Standardized Average Squared Error (SASE)
  
  
  N <- 10000
  testTime <- seq(0, 1, length.out = N)
  
  # true mean
  
  true <- muf1(testTime)
  
  # GP posterior means
  
  source('code/RBF.R')
  
  mu1_rbf <- postMu(testTime, fet.muf1.rbf)
  
  source('code/laplace.R')
  
  mu1_lap <- postMu(testTime, fet.muf1.lap)
  
  source('code/matern52.R')
  
  mu1_m52 <- postMu(testTime, fet.muf1.m52)
  
  source('code/matern32.R')
  
  mu1_m32 <- postMu(testTime, fet.muf1.m32)
  
  
  # pace and sme
  sme_mu1 <- sme_mu(testTime, sme1)
  pace_mu1 <- fpcamu(testTime, fpca1)
  
  
  ase_muf1 <- cbind(residfunc(true, mu1_rbf),
                    residfunc(true, mu1_lap),
                    residfunc(true, mu1_m52),
                    residfunc(true, mu1_m32),
                    residfunc(true, sme_mu1),
                    residfunc(true, pace_mu1))

  
  # MC method
  
  testTime <- sort(runif(N, 0, 1))
  
  # true mean
  
  true <- muf1(testTime)
  
  # GP posterior means
  
  source('code/RBF.R')
  
  mc_mu1_rbf <- postMu(testTime, fet.muf1.rbf)
  
  source('code/laplace.R')
  
  mc_mu1_lap <- postMu(testTime, fet.muf1.lap)
  
  source('code/matern52.R')
  
  mc_mu1_m52 <- postMu(testTime, fet.muf1.m52)
  
  source('code/matern32.R')
  
  mc_mu1_m32 <- postMu(testTime, fet.muf1.m32)
  
  
  # pace and sme
  mc_sme_mu1 <- sme_mu(testTime, sme1)
  mc_pace_mu1 <- fpcamu(testTime, fpca1)
  
  
  mc_ase_muf1 <- cbind(residfunc(true, mc_mu1_rbf),
                    residfunc(true, mc_mu1_lap),
                    residfunc(true, mc_mu1_m52),
                    residfunc(true, mc_mu1_m32),
                    residfunc(true, mc_sme_mu1),
                    residfunc(true, mc_pace_mu1))
  
  return(list(SASE= ase_muf1, MC = mc_ase_muf1))
  
}

sampleSize = 10
muf = muf1
maxt = 10

result <- c()
for(i in 1:2) {
  err <- errorEstMean(sampleSize, muf, maxt)
  result$SASE <- rbind(result$SASE, err$SASE) 
  result$MC <- rbind(result$MC, err$MC)
}


#save(list = ls(all.names = TRUE), file = "result/mse_muf2_n20_maxt10.RData", envir = .GlobalEnv)
