#loading packages
rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(fdapace)
library(sme)
library(dplyr)
library(cubature)
library(smoothmest)
library(kableExtra)
library(colorspace)
library(purrr)

source('code/RBF.R')
source('code/functions.R')
#source('code/functionsNew.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')



#----


# generate data

sampleSize = 10
mint = 2
maxt = 10

source('code/RBF.R')
fdata1 <- fdagen_fs(n = sampleSize, mint= mint, maxt = maxt, muf = muf1, theta = c(1, 5, 2))


#plot
plotFdata(fdata1$id, fdata1$t, fdata1$y)



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

# posterior distribution

pd_muf1_rbf <- postMu(testTime, fet.muf1.rbf)
pd_muf1_lap <- postMu(testTime, fet.muf1.lap)
pd_muf1_m52 <- postMu(testTime, fet.muf1.m52)
pd_muf1_m32 <- postMu(testTime, fet.muf1.m32)




# PACE

fpca.dt1 <- MakeFPCAInputs(IDs = fdata1$id, fdata1$t, fdata1$y, sort = FALSE) 
fpca1 <- FPCA(fpca.dt1$Ly, fpca.dt1$Lt, optns = list(dataType = 'Sparse')) #fit model
pace_mu1 <- fpcamu(testTime, fpca1)





# SME

sme1 <- sme(data.frame(tme = fdata1$t, ind = as.factor(fdata1$id), y = fdata1$y), 
               criteria = 'AIC', maxIter = 10000, initial.lambda.mu = 0, initial.lambda.v = 0)
sme_mu1 <- sme_mu(testTime, sme1)



# Check estimates are working

plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime, pd_muf1_rbf, col = 2, lwd = 4)
lines(testTime, pd_muf1_lap, col = 2, lwd = 4, lty = 2)
lines(testTime, pd_muf1_m52, col = 2, lwd = 4, lty = 3)
lines(testTime, pd_muf1_m32, col = 2, lwd = 4, lty = 4)
lines(testTime, pace_mu1, col = 3, lwd = 4)
lines(testTime, sme_mu1, col = 4, lwd = 4)


# Standardized Average Squared Error (SASE)

N <- 10000
testTime <- seq(0, 1, length.out = N)

# true mean

true <- muf1(testTime)

source('code/RBF.R')

mu1_rbf <- postMu(testTime, fet.muf1.rbf)

# pace and sme
sme_mu1 <- sme_mu(testTime, sme1)
pace_mu1 <- fpcamu(testTime, fpca1)


ase_muf1 <- cbind(residfunc(true, mu1_rbf),
                  residfunc(true, sme_mu1),
                  residfunc(true, pace_mu1))
ase_muf1




# MC method

testTime <- sort(runif(N, 0, 1))

# true mean

true <- muf1(testTime)

source('code/RBF.R')

mu1_rbf <- postMu(testTime, fet.muf1.rbf)

# pace and sme
sme_mu1 <- sme_mu(testTime, sme1)
pace_mu1 <- fpcamu(testTime, fpca1)


mc_ase_muf1 <- cbind(residfunc(true, mu1_rbf),
                  residfunc(true, sme_mu1),
                  residfunc(true, pace_mu1))
mc_ase_muf1










