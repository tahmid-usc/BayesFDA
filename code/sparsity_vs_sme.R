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


# check hyperparameter value

fdata <- fdagen(n = 10, gridSize = 100, sparsity = 1, muf = muf1, theta = c(1,1,1))
matplot(fdata$grid, t(fdata$subData), type = 'l')


#----


# generate data

source('code/RBF.R')
fdata1 <- fdagen_fs(n = 20, mint= 2, maxt = 10, muf = muf1, theta = c(1, 1, 1))

fdata2 <- fdagen(n = 20, gridSize = 100, sparsity = .06, muf = muf1, theta = c(1,1,1))
fdata2 <- fdata2$sparseData



#plot
plotFdata(fdata1$id, fdata1$t, fdata1$y)

plotFdata(fdata2$id, fdata2$t, fdata2$y)



# GP fit for muf1

testTime <- seq(0, 1, length.out = 100)

#Gaussian kernel & posterior distribution
source('code/RBF.R')
fet1.muf1.rbf <- feature(fdata1$id, fdata1$t, fdata1$y)
pd1_muf1_rbf <- postDist(testTime, fet1.muf1.rbf)

fet2.muf1.rbf <- feature(fdata2$id, fdata2$t, fdata2$y)
pd2_muf1_rbf <- postDist(testTime, fet2.muf1.rbf)



# PACE

fpca.dt1 <- MakeFPCAInputs(IDs = fdata1$id, fdata1$t, fdata1$y, sort = FALSE) 
fpca1 <- FPCA(fpca.dt1$Ly, fpca.dt1$Lt, optns = list(dataType = 'Sparse')) #fit model
pace_mu1 <- fpcamu(testTime, fpca1)

fpca.dt2 <- MakeFPCAInputs(IDs = fdata2$id, fdata2$t, fdata2$y, sort = FALSE) 
fpca2 <- FPCA(fpca.dt2$Ly, fpca.dt2$Lt, optns = list(dataType = 'Sparse')) #fit model
pace_mu2 <- fpcamu(testTime, fpca2)



# SME

sme1 <- sme(data.frame(tme = fdata1$t, ind = as.factor(fdata1$id), y = fdata1$y), 
            criteria = 'AICc', maxIter = 10000)
sme_mu1 <- sme_mu(testTime, sme1)


sme2 <- sme(data.frame(tme = fdata2$t, ind = as.factor(fdata2$id), y = fdata2$y), 
            criteria = 'AIC', maxIter = 10000, initial.lambda.mu = 0, initial.lambda.v = 0)
sme_mu2 <- sme_mu(testTime, sme2)


# sme3 <- spline(x = as.numeric(colnames(sme1$coefficients)), y =  sme1$coefficients[1,], 
#                xout = testTime, method = "natural")
# 
# sme4 <- spline(x = as.numeric(colnames(sme2$coefficients)), y =  sme2$coefficients[1,], 
#                xout = testTime, method = "natural")


# plot

plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime, muf1(testTime), lwd = 4)
lines(testTime, pd1_muf1_rbf$mu, col = 2, lwd = 4)
lines(testTime, pace_mu1, col = 3, lwd = 4)
lines(testTime, sme_mu1, col = 4, lwd = 4)
#lines(testTime, sme3$y, col = 5, lwd = 4)
legend('topright', c('True', 'GP', 'LP', 'SME'), col = 1:4, bty = 'n', lty = 1, lwd = 2)



plotFdata(fdata2$id, fdata2$t, fdata2$y)
lines(testTime, muf1(testTime), lwd = 4)
lines(testTime, pd2_muf1_rbf$mu, col = 2, lwd = 4)
lines(testTime, pace_mu2, col = 3, lwd = 4)
lines(testTime, sme_mu2, col = 4, lwd = 4)
#lines(testTime, sme4$y, col = 5, lwd = 4)
legend('topright', c('True', 'GP', 'LP', 'SME'), col = 1:4, bty = 'n', lty = 1, lwd = 2)


