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
#source('code/functions.R')
#source('code/functionsNew.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')

# Sparsified
#mu1

# load

dname <- 'n20_sp10_mu1'
fpath <- paste0("E:/R Project/BayesFDA/data/sim/sparsified/", dname, ".RData")

load(fpath)


testTime <- seq(0, 1, length.out = 100)
source('code/functions.R')

# PACE

fpca.dt1 <- MakeFPCAInputs(IDs = fdata1$id, fdata1$t, fdata1$y, sort = FALSE) 
fpca1 <- FPCA(fpca.dt1$Ly, fpca.dt1$Lt, optns = list(dataType = 'Sparse')) #fit model
pace_mu1 <- fpcamu(testTime, fpca1)





# SME

sme1 <- sme(data.frame(tme = fdata1$t, ind = as.factor(fdata1$id), y = fdata1$y), 
            lambda.mu=1e-6, lambda.v=1, knots=seq(0.05, 0.95, by=0.05))
sme_mu1 <- sme_mu(testTime, sme1)


# plot

plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime, pd_muf1_rbf, col = 2, lwd = 4)
lines(testTime, pd_muf1_lap, col = 2, lwd = 4, lty = 2)
lines(testTime, pd_muf1_m52, col = 2, lwd = 4, lty = 3)
lines(testTime, pd_muf1_m32, col = 2, lwd = 4, lty = 4)
lines(testTime, pace_mu1, col = 3, lwd = 4)
lines(testTime, sme_mu1, col = 4, lwd = 4)




#----


pdf(file = paste0('plot/sim/sparsified/', dname, '_rbf.pdf'), width = 14, height = 6)
plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime, muf1(testTime), lwd = 4)
lines(testTime, pd_muf1_rbf, col = 2, lwd = 4)
lines(testTime, pace_mu1, col = 3, lwd = 4)
lines(testTime, sme_mu1, col = 4, lwd = 4)
legend('topleft', c('True', 'GP', 'LP', 'SME'), lty = 1, lwd =2, col = 1:4, bty = 'n')
dev.off()


pdf(file = paste0('plot/sim/sparsified/', dname, '_lap.pdf'), width = 14, height = 6)
plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime, muf1(testTime), lwd = 4)
lines(testTime, pd_muf1_lap, col = 2, lwd = 4)
lines(testTime, pace_mu1, col = 3, lwd = 4)
lines(testTime, sme_mu1, col = 4, lwd = 4)
legend('topleft', c('True', 'GP', 'LP', 'SME'), lty = 1, lwd =2, col = 1:4, bty = 'n')
dev.off()

pdf(file = paste0('plot/sim/sparsified/', dname, '_m52.pdf'), width = 14, height = 6)
plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime, muf1(testTime), lwd = 4)
lines(testTime, pd_muf1_m52, col = 2, lwd = 4)
lines(testTime, pace_mu1, col = 3, lwd = 4)
lines(testTime, sme_mu1, col = 4, lwd = 4)
legend('topleft', c('True', 'GP', 'LP', 'SME'), lty = 1, lwd =2, col = 1:4, bty = 'n')
dev.off()


pdf(file = paste0('plot/sim/sparsified/', dname, '_m32.pdf'), width = 14, height = 6)
plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime, muf1(testTime), lwd = 4)
lines(testTime, pd_muf1_m32, col = 2, lwd = 4)
lines(testTime, pace_mu1, col = 3, lwd = 4)
lines(testTime, sme_mu1, col = 4, lwd = 4)
legend('topleft', c('True', 'GP', 'LP', 'SME'), lty = 1, lwd =2, col = 1:4, bty = 'n')
dev.off()
