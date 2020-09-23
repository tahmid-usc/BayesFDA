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


# Confidence band 
# mu2
#fully sparse
# load

dname <- 'n20_min6_max14_mu2'
fpath <- paste0("E:/R Project/BayesFDA/data/sim/fully_sparse/", dname, ".RData")

load(fpath)

source('code/functions.R')
testTime <- seq(0, 1, length.out = 100)


#Gaussian kernel & posterior distribution
source('code/RBF.R')
pd_muf2_rbf <- postDist(testTime, fet.muf1.rbf)

source('code/laplace.R')
pd_muf2_lap <- postDist(testTime, fet.muf1.lap)

source('code/matern52.R')
pd_muf2_m52 <- postDist(testTime, fet.muf1.m52)

source('code/matern32.R')
pd_muf2_m32 <- postDist(testTime, fet.muf1.m32)


# plot

plotFdata(fdata1$id, fdata1$t, fdata1$y)
polygon(c(testTime,rev(testTime)),c(pd_muf2_rbf$ll,rev(pd_muf2_rbf$ul)), 
        col= rgb(0,0,.5,.25),border=NA)
lines(testTime, muf2(testTime), lwd = 4)
lines(testTime, pd_muf2_rbf$mu, col = 2, lwd = 4)
legend('topleft', c('True', 'GP'), lty = 1, lwd =2, col = 1:2, bty = 'n')





#----


pdf(file = paste0('plot/sim/fully_sparse/', dname, '_conf_rbf.pdf'), width = 14, height = 6)
plotFdata(fdata1$id, fdata1$t, fdata1$y)
polygon(c(testTime,rev(testTime)), c(pd_muf2_rbf$ll,rev(pd_muf2_rbf$ul)), 
        col= rgb(0,0,.5,.25),border=NA)
lines(testTime, muf2(testTime), lwd = 4)
lines(testTime, pd_muf2_rbf$mu, col = 2, lwd = 4)
legend('topleft', c('True', 'GP'), lty = 1, lwd =2, col = 1:2, bty = 'n')
dev.off()


pdf(file = paste0('plot/sim/fully_sparse/', dname, '_conf_lap.pdf'), width = 14, height = 6)
plotFdata(fdata1$id, fdata1$t, fdata1$y)
polygon(c(testTime,rev(testTime)), c(pd_muf2_lap$ll,rev(pd_muf2_lap$ul)), 
        col= rgb(0,0,.5,.25),border=NA)
lines(testTime, muf2(testTime), lwd = 4)
lines(testTime, pd_muf2_lap$mu, col = 2, lwd = 4)
legend('topleft', c('True', 'GP'), lty = 1, lwd =2, col = 1:2, bty = 'n')
dev.off()

pdf(file = paste0('plot/sim/fully_sparse/', dname, '_conf_m52.pdf'), width = 14, height = 6)
plotFdata(fdata1$id, fdata1$t, fdata1$y)
polygon(c(testTime,rev(testTime)), c(pd_muf2_m52$ll,rev(pd_muf2_m52$ul)), 
        col= rgb(0,0,.5,.25),border=NA)
lines(testTime, muf2(testTime), lwd = 4)
lines(testTime, pd_muf2_m52$mu, col = 2, lwd = 4)
legend('topleft', c('True', 'GP'), lty = 1, lwd =2, col = 1:2, bty = 'n')
dev.off()


pdf(file = paste0('plot/sim/fully_sparse/', dname, '_conf_m32.pdf'), width = 14, height = 6)
plotFdata(fdata1$id, fdata1$t, fdata1$y)
polygon(c(testTime,rev(testTime)), c(pd_muf2_m32$ll,rev(pd_muf2_m32$ul)), 
        col= rgb(0,0,.5,.25),border=NA)
lines(testTime, muf2(testTime), lwd = 4)
lines(testTime, pd_muf2_m32$mu, col = 2, lwd = 4)
legend('topleft', c('True', 'GP'), lty = 1, lwd =2, col = 1:2, bty = 'n')
dev.off()
