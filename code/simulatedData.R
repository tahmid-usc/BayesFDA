rm(list=ls())

library(fdapace)

source("code/functions.R")
#choose kernel
source("code/RBF.R")
#source("code/laplace.R")
#source("code/matern52.R")
#source("code/matern32.R")


# data

fdata_full <- fdagen(n = 20, gridSize = 30, sparsity = .2, muf = muf1)
fdata1 <- fdata_full$sparseData
#fdata1 <- unnest(fdata1)


# Fit the GP model
gpFit1 <- feature(id = fdata1$id, x = fdata1$t, y = fdata1$y)
testTime <- seq(0, 1, length.out = 100)
postDist1 <- postDist(testTime, gpFit1)



# PACE
fpca.dt1 <- MakeFPCAInputs(IDs = fdata1$id, fdata1$t, fdata1$y, sort = FALSE)
fpca1 <- FPCA(fpca.dt1$Ly, fpca.dt1$Lt, optns = list(dataType = 'Sparse'))

plot(fpca1)

fm <- GetMeanCurve(Ly = fpca.dt1$Ly, Lt =fpca.dt1$Lt, optns = list(plot = TRUE))

fittedY <- fitted(fpca1, ciOptns = list(alpha=0.05))


# face package: P-spline

library(face)

argvals <- fdata1$t
subj <- fdata1$id
y <- fdata1$y

data <- data.frame(argvals, subj, y)
fit_face <- face.sparse(data, argvals.new = testTime)

plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime, fit_face$mu.new, lwd = 4)

# Refund package: B spline
library(refund)

ydata <- data.frame(.id = fdata1$id, .index = fdata1$t, .value = fdata1$y)
Fit.MM <- fpca.sc(ydata=ydata, var = TRUE)

library(sme)

fit_sme <- sme(data.frame(tme = fdata1$t, ind = as.factor(fdata1$id), y = fdata1$y), 
               criteria = 'AIC', maxIter = 10000, initial.lambda.mu = 0, initial.lambda.v = 0)
sme_mu <- spline(x = as.numeric(colnames(fit_sme$coefficients)), y = fit_sme$coefficients[1,], 
                       xout = testTime, method = "natural")

#plot
par(cex.main = 1.5, cex.lab=1.5, cex.axis=1.5, oma = c(1,1,1,1), mar = c(6,6,3,3)) 
plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime,  postDist1$mu, lwd = 5, col = transparent(2, trans.val = .4))
lines(testTime, muf1(testTime), lwd = 5, col = transparent(1, trans.val = .4))
lines(fm$workGrid, fm$mu, lwd = 5, col = transparent(6, trans.val = .4))
lines(testTime, fit_face$mu.new, lwd = 5, col = transparent(4, trans.val = .4))
lines(sort(unique(fdata1$t)), Fit.MM$mu, lwd = 5, col = transparent(5, trans.val = .4))
lines(sme_mu$x, sme_mu$y, lwd = 5, col = transparent(7, trans.val = .4))
legend('topleft', c('True','GP', 'P-spline', 'B-Spline', 'Local polynomial (PACE)', 'SME'), lty = 1, lwd = 2, 
       bty = 'n', col = 1:7)


# Subject Specific

sub1_mu <- postDist_sub(testTime, gpFit1, sub_id = 5)  


#plot
par(cex.main = 1.5, cex.lab=1.5, cex.axis=1.5, oma = c(1,1,1,1), mar = c(6,6,3,3)) 
plotFdata(fdata1$id, fdata1$t, fdata1$y)
lines(testTime, sub1_mu$mu, lwd = 3, col = 6)


# plot
#par(mfrow = c(2,2))
#for(j in 1:4){
  
j <- 4
  sub_mu <- postDist_sub(testTime, gpFit1, sub_id = j)  
  
  plot(fdata1$t, fdata1$y, type = 'n' , xlab ="Time", ylab="y")
  lines(testTime, sub_mu$mu, lwd = 4, col = 2)
  lines(pull(fdata1[fdata1$id == j,1]), pull(fdata1[fdata1$id == j,2]), type = 'b', lwd = 4, col = rgb(0, 0, 0, .2))
  lines(fittedY$workGrid, fittedY$fitted[j,], lwd = 4, col = 6)
  lines(fdata_full$grid, fdata_full$subData[j,], lwd = 4, col = 1)
  legend('topleft', c('True subject curve', 'GP', 'PACE'), lty = 1, lwd = 2, bty = 'n', col = c(1,2,6))
#}  
