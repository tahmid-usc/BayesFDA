

library(readr)
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(readxl)


# Extract relevant feature of the data

feature <- function(id, x, y) {
  
  n <- length(unique(id))
  N <- length(id)
  # standardize training data
  ymean <- mean(y)
  ysd <- sd(y)
  xmax <- max(x)
  x <- x/xmax
  y <- (y - ymean)/ysd
  
  rep <- as.numeric(table(id))
# hyper <- Hyper(x = x, y = y, repnum = rep, N = N)
  hyper <- Hyper.ms(x = x, y = y, repnum = rep, N = N)
  kinv <- chol2inv(chol(covmat(x, rep , theta = hyper[1:4]) + hyper[5] * diag(N)))
  return(list(n = n,N = N,rep = rep, trainx = x, trainy = y, hyper = hyper, 
              kinv = kinv, ymean = ymean, ysd = ysd, xmax = xmax))
  
}


# Estimate hypoerparameter

Hyper <- function(x, y, repnum, N) {
  
  marlik <- function(theta) {
    theta <- theta^2
    kxx <- covmat(trainx = x, repnum = repnum, theta = theta[1:4])
    -dmvnorm(x = y, sigma = kxx + theta[5] * diag(N), log = T)
  }
  
  hyp <- optim(par=rep(1, 5), fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 100000))
  print(hyp)
  return(hyp$par^2)
  
}

# multistart
Hyper.ms <- function(x, y, repnum, N) {
  
  marlik <- function(theta) {
    theta <- theta^2
    kxx <- covmat(trainx = x, repnum = repnum, theta = theta[1:4])
    -dmvnorm(x = y, sigma = kxx + theta[5] * diag(N), log = T)
  }
  
  #  parval <- c(.01, .05, .1, .5, 1, 2, 3, 4, 5)
  parmat <- matrix(rexp(50, 1), ncol = 5)
  hyp <- multistart(parmat, fn = marlik, method='Nelder-Mead', control = list(maxit = 10000))
  hyp <- as.numeric(hyp[which(hyp$value == min(hyp$value)), 1:5])^2
  return(hyp)
  
}


# Posterior distribution ----------------------------

postDist <- function(x, train) {
  
  n <- length(x)
  x <- x/train$xmax
  kxx <- covmat(trainx = train$trainx, repnum = train$rep, 
                theta = train$hyper[1:4])
  kx <-  testcov(x = x, y = train$trainx, theta = train$hyper[1:4])
  #kinv <-  chol2inv(chol(kxx + trainlist$hyper[5] * diag(trainlist$N)))
  k <- kx %*% train$kinv
  pred <- k %*% as.matrix(train$trainy)

  sigma <- testmat(x = x, theta = train$hyper[1:4]) + (train$hyper[5] * diag(n)) - (k %*% t(kx))
  sigma <- as.matrix(forceSymmetric(sigma))

  ul <- pred + 1.96 * sqrt(diag(sigma))
  ll <- pred - 1.96 * sqrt(diag(sigma))
  
  pred <-  train$ymean + train$ysd * pred
  ul <- train$ymean + train$ysd * ul
  ll <- train$ymean + train$ysd * ll
  
  return(list(mu = pred, sigma = sigma, ul = ul, ll = ll))
  
}




# Local linear smoother (estimate of FPCA mean)
fpcamu <- function(t, fpca) {
  
  to <- order(t)
  t <- sort(t)
  xt <- unlist(fpca$inputData$Lt)
  yt <- unlist(fpca$inputData$Ly)
  dt <- data.frame(xt, yt)
  dt <- arrange(dt, xt)
  
  mu <- Lwls1D(bw = fpca$bwMu, kernel_type = 'gauss', xin = dt$xt, 
               yin = dt$yt, xout = t)
  mu <- mu[to]
  return(mu)
}


