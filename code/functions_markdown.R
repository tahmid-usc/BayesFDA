library(readr)
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(readxl)
library(tidyverse)
library(smoothmest)

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
  return(list(n = n,N = N,rep = rep, trainx = x, trainy = y, id = id, hyper = hyper, 
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




postDist_sub <- function(x, train, sub_id) {
  
  n <- length(x)
  x <- x/train$xmax
#  kxx <- covmat(trainx = train$trainx, repnum = train$rep, theta = train$hyper[1:4])
  kx <-  testcov(x = x, y = train$trainx, theta = train$hyper[1:4])
  
  time_sub <- train$trainx
  time_sub <- time_sub[train$id == sub_id]
  kx_sub <- testcov_sub(x = x, y = time_sub, theta = train$hyper[1:4])
  sub_ind <- which(train$id == sub_id)
  kx[1:length(x), min(sub_ind):max(sub_ind)] <- kx_sub
  
  k <- kx %*% train$kinv
  pred <- k %*% as.matrix(train$trainy)
  
  sigma <- testmat(x = x, theta = train$hyper[1:4]) + (train$hyper[5] * diag(n)) - (k %*% t(kx))
  sigma <- as.matrix(forceSymmetric(sigma))
  diag(sigma) <- ifelse(diag(sigma)<0, 0, diag(sigma))
  
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


#-----------------------------------

# Simulation functions-----------------------


# Mean functions 

# Mix of normal densities over (0,1)

muf1 <- function(x) {
  fx <- #.25 * dnorm(x, mean = .2, sd = .05)
    dnorm(x, mean = .2, sd = .05) + 
    dnorm(x, mean = .5, sd = .03) +   
    dnorm(x, mean = .8,sd = .02) 
  return(fx)
}


# Mix of laplace densities over (0,1)

muf2 <- function(x) {
  fx <- 
    ddoublex(x, .2, .05) + 
    ddoublex(x, .5,  .03) +   
    ddoublex(x, .8, .02) 
  return(fx)
}





fdagen <- function(n = 10, gridSize = 100, sparsity = .5, muf, theta = rep(1,3)) {
  
  source("../code/RBF.R")
  # n = number of functions in the sample data
  
  grid <- seq(0, 1, length.out = gridSize)
  
  denseData <- data.frame()
  subData <- data.frame()
  for(i in 1:n){
    t <- grid
    mu <- muf(t)
    gt <- mvrnorm(1, rep(0, gridSize), ker(t, l = theta[1], sigf = theta[2])) 
    sub_y <- mu + gt
    y <- sub_y + rnorm(gridSize, 0, theta[3])
    denseData <- rbind(denseData, t(y))
    subData <- rbind(subData, t(sub_y))
  }
  colnames(denseData) <- as.character(1:gridSize)
  colnames(subData) <- as.character(1:gridSize)
  
  
  id <- c()
  sparseData <- data.frame()
  for(i in 1:n) {
    sel <- rbinom(n = gridSize, size = 1, prob = sparsity)
    sel <- as.logical(sel)
    if(sum(sel) < 2) break 
    t <- grid[sel]
    y <- denseData[i,]
    y <- y[sel]
    id <- rep(i, times = sum(sel))
    t <- as.vector(t) 
    y <- matrix(y, ncol = 1) 
    id <- as.vector(id)
    sparseData <- rbind(sparseData, cbind(t, y, id))
  }
  colnames(sparseData) <- c('t', 'y', 'id')
  
  denseData <- unnest(denseData)
  sparseData <- unnest(sparseData)
  subData <- unnest(subData)
  
 return(list(denseData = denseData, sparseData = sparseData, subData = subData, grid = grid))

}



# Generate truly sparse data directly

fdagen_fs <- function(n = 10, mint = 2, maxt = 20, muf, theta = rep(1,3)) {
  
  train <- data.frame()
  #n number of functions in the sample data
  n.time <- sample(mint:maxt, size=n, replace=T) 
  
  for(i in 1:n) {
    id <- rep(i, n.time[i])
    t <- sort(runif(n.time[i], 0 , 1))
    mu <- muf(t)
    gt <- mvrnorm(1, rep(0, n.time[i]), ker(t, l = theta[1], sigf = theta[2])) 
    #y <- mu + gt
    y <- mu + gt + rnorm(n.time[i], 0, theta[3])
    train <- rbind(train, cbind(t, y, id))
  }
  return(train)
}



# plot functional data
plotFdata <- function(id, x, y) {
  plot(x, y, type='n', xlab ="Time", ylab="y")
  for(i in unique(id)){
    lines(x[id == i], y[id == i], col = 'grey', type = 'b', lwd = 2)
  }
  
}


# subject specific prediction

testcov_sub <- function(x, y, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  k1 <- theta[2] * kernelMatrix(rbf1, x = x, y = y)
  rbf2 <- rbfdot(sigma = abs(theta[3]))
  k2 <- abs(theta[4]) * kernelMatrix(rbf2, x = x, y = y)
  return(k1 + k2)
}

gpsmooth_sub <-  function(x, trainlist, sub_id) {
  
  kx <-  testcov(x = x, y = trainlist$trainx, theta = trainlist$hyper[1:4])
  
  time_sub <- trainlist$trainx
  time_sub <- time_sub[trainlist$id == sub_id]
  kx_sub <- testcov_sub(x = x, y = time_sub, theta = trainlist$hyper[1:4])
  sub_ind <- which(trainlist$id == sub_id)
  kx[1:length(x), min(sub_ind):max(sub_ind)] <- kx_sub
  
  k <- kx %*% trainlist$kinv
  pred <- k %*% as.matrix(trainlist$trainy)
  pred <-  trainlist$ymean + trainlist$ysd * pred
  return(pred)
  
}


