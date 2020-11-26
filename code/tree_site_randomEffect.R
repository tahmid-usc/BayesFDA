#Fit tree data assuming 3 sites as random effects

library(readr)
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(readxl)
library(tidyverse)

# RBF kernel

ker <- function(x, l, sigf) {
  rbf <- rbfdot(sigma = 1/l)
  return(sigf * kernelMatrix(rbf, x = x))
}

Jmat <- function(m) {
  return(matrix(rep(1, m^2), ncol = m))
}

BDmat <- function(repnum) {
  mat <- lapply(repnum, Jmat)
  as.matrix(bdiag(mat))
}

covmat <- function(trainx, repnum1, repnum2, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  rbf2 <- rbfdot(sigma = 1/theta[3])
  rbf3 <- rbfdot(sigma = 1/theta[5])
  k1 <- theta[2] * kernelMatrix(rbf1, x = trainx)
  k2 <- theta[4] * kernelMatrix(rbf2, x = trainx)
  k3 <- theta[6] * kernelMatrix(rbf3, x = trainx)
  k2 <- k2 * BDmat(repnum1)
  k3 <- k3 * BDmat(repnum2)
  return(k1 + k2 + k3)
}

testcov <- function(x, y, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  #rbf2 <- rbfdot(sigma = abs(theta[3]))
  k1 <- theta[2] * kernelMatrix(rbf1, x = x, y = y)
  #k2 <- abs(theta[4]) * kernelMatrix(rbf2, x = x, y = y)
  return(k1)
}

testmat <- function(x, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  rbf2 <- rbfdot(sigma = 1/theta[3])
  rbf3 <- rbfdot(sigma = 1/theta[5])
  k1 <- theta[2] * kernelMatrix(rbf1, x = x)
  k2 <- theta[4] * kernelMatrix(rbf2, x = x)
  k3 <- theta[6] * kernelMatrix(rbf3, x = x)
  return(k1 + k2 + k3)
}



#------------------------

# Extract relevant feature of the data

feature <- function(id, site, x, y) {
  
  n <- length(unique(id))
  N <- length(id)
  # standardize training data
  ymean <- mean(y)
  ysd <- sd(y)
  xmax <- max(x)
  x <- x/xmax
  y <- (y - ymean)/ysd
  
  rep1 <- as.numeric(table(id))
  rep2 <- as.numeric(table(site))
  #hyper <- Hyper(x = x, y = y, repnum1 = rep1, repnum2 = rep2, N = N)
  hyper <- Hyper.ms(x = x, y = y, repnum1 = rep1, repnum2 = rep2, N = N)
  kinv <- chol2inv(chol(covmat(x, rep1, rep2, theta = hyper[1:6]) + hyper[7] * diag(N)))
  return(list(n = n,N = N,rep1 = rep1, rep2 = rep2, trainx = x, trainy = y, id = id, hyper = hyper, 
              kinv = kinv, ymean = ymean, ysd = ysd, xmax = xmax, site = site))
  
}




Hyper <- function(x, y, repnum1, repnum2, N) {
  
  marlik <- function(theta) {
    theta <- theta^2
    kxx <- covmat(trainx = x, repnum1 = repnum1, repnum2 = repnum2, theta = theta[1:6])
    -dmvnorm(x = y, sigma = kxx + theta[7] * diag(N), log = T)
  }
  
  hyp <- optim(par=rep(1, 7), fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 100000))
  print(hyp)
  return(hyp$par^2)
  
}


# multistart
Hyper.ms <- function(x, y, repnum1, repnum2, N) {
  
  marlik <- function(theta) {
    theta <- theta^2
    kxx <- covmat(trainx = x, repnum1 = repnum1, repnum2 = repnum2, theta = theta[1:6])
    -dmvnorm(x = y, sigma = kxx + theta[7] * diag(N), log = T)
  }
  
  #  parval <- c(.01, .05, .1, .5, 1, 2, 3, 4, 5)
  parmat <- matrix(rexp(70, 1), ncol = 7)
  hyp <- multistart(parmat, fn = marlik, method='Nelder-Mead', control = list(maxit = 10000))
  hyp <- as.numeric(hyp[which(hyp$value == min(hyp$value)), 1:7])^2
  return(hyp)
  
}


# Posterior distribution ----------------------------

postDist <- function(x, train) {
  
  n <- length(x)
  xstar <- x
  x <- x/train$xmax
  kxx <- covmat(trainx = train$trainx, repnum1 = train$rep1, repnum2 = train$rep2,
                theta = train$hyper[1:6])
  kx <-  testcov(x = x, y = train$trainx, theta = train$hyper[1:6])
  #kinv <-  chol2inv(chol(kxx + trainlist$hyper[5] * diag(trainlist$N)))
  k <- kx %*% train$kinv
  pred <- k %*% as.matrix(train$trainy)
  
  sigma <- testmat(x = x, theta = train$hyper[1:6]) + (train$hyper[7] * diag(n)) - (k %*% t(kx))
  sigma <- as.matrix(forceSymmetric(sigma))
  
  ul <- pred + 1.96 * sqrt(diag(sigma))
  ll <- pred - 1.96 * sqrt(diag(sigma))
  
  pred <-  train$ymean + train$ysd * pred
  ul <- train$ymean + train$ysd * ul
  ll <- train$ymean + train$ysd * ll
  
  return(list(mu = pred, sigma = sigma, ul = ul, ll = ll, xstar = xstar))
  
}



# site  specific prediction

testcov_sub <- function(x, y, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  k1 <- theta[2] * kernelMatrix(rbf1, x = x, y = y)
  rbf2 <- rbfdot(sigma = theta[5])
  k2 <- theta[6] * kernelMatrix(rbf2, x = x, y = y)
  return(k1 + k2)
}



#site specific
postDist_site <- function(x, train, site) {
  
  
  n <- length(x)
  x <- x/train$xmax
  #  kxx <- covmat(trainx = train$trainx, repnum = train$rep, theta = train$hyper[1:4])
  kx <-  testcov(x = x, y = train$trainx, theta = train$hyper[1:6])
  
  time_sub <- train$trainx
  time_sub <- time_sub[train$site == site]
  kx_sub <- testcov_sub(x = x, y = time_sub, theta = train$hyper[1:6])
  sub_ind <- which(train$site == site)
  kx[1:length(x), min(sub_ind):max(sub_ind)] <- kx_sub
  
  k <- kx %*% train$kinv
  pred <- k %*% as.matrix(train$trainy)
  
  sigma <- testmat(x = x, theta = train$hyper[1:6]) + (train$hyper[7] * diag(n)) - (k %*% t(kx))
  sigma <- as.matrix(forceSymmetric(sigma))
  diag(sigma) <- ifelse(diag(sigma)<0, 0, diag(sigma))
  
  ul <- pred + 1.96 * sqrt(diag(sigma))
  ll <- pred - 1.96 * sqrt(diag(sigma))
  
  pred <-  train$ymean + train$ysd * pred
  ul <- train$ymean + train$ysd * ul
  ll <- train$ymean + train$ysd * ll
  
  return(list(mu = pred, sigma = sigma, ul = ul, ll = ll))
  
}



#----------------

tree <- read.csv("E:/R Project/BayesFDA/data/tree/tree.csv")
tree <- tree %>% dplyr::select(Site, Tree, Rep, Sp, Year, Height) %>% filter(!is.na(Height)) %>% 
  mutate(id = paste0(Tree, Rep), logHeight = log(Height)) %>% arrange(Site, id)
splitData <- split(tree, tree$Sp)



testTime <- seq(0, 30, length.out = 100)
gp_Sp <- lapply(splitData, function(L) {feature(id = L$id,site = L$Site, x = L$Year, y = L$Height)})

postDist_Sp <- lapply(gp_Sp, function(L){postDist(testTime, L)})



plotGP(gp_Sp$BS, postDist_Sp$BS, trans = .05)
plotGP(gp_Sp$EL, postDist_Sp$EL, trans = .05)
plotGP(gp_Sp$HL, postDist_Sp$HL, trans = .05)
plotGP(gp_Sp$JL, postDist_Sp$JL, trans = .05)
plotGP(gp_Sp$JP, postDist_Sp$JP, trans = .05)
plotGP(gp_Sp$RP, postDist_Sp$RP, trans = .05)
plotGP(gp_Sp$TL, postDist_Sp$TL, trans = .05)
plotGP(gp_Sp$WS, postDist_Sp$WS, trans = .05)



#------------- Site specific for all sites


pred_bri <- postDist_site(x = testTime, train = gp_Sp$BS, site = "Brighton")
pred_chSt <- postDist_site(x = testTime, train = gp_Sp$BS, site = "Chase Stream")
pred_liBa <- postDist_site(x = testTime, train = gp_Sp$BS, site = "Lily Bay")


df <- data.frame(t = testTime, mu = postDist_Sp$BS$mu, ll= postDist_Sp$BS$ll, ul = postDist_Sp$BS$ul)
df_bri <- data.frame(t = testTime, mu = pred_bri$mu, ll= pred_bri$ll, ul = pred_bri$ul)
df_chSt <- data.frame(t = testTime, mu = pred_chSt$mu, ll= pred_chSt$ll, ul = pred_chSt$ul)
df_liBa <- data.frame(t = testTime, mu = pred_liBa$mu, ll= pred_liBa$ll, ul = pred_liBa$ul)

tree %>% filter(Sp == 'BS') %>% 
  ggplot() + ggtitle('BS') +
  geom_line(aes(x = Year, y = Height, color = Site, group = id)) +
  geom_line(data = df_bri, aes(x = t, y = mu), size = 1.5, col = 'red') +
#  geom_ribbon(data = df_bri, aes(x = t, ymin = ll, ymax = ul), alpha = .5) +
  geom_line(data = df_chSt, aes(x = t, y = mu), size = 1.5, col = 'green') +
#  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
  geom_line(data = df_liBa, aes(x = t, y = mu), size = 1.5, col = 'blue') +
  geom_line(data = df, aes(x = t, y = mu), size = 1.5, col = 'black')
  #  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
  
#-------------------


pred_bri <- postDist_site(x = testTime, train = gp_Sp$EL, site = "Brighton")
pred_chSt <- postDist_site(x = testTime, train = gp_Sp$EL, site = "Chase Stream")
pred_liBa <- postDist_site(x = testTime, train = gp_Sp$EL, site = "Lily Bay")

df <- data.frame(t = testTime, mu = postDist_Sp$EL$mu, ll= postDist_Sp$EL$ll, ul = postDist_Sp$EL$ul)
df_bri <- data.frame(t = testTime, mu = pred_bri$mu, ll= pred_bri$ll, ul = pred_bri$ul)
df_chSt <- data.frame(t = testTime, mu = pred_chSt$mu, ll= pred_chSt$ll, ul = pred_chSt$ul)
df_liBa <- data.frame(t = testTime, mu = pred_liBa$mu, ll= pred_liBa$ll, ul = pred_liBa$ul)

tree %>% filter(Sp == 'EL') %>% 
  ggplot() +  ggtitle('EL') +
  geom_line(aes(x = Year, y = Height, color = Site, group = id)) +
  geom_line(data = df_bri, aes(x = t, y = mu), size = 1.5, col = 'red') +
  #  geom_ribbon(data = df_bri, aes(x = t, ymin = ll, ymax = ul), alpha = .5) +
  geom_line(data = df_chSt, aes(x = t, y = mu), size = 1.5, col = 'green') +
  #  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
  geom_line(data = df_liBa, aes(x = t, y = mu), size = 1.5, col = 'blue') +
  geom_line(data = df, aes(x = t, y = mu), size = 1.5, col = 'black')
#  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)

#------------


pred_bri <- postDist_site(x = testTime, train = gp_Sp$HL, site = "Brighton")
pred_chSt <- postDist_site(x = testTime, train = gp_Sp$HL, site = "Chase Stream")
pred_liBa <- postDist_site(x = testTime, train = gp_Sp$HL, site = "Lily Bay")

df <- data.frame(t = testTime, mu = postDist_Sp$HL$mu, ll= postDist_Sp$HL$ll, ul = postDist_Sp$HL$ul)
df_bri <- data.frame(t = testTime, mu = pred_bri$mu, ll= pred_bri$ll, ul = pred_bri$ul)
df_chSt <- data.frame(t = testTime, mu = pred_chSt$mu, ll= pred_chSt$ll, ul = pred_chSt$ul)
df_liBa <- data.frame(t = testTime, mu = pred_liBa$mu, ll= pred_liBa$ll, ul = pred_liBa$ul)

tree %>% filter(Sp == 'HL') %>% 
  ggplot() + ggtitle('HL') +
  geom_line(aes(x = Year, y = Height, color = Site, group = id)) +
  geom_line(data = df_bri, aes(x = t, y = mu), size = 1.5, col = 'red') +
  #  geom_ribbon(data = df_bri, aes(x = t, ymin = ll, ymax = ul), alpha = .5) +
  geom_line(data = df_chSt, aes(x = t, y = mu), size = 1.5, col = 'green') +
  #  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
  geom_line(data = df_liBa, aes(x = t, y = mu), size = 1.5, col = 'blue') +
  geom_line(data = df, aes(x = t, y = mu), size = 1.5, col = 'black')
#  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)

#------------


pred_bri <- postDist_site(x = testTime, train = gp_Sp$JL, site = "Brighton")
pred_chSt <- postDist_site(x = testTime, train = gp_Sp$JL, site = "Chase Stream")
pred_liBa <- postDist_site(x = testTime, train = gp_Sp$JL, site = "Lily Bay")


df <- data.frame(t = testTime, mu = postDist_Sp$JL$mu, ll= postDist_Sp$JL$ll, ul = postDist_Sp$JL$ul)
df_bri <- data.frame(t = testTime, mu = pred_bri$mu, ll= pred_bri$ll, ul = pred_bri$ul)
df_chSt <- data.frame(t = testTime, mu = pred_chSt$mu, ll= pred_chSt$ll, ul = pred_chSt$ul)
df_liBa <- data.frame(t = testTime, mu = pred_liBa$mu, ll= pred_liBa$ll, ul = pred_liBa$ul)

tree %>% filter(Sp == 'JL') %>% 
  ggplot() + ggtitle('JL') +
  geom_line(aes(x = Year, y = Height, color = Site, group = id)) +
  geom_line(data = df_bri, aes(x = t, y = mu), size = 1.5, col = 'red') +
  #  geom_ribbon(data = df_bri, aes(x = t, ymin = ll, ymax = ul), alpha = .5) +
  geom_line(data = df_chSt, aes(x = t, y = mu), size = 1.5, col = 'green') +
  #  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
  geom_line(data = df_liBa, aes(x = t, y = mu), size = 1.5, col = 'blue') +
  geom_line(data = df, aes(x = t, y = mu), size = 1.5, col = 'black')
#  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)


#------------


pred_bri <- postDist_site(x = testTime, train = gp_Sp$JP, site = "Brighton")
pred_chSt <- postDist_site(x = testTime, train = gp_Sp$JP, site = "Chase Stream")
pred_liBa <- postDist_site(x = testTime, train = gp_Sp$JP, site = "Lily Bay")


df <- data.frame(t = testTime, mu = postDist_Sp$JP$mu, ll= postDist_Sp$JP$ll, ul = postDist_Sp$JP$ul)
df_bri <- data.frame(t = testTime, mu = pred_bri$mu, ll= pred_bri$ll, ul = pred_bri$ul)
df_chSt <- data.frame(t = testTime, mu = pred_chSt$mu, ll= pred_chSt$ll, ul = pred_chSt$ul)
df_liBa <- data.frame(t = testTime, mu = pred_liBa$mu, ll= pred_liBa$ll, ul = pred_liBa$ul)

tree %>% filter(Sp == 'JP') %>% 
  ggplot() + ggtitle('JP') +
  geom_line(aes(x = Year, y = Height, color = Site, group = id)) +
  geom_line(data = df_bri, aes(x = t, y = mu), size = 1.5, col = 'red') +
  #  geom_ribbon(data = df_bri, aes(x = t, ymin = ll, ymax = ul), alpha = .5) +
  geom_line(data = df_chSt, aes(x = t, y = mu), size = 1.5, col = 'green') +
  #  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
  geom_line(data = df_liBa, aes(x = t, y = mu), size = 1.5, col = 'blue') +
  geom_line(data = df, aes(x = t, y = mu), size = 1.5, col = 'black')
#  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)

#------------


pred_bri <- postDist_site(x = testTime, train = gp_Sp$RP, site = "Brighton")
pred_chSt <- postDist_site(x = testTime, train = gp_Sp$RP, site = "Chase Stream")
pred_liBa <- postDist_site(x = testTime, train = gp_Sp$RP, site = "Lily Bay")

df <- data.frame(t = testTime, mu = postDist_Sp$RP$mu, ll= postDist_Sp$RP$ll, ul = postDist_Sp$RP$ul)
df_bri <- data.frame(t = testTime, mu = pred_bri$mu, ll= pred_bri$ll, ul = pred_bri$ul)
df_chSt <- data.frame(t = testTime, mu = pred_chSt$mu, ll= pred_chSt$ll, ul = pred_chSt$ul)
df_liBa <- data.frame(t = testTime, mu = pred_liBa$mu, ll= pred_liBa$ll, ul = pred_liBa$ul)

tree %>% filter(Sp == 'RP') %>% 
  ggplot() + ggtitle('RP') +
  geom_line(aes(x = Year, y = Height, color = Site, group = id)) +
  geom_line(data = df_bri, aes(x = t, y = mu), size = 1.5, col = 'red') +
  #  geom_ribbon(data = df_bri, aes(x = t, ymin = ll, ymax = ul), alpha = .5) +
  geom_line(data = df_chSt, aes(x = t, y = mu), size = 1.5, col = 'green') +
  #  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
  geom_line(data = df_liBa, aes(x = t, y = mu), size = 1.5, col = 'blue') +
  geom_line(data = df, aes(x = t, y = mu), size = 1.5, col = 'black')
#  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)


#------------


pred_bri <- postDist_site(x = testTime, train = gp_Sp$TL, site = "Brighton")
pred_chSt <- postDist_site(x = testTime, train = gp_Sp$TL, site = "Chase Stream")
pred_liBa <- postDist_site(x = testTime, train = gp_Sp$TL, site = "Lily Bay")


df <- data.frame(t = testTime, mu = postDist_Sp$TL$mu, ll= postDist_Sp$TL$ll, ul = postDist_Sp$TL$ul)
df_bri <- data.frame(t = testTime, mu = pred_bri$mu, ll= pred_bri$ll, ul = pred_bri$ul)
df_chSt <- data.frame(t = testTime, mu = pred_chSt$mu, ll= pred_chSt$ll, ul = pred_chSt$ul)
df_liBa <- data.frame(t = testTime, mu = pred_liBa$mu, ll= pred_liBa$ll, ul = pred_liBa$ul)

tree %>% filter(Sp == 'TL') %>% 
  ggplot() + ggtitle('TL') +
  geom_line(aes(x = Year, y = Height, color = Site, group = id)) +
  geom_line(data = df_bri, aes(x = t, y = mu), size = 1.5, col = 'red') +
  #  geom_ribbon(data = df_bri, aes(x = t, ymin = ll, ymax = ul), alpha = .5) +
  geom_line(data = df_chSt, aes(x = t, y = mu), size = 1.5, col = 'green') +
  #  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
  geom_line(data = df_liBa, aes(x = t, y = mu), size = 1.5, col = 'blue') +
  geom_line(data = df, aes(x = t, y = mu), size = 1.5, col = 'black')
#  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)



#------------

pred_bri <- postDist_site(x = testTime, train = gp_Sp$WS, site = "Brighton")
pred_chSt <- postDist_site(x = testTime, train = gp_Sp$WS, site = "Chase Stream")
pred_liBa <- postDist_site(x = testTime, train = gp_Sp$WS, site = "Lily Bay")


df <- data.frame(t = testTime, mu = postDist_Sp$WS$mu, ll= postDist_Sp$WS$ll, ul = postDist_Sp$WS$ul)
df_bri <- data.frame(t = testTime, mu = pred_bri$mu, ll= pred_bri$ll, ul = pred_bri$ul)
df_chSt <- data.frame(t = testTime, mu = pred_chSt$mu, ll= pred_chSt$ll, ul = pred_chSt$ul)
df_liBa <- data.frame(t = testTime, mu = pred_liBa$mu, ll= pred_liBa$ll, ul = pred_liBa$ul)

tree %>% filter(Sp == 'WS') %>% 
  ggplot() + ggtitle('WS') +
  geom_line(aes(x = Year, y = Height, color = Site, group = id)) +
  geom_line(data = df_bri, aes(x = t, y = mu), size = 1.5, col = 'red') +
  #  geom_ribbon(data = df_bri, aes(x = t, ymin = ll, ymax = ul), alpha = .5) +
  geom_line(data = df_chSt, aes(x = t, y = mu), size = 1.5, col = 'green') +
  #  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
  geom_line(data = df_liBa, aes(x = t, y = mu), size = 1.5, col = 'blue') +
  geom_line(data = df, aes(x = t, y = mu), size = 1.5, col = 'black')
#  geom_ribbon(data = df_chSt, aes(x = t, ymin = ll, ymax = ul), alpha = .5)
