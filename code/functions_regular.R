# Extract relevant feature of the data

feature <- function(train) {
  
  n <- dim(train)[1]
  m <- dim(train)[2]
  x <- seq(0,1, length.out = m)
  
  #standardization
  ymean <- mean(as.matrix(train))
  ysd <- sd(as.matrix(train))
  
  train <- (train - ymean) / ysd
  #theta <- Hyper(train[,2:19])
  theta <- Hyper.ms(train)
  
  k1 <- ker1(x, theta)
  k2 <- ker2(x, theta)
  
  k1inv <- chol2inv(chol(k1))
  D <- chol2inv(chol(k2 + theta[5] * diag(m)))
  C <- chol2inv(chol(k1inv + n * D))
  
  ySum <- as.matrix(apply(train, 2, sum), ncol = 1)
  
  mumat <- ( diag(m) - n * D %*% C) %*% D #%*% ySum
  
  return(list(n = n,m = m, theta = theta, mumat = mumat, ySum = ySum, ymean = ymean, ysd = ysd))
  
}


Hyper <- function(dtf) {
  
  marlik <- function(theta) {
    
    n <- dim(dtf)[1]
    m <- dim(dtf)[2]
    x <- seq(0,1, length.out = m)
    theta <- theta^2
    
    k1 <- ker1(x, theta)
    k2 <- ker2(x, theta)
    
    k1inv <- chol2inv(chol(k1))
    G <- k2 + theta[5] * diag(m)
    D <- chol2inv(chol(G))
    C <- chol2inv(chol(k1inv + n * D))
    ySum <- as.matrix(apply(dtf, 2, sum), ncol = 1)
    P <- D %*% ySum
    YD <- sum(apply(dtf, 1, function(x) { x <- matrix(x, ncol = 1); t(x) %*% D %*% x}))
    logl <- YD - (t(P) %*% C %*% P) + determinant((k1inv + n * D), logarithm = T)$modulus + determinant(k1, logarithm = T)$modulus + 
      n * determinant(G, logarithm = T)$modulus
    
    return(logl)
  }
  
  hyp <- optim(par=rep(.01, 5), fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 10000))
  print(hyp)
  return(hyp$par^2)
  
}


Hyper.ms <- function(dtf) {
  
  marlik <- function(theta) {
    
    n <- dim(dtf)[1]
    m <- dim(dtf)[2]
    x <- seq(0,1, length.out = m)
    theta <- theta^2
    
    k1 <- ker1(x, theta)
    k2 <- ker2(x, theta)
    
    k1inv <- chol2inv(chol(k1))
    G <- k2 + theta[5] * diag(m)
    D <- chol2inv(chol(G))
    C <- chol2inv(chol(k1inv + n * D))
    ySum <- as.matrix(apply(dtf, 2, sum), ncol = 1)
    P <- D %*% ySum
    YD <- sum(apply(dtf, 1, function(x) { x <- matrix(x, ncol = 1); t(x) %*% D %*% x}))
    logl <- YD - (t(P) %*% C %*% P) + determinant((k1inv + n * D), logarithm = T)$modulus + determinant(k1, logarithm = T)$modulus + 
      n * determinant(G, logarithm = T)$modulus
    
    return(logl)
  }
  
  parmat <- matrix(rexp(50, 1), ncol = 5)
  hyp <- multistart(parmat, fn = marlik, method='Nelder-Mead', control = list(maxit = 10000))
  hyp <- as.numeric(hyp[which(hyp$value == min(hyp$value)), 1:5])^2
  return(hyp)
  
}

#posterior distribution

postDist <- function(y, train) {
  
  ns <- length(y)
  x <- seq(0,1, length.out = train$m)
  k <- covker(y, x, train$theta)
  kxx <- ker1(y, train$theta)
  pred <- k %*% train$mumat %*% train$ySum
  sigmastar <- kxx + train$theta[5] * diag(ns) - train$n * (k %*% train$mumat %*% t(k))
  sigmastar <- forceSymmetric(sigmastar)
  sigmastar <- as.matrix(sigmastar)
  ll <- pred - 1.96 * sqrt(diag(sigmastar))
  ul <- pred + 1.96 * sqrt(diag(sigmastar))
  
  
  mu <- pred * train$ysd + train$ymean
  ul <- ul * train$ysd + train$ymean
  ll <- ll * train$ysd + train$ymean
  
  return(list(mu = mu, sigma = sigmastar, ll = ll, ul = ul))
  
}

