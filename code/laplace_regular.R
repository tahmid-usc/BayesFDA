## laplace kernel 

ker1 <- function(x, theta, jitter = .00000001) {
  m <- length(x)
  k1ker <- laplacedot(sigma = 1/theta[1])
  k1 <- theta[2] * kernelMatrix(k1ker, x = x)
  k1 <- k1 + jitter * diag(m)
  return(k1)
}


ker2 <- function(x, theta) {
  k2ker <- laplacedot(sigma = 1/theta[3]) 
  k2 <- theta[4] * kernelMatrix(k2ker, x = x) 
  return(k2)
}


covker <- function(x, y, theta) {
  rbf <- laplacedot(sigma = 1/theta[1])
  covker <- theta[2] * kernelMatrix(rbf, x = x, y = y)
  return(covker)
}