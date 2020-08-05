## Matern 3/2 


ker1 <- function(x, theta, jitter = .00000001) {
  m <- length(x)
  r <- as.matrix(dist(x))
  k1 <- theta[2] * (1 + (sqrt(3) * r / theta[1])) * exp(- sqrt(3) * r /theta[1])
  k1 <- k1 + jitter * diag(m)
  return(k1)
}


ker2 <- function(x, theta) {
  m <- length(x)
  r <- as.matrix(dist(x))
  k2 <- theta[4] * (1 + (sqrt(3) * r / theta[3])) * exp(- sqrt(3) * r /theta[3])
  return(k2)
}


covker <- function(x, y, theta) {
  r <- abs(outer(x, y, "-"))
  covker <- theta[2] * (1 + (sqrt(3) * r / theta[1])) * exp(- sqrt(3) * r /theta[1])
  return(covker)
}