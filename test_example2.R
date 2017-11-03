x
data <- IV.model.generate()
gmm <- function(x, a0) {
  obj <- function(x, theta) {
    return(mean(x^2) - theta)
  }
  grad <- function(x, theta) {
    return(-1)
  }
  
  f <- Inf
  df <- Inf
  tol.err <- 1e-9
  theta <- theta0
  while (abs(f) > tol.err) {
    f <- obj(x, theta)
    df <- grad(x, theta)
    theta <- theta - f/df
  }
  return(theta)
}
trials <- 4000
out <- matrix(0, nrow = trials, ncol = 10)
for (ind in 1:10) {
  samplesize <- ind * 10
  for (tmp in 1:trials) {
    x <- rnorm(samplesize, 0, 1)
    theta1 <- gmm(x, 0)
    out[tmp, ind] <- theta1
  }
}
mae <- (colMeans(out) - 1)
plot(mae)


