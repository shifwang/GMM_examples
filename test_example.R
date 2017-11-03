gmm <- function(x, theta0) {
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
el <- function(x, theta0) {
  f <- function(x, theta, lambda) mean(log(1 - lambda * (x^2 - theta)))
  df_theta <- function(x, theta, lambda) mean(1/(1 - lambda * (x^2 - theta)))*lambda 
  df_lambda <- function(x, theta, lambda) mean(-(x^2 - theta)/(1 - lambda * (x^2 - theta)))
  ddf_theta <- function(x, theta, lambda) mean(1/(1 - lambda * (x^2 - theta))^2)*lambda^2 
  ddf_lambda <- function(x, theta, lambda) mean(-(x^2 - theta)^2/(1 - lambda * (x^2 - theta))^2)
  tol.err <- 1e-9
  theta <- theta0
  lambda <- 0
  while (abs(df_lambda(x, theta, lambda)) > tol.err) {
   lambda <- lambda - df_lambda(x, theta, lambda)/ddf_lambda(x, theta, lambda) 
  }
  obj <- Inf
  obj_next <- f(x, theta, lambda)
  while(obj_next < obj - tol.err) {
    obj <- obj_next
    while (abs(df_theta(x, theta, lambda)) > tol.err) {
      theta <- theta - df_theta(x, theta, lambda)/ddf_theta(x, theta, lambda) 
    }
    while (abs(df_lambda(x, theta, lambda)) > tol.err) {
      lambda <- lambda - df_lambda(x, theta, lambda)/ddf_lambda(x, theta, lambda) 
    }
    obj_next <- f(x, theta, lambda)
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


