library(plot3D)

mu.grid <- seq(-5, 5, 0.01)
sigma.grid <- seq(1, 20, 0.1)
mu.sigma.grid <- expand.grid(x = mu.grid, y = sigma.grid)
mu.sigma.grid <- cbind(mu.sigma.grid$x, mu.sigma.grid$y)
out <- apply(mu.sigma.grid, 1, function(mu.sigma){
  integrate(function(x){dnorm(x, mean = mu.sigma[1], sd = mu.sigma[2]) *
      (1 / (1 + exp(-x)))}, upper = Inf, lower = -Inf)$value})

library(plot3D)

scatter2D(x = mu.sigma.grid[ , 1], y = mu.sigma.grid[ , 2], colvar = out,
          pch = 15)

out[mu.sigma.grid[ , 2] == 1]

keep1 <- which((mu.sigma.grid[ , 2] > 0.999)&(mu.sigma.grid[ , 2] < 1.001))
keep1.5 <- which((mu.sigma.grid[ , 2] > 1.499)&(mu.sigma.grid[ , 2] <
                                                  1.501))

plot(mu.sigma.grid[keep1, 1], out[keep1], type = "l")
lines(mu.sigma.grid[keep1.5, 1], out[keep1.5], col = 2)

integrate(function(x){x/(1 + exp(-(sqrt(10) * x))) * dnorm(x)}, lower =
            -Inf, upper = Inf)

out2 <- sapply(seq(0.1, 100, 0.1), function(sig22){
  integrate(function(x){x/(1 + exp(-(sqrt(sig22) * x))) * dnorm(x)},
            lower = -Inf, upper = Inf)$value})
plot(seq(0.1, 100, 0.1), out2, type = "l")


###########

out2 <- sapply(seq(0.1, 100, 0.1), function(sig22){
  integrate(function(x){x/(1 + exp(-(sqrt(sig22) * x))) * dnorm(x)},
            lower = -Inf, upper = Inf)$value})
plot(seq(0.1, 100, 0.1), out2, type = "l")

p <- 3
q <- 2
true.Sigma <- matrix(c(2,1,1,1), nrow = q, ncol = q, byrow = TRUE)

true.beta <- matrix(c(10, 0,  0, 10, 0, 0), nrow = p, ncol = q)

W <- matrix(rnorm(2 * 1e4), 1e4, 2) %*% chol(true.Sigma) - 0

Y1 <- W[ , 1]
Y2 <- rbinom(1e4, 1, prob = 1 / (1 + exp(-W[,2])))

cor(Y1, Y2)

W <- matrix(rnorm(2 * 1e4), 1e4, 2) %*% chol(true.Sigma) + 5

Y1 <- rpois(1e4, lambda = exp(W[ , 1]))
Y2 <- rbinom(1e4, 1, prob = 1 / (1 + exp(-W[,2])))

cor(Y1, Y2)

true.Sigma <- matrix(c(2, 0, 0, sqrt(10)), 2, 2) %*%
  matrix(c(1, 0.9, 0.9, 1), 2, 2) %*% matrix(c(2, 0, 0, sqrt(10)), 2, 2)

true.Sigma <- matrix(c(2,1,1,1), nrow = q, ncol = q, byrow = TRUE)

out3 <- sapply(seq(-5, 5, 0.1), function(mu){
  W <- matrix(rnorm(2 * 1e4), 1e4, 2) %*% chol(true.Sigma) +
    t(replicate(1e4, c(0, mu)))
  
  Y1 <- rpois(1e4, lambda = exp(W[ , 1]))
  Y2 <- rbinom(1e4, 1, prob = 1 / (1 + exp(-W[,2])))
  
  cor(Y1, Y2)})

plot(lowess(seq(-5, 5, 0.1), out3), type = "l")


mu.grid <- seq(-5, 5, 0.01)
sigma.grid <- seq(1, 20, 0.1)
mu.sigma.grid <- expand.grid(x = mu.grid, y = sigma.grid)
mu.sigma.grid <- cbind(mu.sigma.grid$x, mu.sigma.grid$y)
out <- apply(mu.sigma.grid, 1, function(mu.sigma){
  integrate(function(x){dnorm(x, mean = mu.sigma[1], sd = mu.sigma[2]) *
      (1 / (1 + exp(-x)))}, upper = Inf, lower = -Inf)$value})

library(plot3D)

scatter2D(x = mu.sigma.grid[ , 1], y = mu.sigma.grid[ , 2], colvar = out,
          pch = 15)

out[mu.sigma.grid[ , 2] == 1]

keep1 <- which((mu.sigma.grid[ , 2] > 0.999)&(mu.sigma.grid[ , 2] < 1.001))
keep1.5 <- which((mu.sigma.grid[ , 2] > 1.499)&(mu.sigma.grid[ , 2] <
                                                  1.501))

plot(mu.sigma.grid[keep1, 1], out[keep1], type = "l")
lines(mu.sigma.grid[keep1.5, 1], out[keep1.5], col = 2)

integrate(function(x){x/(1 + exp(-(sqrt(10) * x))) * dnorm(x)}, lower =
            -Inf, upper = Inf)

out2 <- sapply(seq(0.1, 100, 0.1), function(sig22){
  integrate(function(x){x/(1 + exp(-(sqrt(sig22) * x))) * dnorm(x)},
            lower = -Inf, upper = Inf)$value})
plot(seq(0.1, 100, 0.1), out2, type = "l")