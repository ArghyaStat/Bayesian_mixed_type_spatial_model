
sim.data <- function(q, N, family, true.beta, true.Sigma, 
                                    true.phi, true.nu, pred.prop){
  
  
  # Validate inputs
  if (length(family) != q) {
    stop("The length of 'family' must be equal to q.")
  }
  
  locations <- matrix(0, nrow = N, ncol = 2)
  
  locations[, 1] <- runif(N)
  locations[, 2] <- runif(N)
  
  X = cbind(1, locations)
  
  p <- ncol(X)
  # Fixed effect 
  mu.true <- X %*% true.beta
  
  # Calculates the Euclidean distance matrix
  
  distmat <- rdist(locations)
  
  # Calculates the correlation matrix from matern kernel
  
  K.true <- fields::Matern(distmat, range = true.phi, smoothness = true.nu)
  
  # Calculates the separable covariance matrix with nugget effect
  
  chol.true.Sigma <- chol(true.Sigma)
  chol.K.true <- chol(K.true)
  
  #Omega <- true.Sigma %x% K.true
  
  # Generating the response vector of dimension Nq*1
  
  W.error <- matrix(rnorm(N*q, 0, 1), nrow = N , ncol = q)
  
  true.W <- mu.true + t(chol.K.true) %*% W.error %*% chol.true.Sigma
  
  Y <- matrix(NA, N, q)
  
  # Generate responses based on family
  for (j in 1:q) {
    if (family[j] == "Gaussian") {
      Y[, j] <- rnorm(N, mean = true.W[, j], sd = 1)
    } else if (family[j] == "Poisson") {
      Y[, j] <- rpois(N, lambda = exp(true.W[, j]))
    } else if (family[j] == "Binomial") {
      Y[, j] <- rbinom(N, size = 1, prob = (exp(true.W[, j]) / (1 + exp(true.W[, j]))))
    } else if (family[j] == "Gamma") {
      Y[, j] <- rgamma(N, shape = exp(true.W[, j]), rate = 1)
    } else if (family[j] == "Negative-Binomial") {
      Y[, j] <- MASS::rnegbin(N, mu = exp(true.W[, j]), theta = 1) # Theta as dispersion parameter
    } else {
      stop(paste("Unsupported family:", family[j]))
    }
  }
  
  # Partition data into observed and predicted
  N.pred <- ceiling(pred.prop * N)

  N.obs <- N - N.pred
  
  pred.indices <- sample(1:N, N.pred)
  pred.locs <- locations[pred.indices, ]
  obs.locs <- locations[-pred.indices, ]
  
  joint.locs <- rbind(pred.locs, obs.locs)  
  
  distmat.obs <- rdist(obs.locs)
  diameter <- max(distmat.obs)
  
  X.obs <- X[-pred.indices,]
  X.pred <- X[pred.indices,]
  
  Y.pred.true <- Y[pred.indices,]
  Y.obs <- Y[-pred.indices,]
  
  true.W.obs <- true.W[-pred.indices,]
  
  data.details <- list(N = N, N.obs = N.obs, N.pred = N.pred, p = p, q = q, 
                       joint.locs = joint.locs, obs.locs = obs.locs, pred.locs = pred.locs,
                       X.obs = X.obs, X.pred = X.pred, diameter = diameter,
                       Y.obs = Y.obs, Y.pred.true = Y.pred.true, true.W.obs = true.W.obs,
                       true.beta = true.beta, true.Sigma = true.Sigma, 
                       true.phi = true.phi, true.nu = true.nu)

  
}
