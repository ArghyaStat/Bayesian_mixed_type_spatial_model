W.pred.ord <- function(W.obs.ord, beta, chol.Sigma, phi, nu, m,
                       X.obs.ord, X.pred.ord, U.joint, U.pred.col, U.obs.col,
                       chol.W.pred.var, N.obs, N.pred, q, p){
  
  W.pred.mean.part <- U.obs.col %*% (W.obs.ord - X.obs.ord %*% beta)
  W.pred.mean.tilde <- crossprod.spam(U.pred.col,  W.pred.mean.part)
  
  W.pred.samp <- (matrix(rnorm(N.pred*q, 0, 1), nrow = N.pred , ncol = q) %*% chol.Sigma)
  
  W.pred.temp <- chol.W.pred.var %*% W.pred.mean.tilde +  W.pred.samp
  
  W.pred.ord <- X.pred.ord %*% beta -  crossprod(chol.W.pred.var,  W.pred.temp)
  
  W.pred.ord
  
}

Y.pred.ord <- function(W.pred.ord, N.pred, q, family){
  
  # Validate inputs
  if (length(family) != q) {
    stop("The length of 'family' must be equal to q.")
  }
  
  Y.mat <- matrix(NA, nrow = N.pred, ncol = q)
  
  # Validate inputs
  if (length(family) != q) {
    stop("The length of 'family' must be equal to q.")
  }
  
  # Generate responses based on family
  for (j in 1:q) {
    if (family[j] == "Gaussian") {
      Y.mat[, j] <- rnorm(N.pred, mean = W.pred.ord[, j], sd = 1)
    } else if (family[j] == "Poisson") {
      Y.mat[, j] <- rpois(N.pred, lambda = exp(W.pred.ord[, j]))
    } else if (family[j] == "Binomial") {
      Y.mat[, j] <- rbinom(N.pred, size = 1, prob = (exp(W.pred.ord[, j]) / (1 + exp(W.pred.ord[, j]))))
    } else if (family[j] == "Gamma") {
      Y.mat[, j] <- rgamma(N.pred, shape = exp(W.pred.ord[, j]), rate = 1)
    } else if (family[j] == "Negative-Binomial") {
      Y.mat[, j] <- MASS::rnegbin(N.pred, mu = exp(W.pred.ord[, j]), theta = 1) # Theta as dispersion parameter
    } else {
      stop(paste("Unsupported family:", family[j]))
    }
  }
  
  Y.mat
}

predictive.output <- function(W.obs.ord.samples, beta.samples, chol.Sigma.samples,
                               phi.samples, nu, m,
                               X.obs.ord, X.pred.ord, distlocs.nn, NNarray.all,
                               N.obs, N.pred, pred.iters, q, p){
  
  Y.pred.ord <- replicate(pred.iters, matrix(0, N.pred, q), simplify = F)
  
  log.pred.joint <- numeric(pred.iters)
  
  start_time <- Sys.time()
  
  for(i in 1:pred.iters){
    
    ## Recompute Vecchia factor for current phi
    U.joint <- U.sgv(distlocs.nn, NNarray.all,
                     phi.post[[i]],
                     nu, m)
    
    U.obs.col  <- U.joint[, 1:N.obs]
    U.pred.col <- U.joint[, (N.obs+1):(N.obs+N.pred)]
    
    ## Predictive covariance of W*
    W.pred.prec <- crossprod.spam(U.pred.col)
    W.pred.var  <- solve.spam(as.spam(W.pred.prec))
    chol.W.pred.var <- chol(W.pred.var)
    
    W.pred.ord.temp <-  W.pred.ord(W.obs.ord = W.obs.ord.samples[[i]],
                                   beta =  beta.samples[[i]],
                                   chol.Sigma = chol.Sigma.samples[[i]],
                                   phi = phi.samples[[i]], nu, m,
                                   X.obs.ord, X.pred.ord, U.joint, U.pred.col, U.obs.col,
                                   chol.W.pred.var, N.obs, N.pred, q, p)
    
    log.pred.joint[i] <- log.likelihood(W = W.pred.ord.temp,
                                        Y = Y.pred.ord.true,
                                        family = family)
    
    Y.pred.ord[[i]] <- Y.pred.ord(W.pred.ord.temp, N.pred, q, family)
    
    if (i %% (pred.iters / 10) == 0) {
      elapsed_time <- as.numeric(
        difftime(Sys.time(), start_time, units = "mins")
      )
      cat(sprintf("ELPD progress: %.0f%% | %.2f minutes\n",
                  100 * i / pred.iters,
                  elapsed_time))
    }
  }
  
  end_time <- Sys.time()
  
  pred.time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("Total ELPD computation time (minutes):",
      pred.time,
      "\n")
  
  results <- list(Y.pred.ord =  Y.pred.ord,
                  log.pred.joint =  log.pred.joint,
                  pred.time = pred.time)
  
}
