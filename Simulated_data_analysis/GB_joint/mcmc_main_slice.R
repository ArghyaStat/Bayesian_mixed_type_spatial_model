###### MCMC main function ######

mcmc.main <- function(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                      distobs.nn, chol.K.obs.ord.inv, NNarray.obs, 
                      family, nu, W.obs.ord, beta, Sigma, phi,
                      M.prior, V.prior, S.prior, df.prior, b_phi,
                      niters, tuning.phi){
  
  
  W.obs.ord.chain <- replicate(niters, matrix(NA, N.obs, q), simplify = F)
  beta.chain <- replicate(niters, matrix(NA, p, q), simplify = F)
  Sigma.chain <- replicate(niters, matrix(NA, q, q), simplify = F)
  chol.Sigma.chain <- replicate(niters, matrix(NA, q, q), simplify = F)
  phi.chain <- rep(NA, niters)
  
  chol.Sigma <- chol(Sigma)
  Sigma.inv <- chol2inv(chol.Sigma)
  chol.Sigma.inv <- chol(Sigma.inv)
  
  
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  chol.prec.V.prior <- chol(prec.V.prior)
  
  # acceptance rate for each metropolis update parameters
  
  acc.phi <- 0
  log.like <- 0
  
  start_time <- Sys.time()  # Start time using Sys.time()
  
  for (iter in 1:niters) {
    
    # Update phi
    phi.update.details <- update.phi(phi, beta, chol.Sigma.inv, nu,
                                     W.obs.ord, distobs.nn, NNarray.obs,
                                     chol.K.obs.ord.inv, acc.phi, tuning.phi, b_phi)
    
    phi <- phi.update.details$phi
    chol.K.obs.ord.inv <- phi.update.details$chol.cormat.obs.inv
    acc.phi <- phi.update.details$acc.phi
    
    
    K_x <- chol.K.obs.ord.inv %*% X.obs.ord
    
    K_x.proj <- crossprod(K_x)
    
    prec.V.tilde <- (K_x.proj + prec.V.prior)
    V.tilde <- chol2inv(chol(prec.V.tilde))
    M.tilde.part <-  crossprod(K_x, chol.K.obs.ord.inv) %*% W.obs.ord + (prec.V.prior %*% M.prior)
    
    chol.V.tilde <- chol(V.tilde)
    
    Sigma.update.details <- update.Sigma(W.obs.ord, X.obs.ord, chol.K.obs.ord.inv, NNarray.obs,
                                         M.prior, chol.prec.V.prior, M.tilde.part, chol.V.tilde,
                                         S.prior, df.prior, N.obs)

    Sigma <- Sigma.update.details$Sigma
    chol.Sigma.inv <- Sigma.update.details$chol.Sigma.inv
    chol.Sigma <- Sigma.update.details$chol.Sigma

    beta <- update.beta(M.tilde.part, chol.Sigma, chol.V.tilde)
    
    W.update.details <- update.W(W.obs.ord, beta, chol.Sigma.inv, chol.K.obs.ord.inv,
                                 Y.obs.ord, X.obs.ord, N.obs, log.like, family)
    
    W.obs.ord <- W.update.details$W.obs.ord
    log.like <- W.update.details$log.like
    
    
    if(iter %% (niters / 10) == 0) {
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      cat(sprintf("Main Phase Progress: %.0f%%, Elapsed time: %.2f minutes\n", 
                  100 * (iter / niters), elapsed_time))
    }
    
    # Save chain
    W.obs.ord.chain[[iter]] <- W.obs.ord
    beta.chain[[iter]] <- beta
    Sigma.chain[[iter]] <- Sigma
    phi.chain[iter] <- phi
    chol.Sigma.chain[[iter]] <- chol.Sigma
    
  }
  
  # Adjust Sigma for binomial family after storing chains
  for (iter in 1:niters) {
    for (j in seq_len(q)) {
      if (family[j] == "Binomial") {
        Sigma.chain[[iter]][j, ] <- Sigma.chain[[iter]][j, ] / sqrt(Sigma.chain[[iter]][j, j])
        Sigma.chain[[iter]][, j] <- Sigma.chain[[iter]][j, ]
        Sigma.chain[[iter]][j, j] <- 1
      }
    }
    chol.Sigma.chain[[iter]] <- chol(Sigma.chain[[iter]])
  }
  
  end_time <- Sys.time()
  total_elapsed_time <- difftime(end_time, start_time, units = "mins")[[1]] # in minutes
  
  theta.chain <- list("W.ord.samples" = W.obs.ord.chain,
                      "beta.samples" = beta.chain,
                      "Sigma.samples" = Sigma.chain,
                      "chol.Sigma.samples" =  chol.Sigma.chain,
                      "phi.samples" = phi.chain,
                      "acceptance.phi" = acc.phi/niters,
                      "log.like.mean" = log.like/(q*niters),
                      "total_time" = total_elapsed_time)
  
  # Calculate acceptance rate

  print(paste("Acceptance.phi = ", acc.phi/(niters)))
  
  return(theta.chain)
  
}