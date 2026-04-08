mixed.workflow.joint <- function(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                                 nu, W.obs.ord, beta, Sigma, phi,
                                 M.prior, V.prior, S.prior, df.prior, b_phi,
                                 niters, tuning.phi,
                                 distobs.nn, coeff.obs, NNarray.obs,
                                 family, family.par = NULL, link = NULL){
  
  ## Storage
  W.obs.ord.chain  <- replicate(niters, matrix(NA, N.obs, q), simplify = FALSE)
  beta.chain       <- replicate(niters, matrix(NA, p, q), simplify = FALSE)
  Sigma.chain      <- replicate(niters, matrix(NA, q, q), simplify = FALSE)
  phi.chain        <- rep(NA, niters)
  
  log.post.W   <- numeric(niters)
  log.post.full <- numeric(niters)
  
  ## Initial quantities
  chol.Sigma <- chol(Sigma)
  
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  
  B.ind <- which(family == "Binomial")
  
  if (length(B.ind) > 0) {
    
    fixed_variances <- diag(Sigma)[B.ind]
    
  } else {
    
    fixed_variances <- NULL
  }
  
  acc.phi <- 0
  
  est_time_total  <- 0
  est_block_prev  <- 0

  
  start_total_time <- Sys.time()
  
  for (iter in 1:niters) {
    
    est_start <- Sys.time()
 
    ## --- phi update ---
    phi.update.details <- update.phi(phi, beta, chol.Sigma, nu,
                                     W.obs.ord, distobs.nn, distobs.nn,
                                     NNarray.obs, NNarray.obs,
                                     coeff.obs, acc.phi, tuning.phi, b_phi)
    
    phi <- phi.update.details$phi
    coeff.obs <- phi.update.details$coeff.all   ## IMPORTANT: now directly obs
    acc.phi <- phi.update.details$acc.phi
    
    B <- coeff.obs$B
    r <- coeff.obs$r
    
    ## --- Beta update ---
    vecchia.prod.out <- vecchia_products(W.obs.ord, X.obs.ord, B, r, NNarray.obs)
    UX <- vecchia.prod.out$UX
    UW <- vecchia.prod.out$UW
    
    prec.V.tilde <- crossprod(UX) + prec.V.prior
    chol.prec.V.tilde <- chol(prec.V.tilde)
    V.tilde <- chol2inv(chol.prec.V.tilde)
    chol.V.tilde <- chol(V.tilde)
    
    M.tilde.part <- crossprod(UX, UW) + prec.V.prior %*% M.prior
    
    ## --- Sigma update ---
    Sigma.update.details <- update.Sigma(UW = UW, 
                                         M.prior = M.prior,
                                         chol.prec.V.prior = chol(prec.V.prior),
                                         M.tilde.part = M.tilde.part,
                                         chol.V.tilde = chol.V.tilde,
                                         S.prior = S.prior,
                                         df.prior = df.prior,
                                         N.obs = N.obs,
                                         family = family,
                                         fixed_variances = fixed_variances)
    
    Sigma <- Sigma.update.details$Sigma
    chol.Sigma <- Sigma.update.details$chol.Sigma
    
    ## --- Beta ---
    beta <- update.beta(M.tilde.part, chol.Sigma, chol.V.tilde)
    
    ## --- W update ---
    W.obs.ord <- update.W(W.obs.ord, beta, Sigma, B, r, NNarray.obs,
                          Y.obs.ord, X.obs.ord, N.obs, log.like, family,
                          family.par = NULL,
                          link = NULL)
    
    ## ===============================
    ## LOG POSTERIOR
    ## ===============================
    
    log.like.val <- log.likelihood(W.obs.ord, Y.obs.ord, family,
                                   family.par = NULL,
                                   link = NULL)$like.total
    
    ## Vecchia prior
    vecchia.prod.out <- vecchia_products(W.obs.ord, X.obs.ord, B, r, NNarray.obs)
    UW <- vecchia.prod.out$UW
    
    R.Sigma <- chol(Sigma)
    
    Sigma.inv <- chol2inv(R.Sigma)
    log.det.Sigma <- 2 * sum(log(diag(R.Sigma)))
    
    log.prior.Sigma <- -(df.prior + q + 1)/2 * log.det.Sigma -
      0.5 * sum(diag(S.prior %*% Sigma.inv))
  
    
    log.prior.W <- vecchia_ll(W.ord = W.obs.ord,
                              M.ord = X.obs.ord %*% beta,
                              NNarray = NNarray.obs, 
                              coeff = coeff.obs, 
                              chol.V = chol.Sigma, log = TRUE)
    
    E <- beta - M.prior
    
    tmp <- forwardsolve(chol.V.prior, E)
    
    Sigma.inv <- chol2inv(chol.Sigma)
    
    quad.beta <- sum((tmp %*% Sigma.inv) * tmp)
    
    log.prior.beta <- -0.5 * quad.beta
    
    ## Phi prior
    log.prior.phi <- ifelse(phi > 0 & phi < b_phi, 0, -Inf)
    
    log.post.W[iter] <- log.like.val + log.prior.W
    
    log.post.full[iter] <- log.post.W[iter] + log.prior.beta + log.prior.Sigma + log.prior.phi
    
    ## Save
    W.obs.ord.chain[[iter]]  <- W.obs.ord
    beta.chain[[iter]]       <- beta
    Sigma.chain[[iter]]      <- Sigma
    phi.chain[iter]          <- phi
    
    est_end <- Sys.time()
    
    est_time_total <- est_time_total +
      as.numeric(difftime(est_end, est_start, units = "mins"))
    
    
    if (iter %% ceiling(niters / 10) == 0) {
      
      est_block_time <- est_time_total - est_block_prev
      
      cat(sprintf(
        "Estimation Progress: %3.0f%% | Last 10%% Est time: %.2f mins\n",
        100 * iter / niters,
        est_block_time
      ))
      
      est_block_prev <- est_time_total
    }
  }
  
  est_time_total  <- est_time_total
  cat(sprintf("Estimation time (minutes): %.3f\n", est_time_total))
 
  cat(sprintf("Acceptance.phi = %.3f\n", acc.phi / niters))
  
  beta.ess = compute_ess(beta.chain)
  Sigma.ess = compute_ess(Sigma.chain)
  phi.ess = compute_ess(phi.chain)
  W.obs.ord.ess = compute_ess(W.obs.ord.chain)
  
  acf.log.post.W = acf(log.post.W, plot = FALSE)
  acf.log.post.full = acf(log.post.full, plot = FALSE)
 
  
  ## Diagnostics
  output <- list("beta.ess" = beta.ess,
                 "Sigma.ess" = Sigma.ess,
                 "phi.ess" = phi.ess,
                 "W.obs.ord.ess" = W.obs.ord.ess,
                 "log.post.W" = log.post.W,
                 "log.post.full" = log.post.full,
                 "acf.log.post.W" = acf.log.post.W,
                 "acf.log.post.full" =  acf.log.post.full,
                 "acc.phi" = acc.phi/niters,
                 "est.time" = est_time_total)
  
  return(output)
}