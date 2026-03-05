#### Mixed effect model inference ###

mixed.workflow.joint <- function(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                           distobs.nn, chol.K.obs.ord.inv, NNarray.obs, 
                           nu, W.obs.ord, beta, Sigma, phi,
                           M.prior, V.prior, S.prior, df.prior, b_phi,
                           niters, tuning.phi,
                           Y.pred.ord.true, X.pred.ord, distlocs.nn, NNarray.all,
                           N.pred, family,
                           family.par = NULL,
                           link = NULL){
  
  ## Storage
  W.obs.ord.chain  <- replicate(niters, matrix(NA, N.obs, q), simplify = FALSE)
  beta.chain       <- replicate(niters, matrix(NA, p, q), simplify = FALSE)
  Sigma.chain      <- replicate(niters, matrix(NA, q, q), simplify = FALSE)
  chol.Sigma.chain <- replicate(niters, matrix(NA, q, q), simplify = FALSE)
  phi.chain        <- rep(NA, niters)
  
  ## Storage for prediction
  
  Y.hat.pred.ord <- replicate(niters, matrix(NA, N.pred, q), simplify = F)
  
  log.pred.joint <- numeric(niters)
  log.pred.comp <-  matrix(NA, nrow = niters, ncol = q)
 
  
  ## Initial Sigma quantities
  chol.Sigma <- chol(Sigma)
  Sigma.inv <- chol2inv(chol.Sigma)
  chol.Sigma.inv <- chol(Sigma.inv)
  
  ## Prior quantities
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  chol.prec.V.prior <- chol(prec.V.prior)
  
  B.ind <- which(family == "Binomial")
  
  if (length(B.ind) > 0) {
    
    fixed_variances <- diag(Sigma)[B.ind]
    
  } else {
    
    fixed_variances <- NULL
  }
  
  acc.phi <- 0
 
  
  est_time_total  <- 0
  pred_time_total <- 0
  
  est_block_prev  <- 0
  pred_block_prev <- 0
  
  start_total_time <- Sys.time()
  
  for (iter in 1:niters) {
    
    est_start <- Sys.time()
    
    ## Update phi
  
    phi.update.details <- update.phi(
      phi, beta, chol.Sigma.inv, nu,
      W.obs.ord, distobs.nn, NNarray.obs,
      chol.K.obs.ord.inv, acc.phi, tuning.phi, b_phi
    )
    
    phi <- phi.update.details$phi
    chol.K.obs.ord.inv <- phi.update.details$chol.cormat.obs.inv
    acc.phi <- phi.update.details$acc.phi
    
 
    ## Beta posterior quantities
   
    K_x <- chol.K.obs.ord.inv %*% X.obs.ord
    K_x.proj <- crossprod(K_x)
    
    prec.V.tilde <- K_x.proj + prec.V.prior
    V.tilde <- chol2inv(chol(prec.V.tilde))
    
    M.tilde.part <-
      crossprod(K_x, chol.K.obs.ord.inv) %*% W.obs.ord +
      prec.V.prior %*% M.prior
    
    chol.V.tilde <- chol(V.tilde)
   
    ## Update Sigma (restricted IW)
   
    Sigma.update.details <- update.Sigma(
      W.obs.ord          = W.obs.ord,
      chol.K.obs.ord.inv = chol.K.obs.ord.inv,
      M.prior            = M.prior,
      chol.prec.V.prior  = chol.prec.V.prior,
      M.tilde.part       = M.tilde.part,
      chol.V.tilde       = chol.V.tilde,
      S.prior            = S.prior,
      df.prior           = df.prior,
      N.obs              = N.obs,
      family             = family,
      fixed_variances    = fixed_variances
    )
    
    Sigma <- Sigma.update.details$Sigma
    chol.Sigma <- Sigma.update.details$chol.Sigma
    chol.Sigma.inv <- Sigma.update.details$chol.Sigma.inv
    
    ## Update beta
    beta <- update.beta(M.tilde.part, chol.Sigma, chol.V.tilde)
    
  
    ## Update W
    W.obs.ord <- update.W(W.obs.ord, beta, chol.Sigma.inv, chol.K.obs.ord.inv,
                          Y.obs.ord, X.obs.ord, N.obs, log.like, family,
                          family.par = NULL,
                          link = NULL)
    
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
    
    ## Save chain
    W.obs.ord.chain[[iter]]  <- W.obs.ord
    beta.chain[[iter]]       <- beta
    Sigma.chain[[iter]]      <- Sigma
    chol.Sigma.chain[[iter]] <- chol.Sigma
    phi.chain[iter]          <- phi
    
    
    pred_start <- Sys.time()
    
    ## Recompute Vecchia factor for current phi
    U.joint <- U.sgv(distlocs.nn, NNarray.all,
                     phi.chain[iter],
                     nu, m)
    
    U.obs.col  <- U.joint[, 1:N.obs]
    U.pred.col <- U.joint[, (N.obs+1):(N.obs+N.pred)]
    
    ## Predictive covariance of W*
    W.pred.prec <- crossprod.spam(U.pred.col)
    W.pred.var  <- solve.spam(as.spam(W.pred.prec))
    chol.W.pred.var <- chol(W.pred.var)
    
    W.pred.ord.temp <-  W.pred.ord(W.obs.ord = W.obs.ord.chain[[iter]],
                                   beta =  beta.chain[[iter]],
                                   chol.Sigma = chol.Sigma.chain[[iter]],
                                   phi = phi.chain[iter], nu, m,
                                   X.obs.ord, X.pred.ord, U.obs.col, U.pred.col,
                                   chol.W.pred.var, N.obs, N.pred, q)
    
    rm(U.joint)
    rm(U.obs.col)
    rm(U.pred.col)
    
    log.like.details <- log.likelihood(W = W.pred.ord.temp,
                                       Y = Y.pred.ord.true,
                                       family,
                                       family.par = NULL,
                                       link = NULL)
    
    log.pred.joint[iter] <- log.like.details$like.total
    log.pred.comp[iter,] <- log.like.details$like.comp
    
    Y.hat.pred.ord[[iter]] <- Y.pred.ord(W.pred.ord = W.pred.ord.temp, 
                                         N.pred = N.pred,
                                         q = q,
                                         family = family,
                                         family.par = NULL,
                                         link = NULL)
    
    pred_end <- Sys.time()
    
    pred_time_total <- pred_time_total +
      as.numeric(difftime(pred_end, pred_start, units = "mins"))
    
    if (iter %% ceiling(niters / 10) == 0) {
      
      pred_block_time <- pred_time_total - pred_block_prev
      
      cat(sprintf(
        "Prediction Progress: %3.0f%% | Last 10%% Pred time: %.2f mins\n",
        100 * iter / niters,
        pred_block_time
      ))
      
      pred_block_prev <- pred_time_total
    }
   
  }
  
  
  end_total_time <- Sys.time()
  
  total_time <- as.numeric(
    difftime(end_total_time, start_total_time, units = "mins")
  )
  
  est_time_total  <- est_time_total
  pred_time_total <- pred_time_total
  total_time      <- total_time     
  
  cat(sprintf("Estimation time (minutes): %.3f\n", est_time_total))
  cat(sprintf("Prediction time (minutes): %.3f\n", pred_time_total))
  cat(sprintf("Total time (minutes): %.3f\n", total_time))
  
  
  cat(sprintf("Acceptance.phi = %.3f\n", acc.phi / niters))

  
  #### ESS calculation for component chains ####
  
  beta.ess <- compute_ess(beta.chain)
  Sigma.ess <- compute_ess(Sigma.chain)
  phi.ess <- compute_ess(phi.chain)
  W.obs.ord.ess <- compute_ess(W.obs.ord.chain)
  
  # Joint ELPD
  elpd.joint <- log_mean_exp(log.pred.joint)
  
  # Precompute marginal component ELPDs
  # (each column gives predictive log density per iteration)
  elpd.comp <- numeric(q)
  
  for (j in seq_len(q)) {
    elpd.comp[j] <- log_mean_exp(log.pred.comp[, j])
  }
  
  # Total possible subsets
  q.subsets <- 2^q
  
  # We exclude empty (0) and full (2^q - 1)
  n.valid <- q.subsets - 2
  
  cond.elpd <- numeric(n.valid)
  
  counter <- 1
  
  for (i in 1:(q.subsets - 2)) {
    
    ind <- as.logical(intToBits(i)[1:q])
    
    elpd.B <- sum(elpd.comp[ind])
    
    cond.elpd[counter] <- elpd.joint - elpd.B
    
    counter <- counter + 1
  }
  
  
  pred.summary <- summary_stats( Y.hat.pred.ord, Y.pred.ord.true)
  pred.coverage <- pred_coverage(Y.pred.ord.true,  Y.hat.pred.ord)
  

beta.stats <- summary_stats(beta.chain, data$true.beta)
Sigma.stats <- summary_stats(Sigma.chain, data$true.Sigma)
phi.stats <- summary_stats(phi.chain, data$true.phi)

output <- list("beta.stats" = beta.stats,
               "Sigma.stats" = Sigma.stats, 
               "phi.stats" = phi.stats,
               "beta.ess" = beta.ess,
               "Sigma.ess" = Sigma.ess,
               "phi.ess" = phi.ess,
               "W.obs.ord.ess" = W.obs.ord.ess,
               "elpd" = elpd.joint,
               "cond.elpd" = cond.elpd,
               "pred.summary" = pred.summary,
               "pred.coverage" = pred.coverage,
               "acc.phi" = acc.phi/niters,
               "est.time" = est_time_total,
               "pred.time" = pred_time_total,
               "total.time" = total_time)

return(output)
  
}