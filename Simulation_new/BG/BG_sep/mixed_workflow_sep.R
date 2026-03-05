mixed.workflow.sep <- function(
    Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
    distobs.nn, chol.K.obs.ord.inv, NNarray.obs,
    nu, W.obs.ord, beta, Sigma, phi,
    M.prior, V.prior, S.prior, df.prior, b_phi,
    niters, tuning.phi,
    Y.pred.ord.true, X.pred.ord, distlocs.nn, NNarray.all,
    N.pred, family,
    family.par = NULL,
    link = NULL
) {
  
  #### Storage ####
  beta.stats.list  <- vector("list", q)
  Sigma.stats.list <- vector("list", q)
  phi.stats.list   <- vector("list", q)
  
  beta.ess  <- matrix(NA, p, q)
  Sigma.ess <- rep(NA, q)
  phi.ess   <- rep(NA, q)
  W.obs.ord.ess <- matrix(NA, N.obs, q)
  
  elpd <- 0
  est_time_total  <- 0
  pred_time_total <- 0
  
  start_total_time <- Sys.time()
  
  for (j in seq_len(q)) {
    
    cat(sprintf("\n--- Separate Model: Component %d ---\n", j))
    
    q_sub <- 1
    sub_family <- family[j]
    
    #### Extract component-specific objects ####
    Yj.obs.ord  <- as.matrix(Y.obs.ord[, j])
    Wj.obs.ord  <- as.matrix(W.obs.ord[, j])
    betaj       <- as.matrix(beta[, j])
    Sigmaj      <- matrix(Sigma[j, j], 1, 1)
    Yj.pred.true <- as.matrix(Y.pred.ord.true[, j])
    
    #### Prior (dimension 1) ####
    M.prior.j  <- matrix(0, p, 1)
    V.prior.j  <- V.prior
    S.prior.j  <- matrix(1)
    df.prior.j <- 2
    
    est_start <- Sys.time()
    
    fit <- mixed.workflow.joint(
      Y.obs.ord = Yj.obs.ord,
      X.obs.ord = X.obs.ord,
      obs.locs.ord = obs.locs.ord,
      m = m,
      N.obs = N.obs,
      p = p,
      q = 1,
      distobs.nn = distobs.nn,
      chol.K.obs.ord.inv = chol.K.obs.ord.inv,
      NNarray.obs = NNarray.obs,
      nu = nu,
      W.obs.ord = Wj.obs.ord,
      beta = betaj,
      Sigma = Sigmaj,
      phi = phi,
      M.prior = M.prior.j,
      V.prior = V.prior.j,
      S.prior = S.prior.j,
      df.prior = df.prior.j,
      b_phi = b_phi,
      niters = niters,
      tuning.phi = tuning.phi,
      Y.pred.ord.true = Yj.pred.true,
      X.pred.ord = X.pred.ord,
      distlocs.nn = distlocs.nn,
      NNarray.all = NNarray.all,
      N.pred = N.pred,
      family = sub_family,
      family.par = family.par,
      link = link
    )
    
    est_time_total  <- est_time_total  + fit$est.time
    pred_time_total <- pred_time_total + fit$pred.time
    
    #### Collect summaries ####
    beta.stats.list[[j]]  <- fit$beta.stats
    Sigma.stats.list[[j]] <- fit$Sigma.stats
    phi.stats.list[[j]]   <- fit$phi.stats
    
    beta.ess[, j]      <- fit$beta.ess
    Sigma.ess[j]       <- fit$Sigma.ess
    phi.ess[j]         <- fit$phi.ess
    W.obs.ord.ess[, j] <- fit$W.obs.ord.ess
    
    elpd <- elpd + fit$elpd
    
    cat(sprintf("Component %d ELPD = %.3f\n", j, fit$elpd))
  }
  
  total_time <- as.numeric(
    difftime(Sys.time(), start_total_time, units = "mins")
  )
  
  cat("\n=====================================\n")
  cat(sprintf("Total Estimation Time: %.3f mins\n", est_time_total))
  cat(sprintf("Total Prediction Time: %.3f mins\n", pred_time_total))
  cat(sprintf("Total Workflow Time: %.3f mins\n", total_time))
  cat("=====================================\n")
  
  out <- list(
    beta.stats = combine_stats(beta.stats.list),
    Sigma.stats = combine_stats(Sigma.stats.list),
    phi.stats = combine_stats(phi.stats.list),
    beta.ess = beta.ess,
    Sigma.ess = Sigma.ess,
    phi.ess = phi.ess,
    W.obs.ord.ess = W.obs.ord.ess,
    elpd = elpd,
    est.time = est_time_total,
    pred.time = pred_time_total,
    total.time = total_time
  )
  
  return(out)
}