mixed.workflow.sep <- function(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                               distobs.nn, chol.K.obs.ord.inv, NNarray.obs, 
                               nu, W.obs.ord, beta, Sigma, phi,
                               M.prior, V.prior, S.prior, df.prior, b_phi,
                               niters, tuning.phi,
                               Y.pred.ord.true, X.pred.ord, distlocs.nn, NNarray.all,
                               N.pred, family,
                               family.par = NULL,
                               link = NULL){
  
  ## Storage (same as sep)
  W.obs.ord.chain  <- replicate(niters, matrix(NA, N.obs, q), simplify = FALSE)
  beta.chain       <- replicate(niters, matrix(NA, p, q), simplify = FALSE)
  Sigma.chain      <- replicate(niters, matrix(NA, q, q), simplify = FALSE)
  chol.Sigma.chain <- replicate(niters, matrix(NA, q, q), simplify = FALSE)
  phi.chain        <- rep(NA, niters)
  
  ## Prediction storage
  Y.hat.pred.ord <- replicate(niters, matrix(NA, N.pred, q), simplify = FALSE)
  
  log.pred.sep <- numeric(niters)
  log.pred.comp  <- matrix(NA, niters, q)
  
  ## Initial Sigma quantities
  chol.Sigma <- chol(Sigma)
  Sigma.inv <- chol2inv(chol.Sigma)
  chol.Sigma.inv <- chol(Sigma.inv)
  
  ## Prior quantities
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  chol.prec.V.prior <- chol(prec.V.prior)
  
  B.ind <- which(family == "Binomial")
  fixed_variances <- if(length(B.ind)>0) diag(Sigma)[B.ind] else NULL
  
  acc.phi <- 0
  
  est_time_total  <- 0
  pred_time_total <- 0
  
  est_block_prev  <- 0
  pred_block_prev <- 0
  
  start_total_time <- Sys.time()
  
  for (iter in 1:niters) {
    
    est_start <- Sys.time()
    
    #### Update phi once (same spatial structure) ####
    phi.update.details <- update.phi(
      phi, beta, chol.Sigma.inv, nu,
      W.obs.ord, distobs.nn, NNarray.obs,
      chol.K.obs.ord.inv, acc.phi, tuning.phi, b_phi)
    
    phi <- phi.update.details$phi
    chol.K.obs.ord.inv <- phi.update.details$chol.cormat.obs.inv
    acc.phi <- phi.update.details$acc.phi
    
    #### Loop over components ####
    for (j in 1:q) {
      
      Yj <- matrix(Y.obs.ord[,j], N.obs, 1)
      Wj <- matrix(W.obs.ord[,j], N.obs, 1)
      betaj <- matrix(beta[,j], p, 1)
      
      #### Beta update ####
      K_x <- chol.K.obs.ord.inv %*% X.obs.ord
      prec.V.tilde <- crossprod(K_x) + prec.V.prior
      V.tilde <- chol2inv(chol(prec.V.tilde))
      
      M.tilde.part <- crossprod(K_x, chol.K.obs.ord.inv) %*% Wj +
        prec.V.prior %*% M.prior[,j,drop=FALSE]
      
      chol.V.tilde <- chol(V.tilde)
      
      #### Sigma update (scalar) ####
      Sigma.details <- update.Sigma(
        Wj, chol.K.obs.ord.inv,
        M.prior[,j,drop=FALSE],
        chol.prec.V.prior,
        M.tilde.part,
        chol.V.tilde,
        S.prior[j,j],
        df.prior,
        N.obs,
        family[j],
        fixed_variances[j])
      
      Sigma[j,j] <- Sigma.details$Sigma
      
      chol.Sigma.jj <- chol(Sigma[j,j])
      chol.Sigma.inv.jj <- chol(chol.Sigma.jj)
      
      betaj <- update.beta(
        M.tilde.part,
        chol.Sigma.jj,
        chol.V.tilde
      )
      
      beta[,j] <- betaj
    
      #### W update ####
      Wj <- update.W(Wj, betaj,
                     chol.Sigma.inv.jj,
                     chol.K.obs.ord.inv,
                     Yj, X.obs.ord, N.obs,
                     log.like, family[j],
                     family.par=NULL, link=NULL)
      
      #### assemble ####
      W.obs.ord[,j] <- Wj
     
    }
    
   
    #### Store chains ####
    W.obs.ord.chain[[iter]]  <- W.obs.ord
    beta.chain[[iter]]       <- beta
    Sigma.chain[[iter]]      <- Sigma
    chol.Sigma.chain[[iter]] <- chol(Sigma)
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
    
    #### Prediction (unchanged from sep) ####
    pred_start <- Sys.time()
    
    U.joint <- U.sgv(distlocs.nn, NNarray.all, phi, nu, m)
    
    U.obs.col  <- U.joint[,1:N.obs]
    U.pred.col <- U.joint[,(N.obs+1):(N.obs+N.pred)]
    
    W.pred.prec <- crossprod.spam(U.pred.col)
    W.pred.var  <- solve.spam(as.spam(W.pred.prec))
    chol.W.pred.var <- chol(W.pred.var)
    
    W.pred.ord.temp <- W.pred.ord(
      W.obs.ord, beta, chol.Sigma,
      phi, nu, m,
      X.obs.ord, X.pred.ord,
      U.obs.col, U.pred.col,
      chol.W.pred.var, N.obs, N.pred, q)
    
    log.like.details <- log.likelihood(
      W.pred.ord.temp,
      Y.pred.ord.true,
      family,
      family.par=NULL,
      link=NULL)
    
    log.pred.comp[iter,] <- log.like.details$like.comp
    
    Y.hat.pred.ord[[iter]] <- Y.pred.ord(
      W.pred.ord.temp, N.pred, q, family)
    
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
  
  beta.ess <- compute_ess(beta.chain)
  Sigma.ess <- compute_ess(Sigma.chain)
  phi.ess <- compute_ess(phi.chain)
  W.obs.ord.ess <- compute_ess(W.obs.ord.chain)
  
  
  elpd.comp <- numeric(q)
  for (j in 1:q)
    elpd.comp[j] <- log_mean_exp(log.pred.comp[,j])
  
  elpd.sep <- sum(elpd.comp)
  
  pred.summary <- summary_stats_sep(Y.hat.pred.ord,Y.pred.ord.true)
  pred.coverage <- pred_coverage(Y.pred.ord.true,Y.hat.pred.ord)
  
  beta.stats <- summary_stats_sep(beta.chain, data$true.beta)
  Sigma.stats <- summary_stats_sep(Sigma.chain, data$true.Sigma)
  phi.stats <- summary_stats_sep(phi.chain, data$true.phi)
  
  output <- list("beta.stats" = beta.stats,
                 "Sigma.stats" = Sigma.stats, 
                 "phi.stats" = phi.stats,
                 "beta.ess" = beta.ess,
                 "Sigma.ess" = Sigma.ess,
                 "phi.ess" = phi.ess,
                 "W.obs.ord.ess" = W.obs.ord.ess,
                 "elpd" = elpd.sep,
                 "elpd.comp" = elpd.comp,
                 "pred.summary" = pred.summary,
                 "pred.coverage" = pred.coverage,
                 "acc.phi" = acc.phi/niters,
                 "est.time" = est_time_total,
                 "pred.time" = pred_time_total,
                 "total.time" = total_time)
  
  return(output)
}