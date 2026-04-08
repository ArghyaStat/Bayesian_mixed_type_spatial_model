mixed.workflow.sep <- function(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                               nu, W.obs.ord, beta, Sigma, phi,
                               M.prior, V.prior, a.prior, b.prior, b_phi,
                               niters, tuning.phi,
                               Y.pred.ord.true, X.pred.ord, N.pred,
                               distall.nn, coeff.all, NNarray.all,
                               family, family.par = NULL, link = NULL, burnin){
  
  ## =========================
  ## Storage (POST-BURNIN ONLY)
  ## =========================
  
  W.obs.ord.chain  <- replicate(niters, matrix(NA, N.obs, q), simplify = FALSE)
  beta.chain       <- replicate(niters, matrix(NA, p, q), simplify = FALSE)
  Sigma.chain      <- replicate(niters, matrix(NA, q, q), simplify = FALSE)
  chol.Sigma.chain <- replicate(niters, matrix(NA, q, q), simplify = FALSE)
  phi.chain        <- rep(NA, niters)
  
  ## Prediction storage
  Y.hat.pred.ord <- replicate(niters, matrix(NA, N.pred, q), simplify = FALSE)
  log.pred.sep   <- numeric(niters)
  log.pred.comp  <- matrix(NA, nrow = niters, ncol = q)
  
  ## =========================
  ## Initial quantities
  ## =========================
  
  chol.Sigma <- chol(Sigma)
  
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  chol.prec.V.prior <- chol(prec.V.prior)
  
  NNarray.obs <- NNarray.all[1:N.obs,]
  distobs.nn  <- distall.nn[1:N.obs]
  
  obs.id  <- which(ord <= N.obs)
  pred.id <- which(ord > N.obs)
  
  B.ind <- which(family == "Binomial")
  fixed_variances <- if(length(B.ind) > 0) diag(Sigma)[B.ind] else NULL
  
  acc.phi <- 0
  
  est_time_total  <- 0
  pred_time_total <- 0
  
  start_total_time <- Sys.time()
  
  total_iters <- burnin + niters
  
  ## =========================
  ## MCMC loop
  ## =========================
  
  for (iter in 1:total_iters) {
    
    est_start <- Sys.time()
    
    ## =========================
    ## Update phi
    ## =========================
    
    phi.update.details <- update.phi(phi, beta, chol.Sigma, nu,
                                     W.obs.ord, distobs.nn, distall.nn,
                                     NNarray.obs, NNarray.all,
                                     coeff.all, acc.phi, tuning.phi, b_phi)
    
    phi <- phi.update.details$phi
    coeff.all <- phi.update.details$coeff.all
    acc.phi <- phi.update.details$acc.phi
    
    coeff.obs <- list(
      B = coeff.all$B[1:N.obs],
      r = coeff.all$r[1:N.obs]
    )
    
    B <- coeff.obs$B
    r <- coeff.obs$r
    
    ## =========================
    ## Beta posterior quantities
    ## =========================
    
    vecchia.prod.out <- vecchia_products(W.obs.ord, X.obs.ord, B, r, NNarray.obs)
    
    UX <- vecchia.prod.out$UX
    UW <- vecchia.prod.out$UW
    
    prec.V.tilde <- crossprod(UX) + prec.V.prior
    chol.prec.V.tilde <- chol(prec.V.tilde)
    V.tilde <- chol2inv(chol.prec.V.tilde)
    chol.V.tilde <- chol(V.tilde)
    
    M.tilde.part <- crossprod(UX, UW) + prec.V.prior %*% M.prior
    
    ## =========================
    ## Component-wise updates
    ## =========================
    
    for (j in 1:q) {
      
      Yj <- matrix(Y.obs.ord[,j], N.obs, 1)
      Wj <- matrix(W.obs.ord[,j], N.obs, 1)
      betaj <- matrix(beta[,j], p, 1)
      
      Sigma.details <- update.Sigma(
        UW[,j],
        M.prior[,j,drop=FALSE],
        chol.prec.V.prior,
        M.tilde.part[,j],
        chol.V.tilde,
        a.prior,
        b.prior,
        N.obs,
        family[j],
        fixed_variances[j]
      )
      
      Sigma[j,j] <- Sigma.details$Sigma
      Sigma_jj <- as.matrix(Sigma[j,j])
      chol.Sigma.jj <- matrix(sqrt(Sigma[j,j]), 1, 1)
      
      betaj <- update.beta(
        M.tilde.part[,j],
        chol.Sigma.jj,
        chol.V.tilde
      )
      
      beta[,j] <- betaj
      
      ## W update
      Wj <- update.W(Wj, betaj,
                     Sigma_jj, B, r, NNarray.obs,
                     Yj, X.obs.ord, N.obs,
                     log.like, family[j],
                     family.par=NULL, link=NULL)
      
      W.obs.ord[,j] <- Wj
    }
    
    est_time_total <- est_time_total +
      as.numeric(difftime(Sys.time(), est_start, units = "mins"))
    
    ## =========================
    ## Save only AFTER burn-in
    ## =========================
    
    if (iter > burnin) {
      
      save_idx <- iter - burnin
      
      W.obs.ord.chain[[save_idx]]  <- W.obs.ord
      beta.chain[[save_idx]]       <- beta
      Sigma.chain[[save_idx]]      <- Sigma
      chol.Sigma.chain[[save_idx]] <- chol(Sigma)
      phi.chain[save_idx]          <- phi
      
      ## =========================
      ## Prediction (post-burnin only)
      ## =========================
      
      pred_start <- Sys.time()
      
      W.pred.ord.temp <- vecchia_pred_draw(
        W.obs.ord = W.obs.ord,
        beta = beta,
        Sigma = Sigma,
        X.obs.ord,
        X.pred.ord,
        NNarray.all,
        coeff.all,
        obs.id,
        pred.id
      )
      
      log.like.details <- log.likelihood(
        W.pred.ord.temp,
        Y.pred.ord.true,
        family,
        family.par=NULL,
        link=NULL
      )
      
      log.pred.comp[save_idx,] <- log.like.details$like.comp
      log.pred.sep[save_idx]   <- log.like.details$like.total
      
      Y.hat.pred.ord[[save_idx]] <- Y.pred.ord(
        W.pred.ord.temp, N.pred, q, family
      )
      
      pred_time_total <- pred_time_total +
        as.numeric(difftime(Sys.time(), pred_start, units = "mins"))
    }
    
    ## =========================
    ## Progress
    ## =========================
    
    if (iter %% ceiling(total_iters / 10) == 0) {
      cat(sprintf("Progress: %3.0f%%\n", 100 * iter / total_iters))
    }
  }
  
  ## =========================
  ## Timing summary
  ## =========================
  
  total_time <- as.numeric(difftime(Sys.time(), start_total_time, units = "mins"))
  
  cat(sprintf("Estimation time: %.3f mins\n", est_time_total))
  cat(sprintf("Prediction time: %.3f mins\n", pred_time_total))
  cat(sprintf("Total time: %.3f mins\n", total_time))
  cat(sprintf("Acceptance.phi = %.3f\n", acc.phi / total_iters))
  
  ## =========================
  ## ELPD
  ## =========================
  
  elpd.sep <- log_mean_exp(log.pred.sep)
  
  elpd.comp <- numeric(q)
  for (j in seq_len(q)) {
    elpd.comp[j] <- log_mean_exp(log.pred.comp[, j])
  }
  
  q.subsets <- 2^q
  cond.elpd <- numeric(q.subsets - 2)
  
  counter <- 1
  for (i in 1:(q.subsets - 2)) {
    ind <- as.logical(intToBits(i)[1:q])
    cond.elpd[counter] <- elpd.sep - sum(elpd.comp[ind])
    counter <- counter + 1
  }
  
  ## =========================
  ## Prediction summaries
  ## =========================
  
  pred.summary  <- summary_stats(Y.hat.pred.ord, Y.pred.ord.true)
  pred.coverage <- pred_coverage(Y.pred.ord.true, Y.hat.pred.ord)
  
  beta.stats  <- summary_stats(beta.chain, data$true.beta)
  Sigma.stats <- summary_stats(Sigma.chain, data$true.Sigma)
  phi.stats   <- summary_stats(phi.chain, data$true.phi)
  
  output <- list(
    beta.stats = beta.stats,
    Sigma.stats = Sigma.stats,
    phi.stats = phi.stats,
    elpd = elpd.sep,
    cond.elpd = cond.elpd,
    pred.summary = pred.summary,
    pred.coverage = pred.coverage,
    acc.phi = acc.phi / total_iters,
    est.time = est_time_total,
    pred.time = pred_time_total,
    total.time = total_time
  )
  
  return(output)
}