mixed.workflow.sep <- function(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                               nu, W.obs.ord, beta, Sigma, phi,
                               M.prior, V.prior, a.prior, b.prior, b_phi,
                               niters, tuning.phi,
                               Y.pred.ord.true, X.pred.ord, N.pred,
                               distall.nn, coeff.all, NNarray.all,
                               family, family.par = NULL, link = NULL){
  
  ## Storage
  W.obs.ord.chain  <- replicate(niters, matrix(NA, N.obs, q), simplify = FALSE)
  beta.chain       <- replicate(niters, matrix(NA, p, q), simplify = FALSE)
  Sigma.chain      <- replicate(niters, matrix(NA, q, q), simplify = FALSE)
  chol.Sigma.chain <- replicate(niters, matrix(NA, q, q), simplify = FALSE)
  phi.chain        <- rep(NA, niters)
  
  
  ## Storage for prediction
  
  Y.hat.pred.ord <- replicate(niters, matrix(NA, N.pred, q), simplify = F)
  
  log.pred.sep <- numeric(niters)
  log.pred.comp  <- matrix(NA, niters, q)
  
  
  ## Initial Sigma quantities
  chol.Sigma <- chol(Sigma)
  
  ## Prior quantities
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  chol.prec.V.prior <- chol(prec.V.prior)
  
  NNarray.obs <- NNarray.all[1:N.obs,]
  distobs.nn <- distall.nn[1:N.obs]
  coeff.obs <- vecchia_coeff(distobs.nn, NNarray.obs, phi, nu)
  
  obs.id  <- which(ord <= N.obs)
  pred.id <- which(ord > N.obs)
  
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
    
    ## Beta posterior quantities
    
    vecchia.prod.out <- vecchia_products(W.obs.ord, X.obs.ord, B, r,  NNarray.obs)
    
    UX <- vecchia.prod.out$UX
    UW <- vecchia.prod.out$UW
    
    prec.V.tilde <- crossprod(UX) + prec.V.prior
    
    chol.prec.V.tilde <- chol(prec.V.tilde)
    V.tilde <- chol2inv(chol(prec.V.tilde))
    
    chol.V.tilde <- chol(V.tilde)
    
    M.tilde.part <-
      crossprod(UX, UW) +
      prec.V.prior %*% M.prior
    
    #### Loop over components ####
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
        fixed_variances[j])
      
      Sigma[j,j] <- Sigma.details$Sigma
      
      Sigma_jj <- as.matrix(Sigma[j,j])
      
      chol.Sigma.jj <- matrix(sqrt(Sigma[j,j]), 1, 1)
      
      betaj <- update.beta(
        M.tilde.part[,j],
        chol.Sigma.jj,
        chol.V.tilde
      )
      
      beta[,j] <- betaj
    
      #### W update ####
      Wj <- update.W(Wj, betaj,
                     Sigma_jj, B, r, NNarray.obs,
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
    
    W.pred.ord.temp <-  vecchia_pred_draw(W.obs.ord = W.obs.ord.chain[[iter]],
                                          beta =  beta.chain[[iter]],
                                          Sigma = Sigma.chain[[iter]],
                                          X.obs.ord,
                                          X.pred.ord,
                                          NNarray.all,
                                          coeff.all,
                                          obs.id,
                                          pred.id)
    
    
    log.like.details <- log.likelihood(
      W.pred.ord.temp,
      Y.pred.ord.true,
      family,
      family.par=NULL,
      link=NULL)
    
    log.pred.comp[iter,] <- log.like.details$like.comp
    log.pred.sep[iter] <- log.like.details$like.total
    
    elpd.sep <- log_mean_exp(log.pred.sep)
    
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
  
  # Total possible subsets
  q.subsets <- 2^q
  
  # We exclude empty (0) and full (2^q - 1)
  n.valid <- q.subsets - 2
  
  cond.elpd <- numeric(n.valid)
  
  counter <- 1
  
  for (i in 1:(q.subsets - 2)) {
    
    ind <- as.logical(intToBits(i)[1:q])
    
    elpd.B <- sum(elpd.comp[ind])
    
    cond.elpd[counter] <- elpd.sep - elpd.B
    
    counter <- counter + 1
  }
  
  
  pred.summary <- summary_stats(Y.hat.pred.ord,Y.pred.ord.true)
  pred.coverage <- pred_coverage(Y.pred.ord.true,Y.hat.pred.ord)
  
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
                 "elpd" = elpd.sep,
                 "elpd.comp" = elpd.comp,
                 "cond.elpd" = cond.elpd,
                 "pred.summary" = pred.summary,
                 "pred.coverage" = pred.coverage,
                 "acc.phi" = acc.phi/niters,
                 "est.time" = est_time_total,
                 "pred.time" = pred_time_total,
                 "total.time" = total_time)
  
  return(output)
}