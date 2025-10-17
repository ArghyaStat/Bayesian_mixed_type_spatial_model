## GP weak dependent 100 ##

rm(list = ls())

libraries <- c("spam", "this.path", "fields", "fBasics", "MCMCpack", "truncnorm", 
               "rlist", "foreach","doParallel", "FNN", "mcmcse")


# Load libraries
invisible(lapply(libraries, library, character.only = TRUE))

mydir <- this.path::here()
setwd(mydir)


source("data_simulation.R")
source("vecchia_slice.R")
source("aux_functions_slice.R")
source("update_parameters_slice-sep.R")
source("mcmc_main_slice.R")
source("pred_vecchia_slice-sep.R")
source("predictive_scores_slice.R")

functions_list <- c("dmatnorm.sgv", "log.likelihood", "update.W", 
                    "update.beta", "update.Sigma", "update.phi", "max_min",
                    "neighbor_matrix", "comb", "dist.nn", "U.sgv", 
                    "tuning.update","compute_ess", "summary_stats", "score_function",
                    "RMSPE", "pred_coverage", "energy_score", "compute_logs", 
                    "compute_dss", "compute_crps")
functions_pred <- c('Y.pred.ord')

reps <- 50
n.cores <- 50
cl <- makeCluster(n.cores)
registerDoParallel(cl)

results <- foreach(r = 1:reps, .packages = libraries, .export = functions_pred) %dopar% {
  
  ### data generation
                     
  set.seed(r + 5)

  p <- 3
  q <- 2
  
  true.beta <- matrix(c(1.0, -0.5,  3,  1.5, -1.2,  0.0), nrow = p, ncol = q, byrow = TRUE)
  true.Sigma <- matrix(c(2, 0, 0, 1), nrow = q, ncol = q, byrow = TRUE)
  true.phi <- 0.5
  true.nu <- 0.5
  pred.prop <- 0.2
  
  family <- c("Gaussian", "Poisson")
  
  data <- sim.data(q = 2, N = 1e2, 
                   family = family,
                   true.beta = true.beta,
                   true.Sigma = true.Sigma, 
                   true.phi = true.phi,
                   true.nu = true.nu,
                   pred.prop = 0.2)
  
  
  # Number of iterations
  niters <- 5e4
  pred.iters <- niters
  N.obs <- data$N.obs
  N.pred <- data$N.pred
  
  Y.pred.sep <- replicate(pred.iters, matrix(NA, N.pred, q), simplify = F)
  W.obs.ord.sep <- replicate(niters, matrix(NA, N.obs, q), simplify = F)
  
  beta.stats.list <- vector("list", q)
  Sigma.stats.list <- vector("list", q)
  phi.stats.list <- vector("list", q)
  
  
  log.like <- rep(NA, q)
  W.obs.ord.mean <- matrix(NA, nrow = N.obs, ncol = q)
  
  
  beta.ess <- matrix(NA, nrow = p, ncol = q)
  Sigma.ess <- matrix(NA, nrow = q, ncol = q)
  W.obs.ord.ess <- matrix(NA, nrow = N.obs, ncol = q)
  phi.ess <- rep(NA, q)
  
  m <- 20
  obs.locs <- data$obs.locs
  Y.obs <- data$Y.obs
  X.obs <- data$X.obs
  
  
  obs.ord <- max_min(obs.locs) 
  # Assume max_min returns max-min order indices
  obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]  # Reorder locations based on max-min ordering
  distobs.ord <- rdist(obs.locs.ord)
  diameter <-  max(distobs.ord)
  b_phi <- diameter/log(10)
  
  
  NNarray.obs <- neighbor_matrix(obs.locs.ord, m)
  distobs.nn <- dist.nn(obs.locs.ord, neighbor_matrix = NNarray.obs)
  
  X.obs.ord <- X.obs[obs.ord, , drop = FALSE]
  
  #### Prediction set up ####
  
  pred.locs <- data$pred.locs
  pred.ord <- max_min(pred.locs)
  ord <- c(obs.ord, pred.ord + nrow(obs.locs))
  locs.ord <- rbind(obs.locs, pred.locs)[ord, , drop = FALSE]
  distlocs.ord <- rdist(locs.ord)
  NNarray.all <- neighbor_matrix(locs.ord, m)
  distlocs.nn <- dist.nn(locs.ord, NNarray.all)
  NNarray.pred=NNarray.all[N.obs + (1:N.pred),]
  
  
  
  X.pred <- data$X.pred
  X.pred.ord <- X.pred[pred.ord, , drop = FALSE]
  
  phi <- data$true.phi
  nu <- data$true.nu
  
  U.joint <- U.sgv(distlocs.nn, NNarray.all, phi, 
                   nu, m)
  
  
  U.obs.col <- U.joint[ , 1:N.obs]
  U.pred.col <- U.joint[ , (N.obs+1):(N.pred+N.obs)]
  W.pred.prec <- crossprod.spam(U.pred.col)
  W.pred.var <- solve.spam(as.spam(W.pred.prec))
  chol.W.pred.var <- chol(W.pred.var)
  
  chol.K.obs.ord.inv <- U.sgv(distobs.nn, NNarray.obs, 
                              phi, nu, m)
  
  Y.pred.true <- data$Y.pred.true
  Y.pred.ord.true <- Y.pred.true[pred.ord, , drop = FALSE]
  
  
  for(j in 1:q){
    
    sub_family <- family[j]
    q_sub <- length(sub_family)
  
  
  #### Prior specifications ####
  
    M.prior <- matrix(0, p, q_sub)
    V.prior <-  1e2*diag(p) 
    S.prior <-  diag(q_sub)
    df.prior <- q_sub + 1

  
    # Reorder Y.obs, X.obs, and W.obs accordingly
    Y.obs.ord <- as.matrix(Y.obs[obs.ord, , drop = FALSE][,j])
    W.obs.ord <- as.matrix(data$true.W.obs[obs.ord, , drop = FALSE][,j])
    beta <- data$true.beta[,j]
    Sigma <- as.matrix(data$true.Sigma[j,j])
    
   
    tuning.phi <- 3e-2
  
    
    fit.main <- mcmc.main(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, 
                          q = q_sub, distobs.nn, chol.K.obs.ord.inv, NNarray.obs, 
                          family = sub_family, nu, W.obs.ord, beta, Sigma, phi,
                          M.prior, V.prior, S.prior, df.prior, b_phi,
                          niters, tuning.phi)
    
   
    
    W.ord.post <- fit.main$W.ord.samples
    beta.post <- fit.main$beta.samples
    Sigma.post <- fit.main$Sigma.samples
    chol.Sigma.post <- fit.main$chol.Sigma.samples
    phi.post <- fit.main$phi.samples
    acc.phi <- fit.main$acceptance.phi
    total_time <- fit.main$total_time
    
    
    log.like[j] <- fit.main$log.like.mean
    W.obs.ord.mean[,j] <- Reduce("+", W.ord.post) / niters
    
    rm(fit.main)
    
    #### Summary of estimation ###
    
    beta.stats.list[[j]] <- summary_stats(beta.post, as.matrix(data$true.beta[,j]))
    Sigma.stats.list[[j]] <- summary_stats(Sigma.post, as.matrix(data$true.Sigma[j,j]))
    phi.stats.list[[j]] <- summary_stats(phi.post, data$true.phi)
    
    
    #### ESS calculation for component chains ####
    
    beta.ess[,j] <- compute_ess(beta.post)
    Sigma.ess[j,j] <- compute_ess(Sigma.post)
    phi.ess[j] <- compute_ess(phi.post)
    W.obs.ord.ess[,j] <- compute_ess(W.ord.post)
  
    
    pred.samples <- predictive.samples(W.ord.post, beta.post, chol.Sigma.post, phi.post, nu, m, 
                                         X.obs.ord, X.pred.ord, U.joint, U.pred.col, U.obs.col, 
                                         chol.W.pred.var, N.obs, N.pred, pred.iters, p, family = sub_family)
    
    Y.pred.sep <- Map(function(mat, col) {
      mat[, j] <- col
      mat
    }, Y.pred.sep, pred.samples)
    
   W.obs.ord.sep <- Map(function(mat, col) {
      mat[, j] <- col
      mat
    }, W.obs.ord.sep, W.ord.post)
    
  
  }
  
 
  log_like_mean <- mean(log.like)

  Y.obs.ord <- as.matrix(Y.obs[obs.ord, , drop = FALSE])
  log_like_post <- log.likelihood(W.obs.ord.mean, Y.obs.ord, family) 
    
  beta.stats <- combine_stats(beta.stats.list)
  Sigma.stats <- combine_stats(Sigma.stats.list)
  phi.stats <- combine_stats(phi.stats.list)
  
  #### Prediction quality assesment ####
  
  logs <- compute_logs(Y_true = Y.pred.ord.true, Y_pred_samples = Y.pred.sep, family = family)
  dss <- compute_dss(Y_true = Y.pred.ord.true, Y_pred_samples = Y.pred.sep, family = family)
  crps <- compute_crps(Y_true = Y.pred.ord.true, Y_pred_samples = Y.pred.sep, family = family)
  es <- energy_score(Y_true = Y.pred.ord.true, Y_pred_samples = Y.pred.sep)
  rmspe <- RMSPE(Y.pred.samples = Y.pred.sep, Y.true = Y.pred.ord.true, pred.iters)
  pred.coverage <- pred_coverage(Y.pred.ord.true, Y.pred.sep)
  pred.summary <- summary_stats(Y.pred.sep, Y.pred.ord.true)
  
  mlpd.pred <- mlpd(W.obs.ord.sep, Y.pred.ord.true, family)
  
  mlpd.obs <- mlpd(W.obs.ord.sep, Y.obs.ord, family)
  waic.obs <- waic(W.obs.ord.sep, Y.obs.ord, mlpd.obs, family)
  
  
  out <- list("beta.stats" = beta.stats,
              "Sigma.stats" = Sigma.stats, 
              "phi.stats" = phi.stats,
              "beta.ess" = beta.ess,
              "Sigma.ess" = Sigma.ess,
              "phi.ess" = phi.ess,
              "W.obs.ord.ess" = W.obs.ord.ess,
              "logs" = logs,
              "dss" = dss,
              "crps" = crps,
              "es" = es,
              "rmspe" = rmspe,
              "mlpd" = mlpd.pred,
              "waic.obs" = waic.obs,
              "pred.summary" = pred.summary,
              "pred.coverage" = pred.coverage,
              "acc.phi" = acc.phi,
              "total_time" = total_time)
  
  out
  
}

stopCluster(cl)

# Save results to an .RData file
list.save(results, file = "gp_strong_sep_ind100.Rdata")
