### Gaussian Poisson phi 0.1 dependent

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
source("update_parameters_slice.R")
source("mcmc_main_slice.R")
source("pred_vecchia_slice.R")
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
  true.phi <- 0.1
  true.nu <- 0.5
  pred.prop <- 0.2
  
  family <- c("Gaussian", "Poisson")
  
  data <- sim.data(q = 2, N = 2.5e3, 
                   family = family,
                   true.beta = true.beta,
                   true.Sigma = true.Sigma, 
                   true.phi = true.phi,
                   true.nu = true.nu,
                   pred.prop = 0.2)
  
  
  
  #### Prior specifications ####
  
  M.prior <- matrix(0, p, q)
  V.prior <-  1e2*diag(p) 
  S.prior <-  diag(q)
  df.prior <- q + 1
  
  
  ### Vecchia specifications and pre-computations ###
  
  m <- 20
  obs.locs <- data$obs.locs
  N.obs <- data$N.obs
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
  
  
  #### Prediction set up ####
  
  pred.locs <- data$pred.locs
  N.pred <- data$N.pred
  pred.ord <- max_min(pred.locs)
  ord <- c(obs.ord, pred.ord + nrow(obs.locs))
  locs.ord <- rbind(obs.locs, pred.locs)[ord, , drop = FALSE]
  distlocs.ord <- rdist(locs.ord)
  NNarray.all <- neighbor_matrix(locs.ord, m)
  distlocs.nn <- dist.nn(locs.ord, NNarray.all)
  NNarray.pred=NNarray.all[N.obs + (1:N.pred),]
  
  
  #neighbor_matrix(pred.locs, m)
  
  X.pred <- data$X.pred
  Y.pred.true <- data$Y.pred.true
  
  X.pred.ord <- X.pred[pred.ord, , drop = FALSE]
  Y.pred.ord.true <- Y.pred.true[pred.ord, , drop = FALSE]
  
  phi <- data$true.phi
  nu <- data$true.nu
  
  U.joint <- U.sgv(distlocs.nn, NNarray.all, phi, 
                   nu, m)
  
  
  U.obs.col <- U.joint[ , 1:N.obs]
  U.pred.col <- U.joint[ , (N.obs+1):(N.pred+N.obs)]
  W.pred.prec <- crossprod.spam(U.pred.col)
  W.pred.var <- solve.spam(as.spam(W.pred.prec))
  chol.W.pred.var <- chol(W.pred.var)
  
  # Reorder Y.obs, X.obs, and W.obs accordingly
  Y.obs.ord <- Y.obs[obs.ord, , drop = FALSE]
  X.obs.ord <- X.obs[obs.ord, , drop = FALSE]
  
  # Sample initial parameters from MLE
 
  beta <- data$true.beta
  Sigma <- data$true.Sigma
  W.obs.ord <- data$true.W.obs[obs.ord, , drop = FALSE]
  
 
  W.obs.ord.mean <- matrix(NA, nrow = N.obs, ncol = q)
 
 
  chol.K.obs.ord.inv <- U.sgv(distobs.nn, NNarray.obs, 
                              phi, nu, m)
 
 
  tuning.phi <- 1.8e-2
  
  # Number of iterations
  niters <- 5e4
  
  
  fit.main <- mcmc.main(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                        distobs.nn, chol.K.obs.ord.inv, NNarray.obs, 
                        family, nu, W.obs.ord, beta, Sigma, phi,
                        M.prior, V.prior, S.prior, df.prior, b_phi,
                        niters, tuning.phi)
  
  pred.iters <- niters
  
  W.ord.post <- fit.main$W.ord.samples
  beta.post <- fit.main$beta.samples
  Sigma.post <- fit.main$Sigma.samples
  chol.Sigma.post <- fit.main$chol.Sigma.samples
  phi.post <- fit.main$phi.samples
  acc.phi <- fit.main$acceptance.phi
  total_time <- fit.main$total_time
  log.like <- fit.main$log.like.mean
  
  rm(fit.main)
  
  #### Summary of estimation ###
  
  W.obs.ord.mean <- Reduce("+", W.ord.post) / niters
  
  log_like_mean <- log.like
  
  log_like_post <- log.likelihood(W.obs.ord.mean, Y.obs.ord, family) 
  

  # summary_stats(W.ord.post, data$true.W.obs.ord)
  beta.stats <- summary_stats(beta.post, data$true.beta)
  Sigma.stats <- summary_stats(Sigma.post, data$true.Sigma)
  phi.stats <- summary_stats(phi.post, data$true.phi)
  
  
  
  #### ESS calculation for component chains ####
  
  beta.ess <- compute_ess(beta.post)
  Sigma.ess <- compute_ess(Sigma.post)
  phi.ess <- compute_ess(phi.post)
  W.obs.ord.ess <- compute_ess(W.ord.post)
  
  Y.pred.joint <- predictive.samples(W.ord.post, beta.post, chol.Sigma.post, 
                                     phi.post, nu, m,
                                     X.obs.ord, X.pred.ord, U.joint, U.pred.col, U.obs.col, 
                                     chol.W.pred.var, N.obs, N.pred, pred.iters, q, p)
  
  pred.summary <- summary_stats(Y.pred.joint, Y.pred.ord.true)
  mlpd.pred <- mlpd(W.ord.post, Y.pred.ord.true, family)
 
  mlpd.obs <- mlpd(W.ord.post, Y.obs.ord, family)
  waic.obs <- waic(W.ord.post, Y.obs.ord, mlpd.obs, family)
  
  #### Prediction quality assesment ####
  
  logs <- compute_logs(Y_true = Y.pred.ord.true, Y_pred_samples = Y.pred.joint, family = family)
  dss <- compute_dss(Y_true = Y.pred.ord.true, Y_pred_samples = Y.pred.joint, family = family)
  crps <- compute_crps(Y_true = Y.pred.ord.true, Y_pred_samples = Y.pred.joint, family = family)
  es <- energy_score(Y_true = Y.pred.ord.true, Y_pred_samples = Y.pred.joint)
  rmspe <- RMSPE(Y.pred.samples = Y.pred.joint, Y.true = Y.pred.ord.true, pred.iters)
  pred.coverage <- pred_coverage(Y.pred.ord.true, Y.pred.joint)
 
  
 
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
list.save(results, file = "gp_weak_corr_ind2500.RData")
