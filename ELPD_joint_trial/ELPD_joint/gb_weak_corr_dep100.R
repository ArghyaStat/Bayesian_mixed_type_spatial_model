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

  
  ### data generation
                     
  set.seed(3092)

  p <- 3
  q <- 2
  
  true.beta <- matrix(c(1.0, -0.5,  3,  1.5, -1.2,  0.0), nrow = p, ncol = q, byrow = TRUE)
  true.Sigma <- matrix(c(2, 1, 1, 1), nrow = q, ncol = q, byrow = TRUE)
  true.phi <- 0.1
  true.nu <- 0.5
  pred.prop <- 0.2
  
  family <- c("Gaussian", "Binomial")
  
  data <- sim.data(q = 2, N = 1e2, 
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
  
  N.obs <- data$N.obs
  obs.locs <- data$obs.locs
  Y.obs <- data$Y.obs
  X.obs <- data$X.obs
  m <- N.obs - 1
  
  obs.ord <- max_min(obs.locs)
  # Assume max_min returns max-min order indices
  obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]  # Reorder locations based on max-min ordering
  distobs.ord <- rdist(obs.locs.ord)
  diameter <-  max(distobs.ord)
  b_phi <- diameter/3
  
  
  
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
 
 
  tuning.phi <- 5e-2
  
  # Number of iterations
  niters <- 1e4
  
  
  fit.main <- mcmc.main(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                        distobs.nn, chol.K.obs.ord.inv, NNarray.obs, 
                        family, nu, W.obs.ord, beta, Sigma, phi,
                        M.prior, V.prior, S.prior, df.prior, b_phi,
                        niters, tuning.phi)
  
  W.ord.post <- fit.main$W.ord.samples
  beta.post <- fit.main$beta.samples
  Sigma.post <- fit.main$Sigma.samples
  chol.Sigma.post <- fit.main$chol.Sigma.samples
  phi.post <- fit.main$phi.samples
  acc.phi <- fit.main$acceptance.phi
  total_time <- fit.main$total_time
  log.like <- fit.main$log.like.mean
  
  rm(fit.main)
  

  #  ELPD COMPUTATION — JOINT   #
 
  pred.iters <- niters   # S posterior draws
  
  ## This computes:
  ## sum_{i=1}^u sum_{j=1}^q log pi_j( y_j(s_i*) | w_j*(s_i*) )
  
  joint_log_pred_density <- function(W.pred.draw, Y.pred.true, family) {
    
    N.pred <- nrow(Y.pred.true)
    q      <- ncol(Y.pred.true)
    
    loglik <- 0
    
    for (j in 1:q) {
      
      if (family[j] == "Gaussian") {
        
        ## Gaussian with unit variance (adjust if variance differs)
        loglik <- loglik +
          sum(dnorm(Y.pred.true[, j],
                    mean = W.pred.draw[, j],
                    sd   = 1,
                    log  = TRUE))
        
      } else if (family[j] == "Binomial") {
        
        ## logistic link assumed
        prob <- 1 / (1 + exp(-W.pred.draw[, j]))
        
        loglik <- loglik +
          sum(dbinom(Y.pred.true[, j],
                     size = 1,
                     prob = prob,
                     log  = TRUE))
      }
    }
    
    return(loglik)
  }
  
  
  log_pred_joint <- numeric(pred.iters)
  
  start_time <- Sys.time()
  
  for (s in 1:pred.iters) {
    
    ## Recompute Vecchia factor for current phi
    U.joint <- U.sgv(distlocs.nn, NNarray.all,
                     phi.post[[s]],
                     nu, m)
    
    U.obs.col  <- U.joint[, 1:N.obs]
    U.pred.col <- U.joint[, (N.obs+1):(N.obs+N.pred)]
    
    ## Predictive covariance of W*
    W.pred.prec <- crossprod.spam(U.pred.col)
    W.pred.var  <- solve.spam(as.spam(W.pred.prec))
    chol.W.pred.var <- chol(W.pred.var)
    
    
  
    ## Draw latent predictive W* | W, theta
   
    
    W.pred.draw <- W.pred.ord(
      W.obs.ord   = W.ord.post[[s]],
      beta        = beta.post[[s]],
      chol.Sigma  = chol.Sigma.post[[s]],
      phi         = phi.post[[s]],
      nu, m,
      X.obs.ord, X.pred.ord,
      U.joint, U.pred.col, U.obs.col,
      chol.W.pred.var,
      N.obs, N.pred, q, p
    )
    
    
  
    log_pred_joint[s] <-
      joint_log_pred_density(W.pred.draw,
                             Y.pred.ord.true,
                             family)
    
    
    if (s %% (pred.iters / 10) == 0) {
      elapsed_time <- as.numeric(
        difftime(Sys.time(), start_time, units = "mins")
      )
      cat(sprintf("ELPD progress: %.0f%% | %.2f minutes\n",
                  100 * s / pred.iters,
                  elapsed_time))
    }
  }
  
  end_time <- Sys.time()
  
  cat("Total ELPD computation time (minutes):",
      as.numeric(difftime(end_time, start_time, units = "mins")),
      "\n")
  
  
 
  m <- max(log_pred_joint)
  
  elpd_joint <-
    m + log(mean(exp(log_pred_joint - m)))
  
  cat("Joint ELPD =", elpd_joint, "\n")
  
  elpd_joint
  