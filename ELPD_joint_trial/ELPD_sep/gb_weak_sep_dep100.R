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
  
  
  # Number of iterations
  niters <- 1e4
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
 
  
  pred.iters <- niters   # T posterior draws
  elpd_sep   <- 0        # accumulate over components
  
  
  for (j in 1:q) {
    
    
    sub_family <- family[j]
    q_sub <- 1
 
    
    M.prior  <- matrix(0, p, q_sub)
    V.prior  <- 1e2 * diag(p)
    S.prior  <- diag(q_sub)
    df.prior <- q_sub + 1
    
    
    Y.obs.ord <- as.matrix(Y.obs[obs.ord, j])
    W.obs.init <- as.matrix(data$true.W.obs[obs.ord, j])
    beta.init  <- data$true.beta[, j]
    Sigma.init <- as.matrix(data$true.Sigma[j, j])
    
    tuning.phi <- 8e-1
    
    
    fit.main <- mcmc.main(
      Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p,
      q = q_sub,
      distobs.nn, chol.K.obs.ord.inv, NNarray.obs,
      family = sub_family, nu,
      W.obs.init, beta.init, Sigma.init, phi,
      M.prior, V.prior, S.prior, df.prior, b_phi,
      niters, tuning.phi
    )
    
    W.post        <- fit.main$W.ord.samples
    beta.post     <- fit.main$beta.samples
    phi.post      <- fit.main$phi.samples
    Sigma.post    <- fit.main$Sigma.samples   # IMPORTANT
    
    
    log_pred_j <- numeric(pred.iters)
    
 
    for (t in 1:pred.iters) {
      
      ## Build Vecchia precision for phi_j
      U.joint <- U.sgv(distlocs.nn, NNarray.all,
                       phi.post[[t]], nu, m)
      
      U.obs.col  <- U.joint[, 1:N.obs]
      U.pred.col <- U.joint[, (N.obs+1):(N.obs+N.pred)]
      
     
      W.pred.prec <- crossprod.spam(U.pred.col)
      W.pred.var  <- solve.spam(as.spam(W.pred.prec))
      chol.W.pred.var <- chol(W.pred.var)
      
      
      
      W.centered <- W.post[[t]] - X.obs.ord %*% beta.post[[t]]
      
      mean.part <- U.obs.col %*% W.centered
      mean.tilde <- crossprod.spam(U.pred.col, mean.part)
      
      W.pred.mean <-
        X.pred.ord %*% beta.post[[t]] -
        crossprod(chol.W.pred.var, mean.tilde)
  
   
      Sigma_jj <- as.numeric(Sigma.post[[t]])
      
   
      z <- matrix(rnorm(N.pred), nrow = N.pred, ncol = 1)
      
     
      W.pred.draw <-
        as.vector(
          W.pred.mean +
            sqrt(Sigma_jj) *
            (chol.W.pred.var %*% z)
        )
      
      
      if (sub_family == "Gaussian") {
        
        log_pred_j[t] <-
          sum(dnorm(Y.pred.ord.true[, j],
                    mean = W.pred.draw,
                    sd = 1,
                    log = TRUE))
        
      } else if (sub_family == "Binomial") {
        
        prob <- plogis(W.pred.draw)
        
        log_pred_j[t] <-
          sum(dbinom(Y.pred.ord.true[, j],
                     size = 1,
                     prob = prob,
                     log = TRUE))
        
      } else {
        stop("Unsupported family.")
      }
    }
    
    
    m <- max(log_pred_j)
    
    elpd_j <-
      m + log(mean(exp(log_pred_j - m)))
    
    cat("ELPD component", j, "=", elpd_j, "\n")
    
    elpd_sep <- elpd_sep + elpd_j
  }
  
elpd_sep
  
