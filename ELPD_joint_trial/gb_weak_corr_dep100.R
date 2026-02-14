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
  niters <- 5e4
  
  
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
  
  
  
  #### ELPD Computations ###
  
  pred.iters <- niters
  
  ##### Sample size for estimating inner expectation in ELPD
  R <- 1e2
  
  ## For joint model we are saving into S * R array
  
  elpd_array <- array(0, dim = c(pred.iters, R))
  
  se.elpd <- numeric(pred.iters)
  
  
  start_time <- Sys.time()  # Start time using Sys.time()
  
  for(i in 1:pred.iters){
    
    
    U.joint <- U.sgv(distlocs.nn, NNarray.all, phi.post[[i]], 
                     nu, m)
    
    
    U.obs.col <- U.joint[ , 1:N.obs]
    U.pred.col <- U.joint[ , (N.obs+1):(N.pred+N.obs)]
    W.pred.prec <- crossprod.spam(U.pred.col)
    W.pred.var <- solve.spam(as.spam(W.pred.prec))
    chol.W.pred.var <- chol(W.pred.var)
    
    for (k in 1:R) {
      
      ### W.pred.ord.temp draws sample from latent predictive process W* | W, B, \Sigma, \phi
      ### W.pred.ord is written in pred_veccia_slice file
      
      W.pred.ord.temp <-  W.pred.ord(W.obs.ord = W.ord.post[[i]], 
                                     beta = beta.post[[i]], 
                                     chol.Sigma = chol.Sigma.post[[i]], 
                                     phi = phi.post[[i]], 
                                     nu, m,
                                     X.obs.ord, X.pred.ord, U.joint, U.pred.col, U.obs.col,
                                     chol.W.pred.var, N.obs, N.pred, q, p)
      
      ### Clog-likelihood is defined in aux_slice.R file
      
      elpd_array[i,k] <- log.likelihood(W.obs.ord =  W.pred.ord.temp, 
                                        Y.obs.ord = Y.pred.ord.true, 
                                        family = family)
      
      
    }
    
    if(i %% (pred.iters / 10) == 0) {
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      cat(sprintf("ELPD computation progress: %.0f%%, Elapsed time: %.2f minutes\n", 
                  100 * (i /pred.iters), elapsed_time))
    }
    
  }
  
  end_time <- Sys.time()
  total_elapsed_time <- difftime(end_time, start_time, units = "mins")[[1]] # in minutes
  
  total_elapsed_time
  
  ### Saving the entire array
  
  save(elpd_array, file = "elpd_array_joint.RData")
  
 
  elpd_joint <- mean(
    apply(elpd_array, 1, function(x) {
      m <- max(x)
      m + log(mean(exp(x - m)))
    })
  )
  
  elpd_joint
  
  ## This is done for numerical stabilization (Suggested by GPT)
  
  mean.inner <- numeric(nrow(elpd_array))
  se.inner   <- numeric(nrow(elpd_array))
  
  for (i in 1:nrow(elpd_array)) {
    
    log_vals <- elpd_array[i, ]
    m <- max(log_vals)
    
    ## stable likelihood values
    Z <- exp(log_vals - m)
    
    ## true mean
    mean_Z <- exp(m) * mean(Z)
    
    ## true sample sd
    sd_Z <- exp(m) * sd(Z)
    
    mean.inner[i] <- mean_Z
    se.inner[i]   <- sd_Z / sqrt(R)
  }
  
  
  rng <- range(log(se.inner), finite = TRUE)
  
  hist(log(se.inner),
       breaks = 30,
       xlim = rng,
       main = "Monte Carlo log(SE) of Inner Expectation",
       xlab = "log(SE)",
       col = "lightblue",
       border = "white",
       probability = TRUE)
  
  log(mean(se.inner))
  
  summary(log(se.inner))