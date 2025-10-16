rm(list = ls())

libraries <- c("spam", "this.path", "fields", "fBasics", "MCMCpack", "truncnorm", 
               "rlist", "foreach","doParallel", "FNN", "mcmcse")

# Load libraries
invisible(lapply(libraries, library, character.only = TRUE))

mydir <- this.path::here()
setwd(mydir)


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


load(file.path("data_full.RData"), eva.data <- new.env())

full.data.DF <- eva.data$data_DF
full.data <- as.list(eva.data)

N.t <- 161
N.months <- 7
N.years <- 23
N.s <- 3503

t.index <- 140

spatial.subset <- seq(N.s * (t.index - 1) + 1, N.s * t.index)
data <- full.data.DF[spatial.subset, ]

lat <- data$lat
lon <- data$lon

locations <- cbind(lon, lat)

# Adding a mean term

X <- cbind(1, locations)

#Number of features in the spatial model
p <- ncol(X)

# Response 

Y_BA <- log(1+data$BA)
Y_CNT <- data$CNT
Y <- cbind(Y_BA, Y_CNT)

N <- nrow(Y)
q <- ncol(Y)

locations <- cbind(lon, lat)

# Adding a mean term

X <- cbind(1, locations)


obs.locs <- locations

X.obs <- X
Y.obs <- Y

m <- 20
obs.ord <- max_min(obs.locs)
obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]
distobs.ord <- rdist(obs.locs.ord)
diameter <- max(distobs.ord)
b_phi <- diameter / log(10)
NNarray.obs <- neighbor_matrix(obs.locs.ord, m)
distobs.nn <- dist.nn(obs.locs.ord, neighbor_matrix = NNarray.obs)

N.obs <- nrow(obs.locs)

phi <- b_phi/2
nu <- 0.3

  
  # Number of iterations
  niters <- 5e4
  N.obs <- nrow(obs.locs)
  
  
  beta.stats.list <- vector("list", q)
  Sigma.stats.list <- vector("list", q)
  phi.stats.list <- vector("list", q)
  
  
  log.like <- rep(NA, q)
  W.obs.ord.mean <- matrix(NA, nrow = N.obs, ncol = q)
  
  
  beta.ess <- matrix(NA, nrow = p, ncol = q)
  Sigma.ess <- matrix(NA, nrow = q, ncol = q)
  W.obs.ord.ess <- matrix(NA, nrow = N.obs, ncol = q)
  phi.ess <- rep(NA, q)
  
  
  # Reorder Y.obs, X.obs, and W.obs accordingly
  Y.obs.ord <- Y.obs[obs.ord, , drop = FALSE]
  X.obs.ord <- X.obs[obs.ord, , drop = FALSE]
  
  family <- c("Gaussian", "Poisson")
  
  chol.K.obs.ord.inv <- U.sgv(distobs.nn, NNarray.obs, phi, nu, m)
  
  Sigma.init <- c(1, 1/3)
  
  W.obs.ord.sep <- replicate(niters, matrix(NA, N.obs, q), simplify = F)
  
  
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
    beta <- M.prior
    Sigma <- as.matrix(Sigma.init[j])
    W.obs.ord <- X.obs.ord %*% beta
 
    tuning.phi <- b_phi/3
    
    
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
    
    rm(fit.main)
    
    W.obs.ord.mean[,j] <- Reduce("+", W.ord.post) / niters
    
    #### Summary of estimation ###
    
    beta.stats.list[[j]] <- summary_stats(beta.post, NULL)
    Sigma.stats.list[[j]] <- summary_stats(Sigma.post, NULL)
    phi.stats.list[[j]] <- summary_stats(phi.post, NULL)
    
    
    #### ESS calculation for component chains ####
    
    beta.ess[,j] <- compute_ess(beta.post)
    Sigma.ess[j,j] <- compute_ess(Sigma.post)
    phi.ess[j] <- compute_ess(phi.post)
    W.obs.ord.ess[,j] <- compute_ess(W.ord.post)
   
    
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
  
  
  mlpd.obs <- mlpd(W.obs.ord.sep, Y.obs.ord, family)
  waic.obs <- waic(W.obs.ord.sep, Y.obs.ord, mlpd.obs, family)
  
 
  out <- list("beta.stats" = beta.stats,
              "Sigma.stats" = Sigma.stats, 
              "phi.stats" = phi.stats,
              "beta.ess" = beta.ess,
              "Sigma.ess" = Sigma.ess,
              "phi.ess" = phi.ess,
              "W.obs.ord.ess" = W.obs.ord.ess,
              "waic" = waic.obs,
              "acc.phi" = acc.phi,
              "total_time" = total_time)

# Save results
save(results, file = "results_separate_est.RData")
