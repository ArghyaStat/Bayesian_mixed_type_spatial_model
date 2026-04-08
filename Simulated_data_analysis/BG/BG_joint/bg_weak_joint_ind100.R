### Binomial Gaussian Poisson phi 0.1 joint ind

rm(list = ls())

libraries <- c("spam", "this.path", "fields", "fBasics", "MCMCpack", "truncnorm", 
               "rlist", "foreach","doParallel", "FNN", "mcmcse")


# Load libraries
invisible(lapply(libraries, library, character.only = TRUE))

mydir <- this.path::here()
setwd(mydir)


source("data_simulation.R")
source("vecchia.R")
source("aux_functions_joint.R")
source("update_parameters_joint.R")
source("mixed_workflow_joint.R")
source("predictive_scores.R")


prime_seeds <- c(
  2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
  31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
  73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229
)


reps <- 50
n.cores <- 50
cl <- makeCluster(n.cores)
registerDoParallel(cl)

results <- foreach(r = 1:reps, .packages = libraries) %dopar% {
  
  ### data generation
  
  set.seed(prime_seeds[r])
  
  p <- 3
  q <- 2
  
  true.beta <- matrix(c(1.0, -0.5,  3,  1.5, -1.2,  0.0), nrow = p, ncol = q, byrow = TRUE)
  true.Sigma <- matrix(c(9,0,0,3), nrow = q, ncol = q, byrow = TRUE)
  true.phi <- 0.1
  true.nu <- 0.5
  pred.prop <- 0.2
  
  family <- c("Binomial", "Gaussian")
  
  data <- sim.data(q = q, N = 1e2, 
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
  
  obs.locs <- data$obs.locs
  N.obs <- data$N.obs
  Y.obs <- data$Y.obs
  X.obs <- data$X.obs
  
  m <- 20
  
  phi <- data$true.phi
  nu <- data$true.nu
  
  
  obs.ord <- max_min(obs.locs) 
  # Assume max_min returns max-min order indices
  obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]  # Reorder locations based on max-min ordering
  distobs.ord <- rdist(obs.locs.ord)
  diameter <-  max(distobs.ord)
  b_phi <- diameter/3
  
  
  
  
  #### Prediction set up ####
  
  pred.locs <- data$pred.locs
  N.pred <- data$N.pred
  pred.ord <- max_min(pred.locs)
  ord <- c(obs.ord, pred.ord + nrow(obs.locs))
  locs.ord <- rbind(obs.locs, pred.locs)[ord, , drop = FALSE]
  distlocs.ord <- rdist(locs.ord)
  NNarray.all <- neighbor_matrix(locs.ord, m)
  distall.nn <- dist.nn(locs.ord, NNarray.all)
  coeff.all <- vecchia_coeff(distall.nn, NNarray.all, phi, nu)
  
  
  X.pred <- data$X.pred
  Y.pred.true <- data$Y.pred.true
  
  X.pred.ord <- X.pred[pred.ord, , drop = FALSE]
  Y.pred.ord.true <- Y.pred.true[pred.ord, , drop = FALSE]
  
  
  
  # Reorder Y.obs, X.obs, and W.obs accordingly
  Y.obs.ord <- Y.obs[obs.ord, , drop = FALSE]
  X.obs.ord <- X.obs[obs.ord, , drop = FALSE]
  
  
  # Sample initial parameters from MLE
  
  beta <- data$true.beta
  Sigma <- data$true.Sigma
  W.obs.ord <- data$true.W.obs[obs.ord, , drop = FALSE]
  
  
  tuning.phi <- 1e-2
  
  # Number of iterations
  niters <- 4e4
  
 
  
  out <- mixed.workflow.joint(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                              nu, W.obs.ord, beta, Sigma, phi,
                              M.prior, V.prior, S.prior, df.prior, b_phi,
                              niters, tuning.phi,
                              Y.pred.ord.true, X.pred.ord, N.pred,
                              distall.nn, coeff.all, NNarray.all,
                              family, family.par = NULL, link = NULL)
  
  out
  
}

stopCluster(cl)

# Save results to an .RData file
list.save(results, file = "bg_weak_joint_ind100.RData")