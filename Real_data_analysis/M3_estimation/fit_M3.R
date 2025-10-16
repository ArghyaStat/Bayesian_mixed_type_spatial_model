# Set cross-validation folds

rm(list = ls())


libraries <- c("spam", "this.path", "fields", "fBasics", "MCMCpack", "truncnorm", 
               "rlist", "foreach","doParallel", "FNN", "mcmcse")


# Load libraries
invisible(lapply(libraries, library, character.only = TRUE))

mydir <- this.path::here()
setwd(mydir)


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
Y_CNT <- log(1+data$CNT)
Y <- cbind(Y_BA, Y_CNT)
Y.sc <- cbind(Y_BA, data$CNT)

N <- nrow(Y)
q <- ncol(Y)


obs.locs <- locations

X.obs <- X
Y.obs <- Y
Y.sc.obs <- Y.sc

m <- 20
obs.ord <- max_min(obs.locs)
obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]
distobs.ord <- rdist(obs.locs.ord)
diameter <- max(distobs.ord)
b_phi <- diameter / log(10)
NNarray.obs <- neighbor_matrix(obs.locs.ord, m)
distobs.nn <- dist.nn(obs.locs.ord, neighbor_matrix = NNarray.obs)


#### Prediction set up ####


N.obs <- nrow(obs.locs)

phi <- b_phi/2
nu <- 0.3


# Reorder Y.obs, X.obs, and W.obs accordingly
Y.obs.ord <- Y.obs[obs.ord, , drop = FALSE]
Y.sc.obs.ord <- Y.sc.obs[obs.ord, , drop = FALSE]
X.obs.ord <- X.obs[obs.ord, , drop = FALSE]

# Prior
M.prior <- matrix(0, p, q)
V.prior <- 1e2 * diag(p)
S.prior <- diag(q)
df.prior <- q + 1
family <- c("Gaussian", "Gaussian")

# Initial params
beta <- M.prior
Sigma <- diag(c(1, 1/3))
W.obs.ord <- X.obs.ord %*% beta
chol.K.obs.ord.inv <- U.sgv(distobs.nn, NNarray.obs, phi, nu, m)


# Vecchia specs
tuning.phi <- b_phi/3
niters <- 5e4

fit.main <- mcmc.main(Y.obs.ord, Y.sc.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
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
log.like.sc.mean <- fit.main$log.like.sc.mean
rm(fit.main)

#### Summary of estimation ###

W.obs.ord.mean <- Reduce("+", W.ord.post) / niters

log_like_post.sc <- log.likelihood(W.obs.ord.mean, Y.sc.obs.ord, 
                                   family = c("Gaussian", "Poisson")) 


# summary_stats(W.ord.post, data$true.W.obs.ord)
beta.stats <- summary_stats(beta.post, data$true.beta)
Sigma.stats <- summary_stats(Sigma.post, data$true.Sigma)
phi.stats <- summary_stats(phi.post, data$true.phi)


#### ESS calculation for component chains ####

beta.ess <- compute_ess(beta.post)
Sigma.ess <- compute_ess(Sigma.post)
phi.ess <- compute_ess(phi.post)
W.obs.ord.ess <- compute_ess(W.ord.post)

mlpd.sc.obs <- mlpd(W.ord.post, Y.sc.obs.ord, 
                    family = c("Gaussian", "Poisson"))
waic.sc.obs <- waic(W.ord.post, Y.sc.obs.ord, mlpd.sc.obs, 
                    family = c("Gaussian", "Poisson"))


out <- list("beta.stats" = beta.stats,
            "Sigma.stats" = Sigma.stats, 
            "phi.stats" = phi.stats,
            "beta.ess" = beta.ess,
            "Sigma.ess" = Sigma.ess,
            "phi.ess" = phi.ess,
            "W.obs.ord.ess" = W.obs.ord.ess,
            "waic.obs" = waic.obs,
            "acc.phi" = acc.phi,
            "total_time" = total_time)



# Save results
save(out, file = "results_joint_Gaussian_est.RData")
